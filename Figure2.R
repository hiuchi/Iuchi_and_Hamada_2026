library(tidyverse)
library(furrr)
library(Seurat)
library(dyno)
library(dyntoy)
library(monocle3)
library(mclust)
library(slingshot)
library(scLS)
library(reticulate)
use_condaenv("scLS-env", required = TRUE)
library(tradeSeq)
library(Lamian)
library(scales)
library(irlba)
library(plotROC)
library(pROC)
library(GGally)
library(edgeR)
library(SingleCellExperiment)

theme_set(theme_bw(base_size = 20))
set.seed(8)

## ユーティリティ ----------------------------------------------------------
num.features <- 2000
cell.range   <- c(250, 500, 1000)
de.rate      <- 0.5
noise        <- 0.4

FQnorm <- function(counts){
  rk          <- apply(counts, 2, rank, ties.method = "min")
  counts.sort <- apply(counts, 2, sort)
  refdist     <- apply(counts.sort, 1, median)
  norm        <- apply(rk, 2, function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  norm
}

dir.create("./models",        showWarnings = FALSE)
dir.create("./plots",         showWarnings = FALSE)
dir.create("./plots/models",  showWarnings = FALSE)
dir.create("./plots/res",     showWarnings = FALSE)
dir.create("./res",           showWarnings = FALSE)

## データ生成 --------------------------------------------------------------
for(model in c("Linear", "Bifurcating")){
  for(num.cells in cell.range){
    set.seed(1)
    dataset <- dyntoy::generate_dataset(
      model                         = tolower(model),
      num_cells                     = num.cells,
      num_features                  = num.features,
      differentially_expressed_rate = de.rate
    ) %>%
      add_root(root_milestone_id = .$prior_information$start_milestones) %>%
      add_pseudotime()
    
    g <- plot_dimred(dataset, color_cells = "pseudotime")
    ggsave(
      g,
      file   = str_c("./plots/models/", model, "_", num.cells, "cells_pseudotime.pdf"),
      width  = 12,
      height = 12,
      dpi    = 300
    )
    save(dataset, file = str_c("./models/", model, "_", num.cells, "cells.RData"))
  }
}

## 解析本体 --------------------------------------------------------------
for(model in c("Linear", "Bifurcating")){
  pt.df <- tibble(
    pt.method  = character(),
    cell.num   = integer(),
    cell       = character(),
    pseudotime = numeric()
  )
  res <- tibble(
    gene      = character(),
    answer    = logical(),
    pt.method = character(),
    method    = character(),
    cell.num  = integer(),
    p         = numeric()
  )
  
  for(num.cells in cell.range){
    load(str_c("./models/", model, "_", num.cells, "cells.RData"))
    
    ## counts / expression
    counts <- t(as.matrix(dataset$counts))
    counts <- counts * matrix(
      rnorm(length(counts), mean = 1, sd = noise),
      nrow = nrow(counts),
      ncol = ncol(counts)
    )
    counts[counts < 0] <- 0
    expression <- edgeR::cpm(counts)
    
    ## 遺伝子ごとの DE ラベル（名前付きベクトル）
    de_vec <- dataset[["tde_overall"]][["differentially_expressed"]]
    names(de_vec) <- dataset$feature_info$feature_id
    
    ## SCE（tradeSeq 用）
    sce <- SingleCellExperiment(
      assays = List(
        counts = counts,
        log2   = log2(counts + 1e-8)
      )
    )
    
    ## Seurat（scLS 用）
    seurat <- CreateSeuratObject(counts = counts)
    seurat <- NormalizeData(seurat)
    VariableFeatures(seurat) <- rownames(seurat)
    
    ## -------- Ground truth --------
    pt_gt_raw <- dataset$pseudotime
    
    pt.df <- pt.df %>%
      add_row(
        pt.method  = "Ground_truth",
        cell.num   = num.cells,
        cell       = colnames(counts),
        pseudotime = pt_gt_raw
      )
    
    seurat <- AddMetaData(
      seurat,
      metadata = scales::rescale(pt_gt_raw, to = c(0, 1)),
      col.name = "pseudotime.gt"
    )
    
    ## scLS (GT)
    scLS.res <- scLS.dynamic(
      object     = seurat,
      time.col   = "pseudotime.gt",
      assay      = "RNA",
      slot       = "data",
      f.min      = 0.1,
      f.max      = 1,
      n.bins     = 500,
      fap.method = "davies",
      n.cores    = 8
    )
    
    ## tradeSeq (GT, Linear と Bifurcatingで同じ扱い)
    cellWeights <- rep(1, num.cells)
    sce_ts <- fitGAM(
      counts      = counts,
      pseudotime  = pt_gt_raw,
      cellWeights = cellWeights,
      nknots      = 6,
      verbose     = FALSE
    )
    sveRes  <- startVsEndTest(sce_ts)
    assoRes <- associationTest(sce_ts)
    
    ## Lamian (GT)
    if(num.cells %% 2 == 0){
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(c("WT1", "WT2"), each = num.cells / 2)
      )
    } else {
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(
          c("WT1", "WT2"),
          c(round(num.cells / 2), round(num.cells / 2) - 1)
        )
      )
    }
    .design <- matrix(
      c(1, 1, 0, 1),
      ncol    = 2,
      dimnames = list(c("WT1", "WT2"), c("intercept", "group"))
    )
    
    lamian.res <- try(
      lamian_test(
        expr       = log2(expression + 1),
        cellanno   = .cellanno,
        pseudotime = pt_gt_raw,
        design     = .design,
        test.type  = "time",
        permuiter  = 100,
        ncores     = 1
      ),
      silent = TRUE
    )
    if(inherits(lamian.res, "try-error")) lamian.res <- NULL
    
    ## Ground truth 用の res ブロック
    genes_ls   <- scLS.res$Feature
    ans_ls     <- de_vec[genes_ls]
    p_ls       <- scLS.res$PeakFAP
    
    genes_sve  <- rownames(sveRes)
    ans_sve    <- de_vec[genes_sve]
    p_sve      <- sveRes$pvalue
    
    genes_asso <- rownames(assoRes)
    ans_asso   <- de_vec[genes_asso]
    p_asso     <- assoRes$pvalue
    
    if(!is.null(lamian.res)){
      genes_lam <- rownames(lamian.res$statistics)
      ans_lam   <- de_vec[genes_lam]
      p_lam     <- lamian.res$statistics$pval.overall
      
      block_gt <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso, genes_lam),
        answer    = c(ans_ls,   ans_sve,   ans_asso,   ans_lam),
        pt.method = "Ground_truth",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso)),
          rep("Lamian",              length(genes_lam))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso, p_lam)
      )
    } else {
      block_gt <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso),
        answer    = c(ans_ls,   ans_sve,   ans_asso),
        pt.method = "Ground_truth",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso)
      )
    }
    
    res <- bind_rows(res, block_gt)
    message(str_c(model, ", ", num.cells, ", Ground truth is done."))
    
    ## -------- monocle3 --------
    gene.metadata <- data.frame(
      id              = dataset$feature_info$feature_id,
      gene_short_name = dataset$feature_info$feature_id
    )
    rownames(gene.metadata) <- gene.metadata$id
    
    monocle3.res <- new_cell_data_set(
      t(dataset$counts),
      gene_metadata = gene.metadata
    ) %>%
      preprocess_cds(num_dim = 50) %>%
      reduce_dimension(
        preprocess_method = "PCA",
        reduction_method  = "UMAP"
      ) %>%
      cluster_cells() %>%
      learn_graph() %>%
      order_cells(root_cells = dataset$prior_information$start_id)
    
    pt_m3 <- monocle3.res@principal_graph_aux@listData$UMAP$pseudotime
    finite_max <- max(pt_m3[is.finite(pt_m3)], na.rm = TRUE)
    pt_m3[!is.finite(pt_m3)] <- finite_max
    
    pt.df <- pt.df %>%
      add_row(
        pt.method  = "monocle3",
        cell.num   = num.cells,
        cell       = colnames(counts),
        pseudotime = pt_m3
      )
    
    seurat <- AddMetaData(
      seurat,
      metadata = scales::rescale(pt_m3, to = c(0, 1)),
      col.name = "pseudotime.monocle3"
    )
    
    scLS.res <- scLS.dynamic(
      object     = seurat,
      time.col   = "pseudotime.monocle3",
      assay      = "RNA",
      slot       = "data",
      f.min      = 0.1,
      f.max      = 1,
      n.bins     = 500,
      fap.method = "davies",
      n.cores    = 8
    )
    
    if(model == "Linear"){
      keep <- is.finite(pt_m3)
      .counts <- counts[, keep, drop = FALSE]
      .pt     <- pt_m3[keep]
      .w      <- rep(1, length(.pt))
      sce_ts <- fitGAM(
        counts      = .counts,
        pseudotime  = .pt,
        cellWeights = .w,
        nknots      = 6,
        verbose     = FALSE
      )
    } else {
      .w <- rep(1, length(pt_m3))
      sce_ts <- fitGAM(
        counts      = counts,
        pseudotime  = pt_m3,
        cellWeights = .w,
        nknots      = 6,
        verbose     = FALSE
      )
    }
    sveRes  <- startVsEndTest(sce_ts)
    assoRes <- associationTest(sce_ts)
    
    if(num.cells %% 2 == 0){
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(c("WT1", "WT2"), each = num.cells / 2)
      )
    } else {
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(
          c("WT1", "WT2"),
          c(round(num.cells / 2), round(num.cells / 2) - 1)
        )
      )
    }
    .design <- matrix(
      c(1, 1, 0, 1),
      ncol    = 2,
      dimnames = list(c("WT1", "WT2"), c("intercept", "group"))
    )
    
    lamian.res <- try(
      lamian_test(
        expr       = log2(expression + 1),
        cellanno   = .cellanno,
        pseudotime = pt_m3,
        design     = .design,
        test.type  = "time",
        permuiter  = 100,
        ncores     = 1
      ),
      silent = TRUE
    )
    if(inherits(lamian.res, "try-error")) lamian.res <- NULL
    
    genes_ls   <- scLS.res$Feature
    ans_ls     <- de_vec[genes_ls]
    p_ls       <- scLS.res$PeakFAP
    
    genes_sve  <- rownames(sveRes)
    ans_sve    <- de_vec[genes_sve]
    p_sve      <- sveRes$pvalue
    
    genes_asso <- rownames(assoRes)
    ans_asso   <- de_vec[genes_asso]
    p_asso     <- assoRes$pvalue
    
    if(!is.null(lamian.res)){
      genes_lam <- rownames(lamian.res$statistics)
      ans_lam   <- de_vec[genes_lam]
      p_lam     <- lamian.res$statistics$pval.overall
      
      block_m3 <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso, genes_lam),
        answer    = c(ans_ls,   ans_sve,   ans_asso,   ans_lam),
        pt.method = "monocle3",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso)),
          rep("Lamian",              length(genes_lam))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso, p_lam)
      )
    } else {
      block_m3 <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso),
        answer    = c(ans_ls,   ans_sve,   ans_asso),
        pt.method = "monocle3",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso)
      )
    }
    
    res <- bind_rows(res, block_m3)
    message(str_c(model, ", ", num.cells, ", monocle3 is done."))
    
    ## -------- slingshot --------
    sce_sl <- SingleCellExperiment(
      assays = List(
        counts = counts,
        log2   = log2(counts + 1e-8)
      )
    )
    assays(sce_sl)$norm <- FQnorm(assays(sce_sl)$counts)
    
    pca <- prcomp(t(log1p(assays(sce_sl)$norm)), scale. = FALSE)
    reducedDims(sce_sl) <- SimpleList(PCA = pca$x[, 1:2])
    
    colData(sce_sl)$GMM <- Mclust(pca$x[, 1:2])$classification
    sce_sl <- slingshot(
      sce_sl,
      clusterLabels = "GMM",
      reducedDim    = "PCA",
      start.clus    = colData(sce_sl)$GMM[
        which(names(colData(sce_sl)$GMM) == dataset$prior_information$start_id)
      ]
    )
    
    if(model == "Linear"){
      pt_sl_raw <- as.numeric(sce_sl$slingPseudotime_1)
      pt_sl     <- pt_sl_raw / max(pt_sl_raw, na.rm = TRUE)
      names(pt_sl) <- colnames(sce_sl)
    } else {
      pt_mat <- cbind(sce_sl$slingPseudotime_1, sce_sl$slingPseudotime_2)
      pt_mat[is.na(pt_mat)] <- 0
      w_mat <- slingCurveWeights(sce_sl)
      best  <- max.col(w_mat, ties.method = "first")
      
      pt_sl_raw <- pt_mat[cbind(seq_len(nrow(pt_mat)), best)]
      pt_sl <- (pt_sl_raw - min(pt_sl_raw, na.rm = TRUE)) /
        (max(pt_sl_raw, na.rm = TRUE) - min(pt_sl_raw, na.rm = TRUE))
      names(pt_sl) <- colnames(sce_sl)
    }
    
    pt.df <- pt.df %>%
      add_row(
        pt.method  = "slingshot",
        cell.num   = num.cells,
        cell       = colnames(counts),
        pseudotime = pt_sl
      )
    
    seurat <- AddMetaData(
      seurat,
      metadata = pt_sl,
      col.name = "pseudotime.slingshot"
    )
    
    scLS.res <- scLS.dynamic(
      object     = seurat,
      time.col   = "pseudotime.slingshot",
      assay      = "RNA",
      slot       = "data",
      f.min      = 0.1,
      f.max      = 1,
      n.bins     = 500,
      fap.method = "davies",
      n.cores    = 8
    )
    
    if(model == "Linear"){
      pseudotime_ts <- matrix(sce_sl$slingPseudotime_1, ncol = 1)
      weights_ts    <- matrix(1, nrow = nrow(pseudotime_ts), ncol = 1)
    } else {
      pseudotime_ts <- cbind(
        sce_sl$slingPseudotime_1,
        sce_sl$slingPseudotime_2
      )
      pseudotime_ts[is.na(pseudotime_ts)] <- 0
      weights_ts <- slingCurveWeights(sce_sl)
    }
    
    sce_ts <- fitGAM(
      counts      = counts,
      pseudotime  = pseudotime_ts,
      cellWeights = weights_ts,
      nknots      = 6,
      verbose     = FALSE
    )
    sveRes  <- startVsEndTest(sce_ts)
    assoRes <- associationTest(sce_ts)
    
    if(num.cells %% 2 == 0){
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(c("WT1", "WT2"), each = num.cells / 2)
      )
    } else {
      .cellanno <- data.frame(
        Cell   = colnames(expression),
        Sample = rep(
          c("WT1", "WT2"),
          c(round(num.cells / 2), round(num.cells / 2) - 1)
        )
      )
    }
    .design <- matrix(
      c(1, 1, 0, 1),
      ncol    = 2,
      dimnames = list(c("WT1", "WT2"), c("intercept", "group"))
    )
    
    lamian.res <- try(
      lamian_test(
        expr       = log2(expression + 1),
        cellanno   = .cellanno,
        pseudotime = pt_sl,
        design     = .design,
        test.type  = "time",
        permuiter  = 100,
        ncores     = 1
      ),
      silent = TRUE
    )
    if(inherits(lamian.res, "try-error")) lamian.res <- NULL
    
    genes_ls   <- scLS.res$Feature
    ans_ls     <- de_vec[genes_ls]
    p_ls       <- scLS.res$PeakFAP
    
    genes_sve  <- rownames(sveRes)
    ans_sve    <- de_vec[genes_sve]
    p_sve      <- sveRes$pvalue
    
    genes_asso <- rownames(assoRes)
    ans_asso   <- de_vec[genes_asso]
    p_asso     <- assoRes$pvalue
    
    if(!is.null(lamian.res)){
      genes_lam <- rownames(lamian.res$statistics)
      ans_lam   <- de_vec[genes_lam]
      p_lam     <- lamian.res$statistics$pval.overall
      
      block_sl <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso, genes_lam),
        answer    = c(ans_ls,   ans_sve,   ans_asso,   ans_lam),
        pt.method = "slingshot",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso)),
          rep("Lamian",              length(genes_lam))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso, p_lam)
      )
    } else {
      block_sl <- tibble(
        gene      = c(genes_ls, genes_sve, genes_asso),
        answer    = c(ans_ls,   ans_sve,   ans_asso),
        pt.method = "slingshot",
        method    = c(
          rep("LS",                  length(genes_ls)),
          rep("tradeSeq_startVsEnd", length(genes_sve)),
          rep("tradeSeq_asso",       length(genes_asso))
        ),
        cell.num  = num.cells,
        p         = c(p_ls, p_sve, p_asso)
      )
    }
    
    res <- bind_rows(res, block_sl)
    message(str_c(model, ", ", num.cells, ", slingshot is done."))
  }
  
  res$method <- factor(
    res$method,
    levels = c("LS", "Lamian", "tradeSeq_startVsEnd", "tradeSeq_asso")
  )
  
  save(res,   file = str_c("./res/res_", model, ".RData"))
  save(pt.df, file = str_c("./res/pt_",  model, ".RData"))
}

## 可視化と AUC（method × cell.num） -------------------------------
for(model in c("Linear", "Bifurcating")){
  load(str_c("./res/res_", model, ".RData"))
  
  res_plot <- res %>%
    dplyr::filter(!is.na(answer), is.finite(p))
  
  g <- ggplot(res_plot, aes(d = answer, m = 1 - p, colour = method)) +
    geom_roc(n.cuts = FALSE, linealpha = 0.8) +
    ggtitle(model) +
    xlab("False positive rate") +
    ylab("True positive rate") +
    facet_grid(pt.method ~ cell.num) +
    theme(legend.position = "none")
  ggsave(
    g,
    file   = str_c("./plots/res/ROC_", model, ".pdf"),
    width  = 12,
    height = 12,
    dpi    = 300
  )
  
  g <- ggplot(res_plot, aes(d = answer, m = 1 - p, colour = method)) +
    geom_roc(n.cuts = FALSE, linealpha = 0.8)
  ggsave(g, file = "./plots/res/legend_.pdf", dpi = 300)
  
  auc <- tibble(method = character(), cell.num = character(), AUC = numeric())
  for(i in c("LS", "tradeSeq_startVsEnd", "tradeSeq_asso", "Lamian")){
    for(j in c(250, 500, 1000)){
      df_ij <- res_plot %>% dplyr::filter(method == i, cell.num == j)
      if(nrow(df_ij) == 0L) next
      
      r_obj <- try(
        pROC::roc(response = df_ij$answer,
                  predictor = -df_ij$p,   # p が小さいほど陽性なので符号を反転
                  quiet = TRUE),
        silent = TRUE
      )
      if(inherits(r_obj, "try-error")) next
      
      this_auc <- as.numeric(pROC::auc(r_obj))
      if(is.finite(this_auc)){
        auc <- auc %>%
          add_row(
            method   = i,
            cell.num = str_c(j, " cells"),
            AUC      = this_auc
          )
      }
    }
  }
  
  auc <- na.omit(auc)
  if(nrow(auc) == 0L){
    message("No valid AUC entries for model = ", model, "; skip barplot/CSV.")
    next
  }
  
  auc$cell.num <- factor(
    auc$cell.num,
    levels = c("250 cells", "500 cells", "1000 cells")
  )
  auc$method <- factor(
    auc$method,
    levels = c("LS", "tradeSeq_startVsEnd", "tradeSeq_asso", "Lamian")
  )
  
  g2 <- ggplot(auc, aes(x = method, y = AUC)) +
    geom_bar(stat = "identity") +
    facet_wrap(. ~ cell.num) +
    coord_cartesian(ylim = c(0.9, 1.0)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )
  ggsave(
    g2,
    file   = str_c("./plots/res/AUC_", model, ".pdf"),
    width  = 12,
    height = 12,
    dpi    = 300
  )
  
  write_csv(auc, file = str_c("./res/AUC_", model, ".csv"))
}


## pt.method × method の AUC テーブル ------------------------------
for(model in c("Linear", "Bifurcating")){
  load(str_c("./res/res_", model, ".RData"))
  
  if(!all(c("answer", "p", "pt.method", "method", "cell.num") %in% colnames(res))){
    message("res for model = ", model, " does not have required columns; skip.")
    next
  }
  
  res_plot <- res %>%
    dplyr::filter(!is.na(answer), is.finite(p))
  
  if(nrow(res_plot) == 0L){
    message("No usable res for model = ", model, "; skip pt.method×method AUC.")
    next
  }
  
  auc_tbl <- res_plot %>%
    dplyr::mutate(pred = -p) %>%
    dplyr::group_by(cell.num, pt.method, method) %>%
    dplyr::summarise(
      AUC = {
        r_obj <- try(
          pROC::roc(response   = answer,
                    predictor  = pred,
                    quiet      = TRUE),
          silent = TRUE
        )
        if(inherits(r_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(r_obj))
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(AUC))
  
  if(nrow(auc_tbl) == 0L){
    message("Empty AUC table for model = ", model, "; no CSV written.")
    next
  }
  
  out_tbl <- auc_tbl %>%
    dplyr::mutate(AUC = round(AUC, 3)) %>%
    tidyr::pivot_wider(names_from = method, values_from = AUC) %>%
    dplyr::arrange(cell.num, pt.method)
  
  write_csv(out_tbl, file = str_c(model, ".csv"))
}

for(model in c("linear", "bifurcating")){
  load(str_c("./res/res_", model, ".RData"))
  load(str_c("res/pt_", model, ".RData"))
  d <- rbind(inner_join(pt.df %>% dplyr::filter(pt.method == "Ground_truth"),
                        pt.df %>% dplyr::filter(pt.method == "monocle3"), by = c("cell.num", "cell")),
             inner_join(pt.df %>% dplyr::filter(pt.method == "Ground_truth"),
                        pt.df %>% dplyr::filter(pt.method == "slingshot"), by = c("cell.num", "cell"))) %>%
    dplyr::mutate(pt.method.y = dplyr::recode(pt.method.y,
                                              monocle3 = "Monocle3")) %>%
    dplyr::mutate(pt.method.y = dplyr::recode(pt.method.y,
                                              slingshot = "Slingshot"))
  
  g <- ggplot(d, aes(x = pseudotime.x, y = pseudotime.y))
  g <- g + geom_point(size = 0.1) + facet_grid(pt.method.y ~ cell.num, scales = "free")
  g <- g + xlab("Pseudotime of groud truth") + ylab("Predicted pseudotime")
  if(model == "linear") {
    g <- g + ggtitle("Linear")
  }else{
    g <- g + ggtitle("Bifurcating")
  }
  g
  ggsave(g, file = str_c("./plots/res/cor_", model, ".pdf"), width = 12, height = 8, dpi = 300)
}
