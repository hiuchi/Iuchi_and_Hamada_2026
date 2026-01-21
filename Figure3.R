library(tidyverse)
library(Seurat)
library(patchwork)
library(dyno)
library(dyngen)
library(igraph)
library(scLS)
library(reticulate)
use_condaenv("scLS-env", required = TRUE)
library(monocle3)
library(mclust)
library(slingshot)
library(tradeSeq)
library(Lamian)
library(plotROC)
library(pROC)
library(SingleCellExperiment)
library(future)
theme_set(theme_bw(base_size = 20))

dir.create("./files", showWarnings = FALSE)
dir.create("./res",   showWarnings = FALSE)
dir.create("./cor",   showWarnings = FALSE)
dir.create("./plots", showWarnings = FALSE)

num_targets <- 250
num_hks     <- 250
ko_gene     <- "C1_TF1"
ko_rate     <- 0
num_cores   <- 10
cell_range  <- c(250, 500, 1000)

result_tbl <- tibble(
  backbone  = character(),
  gene      = character(),
  pt_method = character(),
  method    = character(),
  answer    = integer(),
  num_cells = integer(),
  p         = double()
)

pseudotime_tbl <- tibble(
  backbone   = character(),
  pt_method  = character(),
  sample     = character(),
  cell_num   = integer(),
  cell       = character(),
  pseudotime = double()
)

# parameters and functions ------------------------------------------------

fq_norm <- function(counts) {
  rank_mat       <- apply(counts, 2, rank, ties.method = "min")
  sorted_counts  <- apply(counts, 2, sort)
  reference_dist <- apply(sorted_counts, 1, median)
  normalized     <- apply(rank_mat, 2, function(r) reference_dist[r])
  rownames(normalized) <- rownames(counts)
  normalized
}

# data generation ---------------------------------------------------------

for (model in c("Linear", "Bifurcating")) {
  set.seed(10)
  
  backbone <- if (model == "Linear") {
    backbone_linear()
  } else {
    backbone_bifurcating()
  }
  
  config <- initialise_model(
    backbone        = backbone,
    num_cells       = last(cell_range),
    num_tfs         = nrow(backbone$module_info),
    num_targets     = num_targets,
    num_hks         = num_hks,
    num_cores       = num_cores,
    simulation_params = simulation_default(
      census_interval  = 10,
      ssa_algorithm    = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 100)
    )
  )
  
  ggsave(
    plot_backbone_modulenet(config),
    file = str_c("./files/", model, ".pdf"),
    dpi  = 300
  )
  
  model_common <- config %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard()
  
  model_wt <- model_common %>% generate_cells()
  model_ko <- model_common
  
  model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = 100,
    timepoint       = 0,
    genes           = ko_gene,
    num_genes       = 1,
    multiplier      = ko_rate
  )
  
  model_ko <- model_ko %>% generate_cells()
  
  wt <- model_wt %>%
    generate_experiment() %>%
    as_dyno() %>%
    add_root(root_milestone_id = "sA") %>%
    add_pseudotime(pseudotime = NULL)
  
  ko <- model_ko %>%
    generate_experiment() %>%
    as_dyno() %>%
    add_root(root_milestone_id = "sA") %>%
    add_pseudotime(pseudotime = NULL)
  
  wt_counts_full <- wt$counts %>% as.matrix() %>% t()
  ko_counts_full <- ko$counts %>% as.matrix() %>% t()
  counts_full    <- cbind(wt_counts_full, ko_counts_full)
  
  edges <- model_wt$feature_network %>%
    dplyr::select(from, to) %>%
    dplyr::distinct()
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  downstream_idx   <- subcomponent(g, v = ko_gene, mode = "out")
  downstream_genes <- igraph::V(g)$name[downstream_idx]
  downstream_genes <- setdiff(downstream_genes, ko_gene)
  downstream_genes <- intersect(downstream_genes, wt$feature_ids)
  
  de_genes <- c(ko_gene, downstream_genes)
  
  answer <- rep(0L, nrow(counts_full))
  for (i in seq_len(nrow(counts_full))) {
    if (wt$feature_ids[i] %in% de_genes) {
      answer[i] <- 1L
    }
  }
  
  if (model == "Linear") {
    save(wt,     file = "./files/Linear_wt.RData")
    save(ko,     file = "./files/Linear_ko.RData")
    save(answer, file = "./files/answer_Linear.RData")
  } else {
    save(wt,     file = "./files/Bifurcating_wt.RData")
    save(ko,     file = "./files/Bifurcating_ko.RData")
    save(answer, file = "./files/answer_Bifurcating.RData")
  }
}

# inference / DE along pseudotime ----------------------------------------

for (model in c("Linear", "Bifurcating")) {
  load(str_c("./files/", model, "_wt.RData"))
  load(str_c("./files/", model, "_ko.RData"))
  load(str_c("./files/answer_", model, ".RData"))
  
  for (num_cells in cell_range) {
    message("model = ", model, ", num_cells = ", num_cells)
    
    set.seed(2)
    sampling_indices <- sample(1:1000, num_cells) %>% sort()
    
    wt_expression <- wt$expression[sampling_indices, ] %>%
      as.matrix() %>%
      t()
    wt_counts <- wt$counts[sampling_indices, ] %>%
      as.matrix() %>%
      t()
    wt_pseudotime <- wt$pseudotime[sampling_indices]
    
    ko_expression <- ko$expression[sampling_indices, ] %>%
      as.matrix() %>%
      t()
    ko_counts <- ko$counts[sampling_indices, ] %>%
      as.matrix() %>%
      t()
    ko_pseudotime <- ko$pseudotime[sampling_indices]
    
    counts <- cbind(wt_counts, ko_counts)
    colnames(counts) <- c(
      str_c("WT_", colnames(wt_counts)),
      str_c("KO_", colnames(ko_counts))
    )
    
    expression_mat <- cbind(wt_expression, ko_expression)
    colnames(expression_mat) <- colnames(counts)
    
    genes_all <- rownames(counts)
    
    # ground-truth pseudotime --------------------------------------------
    
    pseudotime_tbl <- pseudotime_tbl %>%
      add_row(
        backbone   = model,
        pt_method  = "Ground_truth",
        sample     = rep(c("WT", "KO"), each = num_cells),
        cell_num   = num_cells,
        cell       = colnames(counts),
        pseudotime = c(wt_pseudotime, ko_pseudotime)
      )
    
    plot_df <- tibble(
      component1 = wt$dimred[sampling_indices, 1],
      component2 = wt$dimred[sampling_indices, 2],
      pseudotime = wt_pseudotime / max(wt_pseudotime)
    )
    
    g_wt <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      ggtitle(str_c("Ground truth, ", num_cells, " cells")) +
      labs(subtitle = "WT") +
      coord_cartesian() +
      theme(legend.position = "none")
    
    plot_df <- tibble(
      component1 = ko$dimred[sampling_indices, 1],
      component2 = ko$dimred[sampling_indices, 2],
      pseudotime = ko_pseudotime / max(ko_pseudotime)
    )
    
    g_ko <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      labs(subtitle = "KO") +
      coord_cartesian() +
      theme(legend.position = "bottom")
    
    combined_plot <- g_wt + g_ko
    
    ggsave(
      combined_plot,
      file = str_c("./plots/Groundtruth_", model, "_", num_cells, "cells.pdf"),
      dpi  = 300
    )
    
    # scLS (ground truth) ----------------------------------------
    
    wt_seu <- CreateSeuratObject(counts = wt_counts) %>%
      AddMetaData(metadata = wt_pseudotime / max(wt_pseudotime), col.name = "Pseudotime") %>%
      NormalizeData() %>%
      FindVariableFeatures()
    ko_seu <- CreateSeuratObject(counts = ko_counts) %>%
      AddMetaData(metadata = ko_pseudotime / max(ko_pseudotime), col.name = "Pseudotime") %>%
      NormalizeData() %>%
      FindVariableFeatures()
    
    res_scLS <- scLS.shift(
      wt_seu, ko_seu,
      time.col1 = "Pseudotime", time.col2 = "Pseudotime",
      n.cores   = num_cores
    )
    
    # Lamian (ground truth) ----------------------------------------------
    
    expression_mat <- cbind(wt_expression, ko_expression)
    colnames(expression_mat) <- c(
      str_c("WT_", colnames(wt_expression)),
      str_c("KO_", colnames(ko_expression))
    )
    pseudotime_vec <- c(wt_pseudotime, ko_pseudotime)
    names(pseudotime_vec) <- colnames(expression_mat)
    cell_anno <- data.frame(
      Cell   = colnames(expression_mat),
      Sample = c(rep("WT", ncol(wt_expression)),
                 rep("KO", ncol(ko_expression)))
    )
    design_mat <- matrix(
      c(1, 1,
        0, 1),
      ncol = 2,
      dimnames = list(
        c("WT", "KO"),
        c("intercept", "group")
      )
    )
    
    res_lamian <- lamian_test(
      expr         = expression_mat,
      cellanno     = cell_anno,
      pseudotime   = pseudotime_vec,
      design       = design_mat,
      test.type    = "variable",
      overall.only = TRUE,
      test.method  = "chisq",
      permuiter    = 100,
      ncores       = num_cores
    )
    
    # tradeSeq (ground truth) --------------------------------------------
    
    colnames(counts) <- colnames(expression_mat)
    
    pseudotime_mat <- matrix(
      c(
        wt_pseudotime,
        rep(0, num_cells),
        rep(0, num_cells),
        ko_pseudotime
      ),
      ncol = 2
    )
    
    cell_weights <- matrix(
      c(
        rep(1, num_cells),
        rep(0, num_cells),
        rep(0, num_cells),
        rep(1, num_cells)
      ),
      ncol = 2
    )
    
    gam <- fitGAM(
      counts      = counts,
      pseudotime  = pseudotime_mat,
      cellWeights = cell_weights,
      nknots      = 6,
      parallel    = TRUE
    )
    
    res_tradeseq_diffEnd <- diffEndTest(gam)
    res_tradeseq_pattern <- patternTest(gam)
    res_tradeseq_earlyDE <- earlyDETest(gam, knots = c(1, 3))
    
    ## ====== Ground truth 用 result_tbl の整形（順序合わせ） ====== ##
    # scLS
    scLS_tbl_gt <- res_scLS %>%
      mutate(
        backbone  = model,
        pt_method = "Ground_truth",
        method    = "scLS",
        num_cells = num_cells,
        answer    = answer[match(gene, genes_all)]
      ) %>%
      select(backbone, gene, pt_method, method, answer, num_cells, p)
    
    # Lamian（match を使わず expr の行順をそのまま利用）
    lamian_tbl_gt <- tibble(
      backbone  = model,
      gene      = genes_all,
      pt_method = "Ground_truth",
      method    = "Lamian",
      answer    = answer,
      num_cells = num_cells,
      p         = as.numeric(res_lamian$statistics$pval.chisq.overall)
    )
    
    # tradeSeq
    tradeseq_diffEnd_tbl_gt <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_diffEnd),
      pt_method = "Ground_truth",
      method    = "tradeSeq (diffEnd)",
      answer    = answer[match(rownames(res_tradeseq_diffEnd), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_diffEnd$pvalue
    )
    
    tradeseq_pattern_tbl_gt <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_pattern),
      pt_method = "Ground_truth",
      method    = "tradeSeq (pattern)",
      answer    = answer[match(rownames(res_tradeseq_pattern), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_pattern$pvalue
    )
    
    tradeseq_earlyDE_tbl_gt <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_earlyDE),
      pt_method = "Ground_truth",
      method    = "tradeSeq (earlyDE)",
      answer    = answer[match(rownames(res_tradeseq_earlyDE), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_earlyDE$pvalue
    )
    
    result_tbl <- bind_rows(
      result_tbl,
      scLS_tbl_gt,
      lamian_tbl_gt,
      tradeseq_diffEnd_tbl_gt,
      tradeseq_pattern_tbl_gt,
      tradeseq_earlyDE_tbl_gt
    )
    ## ====== Ground truth 部分ここまで ====== ##
    
    # Monocle3 ------------------------------------------------------------
    
    gene_metadata <- data.frame(
      id              = wt$feature_info$feature_id,
      gene_short_name = wt$feature_info$feature_id
    )
    rownames(gene_metadata) <- gene_metadata$id
    
    root_cells_wt <- names(wt_pseudotime)[which.min(wt_pseudotime)]
    monocle3_wt <- new_cell_data_set(
      wt_counts,
      gene_metadata = gene_metadata
    ) %>%
      preprocess_cds(num_dim = 50) %>%
      reduce_dimension(preprocess_method = "PCA", reduction_method = "UMAP") %>%
      cluster_cells() %>%
      learn_graph() %>%
      order_cells(root_cells = root_cells_wt)
    
    wt_pseudotime_monocle3 <- monocle3_wt@principal_graph_aux@listData$UMAP$pseudotime
    wt_pseudotime_monocle3[is.infinite(wt_pseudotime_monocle3)] <-
      max(wt_pseudotime_monocle3[!is.infinite(monocle3_wt@principal_graph_aux@listData$UMAP$pseudotime)])
    wt_pseudotime_monocle3 <- wt_pseudotime_monocle3 / max(wt_pseudotime_monocle3)
    
    root_cells_ko <- names(ko_pseudotime)[which.min(ko_pseudotime)]
    monocle3_ko <- new_cell_data_set(
      ko_counts,
      gene_metadata = gene_metadata
    ) %>%
      preprocess_cds(num_dim = 50) %>%
      reduce_dimension(preprocess_method = "PCA", reduction_method = "UMAP") %>%
      cluster_cells() %>%
      learn_graph() %>%
      order_cells(root_cells = names(which.min(ko_pseudotime)))
    
    ko_pseudotime_monocle3 <- monocle3_ko@principal_graph_aux@listData$UMAP$pseudotime
    ko_pseudotime_monocle3[is.infinite(ko_pseudotime_monocle3)] <-
      max(ko_pseudotime_monocle3[!is.infinite(monocle3_ko@principal_graph_aux@listData$UMAP$pseudotime)])
    ko_pseudotime_monocle3 <- ko_pseudotime_monocle3 / max(ko_pseudotime_monocle3)
    
    pseudotime_tbl <- pseudotime_tbl %>%
      add_row(
        backbone   = model,
        pt_method  = "Monocle3",
        sample     = rep(c("WT", "KO"), each = num_cells),
        cell_num   = num_cells,
        cell       = colnames(counts),
        pseudotime = c(wt_pseudotime_monocle3, ko_pseudotime_monocle3)
      )
    
    plot_df <- tibble(
      component1 = wt$dimred[sampling_indices, 1],
      component2 = wt$dimred[sampling_indices, 2],
      pseudotime = wt_pseudotime_monocle3
    )
    
    g_wt <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      ggtitle(str_c("Monocle3, ", num_cells, " cells")) +
      labs(subtitle = "WT") +
      coord_cartesian() +
      theme(legend.position = "none")
    
    plot_df <- tibble(
      component1 = ko$dimred[sampling_indices, 1],
      component2 = ko$dimred[sampling_indices, 2],
      pseudotime = ko_pseudotime_monocle3
    )
    
    g_ko <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      labs(subtitle = "KO") +
      coord_cartesian() +
      theme(legend.position = "bottom")
    
    combined_plot <- g_wt + g_ko
    
    ggsave(
      combined_plot,
      file = str_c("./plots/monocle3_", model, "_", num_cells, "cells.pdf"),
      dpi  = 300
    )
    
    # scLS (Monocle3) --------------------------------------------
    
    wt_seu <- CreateSeuratObject(counts = wt_counts) %>%
      AddMetaData(metadata = wt_pseudotime_monocle3, col.name = "Pseudotime") %>%
      NormalizeData() %>%
      FindVariableFeatures()
    ko_seu <- CreateSeuratObject(counts = ko_counts) %>%
      AddMetaData(metadata = ko_pseudotime_monocle3, col.name = "Pseudotime") %>%
      NormalizeData() %>%
      FindVariableFeatures()
    
    res_scLS <- scLS.shift(
      wt_seu, ko_seu,
      time.col1 = "Pseudotime", time.col2 = "Pseudotime",
      n.cores   = num_cores
    )
    
    # Lamian (Monocle3) ---------------------------------------------------
    
    pseudotime_vec <- c(wt_pseudotime_monocle3, ko_pseudotime_monocle3)
    names(pseudotime_vec) <- colnames(expression_mat)
    
    res_lamian <- lamian_test(
      expr         = expression_mat,
      cellanno     = cell_anno,
      pseudotime   = pseudotime_vec,
      design       = design_mat,
      test.type    = "variable",
      overall.only = TRUE,
      test.method  = "chisq",
      permuiter    = 100,
      ncores       = num_cores
    )
    
    # tradeSeq (Monocle3) -------------------------------------------------
    
    colnames(counts) <- colnames(expression_mat)
    
    pseudotime_mat <- matrix(
      c(
        wt_pseudotime_monocle3,
        rep(0, num_cells),
        rep(0, num_cells),
        ko_pseudotime_monocle3
      ),
      ncol = 2
    )
    
    cell_weights <- matrix(
      c(
        rep(1, num_cells),
        rep(0, num_cells),
        rep(0, num_cells),
        rep(1, num_cells)
      ),
      ncol = 2
    )
    
    gam <- fitGAM(
      counts      = counts,
      pseudotime  = pseudotime_mat,
      cellWeights = cell_weights,
      nknots      = 6,
      parallel    = TRUE
    )
    
    res_tradeseq_diffEnd <- diffEndTest(gam)
    res_tradeseq_pattern <- patternTest(gam)
    res_tradeseq_earlyDE <- earlyDETest(gam, knots = c(1, 3))
    
    ## ====== Monocle3 用 result_tbl 整形 ====== ##
    scLS_tbl_m3 <- res_scLS %>%
      mutate(
        backbone  = model,
        pt_method = "Monocle3",
        method    = "scLS",
        num_cells = num_cells,
        answer    = answer[match(gene, genes_all)]
      ) %>%
      select(backbone, gene, pt_method, method, answer, num_cells, p)
    
    lamian_tbl_m3 <- tibble(
      backbone  = model,
      gene      = genes_all,
      pt_method = "Monocle3",
      method    = "Lamian",
      answer    = answer,
      num_cells = num_cells,
      p         = as.numeric(res_lamian$statistics$pval.chisq.overall)
    )
    
    tradeseq_diffEnd_tbl_m3 <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_diffEnd),
      pt_method = "Monocle3",
      method    = "tradeSeq (diffEnd)",
      answer    = answer[match(rownames(res_tradeseq_diffEnd), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_diffEnd$pvalue
    )
    
    tradeseq_pattern_tbl_m3 <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_pattern),
      pt_method = "Monocle3",
      method    = "tradeSeq (pattern)",
      answer    = answer[match(rownames(res_tradeseq_pattern), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_pattern$pvalue
    )
    
    tradeseq_earlyDE_tbl_m3 <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_earlyDE),
      pt_method = "Monocle3",
      method    = "tradeSeq (earlyDE)",
      answer    = answer[match(rownames(res_tradeseq_earlyDE), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_earlyDE$pvalue
    )
    
    result_tbl <- bind_rows(
      result_tbl,
      scLS_tbl_m3,
      lamian_tbl_m3,
      tradeseq_diffEnd_tbl_m3,
      tradeseq_pattern_tbl_m3,
      tradeseq_earlyDE_tbl_m3
    )
    ## ====== Monocle3 部分ここまで ====== ##
    
    # Slingshot -----------------------------------------------------------
    
    sce <- SingleCellExperiment(
      assays = List(
        counts = wt_counts,
        log2   = log2(wt_counts + 1)
      )
    )
    assays(sce)$norm <- fq_norm(wt_counts)
    
    pca_res <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
    reducedDims(sce) <- SimpleList(PCA = pca_res$x[, 1:2])
    
    set.seed(10)
    colData(sce)$gmm_cluster <- Mclust(pca_res$x[, 1:2])$classification
    
    root_cell  <- names(wt_pseudotime)[which.min(wt_pseudotime)]
    start.clus <- colData(sce)$gmm_cluster[root_cell]
    
    wt_slingshot <- slingshot(
      sce,
      clusterLabels = "gmm_cluster",
      reducedDim    = "PCA",
      start.clus    = start.clus
    )
    
    wt_curve_weights        <- slingCurveWeights(wt_slingshot)
    wt_pseudotime_slingshot <- slingPseudotime(wt_slingshot, na = FALSE)
    
    if (ncol(wt_pseudotime_slingshot) >= 2) {
      normalized_weights <- sweep(
        wt_curve_weights,
        1,
        FUN   = "/",
        STATS = apply(wt_curve_weights, 1, sum)
      )
      
      wt_curve_weights_discrete <- apply(
        normalized_weights,
        1,
        function(prob) stats::rmultinom(n = 1, prob = prob, size = 1)
      ) %>%
        t()
      
      sample_curve_pseudotime <- function(cell_index) {
        sample(
          wt_pseudotime_slingshot[cell_index, ],
          size = 1,
          prob = as.numeric(wt_curve_weights_discrete[cell_index, ])
        )
      }
      
      wt_pseudotime_slingshot <- sapply(
        seq_len(nrow(wt_pseudotime_slingshot)),
        sample_curve_pseudotime
      )
    } else {
      wt_pseudotime_slingshot <- wt_pseudotime_slingshot[, 1]
    }
    names(wt_pseudotime_slingshot) <- colnames(sce)
    
    sce <- SingleCellExperiment(
      assays = List(
        counts = ko_counts,
        log2   = log2(ko_counts + 1)
      )
    )
    assays(sce)$norm <- fq_norm(ko_counts)
    
    pca_res <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
    reducedDims(sce) <- SimpleList(PCA = pca_res$x[, 1:2])
    
    set.seed(10)
    colData(sce)$gmm_cluster <- Mclust(pca_res$x[, 1:2])$classification
    
    root_cell  <- names(ko_pseudotime)[which.min(ko_pseudotime)]
    start.clus <- colData(sce)$gmm_cluster[root_cell]
    
    ko_slingshot <- slingshot(
      sce,
      clusterLabels = "gmm_cluster",
      reducedDim    = "PCA",
      start.clus    = start.clus
    )
    
    ko_curve_weights        <- slingCurveWeights(ko_slingshot)
    ko_pseudotime_slingshot <- slingPseudotime(ko_slingshot, na = FALSE)
    
    if (ncol(ko_pseudotime_slingshot) >= 2) {
      normalized_weights <- sweep(
        ko_curve_weights,
        1,
        FUN   = "/",
        STATS = apply(ko_curve_weights, 1, sum)
      )
      
      ko_curve_weights_discrete <- apply(
        normalized_weights,
        1,
        function(prob) stats::rmultinom(n = 1, prob = prob, size = 1)
      ) %>%
        t()
      
      sample_curve_pseudotime <- function(cell_index) {
        sample(
          ko_pseudotime_slingshot[cell_index, ],
          size = 1,
          prob = as.numeric(ko_curve_weights_discrete[cell_index, ])
        )
      }
      
      ko_pseudotime_slingshot <- sapply(
        seq_len(nrow(ko_pseudotime_slingshot)),
        sample_curve_pseudotime
      )
    } else {
      ko_pseudotime_slingshot <- ko_pseudotime_slingshot[, 1]
    }
    names(ko_pseudotime_slingshot) <- colnames(sce)
    
    pseudotime_tbl <- pseudotime_tbl %>%
      add_row(
        backbone   = model,
        pt_method  = "Slingshot",
        sample     = rep(c("WT", "KO"), each = num_cells),
        cell_num   = num_cells,
        cell       = colnames(counts),
        pseudotime = c(wt_pseudotime_slingshot, ko_pseudotime_slingshot)
      )
    
    plot_df <- tibble(
      component1 = wt$dimred[sampling_indices, 1],
      component2 = wt$dimred[sampling_indices, 2],
      pseudotime = wt_pseudotime_slingshot / max(wt_pseudotime_slingshot)
    )
    
    g_wt <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      ggtitle(str_c("Slingshot, ", num_cells, " cells")) +
      labs(subtitle = "WT") +
      coord_cartesian() +
      theme(legend.position = "none")
    
    plot_df <- tibble(
      component1 = ko$dimred[sampling_indices, 1],
      component2 = ko$dimred[sampling_indices, 2],
      pseudotime = ko_pseudotime_slingshot / max(ko_pseudotime_slingshot)
    )
    
    g_ko <- ggplot(plot_df, aes(x = component1, y = component2, colour = pseudotime)) +
      geom_point() +
      labs(subtitle = "KO") +
      coord_cartesian() +
      theme(legend.position = "bottom")
    
    combined_plot <- g_wt + g_ko
    
    ggsave(
      combined_plot,
      file = str_c("./plots/slingshot_", model, "_", num_cells, "cells.pdf"),
      dpi  = 300
    )
    
    # scLS (Slingshot) -------------------------------------------
    
    wt_seu <- CreateSeuratObject(counts = wt_counts) %>%
      AddMetaData(
        metadata = wt_pseudotime_slingshot / max(wt_pseudotime_slingshot),
        col.name = "Pseudotime"
      ) %>%
      NormalizeData() %>%
      FindVariableFeatures()
    ko_seu <- CreateSeuratObject(counts = ko_counts) %>%
      AddMetaData(
        metadata = ko_pseudotime_slingshot / max(ko_pseudotime_slingshot),
        col.name = "Pseudotime"
      ) %>%
      NormalizeData() %>%
      FindVariableFeatures()
    
    res_scLS <- scLS.shift(
      wt_seu, ko_seu,
      time.col1 = "Pseudotime", time.col2 = "Pseudotime",
      n.cores   = num_cores
    )
    
    # Lamian (Slingshot) --------------------------------------------------
    
    pseudotime_vec <- c(wt_pseudotime_slingshot, ko_pseudotime_slingshot)
    names(pseudotime_vec) <- colnames(expression_mat)
    
    res_lamian <- lamian_test(
      expr         = expression_mat,
      cellanno     = cell_anno,
      pseudotime   = pseudotime_vec,
      design       = design_mat,
      test.type    = "variable",
      overall.only = TRUE,
      test.method  = "chisq",
      permuiter    = 100,
      ncores       = num_cores
    )
    
    # tradeSeq (Slingshot) ------------------------------------------------
    
    colnames(counts) <- colnames(expression_mat)
    
    pseudotime_mat <- matrix(
      c(
        wt_pseudotime_slingshot,
        rep(0, num_cells),
        rep(0, num_cells),
        ko_pseudotime_slingshot
      ),
      ncol = 2
    )
    
    if (model == "Linear") {
      cell_weights <- matrix(
        c(
          rep(1, num_cells),
          rep(0, num_cells),
          rep(0, num_cells),
          rep(1, num_cells)
        ),
        ncol = 2
      )
    } else {
      cell_weights <- rbind(wt_curve_weights, ko_curve_weights)
    }
    
    gam <- fitGAM(
      counts      = counts,
      pseudotime  = pseudotime_mat,
      cellWeights = cell_weights,
      nknots      = 6,
      parallel    = TRUE
    )
    
    res_tradeseq_diffEnd <- diffEndTest(gam)
    res_tradeseq_pattern <- patternTest(gam)
    res_tradeseq_earlyDE <- earlyDETest(gam, knots = c(1, 3))
    
    ## ====== Slingshot 用 result_tbl 整形 ====== ##
    scLS_tbl_sl <- res_scLS %>%
      mutate(
        backbone  = model,
        pt_method = "Slingshot",
        method    = "scLS",
        num_cells = num_cells,
        answer    = answer[match(gene, genes_all)]
      ) %>%
      select(backbone, gene, pt_method, method, answer, num_cells, p)
    
    lamian_tbl_sl <- tibble(
      backbone  = model,
      gene      = genes_all,
      pt_method = "Slingshot",
      method    = "Lamian",
      answer    = answer,
      num_cells = num_cells,
      p         = as.numeric(res_lamian$statistics$pval.chisq.overall)
    )
    
    tradeseq_diffEnd_tbl_sl <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_diffEnd),
      pt_method = "Slingshot",
      method    = "tradeSeq (diffEnd)",
      answer    = answer[match(rownames(res_tradeseq_diffEnd), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_diffEnd$pvalue
    )
    
    tradeseq_pattern_tbl_sl <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_pattern),
      pt_method = "Slingshot",
      method    = "tradeSeq (pattern)",
      answer    = answer[match(rownames(res_tradeseq_pattern), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_pattern$pvalue
    )
    
    tradeseq_earlyDE_tbl_sl <- tibble(
      backbone  = model,
      gene      = rownames(res_tradeseq_earlyDE),
      pt_method = "Slingshot",
      method    = "tradeSeq (earlyDE)",
      answer    = answer[match(rownames(res_tradeseq_earlyDE), genes_all)],
      num_cells = num_cells,
      p         = res_tradeseq_earlyDE$pvalue
    )
    
    result_tbl <- bind_rows(
      result_tbl,
      scLS_tbl_sl,
      lamian_tbl_sl,
      tradeseq_diffEnd_tbl_sl,
      tradeseq_pattern_tbl_sl,
      tradeseq_earlyDE_tbl_sl
    )
    ## ====== Slingshot 部分ここまで ====== ##
  }
}

# save results ------------------------------------------------------------

result_tbl$method <- factor(
  result_tbl$method,
  levels = c(
    "scLS",
    "tradeSeq (earlyDE)",
    "tradeSeq (diffEnd)",
    "tradeSeq (pattern)",
    "Lamian"
  )
)

save(result_tbl, file = "./files/res.RData")

pseudotime_tbl$cell <- pseudotime_tbl$cell %>%
  str_remove_all("WT_") %>%
  str_remove_all("KO_")

save(pseudotime_tbl, file = "./files/pt.RData")

# correlation of pseudotime -----------------------------------------------

for (model in c("Linear", "Bifurcating")) {
  for (sample_label in c("WT", "KO")) {
    cor_df <- rbind(
      inner_join(
        pseudotime_tbl %>%
          dplyr::filter(pt_method == "Ground_truth", backbone == model, sample == sample_label),
        pseudotime_tbl %>%
          dplyr::filter(pt_method == "Monocle3", backbone == model, sample == sample_label),
        by = c("sample", "cell_num", "cell")
      ),
      inner_join(
        pseudotime_tbl %>%
          dplyr::filter(pt_method == "Ground_truth", backbone == model, sample == sample_label),
        pseudotime_tbl %>%
          dplyr::filter(pt_method == "Slingshot", backbone == model, sample == sample_label),
        by = c("sample", "cell_num", "cell")
      )
    )
    
    g <- ggplot(cor_df, aes(x = pseudotime.x, y = pseudotime.y)) +
      geom_point(size = 0.1) +
      facet_grid(pt_method.y ~ cell_num, scales = "free") +
      coord_cartesian() +
      ggtitle(str_c(model, ", ", sample_label)) +
      xlab("Pseudotime of ground truth") +
      ylab("Inferred pseudotime")
    
    ggsave(
      g,
      file = str_c("./cor/cor_pseudotime_", sample_label, "_", model, ".pdf"),
      dpi  = 300
    )
  }
}

# ROC curves and barplots -------------------------------------------------

compute_auc <- function(label, score) {
  if (length(unique(label)) < 2L) return(NA_real_)
  
  ord <- order(score, decreasing = TRUE)
  y   <- label[ord]
  pos <- sum(y == 1)
  neg <- sum(y == 0)
  
  tpr <- cumsum(y == 1) / pos
  fpr <- cumsum(y == 0) / neg
  
  tpr <- c(0, tpr)
  fpr <- c(0, fpr)
  
  dx <- diff(fpr)
  ym <- (head(tpr, -1) + tail(tpr, -1)) / 2
  sum(dx * ym)
}

for (model in c("Linear", "Bifurcating")) {
  df_model <- result_tbl %>%
    dplyr::filter(
      backbone == model,
      is.finite(p),
      !is.na(answer)
    )
  
  g_roc <- ggplot(
    df_model,
    aes(d = answer, m = 1 - p, colour = method)
  ) +
    geom_roc(n.cuts = FALSE, linealpha = 1) +
    xlab("False positive rate") +
    ylab("True positive rate") +
    ggtitle(model) +
    facet_grid(pt_method ~ num_cells) +
    coord_cartesian() +
    theme_bw(base_size = 16) +
    theme(legend.position = "none") +
    coord_fixed()
  
  ggsave(
    g_roc,
    file   = str_c("./res/", model, "_ROC.pdf"),
    dpi    = 300,
    width  = 9,
    height = 7
  )
  
  g_roc <- ggplot(
    df_model,
    aes(d = answer, m = 1 - p, colour = method)
  ) +
    geom_roc(n.cuts = FALSE, linealpha = 1) +
    theme(
      legend.title    = element_blank()
    )
  ggsave(g_roc, file = "./plots/legend.pdf")
  
  roc_df <- df_model %>%
    dplyr::group_by(num_cells, method, pt_method) %>%
    dplyr::summarise(
      AUC = compute_auc(label = answer, score = 1 - p),
      .groups = "drop"
    )
  
  g_bar <- ggplot(roc_df, aes(x = method, y = AUC, fill = method)) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_grid(pt_method ~ num_cells) +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle(model) +
    ylab("AUC") +
    xlab("") +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title    = element_blank()
    )
  
  ggsave(
    g_bar,
    file   = str_c("./res/", model, "_AUC.pdf"),
    dpi    = 300,
    width  = 9,
    height = 7
  )
  
  ## ==== AUC table ==== ##
  method_order <- c(
    "scLS",
    "tradeSeq (diffEnd)",
    "tradeSeq (pattern)",
    "tradeSeq (earlyDE)",
    "Lamian"
  )
  
  auc_table <- roc_df %>%
    tidyr::pivot_wider(names_from = method, values_from = AUC) %>%
    dplyr::arrange(num_cells, pt_method)
  
  present_methods <- intersect(method_order, colnames(auc_table))
  
  auc_table <- auc_table %>%
    dplyr::select(
      num_cells,
      pt_method,
      dplyr::all_of(present_methods)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(present_methods),
        .fns  = ~ round(.x, 3)
      )
    )
  
  write_csv(auc_table, file = str_c("./res/", model, "_AUC.csv"))
}
