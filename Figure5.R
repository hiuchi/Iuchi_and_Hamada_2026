library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scMerge)
library(scLS)
library(reticulate)
library(tradeSeq)
library(phateR)
library(patchwork)
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(methods)
set.seed(8)
theme_set(theme_bw(base_size = 20))

use_python("~/.venvs/phate/bin/python", required = TRUE)
py_discover_config(required_module = "phate")

## PHATE and trajectory analysis
load("./data/mutant_WT_integrated_seurat_scaleall")
WT <- readRDS("./data/WT_sce_merged.rds")
KO <- readRDS("./data/KO_sce.rds")

features <- ref.integrated@assays[["integrated"]]@var.features
data <- ref.integrated@assays[["integrated"]]@data
data_subset <- t(data[features, ])
phate_output <- phate(data_subset)
phate_output <- as.matrix(phate_output)
rownames(phate_output) <- rownames(phate_output) %>% str_replace_all(pattern = "\\.", replacement = "-")
WT$batch <- "WT"
KO$batch <- "KO"
combined <- sce_cbind(list(WT, KO), batch_names = c("WT", "KO"))
combined$Pseudotime <- c(WT$slingPseudotime_1, KO$slingPseudotime_1)
phate_output <- phate_output[colnames(combined), , drop = FALSE]
combined$PHATE1 <- phate_output[,1]
combined$PHATE2 <- phate_output[,2]
save(combined, file = "./res/combined.RData")

data <- tibble(batch = combined$batch,
               PHATE1 = combined$PHATE1,
               PHATE2 = combined$PHATE2,
               Pseudotime = combined$Pseudotime)
g <- ggplot(data, aes(x = PHATE1, y = PHATE2, colour = batch)) +
  geom_point(size = 0.3) +
  theme(legend.title = element_blank()) +
  coord_cartesian()
g
ggsave(g, file = "./plots/PHATE.pdf", width = 8, height = 7, dpi = 300)

g <- ggplot(data, aes(x = PHATE1, y = PHATE2, colour = Pseudotime)) +
  geom_point(size = 0.3) +
  scale_color_gradient2(low  = "grey",
                        high = "royalblue") +
  theme(legend.title = element_blank()) +
  coord_cartesian()
g
ggsave(g, file = "./plots/PHATE_pseudotime.pdf", width = 8, height = 7, dpi = 300)

#tradeSeq
pseudotime <- data.frame(WT = combined$Pseudotime, KO = combined$Pseudotime, row.names = colnames(combined))
cellWeights <- data.frame(WT = rep(c(1, 0), c(ncol(WT), ncol(KO))),
                          KO = rep(c(0, 1), c(ncol(WT), ncol(KO))))
sce <- fitGAM(counts = combined@assays@data@listData[["counts"]],
              pseudotime = pseudotime,
              cellWeights = cellWeights,
              parallel = TRUE)
save(sce, file = "./res/sce.RData")

#reset Rstudio here
load("./res/sce.RData")
load("./res/combined.RData")
tradeSeq_det <- diffEndTest(sce)
tradeSeq_det$pvalue[tradeSeq_det$pvalue == 0] <- min(tradeSeq_det$pvalue[tradeSeq_det$pvalue != 0], na.rm = TRUE)
tradeSeq_pt <- patternTest(sce)
tradeSeq_pt$pvalue[tradeSeq_pt$pvalue == 0] <- min(tradeSeq_pt$pvalue[tradeSeq_pt$pvalue != 0], na.rm = TRUE)
tradeSeq_edt <- earlyDETest(sce, knots = c(1, 2))
tradeSeq_edt$pvalue[tradeSeq_edt$pvalue == 0] <- min(tradeSeq_edt$pvalue[tradeSeq_edt$pvalue != 0], na.rm = TRUE)

# scLS
reticulate::use_condaenv("scLS-env",
                         conda = "/opt/homebrew/bin/conda",
                         required = TRUE)
seu <- as.Seurat(combined, counts = "counts", data = "logcounts")
seu$batch <- combined$batch
seu$Pseudotime <- combined$Pseudotime

DefaultAssay(seu) <- "originalexp"
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures()
seu_WT <- subset(seu, subset = batch == "WT")
seu_KO <- subset(seu, subset = batch == "KO")

seu_WT$PT <- scales::rescale(seu_WT$Pseudotime, to = c(0,1), na.rm = TRUE)
seu_KO$PT <- scales::rescale(seu_KO$Pseudotime, to = c(0,1), na.rm = TRUE)

res_scLS.shift <- scLS.shift(group1.object = seu_WT,
                             group2.object = seu_KO,
                             time.col1 = "PT",
                             time.col2 = "PT",
                             features = rownames(seu_WT),
                             assay = "originalexp",
                             n.cores = 12,
                             seed = 8)
res_scLS.shift$p[res_scLS.shift$p == 0] <- min(res_scLS.shift$p[res_scLS.shift$p != 0], na.rm = TRUE)
save(res_scLS.shift, file = "./res/res_scLS.shift.RData")

#plot
res <- tibble(method = rep(c("scLS", "tradeSeq (diffEndTest)", "tradeSeq (patternTest)", "tradeSeq (earlyDETest)"),
                           c(nrow(res_scLS.shift), nrow(tradeSeq_det), nrow(tradeSeq_pt), nrow(tradeSeq_edt))),
              gene = c(res_scLS.shift$gene, rownames(tradeSeq_det), rownames(tradeSeq_pt), rownames(tradeSeq_edt)),
              p_value = c(res_scLS.shift$p, tradeSeq_det$pvalue, tradeSeq_pt$pvalue, tradeSeq_edt$pvalue))
res.wide <- res %>% pivot_wider(names_from = method, values_from = p_value) %>%
  mutate(gene.name = case_when(gene == "Tbrucei---Tb927.7.2660" ~ "ZC3H20", TRUE ~ "")) %>%
  mutate(isZC3H20 = gene == "Tbrucei---Tb927.7.2660") %>%
  arrange(isZC3H20)

log10_lab <- scales::trans_format("log10", scales::math_format(10^.x))
breaks10  <- scales::log_breaks(base = 10) 
theme_set(theme_bw(base_size = 40))
g1 <- ggplot(res.wide, aes(x = scLS, y = `tradeSeq (diffEndTest)`)) +
  geom_point(aes(colour = isZC3H20, size = isZC3H20)) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "magenta"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 0.1, `TRUE` = 2), guide = "none") +
  ggrepel::geom_text_repel(data = subset(res.wide, isZC3H20),
                           aes(label = "ZC3H20"),
                           color = "magenta",
                           size = 10,
                           max.overlaps = 100000) +
  scale_x_log10(breaks = breaks10, labels = log10_lab) +
  scale_y_log10(breaks = breaks10, labels = log10_lab) +
  xlab("scLS's p") +
  ylab("tradeSeq's p (diffEndTest)") +
  coord_cartesian()
g2 <- ggplot(res.wide, aes(x = scLS, y = `tradeSeq (patternTest)`)) +
  geom_point(aes(colour = isZC3H20, size = isZC3H20)) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "magenta"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 0.1, `TRUE` = 2), guide = "none") +
  ggrepel::geom_text_repel(data = subset(res.wide, isZC3H20),
                           aes(label = "ZC3H20"),
                           color = "magenta",
                           size = 10,
                           max.overlaps = 100000) +
  scale_x_log10(breaks = breaks10, labels = log10_lab) +
  scale_y_log10(breaks = breaks10, labels = log10_lab) +
  xlab("scLS's p") +
  ylab("tradeSeq's p (patternTest)") +
  coord_cartesian()
g3 <- ggplot(res.wide, aes(x = scLS, y = `tradeSeq (earlyDETest)`)) +
  geom_point(aes(colour = isZC3H20, size = isZC3H20)) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "magenta"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 0.1, `TRUE` = 2), guide = "none") +
  ggrepel::geom_text_repel(data = subset(res.wide, isZC3H20),
                           aes(label = "ZC3H20"),
                           color = "magenta",
                           size = 10,
                           max.overlaps = 100000) +
  scale_x_log10(breaks = breaks10, labels = log10_lab) +
  scale_y_log10(breaks = breaks10, labels = log10_lab) +
  xlab("scLS's p") +
  ylab("tradeSeq's p (earlyDETest)") +
  coord_cartesian()
g <- g1 + g2 + g3
g
ggsave(g, file = "./plots/cor.pdf", width = 30, height = 10)

g <- plotSmoothers(sce, counts(combined), gene = "Tbrucei---Tb927.7.2660")
ggsave(g, file = "./plots/ZC3H20.pdf", dpi = 300)

#GO analysis
res_scLS.shift$q <- p.adjust(res_scLS.shift$p, method = "BH")
res_scLS.shift$id <- res_scLS.shift$gene %>% str_split(pattern = "--", simplify = TRUE) %>% .[,2] %>% str_sub(start = 2, end = -1)
sig <- res_scLS.shift %>% dplyr::filter(q < 0.01)
bg <- rownames(seu) %>% str_split(pattern = "--", simplify = TRUE) %>% .[,2] %>% str_sub(start = 2, end = -1)

gost_res <- gost(query = sig$id,
                 organism = "tbrucei",
                 domain_scope = "custom",
                 custom_bg = bg,
                 significant = FALSE)

gp <- gost_res$result %>%
  dplyr::filter(!is.na(p_value)) %>%
  dplyr::filter(term_id != "KEGG:00000") %>%
  mutate(
    ID          = term_id,
    Description = term_name,
    pvalue      = p_value,
    p.adjust    = p_value,
    qvalue      = p_value,
    GeneRatio   = paste0(intersection_size, "/", query_size),
    BgRatio     = paste0(term_size, "/", effective_domain_size),
    Count       = intersection_size,
    geneID      = if ("intersection" %in% names(.)) {
      if (is.list(intersection)) {
        vapply(intersection, \(x) paste(x, collapse = "/"), character(1))
      } else {
        gsub("[,; ]+", "/", intersection)
      }
    } else {
      NA_character_
    }
  ) %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, source)

er <- methods::new("enrichResult")
er@result <- as.data.frame(gp)
er@pvalueCutoff <- 0.05
er@pAdjustMethod <- "g_SCS"
er@qvalueCutoff <- 0.05
er@gene <- sig$id

g <- enrichplot::dotplot(er, showCategory = 20, x = "GeneRatio", color = "p.adjust")  +
  theme_bw(base_size = 20) +
  guides(color = guide_colorbar(title = "Adjusted p-value"))

ggsave(g, file = "plots/GO.pdf", dpi = 300)


