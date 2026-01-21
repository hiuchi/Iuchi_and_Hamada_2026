library(tidyverse)
library(Seurat)
library(reticulate)
library(scLS)
library(SingleCellExperiment)
library(monocle3)
library(tradeSeq)
library(gprofiler2)
library(patchwork)
theme_set(theme_bw(base_size = 30))
use_condaenv("scLS-env", required = TRUE)
set.seed(8)

load("./doi_10_5061_dryad_m63xsj454__v20230113/single_cell_monocle_object.Robj")

root_cells <- colnames(cds_subset)[colData(cds_subset)$cluster_id %in% c("01_GSCs and spermatogonia")]
cds_subset <- order_cells(cds_subset, root_cells = root_cells)
plot_cells(cds_subset,
           reduction_method = "UMAP",
           color_cells_by   = "pseudotime",
           show_trajectory_graph = FALSE)

counts <- monocle3::exprs(cds_subset)
counts <- as(counts, "dgCMatrix")
meta <- as.data.frame(SummarizedExperiment::colData(cds_subset))
meta$pseudotime <- monocle3::pseudotime(cds_subset)
data <- CreateSeuratObject(counts   = counts,
                           meta.data = meta,
                           assay    = "RNA")
umap <- reducedDims(cds_subset)$UMAP
umap <- umap[colnames(data), , drop = FALSE] 
colnames(umap) <- c("UMAP_1", "UMAP_2")
data[["umap"]] <- CreateDimReducObject(embeddings = umap,
                                       key = "UMAP_",
                                       assay = DefaultAssay(data))
data$pseudotime[is.infinite(data$pseudotime)] <- max(data$pseudotime[is.finite(data$pseudotime)])

g <- FeaturePlot(data, reduction = "umap", features = "pseudotime") +
  coord_fixed() +
  theme_bw(base_size = 30)
g
ggsave(g, file = "./plots/pseudotime.pdf", dpi = 300)

data <- data %>%
  NormalizeData() %>%
  FindVariableFeatures()
res_scLS <- scLS.dynamic(data, time.col = "pseudotime") %>%
  dplyr::rename(gene = Feature, scLS_pvalue = PeakFAP) %>%
  dplyr::select(gene, scLS_pvalue)

cw <- matrix(1, nrow = ncol(cds_subset), ncol = 1)
gam <- fitGAM(counts = counts(cds_subset)[VariableFeatures(data),],
              pseudotime = data$pseudotime,
              cellWeights = cw,
              nknots = 6,
              parallel = TRUE)
res_assoc <- associationTest(gam) %>%
  rownames_to_column("gene") %>%
  as_tibble %>%
  dplyr::select(gene, pvalue) %>%
  dplyr::rename(assoc_pvalue = pvalue)
res_sve <- startVsEndTest(gam) %>%
  rownames_to_column("gene") %>%
  as_tibble %>%
  dplyr::select(gene, pvalue) %>%
  dplyr::rename(sve_pvalue = pvalue)

res <- left_join(res_scLS, res_assoc, by = "gene") %>%
  left_join(res_sve, by = "gene") %>%
  mutate(scLS_q = p.adjust(scLS_pvalue,  method = "BH"),
         assoc_q = p.adjust(assoc_pvalue, method = "BH"),
         sve_q = p.adjust(sve_pvalue,   method = "BH")) %>%
  mutate(label = case_when(
    scLS_q > 1e-20 & sve_q < 10^(-12.5) ~ "skyblue",
    scLS_q < 1e-300 & sve_q > 10^(-2.5) ~ "magenta",
    TRUE ~ "black")) %>%
  dplyr::select(gene, scLS_q, assoc_q, sve_q, label)

second_min <- min(res$scLS_q[res$scLS_q > 0], na.rm = TRUE)
res$scLS_q[res$scLS_q == 0] <- second_min

second_min <- min(res$assoc_q[res$assoc_q > 0], na.rm = TRUE)
res$assoc_q[res$assoc_q == 0] <- second_min

second_min <- min(res$sve_q[res$sve_q > 0], na.rm = TRUE)
res$sve_q[res$sve_q == 0] <- second_min

g <- res %>%
  ggplot(aes(x = -log10(scLS_q),
             y = -log10(sve_q),
             color = label,
             size  = label)) +
  geom_point() +
  scale_color_identity() +
  scale_size_manual(
    values = c(
      skyblue = 1,
      magenta = 1,
      black   = 0.1
    )
  ) +
  labs(x = "-log10(scLS q-value)",
       y = "-log10(tradeSeq startVsEndTest q-value)") +
  coord_cartesian() +
  theme_bw(base_size = 27) +
  theme(legend.position = "none")

ggsave(g, file = "./plots/scLS_vs_sve.pdf", dpi = 300)

g <- res %>%
  ggplot(aes(x = -log10(scLS_q),
             y = -log10(assoc_q))) +
  geom_point(size = 0.1) +
  labs(x = "-log10(scLS q-value)",
       y = "-log10(tradeSeq association q-value)") +
  coord_cartesian() +
  theme_bw(base_size = 30)
ggsave(g, file = "./plots/scLS_vs_assoc.pdf", dpi = 300)

#GO analysis
topN <- 100
genes_scLS <- res %>%
  dplyr::filter(!is.na(scLS_q), scLS_q < 0.05) %>%
  arrange(scLS_q) %>%
  pull(gene) %>%
  head(topN)

p <- gostplot(gost_res)
p

publish_gostplot(
  gost_res,
  highlight_terms = gost_res$result$term_id[1:10],
  filename = "./GO/scLS_gostplot_pub.pdf",
)

genes_assoc <- res %>%
  dplyr::filter(!is.na(assoc_q), assoc_q < 0.05) %>%
  arrange(assoc_q) %>%
  pull(gene) %>%
  head(topN)

genes_sve <- res %>%
  dplyr::filter(!is.na(sve_q), sve_q < 0.05) %>%
  arrange(sve_q) %>%
  pull(gene) %>%
  head(topN)

gost_res <- gost(query = genes_scLS,
                 organism = "dmelanogaster",
                 user_threshold = 0.05,
                 domain_scope = "custom",
                 custom_bg = VariableFeatures(data))
write_csv(gost_res$result, file = "./GO/scLS.csv")

## 必要パッケージ
library(dplyr)
library(enrichplot)
library(DOSE)     # enrichResult クラス
library(methods)  # new()

gp <- gost_res$result %>%
  filter(!is.na(p_value)) %>%
  filter(term_id != "KEGG:00000") %>%
  mutate(
    ID          = term_id,
    Description = term_name,
    pvalue      = p_value,
    p.adjust    = p_value,  # g:SCS 補正済み p
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
  select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, source)

er <- methods::new("enrichResult")
er@result <- as.data.frame(gp)
er@pvalueCutoff <- 0.05
er@pAdjustMethod <- "g_SCS"
er@qvalueCutoff <- 0.05
er@gene <- genes_scLS

g <- enrichplot::dotplot(er, showCategory = 20, x = "GeneRatio", color = "p.adjust")  +
  theme_bw(base_size = 20) +
  guides(color = guide_colorbar(title = "Adjusted p-value"))

ggsave(g, file = "plots/GO.pdf", dpi = 300)

res_long <- res %>%
  pivot_longer(-c(gene, label, scLS_q), names_to = "method", values_to = "q")

genes.magenta <- res_long %>%
  dplyr::filter(label == "magenta") %>%
  pull(gene) %>%
  unique()
for(gene in genes.magenta){
  g <- plotSmoothers(models = gam,
                     counts = counts(cds_subset),
                     gene = gene) +
    labs(title = gene) +
    theme(legend.position = "none",
          plot.title      = element_text(size = 30),
          axis.title      = element_text(size = 22),
          axis.text       = element_text(size = 22))
  ggsave(g, file = str_c("./plots/genes/magenta/", gene, ".pdf"), dpi = 300)
}

genes.skyblue <- res_long %>%
  dplyr::filter(label == "skyblue") %>%
  pull(gene) %>%
  unique()
for(gene in genes.skyblue){
  g <- plotSmoothers(models = gam,
                     counts = counts(cds_subset),
                     gene = gene)+
    labs(title = gene) +
    theme(legend.position = "none",
          plot.title      = element_text(size = 30),
          axis.title      = element_text(size = 22),
          axis.text       = element_text(size = 22))
  ggsave(g, file = str_c("./plots/genes/skyblue/", gene, ".pdf"), dpi = 300)
}

plot_list <- lapply(genes.skyblue, function(gene) {
  plotSmoothers(models = gam,
                counts = counts(cds_subset),
                gene   = gene) +
    labs(title = gene) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 30, color = "skyblue"),  # ←ここ
      axis.title      = element_text(size = 22),
      axis.text       = element_text(size = 22)
    )
})

g1 <- wrap_plots(plot_list, ncol = 4)

ggsave(plot = g1,
       filename = "./plots/genes/skyblue.pdf",
       dpi = 300,
       width = 16,
       height = 4 * ceiling(length(plot_list) / 3))

plot_list <- lapply(genes.magenta, function(gene) {
  plotSmoothers(models = gam,
                counts = counts(cds_subset),
                gene   = gene) +
    labs(title = gene) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 30, color = "magenta"),  # ←ここ
      axis.title      = element_text(size = 22),
      axis.text       = element_text(size = 22)
    )
})

g2 <- wrap_plots(plot_list, ncol = 3)

ggsave(plot = g2,
       filename = "./plots/genes/magenta.pdf",
       dpi = 300,
       width = 12,
       height = 4 * ceiling(length(plot_list) / 3))
