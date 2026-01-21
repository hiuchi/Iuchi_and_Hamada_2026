library(tidyverse)
library(Seurat)
library(reticulate)
library(scLS)
library(tradeSeq)
library(Lamian)
library(SingleCellExperiment)
library(slingshot)
library(patchwork)
library(dyntoy)
library(dyngen)
library(edgeR)

theme_set(theme_bw(base_size = 20))
use_condaenv("scLS-env", required = TRUE)

if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

set.seed(1)

num_cells   <- 1000L
num_genes   <- 2000L
min_nnz     <- 10L
min_var     <- 1e-5

## ----------------------------------------------------------------------
## Helpers
## ----------------------------------------------------------------------

run_benchmark <- function(fun, n = 3L) {
  replicate(n, unname(system.time(fun())[["elapsed"]]))
}

filter_genes_for_lamian <- function(expr_log) {
  nz <- rowSums(expr_log > 0)
  v  <- apply(expr_log, 1L, var)
  keep <- nz >= min_nnz & is.finite(v) & v > min_var
  keep[is.na(keep)] <- FALSE
  keep
}

# Compute pseudotime directly from milestone_network + progressions
# (no trajectory inference).
compute_pseudotime_from_progressions <- function(milestone_network,
                                                 progressions,
                                                 scale = TRUE) {
  milnet <- as_tibble(milestone_network)
  prog   <- as_tibble(progressions)
  
  milestones <- union(milnet$from, milnet$to)
  
  root_candidates <- setdiff(milnet$from, milnet$to)
  if (length(root_candidates) == 0L) {
    root_candidates <- milestones[1]
  }
  root <- root_candidates[1]
  
  milestone_pos <- rep(NA_real_, length(milestones))
  names(milestone_pos) <- milestones
  milestone_pos[root] <- 0
  
  # simple DP along linear / acyclic backbone
  while (any(is.na(milestone_pos))) {
    progressed <- FALSE
    for (i in seq_len(nrow(milnet))) {
      from <- milnet$from[i]
      to   <- milnet$to[i]
      len  <- milnet$length[i]
      if (!is.na(milestone_pos[[from]]) && is.na(milestone_pos[[to]])) {
        milestone_pos[[to]] <- milestone_pos[[from]] + len
        progressed <- TRUE
      }
    }
    if (!progressed) break
  }
  
  if (any(is.na(milestone_pos))) {
    warning("Some milestones have undefined positions; check topology.")
  }
  
  prog2 <- prog %>%
    left_join(milnet, by = c("from", "to")) %>%
    mutate(
      pos_from   = milestone_pos[from],
      pos_to     = milestone_pos[to],
      seg_length = ifelse(is.na(length), pos_to - pos_from, length),
      seg_length = ifelse(is.na(seg_length), 0, seg_length),
      pseudo_loc = pos_from + percentage * seg_length
    )
  
  na_loc <- is.na(prog2$pseudo_loc)
  prog2$pseudo_loc[na_loc] <- prog2$pos_from[na_loc]
  
  pt_tbl <- prog2 %>%
    group_by(cell_id) %>%
    summarise(
      total_weight = sum(percentage),
      pt = ifelse(
        total_weight > 0,
        sum(pseudo_loc * percentage) / total_weight,
        mean(pseudo_loc)
      ),
      .groups = "drop"
    )
  
  pt <- pt_tbl$pt
  names(pt) <- pt_tbl$cell_id
  
  if (scale && length(pt) > 1L) {
    rmin <- min(pt, na.rm = TRUE)
    rmax <- max(pt, na.rm = TRUE)
    if (is.finite(rmin) && is.finite(rmax) && rmax > rmin) {
      pt <- (pt - rmin) / (rmax - rmin)
    } else {
      pt <- pt - rmin
    }
  }
  
  pt
}

## ============================================================
## 1. Dynamic test (single trajectory, dyntoy)
## ============================================================

message("Simulating dynamic dataset with dyntoy ...")

dataset_dyn <- dyntoy::generate_dataset(
  model                         = "linear",
  num_cells                     = num_cells,
  num_features                  = num_genes,
  differentially_expressed_rate = 0.5,
  allow_tented_progressions     = FALSE
)

counts_dyn <- t(as.matrix(dataset_dyn$counts))
if (!is.null(dataset_dyn$feature_ids)) {
  rownames(counts_dyn) <- dataset_dyn$feature_ids
}
colnames(counts_dyn) <- dataset_dyn$cell_ids
stopifnot(nrow(counts_dyn) == num_genes)

pseudotime_dyn <- compute_pseudotime_from_progressions(
  milestone_network = dataset_dyn$milestone_network,
  progressions      = dataset_dyn$progressions,
  scale             = TRUE
)
pseudotime_dyn <- pseudotime_dyn[colnames(counts_dyn)]
stopifnot(length(pseudotime_dyn) == ncol(counts_dyn))
names(pseudotime_dyn) <- colnames(counts_dyn)

expr_dyn_log <- edgeR::cpm(counts_dyn, log = TRUE, prior.count = 1)

keep_dyn     <- filter_genes_for_lamian(expr_dyn_log)
counts_dyn   <- counts_dyn[keep_dyn, , drop = FALSE]
expr_dyn_log <- expr_dyn_log[keep_dyn, , drop = FALSE]

message(
  "Dynamic dataset dimensions after filtering (genes x cells): ",
  paste(dim(counts_dyn), collapse = " x ")
)

seurat_dyn <- CreateSeuratObject(
  counts  = counts_dyn,
  assay   = "RNA",
  project = "dynamic"
)
seurat_dyn$pseudotime <- pseudotime_dyn[colnames(seurat_dyn)]

features_dyn <- rownames(seurat_dyn)
VariableFeatures(seurat_dyn) <- features_dyn

cells_dyn <- colnames(counts_dyn)

cellanno_dyn <- data.frame(
  cell   = cells_dyn,
  sample = rep(c("sample1", "sample2"), length.out = length(cells_dyn)),
  stringsAsFactors = FALSE
)
rownames(cellanno_dyn) <- cells_dyn

pt_lamian_dyn <- jitter(pseudotime_dyn[cells_dyn], amount = 1e-5)
names(pt_lamian_dyn) <- cells_dyn

design_dyn <- matrix(
  c(1, 0,
    1, 1),
  nrow  = 2,
  byrow = TRUE,
  dimnames = list(
    c("sample1", "sample2"),
    c("intercept", "group")
  )
)

stopifnot(
  all(sort(unique(cellanno_dyn$sample)) == rownames(design_dyn))
)

times_scLS_dynamic <- run_benchmark(function() {
  scLS.dynamic(
    object      = seurat_dyn,
    time.col    = "pseudotime",
    center      = TRUE,
    window.func = "hanning",
    f.min       = 0,
    f.max       = 2.0,
    n.bins      = 500,
    fap.method  = "baluev"
  )
})

tradeSeq_dynamic <- function(test_fun) {
  run_benchmark(function() {
    sce_ts <- fitGAM(
      counts      = counts_dyn,
      pseudotime  = pseudotime_dyn,
      cellWeights = rep(1, length(pseudotime_dyn)),
      nknots      = 6,
      verbose     = TRUE,
      parallel    = FALSE
    )
    test_fun(sce_ts)
  })
}

times_trade_assoc    <- tradeSeq_dynamic(associationTest)
times_trade_startEnd <- tradeSeq_dynamic(startVsEndTest)

times_lamian_dynamic <- run_benchmark(function() {
  lamian_test(
    expr       = expr_dyn_log,
    cellanno   = cellanno_dyn,
    pseudotime = pt_lamian_dyn,
    design     = design_dyn,
    test.type  = "time",
    permuiter  = 100,
    ncores     = 1
  )
})

time_list1 <- list(
  scLS              = times_scLS_dynamic,
  tradeSeq_assoc    = times_trade_assoc,
  tradeSeq_startEnd = times_trade_startEnd,
  Lamian            = times_lamian_dynamic
)

res_dynamic <- enframe(time_list1, name = "Method", value = "times") %>%
  mutate(
    Mean_sec = purrr::map_dbl(times, mean),
    SD_sec   = purrr::map_dbl(times, sd),
    Mean_min = Mean_sec / 60,
    SD_min   = SD_sec / 60
  ) %>%
  select(-times) %>%
  mutate(
    Method = factor(
      Method,
      levels = c("scLS", "tradeSeq_assoc", "tradeSeq_startEnd", "Lamian")
    )
  )

save(res_dynamic, file = "res_dynamic.RData")

g_dynamic <- ggplot(res_dynamic, aes(x = Method, y = Mean_min)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_min - SD_min, ymax = Mean_min + SD_min),
                width = 0.2) +
  scale_y_continuous(name = "Run time (min)") +
  scale_x_discrete(
    name   = "Method",
    labels = c(
      scLS              = "scLS",
      tradeSeq_assoc    = "tradeSeq\n(associationTest)",
      tradeSeq_startEnd = "tradeSeq\n(startVsEndTest)",
      Lamian            = "Lamian"
    )
  ) +
  theme_bw(base_size = 30) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_blank()
  ) +
  ggtitle("A. Dynamic test")

print(g_dynamic)
ggsave("./plots/dynamic_sim.pdf", g_dynamic, dpi = 300)

## ============================================================
## 2. Shifted test (WT vs KO, same backbone, dyngen)
## ============================================================

message("Simulating WT vs KO with dyngen (simulation_type_knockdown) ...")

num_cells_shift_total <- 2L * num_cells
num_genes_shift       <- num_genes

backbone <- dyngen::backbone_linear()
num_tfs  <- nrow(backbone$module_info)

num_targets <- floor((num_genes_shift - num_tfs) / 2)
num_hks     <- num_genes_shift - num_tfs - num_targets
stopifnot(num_targets > 0, num_hks > 0)

num_cores <- 1L

config_shift <- dyngen::initialise_model(
  backbone    = backbone,
  num_cells   = num_cells_shift_total,
  num_tfs     = num_tfs,
  num_targets = num_targets,
  num_hks     = num_hks,
  num_cores   = num_cores,
  simulation_params = dyngen::simulation_default(
    census_interval = 10,
    ssa_algorithm   = dyngen::ssa_etl(tau = 300 / 3600),
    experiment_params = dyngen::simulation_type_wild_type(
      num_simulations = 100
    )
  )
)

model_common <- config_shift %>%
  dyngen::generate_tf_network() %>%
  dyngen::generate_feature_network() %>%
  dyngen::generate_kinetics() %>%
  dyngen::generate_gold_standard()

model_wt <- model_common %>%
  dyngen::generate_cells()

ko_gene <- "C1_TF1"
ko_rate <- 0

model_ko <- model_common
model_ko$simulation_params$experiment_params <- dyngen::simulation_type_knockdown(
  num_simulations = 100,
  timepoint       = 0,
  genes           = ko_gene,
  num_genes       = 1,
  multiplier      = ko_rate
)
model_ko <- model_ko %>%
  dyngen::generate_cells()

wt_model <- model_wt %>%
  dyngen::generate_experiment()

ko_model <- model_ko %>%
  dyngen::generate_experiment()

wt_ds <- dyngen::as_dyno(wt_model)
ko_ds <- dyngen::as_dyno(ko_model)

wt_counts <- t(as.matrix(wt_ds$counts))
if (!is.null(wt_ds$feature_ids)) {
  rownames(wt_counts) <- wt_ds$feature_ids
}
cells_wt_all <- paste0("WT_", wt_ds$cell_ids)
colnames(wt_counts) <- cells_wt_all

ko_counts <- t(as.matrix(ko_ds$counts))
if (!is.null(ko_ds$feature_ids)) {
  rownames(ko_counts) <- ko_ds$feature_ids
}
cells_ko_all <- paste0("KO_", ko_ds$cell_ids)
colnames(ko_counts) <- cells_ko_all

common_genes_shift <- intersect(rownames(wt_counts), rownames(ko_counts))
wt_counts <- wt_counts[common_genes_shift, , drop = FALSE]
ko_counts <- ko_counts[common_genes_shift, , drop = FALSE]

pt_wt_raw <- compute_pseudotime_from_progressions(
  milestone_network = wt_ds$milestone_network,
  progressions      = wt_ds$progressions,
  scale             = FALSE
)
pt_ko_raw <- compute_pseudotime_from_progressions(
  milestone_network = ko_ds$milestone_network,
  progressions      = ko_ds$progressions,
  scale             = FALSE
)

names(pt_wt_raw) <- paste0("WT_", names(pt_wt_raw))
names(pt_ko_raw) <- paste0("KO_", names(pt_ko_raw))

pt_wt_raw <- pt_wt_raw[cells_wt_all]
pt_ko_raw <- pt_ko_raw[cells_ko_all]

target_cells_total <- num_cells

n_wt_all <- length(cells_wt_all)
n_ko_all <- length(cells_ko_all)
n_all    <- n_wt_all + n_ko_all

if (n_all <= target_cells_total) {
  warning("Total cells from dyngen (", n_all,
          ") is <= target; using all cells.")
  cells_wt_sel <- cells_wt_all
  cells_ko_sel <- cells_ko_all
} else {
  frac_wt   <- n_wt_all / n_all
  target_wt <- max(10L, round(target_cells_total * frac_wt))
  target_ko <- max(10L, target_cells_total - target_wt)
  target_wt <- min(target_wt, n_wt_all)
  target_ko <- min(target_ko, n_ko_all)
  set.seed(2)
  cells_wt_sel <- sample(cells_wt_all, target_wt)
  cells_ko_sel <- sample(cells_ko_all, target_ko)
}

wt_counts_sel <- wt_counts[, cells_wt_sel, drop = FALSE]
ko_counts_sel <- ko_counts[, cells_ko_sel, drop = FALSE]

pt_wt_sel_raw <- pt_wt_raw[cells_wt_sel]
pt_ko_sel_raw <- pt_ko_raw[cells_ko_sel]

counts_shift <- cbind(wt_counts_sel, ko_counts_sel)
stopifnot(identical(colnames(counts_shift), c(cells_wt_sel, cells_ko_sel)))

pseudotime_shift_raw <- c(pt_wt_sel_raw, pt_ko_sel_raw)
pt_min <- min(pseudotime_shift_raw, na.rm = TRUE)
pt_max <- max(pseudotime_shift_raw, na.rm = TRUE)
if (is.finite(pt_min) && is.finite(pt_max) && pt_max > pt_min) {
  pseudotime_shift <- (pseudotime_shift_raw - pt_min) / (pt_max - pt_min)
} else {
  pseudotime_shift <- pseudotime_shift_raw - pt_min
}
pseudotime_shift <- pseudotime_shift[colnames(counts_shift)]
names(pseudotime_shift) <- colnames(counts_shift)

expr_shift_log <- edgeR::cpm(counts_shift, log = TRUE, prior.count = 1)

keep_shift     <- filter_genes_for_lamian(expr_shift_log)
counts_shift   <- counts_shift[keep_shift, , drop = FALSE]
expr_shift_log <- expr_shift_log[keep_shift, , drop = FALSE]

message(
  "Shifted dataset dimensions after filtering (genes x cells): ",
  paste(dim(counts_shift), collapse = " x ")
)

cells_shift <- colnames(counts_shift)
cells_wt_sel <- intersect(cells_wt_sel, cells_shift)
cells_ko_sel <- intersect(cells_ko_sel, cells_shift)

seurat_wt <- CreateSeuratObject(
  counts  = counts_shift[, cells_wt_sel, drop = FALSE],
  assay   = "RNA",
  project = "WT"
)
seurat_wt$pseudotime <- pseudotime_shift[colnames(seurat_wt)]

seurat_ko <- CreateSeuratObject(
  counts  = counts_shift[, cells_ko_sel, drop = FALSE],
  assay   = "RNA",
  project = "KO"
)
seurat_ko$pseudotime <- pseudotime_shift[colnames(seurat_ko)]

# harmonise feature names for scLS.shift
common_features <- intersect(rownames(seurat_wt), rownames(seurat_ko))
common_features <- sort(common_features)

VariableFeatures(seurat_wt) <- common_features
VariableFeatures(seurat_ko) <- common_features
features_shift <- common_features

cells_wt_sel <- colnames(seurat_wt)
cells_ko_sel <- colnames(seurat_ko)

set.seed(3)
wt_half <- ceiling(length(cells_wt_sel) / 2)
ko_half <- ceiling(length(cells_ko_sel) / 2)

wt_group <- c(
  rep("WT1", wt_half),
  rep("WT2", length(cells_wt_sel) - wt_half)
)
ko_group <- c(
  rep("KO1", ko_half),
  rep("KO2", length(cells_ko_sel) - ko_half)
)

cellanno_shift <- data.frame(
  cell   = c(cells_wt_sel, cells_ko_sel),
  sample = c(wt_group,     ko_group),
  stringsAsFactors = FALSE
)
rownames(cellanno_shift) <- cellanno_shift$cell

design_shift <- matrix(
  c(1, 0,
    1, 0,
    1, 1,
    1, 1),
  nrow  = 4,
  byrow = TRUE,
  dimnames = list(
    c("WT1", "WT2", "KO1", "KO2"),
    c("intercept", "group")
  )
)

stopifnot(
  all(sort(unique(cellanno_shift$sample)) == sort(rownames(design_shift)))
)

pt_lamian_shift <- jitter(
  pseudotime_shift[rownames(cellanno_shift)],
  amount = 1e-5
)
names(pt_lamian_shift) <- rownames(cellanno_shift)

seurat_wt <- NormalizeData(seurat_wt, verbose = FALSE)
seurat_ko <- NormalizeData(seurat_ko, verbose = FALSE)

times_scLS_shift <- run_benchmark(function() {
  scLS.shift(
    group1.object = seurat_wt,
    group2.object = seurat_ko,
    time.col1     = "pseudotime",
    time.col2     = "pseudotime",
    center        = TRUE,
    window.func   = "hanning",
    f.min         = 0,
    f.max         = 2.0,
    n.bins        = 500,
    n.cores       = 1,
    features      = features_shift
  )
})

cells_all_shift <- colnames(counts_shift)

pt_mat_shift <- cbind(
  WT = pseudotime_shift[cells_all_shift],
  KO = pseudotime_shift[cells_all_shift]
)

cw_mat_shift <- matrix(
  0,
  nrow = length(cells_all_shift),
  ncol = 2,
  dimnames = list(
    cells_all_shift,
    c("WT", "KO")
  )
)
cw_mat_shift[cells_wt_sel, "WT"] <- 1
cw_mat_shift[cells_ko_sel, "KO"] <- 1

stopifnot(!any(is.na(pt_mat_shift)))
stopifnot(!any(is.na(cw_mat_shift)))

times_trade_diffEnd <- run_benchmark(function() {
  sce_ts <- fitGAM(
    counts      = counts_shift,
    pseudotime  = pt_mat_shift,
    cellWeights = cw_mat_shift,
    nknots      = 6,
    verbose     = FALSE,
    parallel    = FALSE
  )
  diffEndTest(sce_ts)
})

times_trade_pattern <- run_benchmark(function() {
  sce_ts <- fitGAM(
    counts      = counts_shift,
    pseudotime  = pt_mat_shift,
    cellWeights = cw_mat_shift,
    nknots      = 6,
    verbose     = FALSE,
    parallel    = FALSE
  )
  patternTest(sce_ts)
})

times_trade_earlyDE <- run_benchmark(function() {
  sce_ts <- fitGAM(
    counts      = counts_shift,
    pseudotime  = pt_mat_shift,
    cellWeights = cw_mat_shift,
    nknots      = 6,
    verbose     = FALSE,
    parallel    = FALSE
  )
  earlyDETest(sce_ts)
})

times_lamian_shift <- run_benchmark(function() {
  lamian_test(
    expr       = expr_shift_log,
    cellanno   = cellanno_shift,
    pseudotime = pt_lamian_shift,
    design     = design_shift,
    test.type  = "variable",
    testvar    = 2,
    permuiter  = 100,
    ncores     = 1
  )
})

time_list2 <- list(
  scLS_shift       = times_scLS_shift,
  tradeSeq_diffEnd = times_trade_diffEnd,
  tradeSeq_pattern = times_trade_pattern,
  tradeSeq_earlyDE = times_trade_earlyDE,
  Lamian_shift     = times_lamian_shift
)

res_shift <- enframe(time_list2, name = "Method", value = "times") %>%
  mutate(
    Mean_sec = purrr::map_dbl(times, mean),
    SD_sec   = purrr::map_dbl(times, sd),
    Mean_min = Mean_sec / 60,
    SD_min   = SD_sec / 60
  ) %>%
  select(-times) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "scLS_shift",
        "tradeSeq_diffEnd",
        "tradeSeq_pattern",
        "tradeSeq_earlyDE",
        "Lamian_shift"
      )
    )
  )

save(res_shift, file = "res_shift.RData")

g_shift <- ggplot(res_shift, aes(x = Method, y = Mean_min)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_min - SD_min, ymax = Mean_min + SD_min),
                width = 0.2) +
  scale_y_continuous(name = "Run time (min)") +
  scale_x_discrete(
    name   = "Method",
    labels = c(
      scLS_shift       = "scLS",
      tradeSeq_diffEnd = "tradeSeq\n(diffEndTest)",
      tradeSeq_pattern = "tradeSeq\n(patternTest)",
      tradeSeq_earlyDE = "tradeSeq\n(earlyDETest)",
      Lamian_shift     = "Lamian"
    )
  ) +
  theme_bw(base_size = 30) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_blank()
  ) +
  ggtitle("B. Shifted test")

print(g_shift)
ggsave("./plots/shift_sim.pdf", g_shift, dpi = 300, width = 12)

## ============================================================
## 3. Combined figure
## ============================================================

g_combined <- g_dynamic + g_shift
print(g_combined)
ggsave("./plots/fig_sim.pdf", g_combined, dpi = 300, width = 12)
