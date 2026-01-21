This is all the source code and a small example related to Iuchi and Hamada, 2026.  
[Figure2](https://github.com/hiuchi/LS/blob/main/Figure_2)  
[Figure3](https://github.com/hiuchi/LS/blob/main/Figure_3)  
[Figure4](https://github.com/hiuchi/LS/blob/main/Figure_4)  
[Figure5](https://github.com/hiuchi/LS/blob/main/Figure_5)  


```
library(tidyverse)
library(dyno)
library(dyngen)
library(dyntoy)
library(edgeR)
library(spectral)

#Function difinition
LS <- function(g, time1, exp1, time2 = NULL, exp2 = NULL, min.freq = 0, max.freq = 1, length.freq = 101, method.ls = "generalized", method.dist = "canberra", sim.num = 100){
  .quiet <- function(x){ 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  time1 <- time1[!is.infinite(time1)]
  exp1 <- exp1[,!is.infinite(time1)]
  if(is.null(time2)){
    seq1 <- exp1[which(rownames(exp1) == g),]
    ls1 <- spec.lomb(x = time1, y = seq1,
                     f = seq(min.freq, max.freq, length = length.freq), mode = method.ls) %>% .quiet
    return(tibble(gene = g, p.dynamic = min(ls1$p), p.shift = NA))
  }else{
    time2 <- time2[!is.infinite(time2)]
    exp2 <- exp2[,!is.infinite(time2)]
    seq1 <- exp1[which(rownames(exp1) == g),]
    seq2 <- exp2[which(rownames(exp2) == g),]
    ls1 <- spec.lomb(x = time1, y = seq1,
                     f = seq(min.freq, max.freq, length = length.freq), mode = method.ls) %>% .quiet
    ls2 <- spec.lomb(x = time2, y = seq2,
                     f = seq(min.freq, max.freq, length = length.freq), mode = method.ls) %>% .quiet
    
    null.hypothesis <- NULL
    for(i in seq(sim.num)){
      sim.ls1 <- spec.lomb(x = sample(time1, length(time1)),
                           y = seq1,
                           f = seq(min.freq, max.freq, length = length.freq), mode = method.ls) %>% .quiet
      sim.ls2 <- spec.lomb(x = sample(time2, length(time2)),
                           y = seq2,
                           f = seq(min.freq, max.freq, length = length.freq), mode = method.ls) %>% .quiet
      null.hypothesis <- c(null.hypothesis, dist(rbind(sim.ls1$A, sim.ls2$A), method = method.dist))
    }
    p <- sum(dist(rbind(ls1$A, ls2$A), method = method.dist) <  null.hypothesis) / sim.num
    return(tibble(gene = g, p.dynamic = min(c(ls1$p, ls2$p)), p.shift = p))
  }
}

#Dynamic test
#Generation of simulation data
set.seed(1)
model <- "linear"
num.cells <- 1000
num.features <- 2000
de.rate <- 0.5
dataset <- generate_dataset(model = model,
                            num_cells = num.cells,
                            num_features = num.features,
                            differentially_expressed_rate = de.rate) %>%
  add_root(root_milestone_id = .$prior_information$start_milestones) %>%
  add_pseudotime()
expression <- cpm(t(as.matrix(dataset$counts)))
pseudotime <- dataset$pseudotime / max(dataset$pseudotime)
pt.df <- tibble(pt.method = "Ground_truth",
                cell.num = num.cells,
                cell = colnames(expression),
                pseudotime = dataset$pseudotime)

#Perform a dynamic test
res <- map_dfr(rownames(expression), LS, pseudotime, expression)

#Shifted test
#Generation of simulation data
num.targets <- 250
num.hks <- 250
ko_gene <- "C1_TF1"
ko.rate <- 0
num.cores <- 8
backbone <- backbone_linear()
config <- initialise_model(backbone = backbone,
                           num_cells = num.cells,
                           num_tfs = nrow(backbone$module_info),
                           num_targets = num.targets,
                           num_hks = num.hks,
                           num_cores = num.cores,
                           simulation_params = simulation_default(census_interval = 10,
                                                                  ssa_algorithm = ssa_etl(tau = 300 / 3600),
                                                                  experiment_params = simulation_type_wild_type(num_simulations = 100)))
model_common <- config %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()
model_wt <- model_common %>% generate_cells()
model_ko <- model_common
model_ko$simulation_params$experiment_params <- simulation_type_knockdown(num_simulations = 100,
                                                                          timepoint = 0, 
                                                                          genes = ko_gene,
                                                                          num_genes = 1,
                                                                          multiplier = ko.rate)
model_ko <- model_ko %>% generate_cells()
wt <- model_wt %>% generate_experiment() %>% as_dyno() %>% add_root(root_milestone_id = "sA") %>% add_pseudotime(pseudotime = NULL)
ko <- model_ko %>% generate_experiment() %>% as_dyno() %>% add_root(root_milestone_id = "sA") %>% add_pseudotime(pseudotime = NULL)
wt.expression <- wt$expression %>% as.matrix() %>% t()
wt.pseudotime <- wt$pseudotime
ko.expression <- ko$expression %>% as.matrix() %>% t()
ko.pseudotime <- ko$pseudotime

#Perform a shifted test
res <- map_dfr(rownames(wt.expression), LS, wt.pseudotime, wt.expression, ko.pseudotime, ko.expression)

res
# A tibble: 519 × 3
   gene      p.dynamic  p.de
   <chr>         <dbl> <dbl>
 1 Burn1_TF1  1.13e-29  0.24
 2 Burn2_TF1  1.61e-35  0.4 
 3 Burn3_TF1  3.80e-23  0.14
 4 Burn4_TF1  2.56e-35  0.4 
 5 Burn5_TF1  1.29e-23  0.19
 6 A1_TF1     1.89e-28  0.83
 7 A2_TF1     1.40e-52  0.68
 8 A3_TF1     9.26e- 8  0.68
 9 A4_TF1     1.93e-56  0.72
10 B1_TF1     8.38e-64  0.34
# ℹ 509 more rows
# ℹ Use `print(n = ...)` to see more rows
```
