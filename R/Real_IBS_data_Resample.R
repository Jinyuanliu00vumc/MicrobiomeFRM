knitr::opts_chunk$set(echo = TRUE)
# Load libraries
rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## ---- Clean prep: merge, keep complete 0-4, build per-visit OTU ----
library(dplyr)
library(stringr)
library(vegan)
library(GLMMMiRKAT)
library(MASS)

# Import functions
sourceCpp("chisq_stat.cpp")
sourceCpp("dist2mat.cpp")

# UGEE
sourceCpp("ugee3_general_V2.cpp", showOutput = FALSE)

Data1 <- read.delim("genera.counts.tsv", 
                    header = TRUE,      # or FALSE if no header
                    stringsAsFactors = FALSE)
Meta1 <- read.delim("metadata.tsv", 
                    header = TRUE,      # or FALSE if no header
                    stringsAsFactors = FALSE)
Other1 <- read.delim("mtb.tsv", 
                     header = TRUE,      # or FALSE if no header
                     stringsAsFactors = FALSE)

otu_vars <- colnames(Data1)[-1]


# 1) Merge by Sample and derive ID/Timepoint
master_dat <- Meta1 %>%
  inner_join(Data1, by = "Sample") %>%
  mutate(
    Timepoint = str_extract(Sample, "[^\\.]+$"),
    ID        = str_extract(Sample, "^[^\\.]+")
  ) %>%
  filter(Timepoint %in% c("0","1","2","3","4"),
         Timepoint != "Flare")

# 2) Keep subjects (IDs) with all five timepoints
master_dat_sub <- master_dat %>%
  group_by(ID) %>%
  filter(all(c("0","1","2","3","4") %in% Timepoint)) %>%
  ungroup()

stopifnot(all(table(master_dat_sub$ID, master_dat_sub$Timepoint) > 0))

# 3) Per-visit OTU tables (keep order by ID for all visits)
visit_levels <- c("0","1","2","3","4")
make_visit <- function(tp) {
  master_dat_sub %>%
    filter(Timepoint == tp) %>%
    arrange(ID) %>%
    select(Sample, ID, Timepoint, Study.Group, any_of(otu_vars))
}
otu_t_list <- setNames(lapply(visit_levels, make_visit), paste0("t", visit_levels))

# Helper to convert to relative abundance with row guards
to_rel_abund <- function(X) {
  X <- as.matrix(X)
  rs <- rowSums(X, na.rm = TRUE)
  zero <- which(rs == 0 | !is.finite(rs))
  if (length(zero)) {
    # add a tiny pseudocount to avoid NaNs
    X[zero, ] <- X[zero, ] + 1e-8
    rs <- rowSums(X)
  }
  sweep(X, 1, rs, "/")
}

## ---- UGEE inputs: distance matrices per time (Bray-Curtis) ----
fij.mat.all <- lapply(otu_t_list, function(tab) {
  X <- to_rel_abund(tab[, -c(1:4), drop = FALSE])
  d <- vegdist(X, method = "bray")
  dist2mat(d, nrow(X))  # dynamic size
})

## ---- Group and covariates from baseline (timepoint '0') ----
baseline <- otu_t_list[["t0"]] %>%
  select(ID, Study.Group) %>%
  left_join(master_dat_sub %>% filter(Timepoint == "0") %>% select(ID, Age, BMI), by = "ID")

# map groups: C=1, D=2, H=3
grp_map <- c(C = 1, D = 2, H = 3)
group.num <- unname(grp_map[baseline$Study.Group])

age <- scale(baseline$Age)[,1]
bmi <- scale(baseline$BMI)[,1]

## ---- Wald helper (unchanged) ----
getwald <- function(results, terms) {
  betahat <- results$theta[terms]
  Sigma   <- results$Sigma_theta[terms, terms, drop = FALSE]
  Wstat   <- as.numeric(t(betahat) %*% ginv(Sigma) %*% betahat)
  df      <- length(terms)
  pchisq(Wstat, df = df, lower.tail = FALSE)
}

## ========= A) Fit UGEE on observed data (choose ONE contrast) =========
# Choose contrast: "DvsAll" or "CvsAll"
contrast <- "DvsAll"

grp_bin <- if (contrast == "DvsAll") as.integer(group.num == 2) else as.integer(group.num == 1)

fit_ugee <- ugee3_general(
  dij_all = fij.mat.all,
  grp_all = grp_bin,
  z_ = cbind(age, bmi),
  interaction_type = "modality"
)

# Omnibus p-value (adjust term indices if your model parameterization differs)
p_ugee <- getwald(fit_ugee, terms = 2:17)
cat("UGEE omnibus p-value:", format(p_ugee, digits = 6), "\n")

## ========= B) Fit GLMM-MiRKAT on observed data (same contrast) =========
# Build stacked OTU (rows = all visits stacked in 0..4 order, within each visit ordered by ID)
Matrix_OTU <- do.call(rbind, lapply(otu_t_list, function(tab) tab[, -c(1:4), drop = FALSE]))
Matrix_OTU <- to_rel_abund(Matrix_OTU)

# Bray–Curtis distance and kernel; stabilize kernel
D.BC <- as.matrix(vegdist(Matrix_OTU, method = "bray"))
K.BC <- GLMMMiRKAT::D2K(D.BC)
K.BC <- 0.5 * (K.BC + t(K.BC))                      # symmetrize
diag(K.BC) <- diag(K.BC) + 1e-6                     # tiny ridge
nK <- nrow(K.BC)
H  <- diag(nK) - matrix(1/nK, nK, nK)
Kc <- H %*% K.BC %*% H                               # center

# Build meta aligned with stacking: each visit contributes n_subject rows in ID order
n_subj <- nrow(baseline)
stopifnot(nK == n_subj * length(visit_levels))

meta <- data.frame(
  X_t  = rep(group.num, length(visit_levels)),
  Z_1  = rep(age,       length(visit_levels)),
  Z_2  = rep(bmi,       length(visit_levels)),
  id   = rep(seq_len(n_subj), length(visit_levels)),
  visit= rep(seq_along(visit_levels), each = n_subj)
)
meta$X_t <- meta$X_t - 1L
covs <- as.data.frame(cbind(meta$Z_1, meta$Z_2))

# Binary outcome per chosen contrast
y <- if (contrast == "DvsAll") as.integer(meta$X_t + 1L == 2L) else as.integer(meta$X_t + 1L == 1L)

# Guard: both classes present
if (length(unique(y)) < 2L) {
  p_mirkat <- NA_real_
  warning("GLMM-MiRKAT: outcome has a single class under observed data; returning NA.")
} else {
  fit_glmm <- GLMMMiRKAT::GLMMMiRKAT(
    y, cov = covs, id = meta$id, Ks = list(Kc), model = "binomial"
  )
  # First kernel p-value
  p_mirkat <- suppressWarnings(as.numeric(fit_glmm$ItembyItem[1]))
}

cat("GLMM-MiRKAT p-value:", format(p_mirkat, digits = 6), "\n")


## ===================== Stratified resampling driver =====================

# 0) A stable baseline table we’ll subset from each iteration
baseline_df <- baseline %>%
  mutate(
    grp_num = unname(c(C=1, D=2, H=3)[Study.Group]),
    age_z   = as.numeric(scale(Age)),  # or keep original scaling you used
    bmi_z   = as.numeric(scale(BMI))
  ) %>%
  arrange(ID)

ids_by_group <- split(baseline_df$ID, baseline_df$Study.Group)  # C/D/H

# ---- Helpers reused inside the loop ----

to_rel_abund <- function(X) {
  X <- as.matrix(X)
  rs <- rowSums(X, na.rm = TRUE)
  z  <- which(rs == 0 | !is.finite(rs))
  if (length(z)) { X[z, ] <- X[z, ] + 1e-8; rs <- rowSums(X) }
  sweep(X, 1, rs, "/")
}

# Build per-visit distance matrices (list of 5) for a given ID set
make_fij_list_subset <- function(id_set) {
  lapply(otu_t_list, function(tab) {
    sub <- tab %>% semi_join(tibble::tibble(ID = id_set), by = "ID") %>% arrange(ID)
    X   <- to_rel_abund(sub[, -c(1:4), drop = FALSE])
    d   <- vegan::vegdist(X, method = "bray")
    dist2mat(d, nrow(X))
  })
}

# Build centered kernel on stacked visits for a given ID set
make_kernel_subset <- function(id_set, eps = 1e-6) {
  mats <- lapply(otu_t_list, function(tab) {
    tab %>% semi_join(tibble::tibble(ID = id_set), by = "ID") %>% arrange(ID) %>%
      select(-ID, -Study.Group, -Sample, -Timepoint)
  })
  Xall <- do.call(rbind, mats)
  Xall <- to_rel_abund(Xall)
  
  D <- as.matrix(vegan::vegdist(Xall, method = "bray"))
  K <- GLMMMiRKAT::D2K(D)
  K <- 0.5 * (K + t(K))
  diag(K) <- diag(K) + eps
  n <- nrow(K); H <- diag(n) - matrix(1/n, n, n)
  H %*% K %*% H
}

# Build meta & covs aligned to the stacked order for a given ID set
make_meta_subset <- function(id_set) {
  base <- baseline_df %>%
    semi_join(tibble::tibble(ID = id_set), by = "ID") %>%
    arrange(ID)
  
  n_subj <- nrow(base)
  visits <- length(otu_t_list)  # 5
  
  meta <- data.frame(
    X_t   = rep(base$grp_num, visits),
    Z_1   = rep(base$age_z,  visits),
    Z_2   = rep(base$bmi_z,  visits),
    id    = rep(seq_len(n_subj), visits),
    visit = rep(seq_len(visits), each = n_subj)
  )
  meta$X_t <- meta$X_t - 1L
  covs <- as.data.frame(cbind(meta$Z_1, meta$Z_2))
  list(meta = meta, covs = covs, base = base)
}

# One UGEE fit (returns p-value)
fit_UGEE_subset <- function(id_set, contrast = c("DvsAll","CvsAll")) {
  contrast <- match.arg(contrast)
  fij <- make_fij_list_subset(id_set)
  M   <- make_meta_subset(id_set)
  base <- M$base
  
  grp_bin <- if (contrast == "DvsAll") as.integer(base$grp_num == 2L) else as.integer(base$grp_num == 1L)
  
  out <- tryCatch({
    fit <- ugee3_general(
      dij_all = fij,
      grp_all = grp_bin,
      z_      = as.matrix(cbind(base$age_z, base$bmi_z)),
      interaction_type = "modality"
    )
    suppressWarnings(as.numeric(getwald(fit, terms = 2:17)))
  }, error = function(e) NA_real_)
  out
}

# One GLMM-MiRKAT fit (returns p-value)
fit_GLMM_subset <- function(id_set, contrast = c("DvsAll","CvsAll")) {
  contrast <- match.arg(contrast)
  Kc <- tryCatch(make_kernel_subset(id_set), error = function(e) NULL)
  if (is.null(Kc)) return(NA_real_)
  
  M <- make_meta_subset(id_set)
  meta <- M$meta; covs <- M$covs
  
  y <- if (contrast == "DvsAll") as.integer(meta$X_t + 1L == 2L) else as.integer(meta$X_t + 1L == 1L)
  
  # Guard: classes and alignment
  if (length(unique(y)) < 2L) return(NA_real_)
  if (length(y) != nrow(Kc))  return(NA_real_)
  
  out <- tryCatch({
    fit <- GLMMMiRKAT::GLMMMiRKAT(y, cov = covs, id = meta$id, Ks = list(Kc), model = "binomial")
    if (!is.null(fit$ItembyItem)) {
      suppressWarnings(as.numeric(fit$ItembyItem[1]))
    } else if (!is.null(fit$p.values)) {
      suppressWarnings(as.numeric(fit$p.values[1]))
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  out
}

# ---- One stratified resample (bootstrap by default) ----
one_resample <- function(per_group = c(C = NULL, D = NULL, H = NULL),
                         replace = TRUE,
                         contrast = c("DvsAll","CvsAll"),
                         seed = NULL) {
  contrast <- match.arg(contrast)
  if (!is.null(seed)) set.seed(seed)
  
  # Default sizes: keep original group sizes if not provided
  nC <- per_group["C"]; nD <- per_group["D"]; nH <- per_group["H"]
  if (is.na(nC) || is.null(nC)) nC <- length(ids_by_group$C)
  if (is.na(nD) || is.null(nD)) nD <- length(ids_by_group$D)
  if (is.na(nH) || is.null(nH)) nH <- length(ids_by_group$H)
  
  sample_ids <- c(
    sample(ids_by_group$C, nC, replace = replace),
    sample(ids_by_group$D, nD, replace = replace),
    sample(ids_by_group$H, nH, replace = replace)
  )
  
  p_mirkat <- fit_GLMM_subset(sample_ids, contrast)
  p_ugee   <- fit_UGEE_subset(sample_ids, contrast)
  
  c(GLMMMiRKAT = p_mirkat, UGEE = p_ugee)
}

# ---- Run B resamples and collect p-values ----
resample_loop <- function(B = 200,
                          per_group = c(C = NULL, D = NULL, H = NULL),
                          replace = TRUE,
                          contrast = c("DvsAll","CvsAll"),
                          alpha = 0.05,
                          parallel = FALSE, ncores = max(1, parallel::detectCores() - 1),
                          seed = 2025) {
  contrast <- match.arg(contrast)
  set.seed(seed)
  
  work <- function(b) one_resample(per_group, replace, contrast, seed = seed + b)
  
  if (parallel) {
    out_list <- parallel::mclapply(seq_len(B), work, mc.cores = ncores)
  } else {
    out_list <- lapply(seq_len(B), work)
  }
  P <- do.call(rbind, out_list)
  colnames(P) <- c("GLMMMiRKAT","UGEE")
  
  # simple summaries (empirical “power” at observed effect)
  pow <- colMeans(P < alpha, na.rm = TRUE)
  se  <- sqrt(pow * (1 - pow) / rowSums(!is.na(P)))
  
  list(pvals = P, prop_sig = pow, se = se, alpha = alpha, contrast = contrast,
       per_group = per_group, replace = replace, B = B)
}

## ===================== Example calls =====================

### 1) Bootstrap resamples at original group sizes (13 vs. 22 vs. 16)
res1_C <- resample_loop(B = 500, contrast = "CvsAll", replace = TRUE)
res1_C$prop_sig    # proportion significant (GLMMMiRKAT vs UGEE)

res1_D <- resample_loop(B = 500, contrast = "DvsAll", replace = TRUE)
res1_D$prop_sig    # proportion significant (GLMMMiRKAT vs UGEE)

# 2) Subsample balanced n per group (n=20 each)
res2_C <- resample_loop(B = 500, per_group = c(C=20, D=20, H=20), 
                      replace = T, contrast = "CvsAll")
res2_C$prop_sig # proportion significant (GLMMMiRKAT vs UGEE)

res2_D <- resample_loop(B = 500, per_group = c(C=20, D=20, H=20), 
                        replace = T, contrast = "DvsAll")
res2_D$prop_sig # proportion significant (GLMMMiRKAT vs UGEE)
