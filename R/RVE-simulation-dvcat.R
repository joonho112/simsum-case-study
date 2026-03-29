source("simulations/simulation-functions.R")
load(file = "simulations/dv-model-data.Rdata")
source_obj <- ls()

#-------------------------------------------------------------------------------
# Simulation design
#-------------------------------------------------------------------------------

set.seed(20210128)

design_factors <- list(
  n_studies = c(30, 60),
  struct = c("CS","HCS"),
  tau = c(0.05, 0.10, 0.20),
  omega = c(0.00, 0.10),
  rho = c(0.2, 0.4, 0.6, 0.8),
  nu = 9
)

batches <- 34
total_reps <- 4012

lengths(design_factors)

params <- expand.grid(c(design_factors, list(batch = 1:batches)))
params <- subset(params, struct == "CS" | omega == 0.00)

(n_cond <- nrow(params) / batches)
params$iterations <- total_reps / batches
params$seed <- round(runif(1) * 2^30) + 1:nrow(params)
nrow(params)
head(params)

params$batch <- NULL
head(params)


#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, setup = "register")

tm <- system.time(
  results <- plyr::mdply(params, .f = purrr::possibly(simulate_meta, otherwise = data.frame()),
                         covariates = covariates, studyID = studyID, 
                         subgroup = subgroup,
                         sample_sizes = sample_sizes,
                         beta = beta,
                         tau_covariate = tau_covariate, tau_corr = tau_corr, tau_scalar = tau_scalar,
                         min_studies = 4, max_effects = 12,
                         mods = c("robu","CHE","SCE"),
                         r_imputed = 0.6,
                         .parallel = TRUE)
)

tm 
parallel::stopCluster(cluster)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(tm, params, results, session_info, run_date, file = "RVE-simulation-dvcat-results.Rdata")
