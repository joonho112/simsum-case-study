#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------



generate_smd <- function(delta, k, N, Sigma) {
  
  # make sure delta is a vector
  delta_vec <- rep(delta, length.out = k)
  
  # create Sigma matrix assuming equicorrelation
  if (!is.matrix(Sigma)) Sigma <- Sigma + diag(1 - Sigma, nrow = k) # cor matrix for 1 study
  
  # generate numerator of SMD 
  mean_diff <- rmvnorm(n = 1, mean = delta_vec, sigma = (4/N) * Sigma) 
  
  # covariance 
  cov_mat <- as.matrix(rWishart(n = 1, df = N - 2, Sigma = Sigma)[,,1])
  sigma_sq <- diag(cov_mat) / (N - 2)
  
  # SMD
  d <- as.vector(mean_diff / sqrt(sigma_sq))  # cohen's d 
  J <- (1 - (3/((4 * (N - 2)) - 1)))
  g <- d * J # Hedges g
  var_g <- J^2 * (4 / N + d^2 / (2 * (N - 2)))
  
  dat <- tibble(g = g, var_g = var_g)
  
  return(dat)
}



generate_meta <- function(covariates,
                          studyID,
                          subgroup = NULL,
                          sample_sizes = pmin(20 + 2 * rpois(length(unique(studyID)), 20), 200),
                          beta = rep(0, ncol(covariates)),
                          tau = 0, omega = 0, rho = 0.6, nu = 50,
                          tau_covariate = NULL, tau_corr = NULL, tau_scalar = 1,
                          return_study_params = FALSE,
                          seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  studyID_vec <- levels(studyID)
  n_studies <- length(studyID_vec)
  n_ES_total <- nrow(covariates)
  
  # build regression
  
  if (is.null(tau_covariate) | is.null(tau_corr)) {
    
    # random intercept at study level
    v_j <- rnorm(n_studies, 0, tau)[studyID]
    
  } else {
    
    # heteroskedastic compound-symmetric random effects at study level
    tau_covariate <- as.factor(tau_covariate)
    tau_dim <- nlevels(tau_covariate)
    tau_scaled <- tau * rep(tau_scalar, length.out = tau_dim)
    Sigma_tau <- tcrossprod(tau_scaled) * (tau_corr + diag(1 - tau_corr, nrow = tau_dim))
    v_j <- rmvnorm(n_studies, mean = rep(0, tau_dim), sigma = Sigma_tau)
    z_j <- split(tau_covariate, studyID)
    v_list <- list()
    for (j in 1:n_studies) {
      v_list[[j]] <- v_j[j,][z_j[[j]]]
    }
    v_j <- unsplit(v_list, studyID)
    
  }
  
  # within-study random effect
  u_ij <- rnorm(n_ES_total, 0, omega)
  
  # meta-regression
  Xbeta <- as.numeric(covariates %*% beta)
  
  # build study-level parameters
  study_data <- 
    tibble(
      delta = Xbeta + v_j + u_ij,
      studyID = studyID
    ) %>%
    group_by(studyID) %>% 
    summarize(
      delta = list(delta),
      k = n(),
      .groups = "drop"
    ) %>%
    mutate(
      N = sample_sizes[studyID_vec],
      Sigma = rbeta(n_studies, rho * nu, (1 - rho) * nu)
    )

  if (return_study_params) return(study_data)
  
  # Generate full meta data
  
  X_data <- 
    as.data.frame(covariates) %>%
    mutate(
      studyID = studyID,
      esID = 1:n_ES_total,
      subgroup = subgroup
    ) %>%
    group_nest(studyID, .key = "X")
  
  meta_reg_dat <- 
    study_data %>%
    mutate(
      smds = pmap(select(., -studyID), generate_smd)
    ) %>%
    left_join(X_data, by = "studyID") %>%
    select(-delta, -k, -N, -Sigma) %>%
    unnest(cols = c(smds, X))
  
  return(meta_reg_dat)
}


#------------------------------------------------------
# Handmade robu
#------------------------------------------------------

# split a matrix
mat_split <- function(mat, f) {
  lapply(unique(f), function(l) mat[f==l,,drop=FALSE])
}

# trace product 
trace_product <- function(A, B) {
  
  a_vec <- as.vector(t(A))
  b_vec <- as.vector(B)
  sum(a_vec * b_vec)
  
}

robu_handmade <- function(X, y, v, cluster, rho = .8, calc_vcov = NULL) {
  
  k_j <- as.numeric(table(cluster))
  sigma_sq_j <- tapply(v, cluster, mean)
  
  m <- length(k_j)
  
  w_tilde <- 1 / (k_j * sigma_sq_j)
  w_tilde_j <- rep(w_tilde, k_j)
  
  mod_prelim <- lm.wfit(x = X, y = y, w = w_tilde_j)
  
  res <- residuals(mod_prelim) 
  
  # calculate weighted residual sum of squares
  QE <- sum(w_tilde_j * res^2)
  
  # split the design matrix by study
  X_j <- mat_split(X, cluster)
  
  # Create M tilde ----------------------------------------------------------
  M_tilde <- chol2inv(chol(crossprod(X, w_tilde_j * X)))
  p_j <- lapply(X_j, colSums)
  
  # trace products ----------------------------------------------------------
  
  # the first B for the trace product numerator
  w_over_k <- w_tilde / k_j
  B_num_all_1 <- Map(function(w, X) w * crossprod(X), w = w_over_k, X = X_j)
  B_num_1 <- reduce(B_num_all_1, `+`)
  
  # the second B for the trace product numerator
  XJX_j <- lapply(p_j, tcrossprod)
  B_num_all_2 <- Map(function(w, xjx) w * xjx, w = w_over_k, xjx = XJX_j)
  B_num_2 <- reduce(B_num_all_2, `+`)
  
  # the B for the trade product denominator
  B_den_all <- Map(function(w, xjx) w^2 * xjx, w = w_tilde, xjx = XJX_j)
  B_den <- reduce(B_den_all, `+`)
  
  num_minus <- m - (1 - rho) * trace_product(M_tilde, B_num_1) - rho * trace_product(M_tilde, B_num_2)
  den <- sum(k_j * w_tilde) - trace_product(M_tilde, B_den)
  
  tau_sq <- (QE - num_minus) / den
  tau_sq <- ifelse(tau_sq < 0, 0, tau_sq)  # added this like robu
  
  # CE weights
  w_j <- 1 / (k_j * (sigma_sq_j + tau_sq))
  w_ij <- rep(w_j, k_j)
  
  # fit WLS regression
  mod_CE <- lm.wfit(x = X, y = y, w = w_ij)
  
  res <- mod_CE[c("coefficients","residuals","fitted.values","weights")]
  res$tau_sq <- tau_sq
  
  res$X <- X
  res$w_ij <- w_ij
  res$cluster <- cluster
  res$k_j <- k_j
  res$sigma_sq_j <- sigma_sq_j
  res$w_tilde_j <- w_tilde_j
  res$num_minus <- num_minus
  res$den <- den
  res$nobs <- sum(k_j)
  
  class(res) <- "handmade"
  
  if (!is.null(calc_vcov)) {
    res$vcov <- vcovCR(res, cluster = cluster, type = calc_vcov, inverse_var = TRUE)
  }
  
  return(res)
  
}

model.matrix.handmade <- function(object, ...) object$X

bread.handmade <- function(x, ...) {
  x$nobs * chol2inv(chol(crossprod(x$X, x$w_ij * x$X)))
}

#------------------------------------------------------
# Estimation methods
#------------------------------------------------------

estimate_robu <- function(dat, formula_string, r) {

  X <- model.matrix(as.formula(formula_string), data = dat)
  
  mod <- robu_handmade(X = X, y = dat$g, v = dat$var_g, 
                        cluster = dat$studyID, 
                        rho = r)  
          
  varcomp <- tibble(coef = "tau2", est = mod$tau_sq)

  coef_CI <- conf_int(mod, cluster = dat$studyID, vcov = "CR2", inverse_var = TRUE)
  
  tibble(
    coef = rownames(coef_CI),
    est = coef_CI$beta,
    df = coef_CI$df,
    ci_lo = coef_CI$CI_L,
    ci_hi = coef_CI$CI_U
  ) %>%
    bind_rows(varcomp) %>%
    mutate(
      mod = "robu",
      cnvg = TRUE
    )
}

estimate_CHE <- function(dat, formula_string, r, smooth_vi = TRUE) {
  
  V_mat <- impute_covariance_matrix(dat$var_g, 
                                    cluster = dat$studyID, 
                                    r = r, 
                                    smooth_vi = smooth_vi)
  
  mod <- rma.mv(as.formula(formula_string),
                V = V_mat, 
                random = ~ 1 | studyID / esID,
                data = dat, sparse = TRUE)
  
  varcomp <- tibble(coef = c("tau2","omega2"), est = mod$sigma2)
  coef_CI <- conf_int(mod, vcov = "CR2")
  
  tibble(
    coef = rownames(coef_CI),
    est = coef_CI$beta,
    df = coef_CI$df,
    ci_lo = coef_CI$CI_L,
    ci_hi = coef_CI$CI_U
  ) %>%
    bind_rows(varcomp) %>%
    mutate(
      mod = "CHE",
      cnvg = TRUE
    )
  
}

starting_taus <- function(dat, V_mat, formula_string) {
 
  f <- function(sub, frml) suppressWarnings(
    rma.mv(as.formula(frml), V = V_mat, random = ~ 1 | studyID, 
           data = dat, sparse = TRUE, subset = subgroup == sub)$sigma2
  ) 
  
  tau2_starts <- lapply(levels(dat$subgroup), safely(f, NA), frml = formula_string)
  tau2_starts <- sapply(tau2_starts, function(x) x$result)
  
  if (any(is.na(tau2_starts))) { 
    tau2_NULLs <- lapply(levels(dat$subgroup)[is.na(tau2_starts)], safely(f, NA), frml = "g ~ 1")
    tau2_NULLs <- sapply(tau2_NULLs, function(x) x$result)
    
    tau2_starts[is.na(tau2_starts)] <- tau2_NULLs
  }
  
  tau2_starts
  
    
  
}

estimate_MV <- function(dat, formula_string, r, struct = "HCS", smooth_vi = TRUE, pre_estimate = FALSE) {
  
  sub <- if (struct == "DIAG") dat$subgroup else NULL
    
  V_mat <- impute_covariance_matrix(dat$var_g,
                                    cluster = dat$studyID,
                                    r = r,
                                    subgroup = sub,
                                    smooth_vi = smooth_vi)
    

  if (pre_estimate) {
    
    tau2_starts <- starting_taus(dat, V_mat, formula_string)
      
    mod <- safely(rma.mv)(as.formula(formula_string),
                          V = V_mat, 
                          random = ~ subgroup | studyID, struct = struct,
                          data = dat, sparse = TRUE, 
                          control = list(tau2.init = tau2_starts))
    
    if (is.null(mod$result)) {
      mod <- rma.mv(as.formula(formula_string),
                    V = V_mat, 
                    random = ~ subgroup | studyID, struct = struct,
                    data = dat, sparse = TRUE,
                    tau2 = tau2_starts)
    } else {
      mod <- mod$result
    }
    
  } else {
    mod <- rma.mv(as.formula(formula_string),
                  V = V_mat, 
                  random = ~ subgroup | studyID, struct = struct,
                  data = dat, sparse = TRUE)
  }
  
  varcomp <- tibble(
    coef = paste("tau2", names(mod$g.levels.k), sep = "_"), 
    est = mod$tau2
  )
  
  coef_CI <- conf_int(mod, vcov = "CR2")
  
  tibble(
    coef = rownames(coef_CI),
    est = coef_CI$beta,
    df = coef_CI$df,
    ci_lo = coef_CI$CI_L,
    ci_hi = coef_CI$CI_U
  ) %>%
    bind_rows(varcomp) %>%
    mutate(
      mod = paste("MV", struct, sep = "-"),
      cnvg = TRUE
    )
}

modify_res <- function(res, new_mod, subgroup_levels) {
  
  varcomp <- 
    res %>%
    filter(coef == "tau2") %>%
    select(-coef) %>%
    crossing(coef = paste("tau2", subgroup_levels, sep = "_"))
    
  res %>%
    filter(! coef %in% c("tau2","omega2")) %>%
    bind_rows(varcomp) %>%
    mutate(mod = new_mod, cnvg = FALSE)
  
}

estimate_models <- function(dat, formula_string, r = 0.6, smooth_vi = TRUE, 
                            mods = c("robu","CHE","SCE","CMVE")) {
  
  res <- tibble()
  
  if ("robu" %in% mods) {
    res <- 
      estimate_robu(dat, formula_string, r) %>%
      bind_rows(res)
  } 
  
  if ("CHE" %in% mods) {
    res_CHE <- estimate_CHE(dat, formula_string, r, smooth_vi = smooth_vi)
    
    res <- bind_rows(res, res_CHE)
  }
  
  
  if ("SCE" %in% mods) {
    res_backup <- modify_res(res_CHE, new_mod = "MV-DIAG", subgroup_levels = levels(dat$subgroup))
    
    res <- 
      possibly(estimate_MV, res_backup)(dat, formula_string, r, struct = "DIAG", smooth_vi = smooth_vi, pre_estimate = TRUE) %>%
      bind_rows(res)
  }
  
  if ("CMVE" %in% mods) {
    res_backup <- modify_res(res_CHE, new_mod = "MV-HCS", subgroup_levels = levels(dat$subgroup))
    
    res <- 
      possibly(estimate_MV, res_backup)(dat, formula_string, r, struct = "HCS", smooth_vi = smooth_vi, pre_estimate = TRUE) %>%
      bind_rows(res)
  }
  
  res  
}

#------------------------------------------------------
# Performance criteria
#------------------------------------------------------

assess_performance <- function(res, beta, tau, omega, subgroup_levels, tau_scalar) {
  
  true_beta <- tibble(
    coef = names(beta),
    param = beta
  )
  
  true_varcomp <- tibble(
    coef = c("tau2","omega2", paste("tau2", subgroup_levels, sep = "_")),
    param = c(tau^2, omega^2, tau^2 * rep(tau_scalar, length.out = length(subgroup_levels))^2)
  )
  
  bind_rows(true_beta, true_varcomp) %>%
    full_join(res, by = "coef") %>% 
    group_by(mod, coef) %>%
    mutate(
      err = est - param
    ) %>%
    summarise(
      iterations = n(),
      convergence = mean(cnvg),
      bias = mean(err),
      rmse = sqrt(mean(err^2)),
      df_M = mean(df),
      df_SD = sd(df),
      cover = mean(ci_lo < param & param < ci_hi),
      width = mean(ci_hi - ci_lo),
      .groups = "drop"
    )
}

#------------------------------------------------------
# Simulation driver
#------------------------------------------------------

generate_estimate <- function(iteration, 
                              n_studies,
                              tau, omega, rho, nu,
                              covariates, studyID, 
                              subgroup,
                              sample_sizes,
                              beta,
                              tau_covariate, tau_corr, tau_scalar,
                              formula_string, mods, r_imputed,
                              min_studies, max_effects = Inf) {
  
  isPD_subtb <- FALSE
  
  while (!isPD_subtb) {
    # draw sample of studies
    studyID_sample <- sample(levels(studyID), size = n_studies, replace = FALSE)
    
    # sample effect sizes from each study
    es_sample <- tapply(1:length(studyID), studyID, function(x) x %in% if (length(x) <= max_effects) x else sample(x, size = max_effects))
    es_sample <- unsplit(es_sample, studyID)
    
    # calculate indices of included studies and effect sizes
    sample_indices <- (studyID %in% studyID_sample) & es_sample
    
    sub_tb <- rowSums(table(subgroup[sample_indices], studyID[sample_indices]) > 0)
    isPD <- qr(covariates[sample_indices,])$rank == ncol(covariates)
    isPD_subtb <- isPD & all(sub_tb >= min_studies)
  }
  
  dat <- generate_meta(covariates = covariates[sample_indices,],
                       studyID = droplevels(studyID[sample_indices]),
                       subgroup = subgroup[sample_indices],
                       sample_sizes = sample_sizes[levels(studyID) %in% studyID_sample],
                       beta = beta,
                       tau = tau, omega = omega, rho = rho, nu = nu,
                       tau_covariate = tau_covariate[sample_indices], 
                       tau_corr = tau_corr, tau_scalar = tau_scalar)  
  
  estimate_models(dat, formula_string = formula_string, 
                  r = r_imputed, smooth_vi = TRUE, 
                  mods = mods)
}

simulate_meta <- function(iterations, 
                          n_studies,
                          struct = "CS",
                          tau = 0, omega = 0, rho = 0.6, nu = 50,
                          covariates, studyID, 
                          subgroup = NULL,
                          sample_sizes = pmin(20 + 2 * rpois(length(unique(studyID)), 20), 200),
                          beta = rep(0, ncol(covariates)),
                          tau_covariate = NULL, tau_corr = NULL, tau_scalar = 1,
                          mods = c("robu","CHE","SCE","CMVE"),
                          r_imputed = 0.6,
                          min_studies = 3, max_effects = Inf,
                          seed = NULL) {
  
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(mvtnorm, quietly = TRUE, warn.conflicts = FALSE)
  suppressPackageStartupMessages(require(metafor, quietly = TRUE, warn.conflicts = FALSE))
  # require(robumeta, quietly = TRUE, warn.conflicts = FALSE)
  require(clubSandwich, quietly = TRUE, warn.conflicts = FALSE)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (struct == "CS") {
    tau_covariate <- NULL
    tau_corr <- NULL
    tau_scalar <- 1
  }
  
  studyID <- as.factor(studyID)
  subgroup <- as.factor(subgroup)
  
  formula_string <- paste(colnames(covariates), collapse = " + ")
  formula_string <- paste("g ~ 0 +", formula_string)
  
  res <- map_dfr(1:iterations, 
                 .f = possibly(generate_estimate, otherwise = NULL), 
                 n_studies = n_studies,
                 tau = tau, omega = omega, rho = rho, nu = nu,
                 covariates = covariates, studyID = studyID, 
                 subgroup = subgroup,
                 sample_sizes = sample_sizes,
                 beta = beta,
                 tau_covariate = tau_covariate, 
                 tau_corr = tau_corr, tau_scalar = tau_scalar,
                 formula_string = formula_string, 
                 mods = mods, 
                 r_imputed = r_imputed, 
                 min_studies = min_studies, max_effects = max_effects)  
  
  assess_performance(res, beta = beta, tau = tau, omega = omega, 
                     subgroup_levels = levels(subgroup), 
                     tau_scalar = tau_scalar)
  
}
  

  
