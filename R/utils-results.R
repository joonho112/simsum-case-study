# Shared utility functions for analyzing simulation results
# Used by Chapters 11, 12, 13, and 14

#' Aggregate batch-level results into condition-level summaries
#'
#' The P&T simulation stores results per batch (34 batches per condition).
#' This function computes weighted averages across batches, using iteration
#' counts as weights. Note: df_SD is pooled using the law of total variance
#' (within-batch variance + between-batch variance of means), not a simple
#' weighted mean of per-batch SDs.
aggregate_batches <- function(df) {
  # Step 1: Compute the pooled mean of df_M for use in the variance calculation
  df |>
    dplyr::group_by(n_studies, struct, tau, omega, rho, nu, mod, coef) |>
    dplyr::summarize(
      total_iters = sum(iterations),
      convergence = weighted.mean(convergence, iterations),
      bias = weighted.mean(bias, iterations),
      rmse = sqrt(weighted.mean(rmse^2, iterations)),
      # Pool df_SD via law of total variance:
      # Var(df) = E[Var_within] + Var[E_within]
      df_SD = sqrt(
        weighted.mean(df_SD^2, iterations) +
        weighted.mean((df_M - weighted.mean(df_M, iterations))^2, iterations)
      ),
      df_M = weighted.mean(df_M, iterations),
      cover = weighted.mean(cover, iterations),
      width = weighted.mean(width, iterations),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Approximate MCSEs from aggregated summaries
      # These are approximate because we don't have iteration-level data
      sd_est = sqrt(pmax(rmse^2 - bias^2, 0)),
      mcse_bias = sd_est / sqrt(total_iters),
      mcse_cover = sqrt(cover * (1 - cover) / total_iters)
    )
}

#' Add readable labels for DGP and estimation method
add_labels <- function(df) {
  df |>
    dplyr::mutate(
      dgp = dplyr::case_when(
        struct == "CS" & omega == 0.00 ~ "CE",
        struct == "CS" & omega == 0.10 ~ "CHE",
        struct == "HCS" & omega == 0.00 ~ "CMVE"
      ),
      dgp = factor(dgp, levels = c("CE", "CHE", "CMVE")),
      method = dplyr::case_when(
        mod == "robu" ~ "CE",
        mod == "CHE" ~ "CHE",
        mod == "MV-DIAG" ~ "SCE"
      ),
      method = factor(method, levels = c("CE", "CHE", "SCE")),
      n_lab = paste0("J = ", n_studies)
    )
}

#' Compute relative RMSE ratios (CHE/CE and SCE/CE)
compute_ratios <- function(df) {
  df |>
    dplyr::select(n_studies, n_lab, struct, tau, omega, rho, dgp, coef, mod, rmse) |>
    tidyr::pivot_wider(names_from = mod, values_from = rmse, names_prefix = "rmse_") |>
    dplyr::mutate(
      ratio_CHE = rmse_CHE / rmse_robu,
      ratio_SCE = `rmse_MV-DIAG` / rmse_robu
    )
}

#' Load and prepare both result datasets
load_and_prepare_results <- function(dvcat_path = "data/RVE-simulation-dvcat-results.Rdata",
                                      ctype_path = "data/RVE-simulation-Ctype-results.Rdata") {
  load(dvcat_path, envir = e <- new.env())
  dvcat_agg <- aggregate_batches(e$results) |> add_labels()

  load(ctype_path, envir = e2 <- new.env())
  ctype_agg <- aggregate_batches(e2$results) |> add_labels()

  list(dvcat = dvcat_agg, ctype = ctype_agg)
}

#' Consistent color palette for estimation methods
method_colors <- c(CE = "#E69F00", CHE = "#56B4E9", SCE = "#009E73")
