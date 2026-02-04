# Pairwise statistical tests for differential abundance analysis
# These functions match the implementations in peppwR for consistency

# =============================================================================
# Wilcoxon Rank-Sum Test (wrapper around stats::wilcox.test)
# =============================================================================

#' Wilcoxon rank-sum test for two groups
#'
#' Performs a two-sample Wilcoxon rank-sum test (Mann-Whitney U test)
#' to compare abundance values between control and treatment groups.
#'
#' @param control Numeric vector of control group values
#' @param treatment Numeric vector of treatment group values
#' @param ... Additional arguments passed to [stats::wilcox.test()]
#'
#' @return A list with components:
#'   \item{p_value}{The p-value from the test}
#'   \item{statistic}{The test statistic W}
#'   \item{method}{"wilcoxon"}
#'
#' @export
#' @examples
#' ctrl <- c(100, 120, 110, 105)
#' trt <- c(200, 220, 180, 210)
#' test_wilcoxon(ctrl, trt)
test_wilcoxon <- function(control, treatment, ...) {
  # Remove NAs
  control <- control[!is.na(control)]
  treatment <- treatment[!is.na(treatment)]

  # Need at least 2 observations per group
  if (length(control) < 2 || length(treatment) < 2) {
    return(list(
      p_value = NA_real_,
      statistic = NA_real_,
      method = "wilcoxon"
    ))
  }

  result <- tryCatch(
    {
      wt <- stats::wilcox.test(treatment, control, ...)
      list(
        p_value = wt$p.value,
        statistic = unname(wt$statistic),
        method = "wilcoxon"
      )
    },
    warning = function(w) {
      # Wilcoxon can warn about ties - still return result
      wt <- suppressWarnings(stats::wilcox.test(treatment, control, ...))
      list(
        p_value = wt$p.value,
        statistic = unname(wt$statistic),
        method = "wilcoxon"
      )
    },
    error = function(e) {
      list(
        p_value = NA_real_,
        statistic = NA_real_,
        method = "wilcoxon"
      )
    }
  )

  result
}


# =============================================================================
# Bootstrap t-test
# =============================================================================

#' Bootstrap t-test for two groups
#'
#' Performs a bootstrap-based t-test comparing two groups. This is more
#' robust than a standard t-test when assumptions of normality may not hold.
#'
#' @param control Numeric vector of control group values
#' @param treatment Numeric vector of treatment group values
#' @param n_boot Number of bootstrap iterations (default 1000)
#' @param seed Optional random seed for reproducibility
#'
#' @return A list with components:
#'   \item{p_value}{The two-sided p-value}
#'   \item{t_obs}{The observed t-statistic}
#'   \item{method}{"bootstrap_t"}
#'
#' @export
#' @examples
#' ctrl <- c(100, 120, 110, 105)
#' trt <- c(200, 220, 180, 210)
#' test_bootstrap_t(ctrl, trt, n_boot = 500)
test_bootstrap_t <- function(control, treatment, n_boot = 1000, seed = NULL) {
  # Remove NAs
  control <- control[!is.na(control)]
  treatment <- treatment[!is.na(treatment)]

  n_ctrl <- length(control)
  n_trt <- length(treatment)

  # Need at least 2 observations per group
  if (n_ctrl < 2 || n_trt < 2) {
    return(list(
      p_value = NA_real_,
      t_obs = NA_real_,
      method = "bootstrap_t"
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Calculate observed t-statistic
  t_obs <- calc_t_statistic(treatment, control)

  if (is.na(t_obs)) {
    return(list(
      p_value = NA_real_,
      t_obs = NA_real_,
      method = "bootstrap_t"
    ))
  }

  # Pool the data
  pooled <- c(control, treatment)
  n_total <- length(pooled)

  # Bootstrap under null hypothesis (no difference)
  t_boot <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    # Resample from pooled data
    boot_idx <- sample.int(n_total, n_total, replace = TRUE)
    boot_ctrl <- pooled[boot_idx[seq_len(n_ctrl)]]
    boot_trt <- pooled[boot_idx[(n_ctrl + 1):n_total]]
    t_boot[i] <- calc_t_statistic(boot_trt, boot_ctrl)
  }

  # Remove NAs from bootstrap samples
  t_boot <- t_boot[!is.na(t_boot)]

  if (length(t_boot) == 0) {
    return(list(
      p_value = NA_real_,
      t_obs = t_obs,
      method = "bootstrap_t"
    ))
  }

  # Two-sided p-value: proportion of bootstrap t-stats as extreme as observed
  p_value <- mean(abs(t_boot) >= abs(t_obs))

  list(
    p_value = p_value,
    t_obs = t_obs,
    method = "bootstrap_t"
  )
}


#' Calculate t-statistic for two groups
#'
#' @param x First group values
#' @param y Second group values
#' @return The t-statistic
#' @keywords internal
calc_t_statistic <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)

  if (n_x < 2 || n_y < 2) return(NA_real_)

  mean_x <- mean(x)
  mean_y <- mean(y)
  var_x <- stats::var(x)
  var_y <- stats::var(y)

  # Pooled standard error
  se <- sqrt(var_x / n_x + var_y / n_y)

  if (se == 0) return(NA_real_)

  (mean_x - mean_y) / se
}


# =============================================================================
# Bayes Factor t-test (JZS approximation)
# =============================================================================

#' Bayes factor t-test for two groups
#'
#' Computes a Bayes factor comparing the alternative hypothesis (group difference)
#' to the null hypothesis (no difference) using the JZS (Jeffreys-Zellner-Siow)
#' prior. Uses an analytical approximation for computational efficiency.
#'
#' @param control Numeric vector of control group values
#' @param treatment Numeric vector of treatment group values
#' @param r_scale Scale parameter for the Cauchy prior on effect size (default 0.707)
#'
#' @return A list with components:
#'   \item{bf}{Bayes factor (BF10) - evidence for alternative vs null}
#'   \item{effect_size}{Cohen's d effect size}
#'   \item{method}{"bayes_t"}
#'
#' @details
#' The Bayes factor is interpreted as:
#' - BF10 > 10: Strong evidence for difference
#' - BF10 > 3: Moderate evidence for difference
#' - BF10 0.33-3: Inconclusive
#' - BF10 < 0.33: Moderate evidence for no difference
#' - BF10 < 0.1: Strong evidence for no difference
#'
#' Unlike p-values, Bayes factors are NOT converted to pseudo-p-values.
#' Use [classify_bf_evidence()] to interpret BF values categorically.
#'
#' @export
#' @examples
#' ctrl <- c(100, 120, 110, 105)
#' trt <- c(200, 220, 180, 210)
#' test_bayes_t(ctrl, trt)
test_bayes_t <- function(control, treatment, r_scale = 0.707) {
  # Remove NAs
  control <- control[!is.na(control)]
  treatment <- treatment[!is.na(treatment)]

  n_ctrl <- length(control)
  n_trt <- length(treatment)

  # Need at least 2 observations per group
  if (n_ctrl < 2 || n_trt < 2) {
    return(list(
      bf = NA_real_,
      effect_size = NA_real_,
      method = "bayes_t"
    ))
  }

  # Calculate t-statistic and degrees of freedom
  t_stat <- calc_t_statistic(treatment, control)
  df <- n_ctrl + n_trt - 2

  if (is.na(t_stat)) {
    return(list(
      bf = NA_real_,
      effect_size = NA_real_,
      method = "bayes_t"
    ))
  }

  # Calculate effect size (Cohen's d)
  pooled_var <- ((n_ctrl - 1) * stats::var(control) + (n_trt - 1) * stats::var(treatment)) / df
  effect_size <- (mean(treatment) - mean(control)) / sqrt(pooled_var)

  # JZS Bayes factor approximation (Rouder et al., 2009)
  # Using the BIC approximation for computational efficiency
  n_eff <- (n_ctrl * n_trt) / (n_ctrl + n_trt)  # Effective sample size

  # BF10 approximation based on t-statistic
  bf <- jzs_bf_approx(t_stat, n_eff, r_scale)

  list(
    bf = bf,
    effect_size = effect_size,
    method = "bayes_t"
  )
}


#' JZS Bayes factor approximation
#'
#' Computes BF10 using the Savage-Dickey density ratio approximation.
#'
#' @param t_stat T-statistic
#' @param n_eff Effective sample size
#' @param r_scale Scale parameter for Cauchy prior
#' @return Bayes factor (BF10)
#' @keywords internal
jzs_bf_approx <- function(t_stat, n_eff, r_scale = 0.707) {
  # Savage-Dickey approximation for JZS Bayes factor
  # Based on Wetzels et al. (2012) approximation

  # Effect size estimate
  d <- t_stat / sqrt(n_eff)

  # Prior density at d = 0 (Cauchy)
  prior_at_0 <- 1 / (pi * r_scale)

  # Posterior density at d = 0 (approximately normal)
  se <- 1 / sqrt(n_eff)
  posterior_at_0 <- stats::dnorm(0, mean = d, sd = se)

  # BF01 = posterior at 0 / prior at 0
  bf01 <- posterior_at_0 / prior_at_0

  # BF10 = 1 / BF01
  bf10 <- 1 / bf01

  # Bound to reasonable range
  bf10 <- max(min(bf10, 1e10), 1e-10)

  bf10
}


#' Classify Bayes factor into evidence categories
#'
#' Converts numeric Bayes factors into categorical evidence levels following
#' conventional thresholds (Jeffreys, 1961; Lee & Wagenmakers, 2013).
#'
#' @param bf Numeric vector of Bayes factors (BF10)
#'
#' @return An ordered factor with levels:
#'   \item{strong_null}{BF < 0.1 - Strong evidence for null hypothesis}
#'   \item{moderate_null}{BF 0.1-0.33 - Moderate evidence for null}
#'   \item{inconclusive}{BF 0.33-3 - Evidence is inconclusive}
#'   \item{moderate_alt}{BF 3-10 - Moderate evidence for alternative}
#'   \item{strong_alt}{BF > 10 - Strong evidence for alternative}
#'
#' @export
#' @examples
#' classify_bf_evidence(c(0.05, 0.2, 1, 5, 20))
classify_bf_evidence <- function(bf) {
  evidence_levels <- c("strong_null", "moderate_null", "inconclusive",
                       "moderate_alt", "strong_alt")

  evidence <- dplyr::case_when(
    bf < 0.1 ~ "strong_null",
    bf < 1/3 ~ "moderate_null",
    bf < 3 ~ "inconclusive",
    bf < 10 ~ "moderate_alt",
    TRUE ~ "strong_alt"
  )

  factor(evidence, levels = evidence_levels, ordered = TRUE)
}


# =============================================================================
# Rank Products Test
# =============================================================================

#' Rank products test for two groups
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated because Rank Products requires ranking across
#' ALL peptides, not within single peptides. The per-peptide permutation
#' approach produces unreliable p-values.
#'
#' Use `compare()` with `test = "rankprod"` instead, which properly uses the

#' RankProd package to rank across all peptides.
#'
#' @param control Numeric vector of control group values
#' @param treatment Numeric vector of treatment group values
#' @param n_perm Number of permutations for p-value estimation (default 1000)
#' @param seed Optional random seed for reproducibility
#'
#' @return A list with components:
#'   \item{p_value_up}{P-value for upregulation (treatment > control)}
#'   \item{p_value_down}{P-value for downregulation (treatment < control)}
#'   \item{p_value}{Combined two-sided p-value (minimum of up/down)}
#'   \item{rp_up}{Rank product for upregulation}
#'   \item{rp_down}{Rank product for downregulation}
#'   \item{method}{"rankprod"}
#'
#' @export
#' @examples
#' ctrl <- c(100, 120, 110, 105)
#' trt <- c(200, 220, 180, 210)
#' \dontrun{
#' test_rankprod(ctrl, trt, n_perm = 100)  # Deprecated
#' }
test_rankprod <- function(control, treatment, n_perm = 1000, seed = NULL) {
  .Deprecated(
    msg = "test_rankprod() is deprecated. Use compare() with test='rankprod' instead, which properly uses the RankProd package to rank across all peptides."
  )
  # Remove NAs
  control <- control[!is.na(control)]
  treatment <- treatment[!is.na(treatment)]

  n_ctrl <- length(control)
  n_trt <- length(treatment)

  # Need at least 2 observations per group
  if (n_ctrl < 2 || n_trt < 2) {
    return(list(
      p_value_up = NA_real_,
      p_value_down = NA_real_,
      p_value = NA_real_,
      rp_up = NA_real_,
      rp_down = NA_real_,
      method = "rankprod"
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Calculate observed rank products
  rp_obs <- calc_rank_product(control, treatment)

  # Permutation test
  pooled <- c(control, treatment)
  n_total <- length(pooled)

  rp_up_perm <- numeric(n_perm)
  rp_down_perm <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    perm_idx <- sample.int(n_total, n_total)
    perm_ctrl <- pooled[perm_idx[seq_len(n_ctrl)]]
    perm_trt <- pooled[perm_idx[(n_ctrl + 1):n_total]]
    rp_perm <- calc_rank_product(perm_ctrl, perm_trt)
    rp_up_perm[i] <- rp_perm$rp_up
    rp_down_perm[i] <- rp_perm$rp_down
  }

  # P-values: proportion of permutations with RP as small as observed
  # (smaller RP = stronger evidence)
  p_value_up <- mean(rp_up_perm <= rp_obs$rp_up)
  p_value_down <- mean(rp_down_perm <= rp_obs$rp_down)

  # Combined p-value (two-sided)
  p_value <- 2 * min(p_value_up, p_value_down)
  p_value <- min(p_value, 1)  # Cap at 1

  list(
    p_value_up = p_value_up,
    p_value_down = p_value_down,
    p_value = p_value,
    rp_up = rp_obs$rp_up,
    rp_down = rp_obs$rp_down,
    method = "rankprod"
  )
}


#' Calculate rank products for two groups
#'
#' For a single peptide, calculates the rank product statistic based on
#' all pairwise fold changes between treatment and control replicates.
#'
#' @param control Control group values
#' @param treatment Treatment group values
#' @return List with rp_up and rp_down
#' @keywords internal
calc_rank_product <- function(control, treatment) {
  n_ctrl <- length(control)
  n_trt <- length(treatment)

  # Calculate all pairwise fold changes (treatment / control)
  fc_values <- as.vector(outer(treatment, control, `/`))
  n_fc <- length(fc_values)

  # Rank the fold changes
  # For upregulation: higher FC should get lower (better) rank
  ranks_up <- rank(-fc_values, ties.method = "average")
  # Geometric mean of ranks
  rp_up <- exp(mean(log(ranks_up)))

  # For downregulation: lower FC should get lower (better) rank
  ranks_down <- rank(fc_values, ties.method = "average")
  rp_down <- exp(mean(log(ranks_down)))

  list(rp_up = rp_up, rp_down = rp_down)
}


# =============================================================================
# Legacy Wrapper Functions (for backwards compatibility)
# =============================================================================

#' Legacy: Get bootstrap t-test p-values for matrix data
#'
#' @param treatment Matrix of treatment data (rows = peptides, cols = replicates)
#' @param control Matrix of control data
#' @param iters Number of bootstrap iterations
#' @return Data frame with bootstrap_t_p_val and bootstrap_t_fdr columns
#' @export
get_bootstrap_percentile <- function(treatment, control, iters = 1000) {
  peptide_count <- nrow(control)
  result <- numeric(peptide_count)

  for (i in seq_len(peptide_count)) {
    test_result <- test_bootstrap_t(control[i, ], treatment[i, ], n_boot = iters)
    result[i] <- test_result$p_value
  }

  data.frame(
    bootstrap_t_p_val = result,
    bootstrap_t_fdr = stats::p.adjust(result, method = "bonferroni")
  )
}


#' Legacy: Get Wilcoxon test p-values for matrix data
#'
#' @param treatment Matrix of treatment data
#' @param control Matrix of control data
#' @return Data frame with wilcoxon_p_val and wilcoxon_fdr columns
#' @export
get_wilcoxon_percentile <- function(treatment, control) {
  peptide_count <- nrow(control)
  result <- numeric(peptide_count)

  for (i in seq_len(peptide_count)) {
    test_result <- test_wilcoxon(control[i, ], treatment[i, ])
    result[i] <- test_result$p_value
  }

  data.frame(
    wilcoxon_p_val = result,
    wilcoxon_fdr = stats::p.adjust(result, method = "bonferroni")
  )
}


#' Legacy: Get Kruskal-Wallis test p-values for matrix data
#'
#' @param treatment Matrix of treatment data
#' @param control Matrix of control data
#' @return Data frame with kruskal_p_val and kruskal_fdr columns
#' @export
get_kruskal_percentile <- function(treatment, control) {
  peptide_count <- nrow(control)
  result <- numeric(peptide_count)

  for (i in seq_len(peptide_count)) {
    ctrl_vals <- control[i, ]
    trt_vals <- treatment[i, ]
    ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
    trt_vals <- trt_vals[!is.na(trt_vals)]

    if (length(ctrl_vals) >= 2 && length(trt_vals) >= 2) {
      kw_result <- tryCatch(
        {
          d <- data.frame(
            value = c(ctrl_vals, trt_vals),
            group = factor(c(rep("ctrl", length(ctrl_vals)), rep("trt", length(trt_vals))))
          )
          stats::kruskal.test(value ~ group, data = d)$p.value
        },
        error = function(e) NA_real_
      )
      result[i] <- kw_result
    } else {
      result[i] <- NA_real_
    }
  }

  data.frame(
    kruskal_p_val = result,
    kruskal_fdr = stats::p.adjust(result, method = "bonferroni")
  )
}


#' Legacy: Get Rank Products test p-values for matrix data
#'
#' @param treatment Matrix of treatment data
#' @param control Matrix of control data
#' @return Data frame with rank product p-values and FDR
#' @export
get_rp_percentile <- function(treatment, control) {
  # Use RankProd package if available, otherwise fall back to our implementation
  if (requireNamespace("RankProd", quietly = TRUE)) {
    nt <- ncol(treatment)
    nc <- ncol(control)
    cl <- c(rep(1, nt), rep(0, nc))
    d <- cbind(treatment, control)
    r <- RankProd::RankProducts(d, cl, logged = FALSE, plot = FALSE, na.rm = TRUE)
    return(data.frame(
      rank_prod_p2_p_val = r$pval[, 1],
      rank_prod_p1_p_val = r$pval[, 2],
      rank_prod_p2_fdr = r$pfp[, 1],
      rank_prod_p1_fdr = r$pfp[, 2]
    ))
  }

  # Fallback implementation
  peptide_count <- nrow(control)
  p_up <- numeric(peptide_count)
  p_down <- numeric(peptide_count)

  for (i in seq_len(peptide_count)) {
    test_result <- test_rankprod(control[i, ], treatment[i, ], n_perm = 1000)
    p_up[i] <- test_result$p_value_up
    p_down[i] <- test_result$p_value_down
  }

  data.frame(
    rank_prod_p2_p_val = p_down,
    rank_prod_p1_p_val = p_up,
    rank_prod_p2_fdr = stats::p.adjust(p_down, method = "BH"),
    rank_prod_p1_fdr = stats::p.adjust(p_up, method = "BH")
  )
}
