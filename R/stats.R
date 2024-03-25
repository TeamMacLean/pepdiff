

#' do a power analysis for each peptide
#'
#' takes the treatment and control values for a peptide and
#' assesses the probability of finding a difference if one exists at the
#' calculated Cohen's d effect size, the standard deviation and sample size for
#' these data. The power required for a two-tailed t-test, Cohen's d and an
#' an estimate of the number of replicates required to reach a given proabability
#' of detecting a difference at that effect size are returned
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @param sig_level p-value expected to be used in statistical tests
#' @param b minimal power of test expected
#' @param max_n maximum number of reps to try when searching for replicates needed to reach power of b
#' @export
get_power <- function(treatment, control, sig_level=0.05, b=0.8, max_n=50) {

  n <- min( c(ncol(treatment), ncol(control)) )
  ncols <- dim(treatment)[1]
  ds <- rep(NA, ncols)

  powers <- rep(NA, ncols)
  target_reps <- rep(NA,ncols)


  for (i in 1:ncols){
    ds[i] <- lsr::cohensD(treatment[i,], control[i,])
    powers[i] <- pwr::pwr.t.test(n = n, d = ds[i], sig.level = sig_level)$power

    target_reps[i] <-
      tryCatch(
        {pwr::pwr.t.test(d=ds[i], sig.level = sig_level, power=b, type="one.sample",alternative="two.sided")$n},
        warning = function(w){ return(NA) },
        error = function(e){
          return( NA )
        },
        finally = {}
      )
      target_reps[target_reps > max_n] <- max_n

  }

  return(list(min_reps = target_reps,
              d = ds,
              power = powers
              ))
}




#' Calculate p-values based on iterative resampling using the lowest observed value replacement.
#'
#' This function computes p-values for each peptide based on iterative resampling with replacement. For each peptide, the function generates synthetic datasets by resampling the control group values while replacing the lowest observed value. It then compares the treatment group's mean to the distribution of means from the synthetic datasets to estimate p-values.
#'
#' @param treatment A dataframe representing the treatment group with quantitative values for each peptide.
#' @param control A dataframe representing the control group with quantitative values for each peptide.
#' @param iters The number of iterations to perform for each peptide.
#'
#' @return A dataframe containing p-values based on the lowest observed value replacement approach.
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' # Example using two dataframes for treatment and control with 1000 iterations:
#' treatment_data <- data.frame(peptide1 = c(10, 20, 30), peptide2 = c(15, 25, 35))
#' control_data <- data.frame(peptide1 = c(12, 22, 32), peptide2 = c(18, 28, 38))
#' result <- get_percentile_lowest_observed_value_iterative(treatment_data, control_data, iters = 1000)
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname get_percentile_lowest_observed_value_iterative
get_percentile_lowest_observed_value_iterative <- function(treatment, control, iters = 1000){

  peptide_means <- apply(control, MARGIN = 1, mean, na.rm = TRUE)
  peptide_sds <- apply(control, MARGIN = 1, sd, na.rm = TRUE)

  ## control samples
  nobs <- dim(control)[1]
  m <- matrix(NA, nrow = nobs, ncol = iters)
  for (i in 1:iters) {
    m[,i] <-  rnorm(nobs, peptide_means, peptide_sds)
  }

  #vector of treatment means
  t_means <- apply(treatment, MARGIN = 1, mean, na.rm = TRUE)

  find_p <- function(x, perc){ecdf(x)(perc)}
  r <- rep(NA, nobs)
  for (i in 1:nobs ){
    r[i] <- find_p(m[i,], t_means[i])
  }
  return(data.frame(norm_quantile_p_val = r ))
}

#' Calculate Bootstrap t-test p-values and adjusted false discovery rates for peptide comparisons.
#'
#' This function computes p-values and adjusted false discovery rates (FDR) using the Bootstrap t-test for each peptide comparison between treatment and control groups. The Bootstrap t-test is used to estimate the sampling distribution of the t-statistic by resampling the data, making it a robust method for identifying differences between groups.
#'
#' @param treatment A dataframe representing the treatment group with quantitative values for each peptide.
#' @param control A dataframe representing the control group with quantitative values for each peptide.
#' @param iters The number of bootstrap iterations to perform for each peptide comparison.
#'
#' @return A dataframe containing Bootstrap t-test p-values and corresponding Bonferroni-adjusted FDR values for each peptide comparison.
#'
#' @import MKinfer
#'
#' @examples
#' \dontrun{
#' # Example using two dataframes for treatment and control with 1000 bootstrap iterations:
#' treatment_data <- data.frame(peptide1 = c(10, 20, 30),
#' peptide2 = c(15, 25, 35))
#' control_data <- data.frame(peptide1 = c(12, 22, 32),
#' peptide2 = c(18, 28, 38))
#' result <- get_bootstrap_percentile(treatment_data, control_data,
#' iters = 1000)
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname get_bootstrap_percentile
get_bootstrap_percentile <- function(treatment, control, iters){


  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    result[i] <- tryCatch(
      {MKinfer::boot.t.test(treatment[i,], control[i,], R = iters)$p.value},
      warning = function(w){ return(NA) },
      error = function(e){
        return( NA )
      },
      finally = {})

  }
  return(data.frame(bootstrap_t_p_val = result, bootstrap_t_fdr = p.adjust(result, method = 'bonferroni')))
}

#' Calculate Wilcoxon rank-sum test p-values and adjusted false discovery rates for peptide comparisons.
#'
#' This function computes Wilcoxon rank-sum test (Mann-Whitney U test) p-values for each peptide comparison between treatment and control groups. It also calculates adjusted false discovery rates (FDR) using the Bonferroni method. The Wilcoxon test assesses whether there are significant differences in the distributions of quantitative values between two independent groups.
#'
#' @param treatment A dataframe representing the treatment group with quantitative values for each peptide.
#' @param control A dataframe representing the control group with quantitative values for each peptide.
#'
#' @return A dataframe containing Wilcoxon rank-sum test p-values and corresponding Bonferroni-adjusted FDR values for each peptide comparison.
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' # Example using two dataframes for treatment and control:
#' treatment_data <- data.frame(peptide1 = c(10, 20, 30),
#' peptide2 = c(15, 25, 35))
#' control_data <- data.frame(peptide1 = c(12, 22, 32),
#' peptide2 = c(18, 28, 38))
#' result <- get_wilcoxon_percentile(treatment_data, control_data)
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname get_wilcoxon_percentile
get_wilcoxon_percentile <- function(treatment, control){

  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      {wcox <- wilcox.test(treatment[i,], control[i,])
      result[i] <- wcox$p.value
      },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

  }
  return(data.frame(wilcoxon_p_val = result, wilcoxon_fdr = p.adjust(result, method = 'bonferroni')))
}

#' Calculate Kruskal-Wallis test p-values and adjusted false discovery rates for 2 peptide comparisons.
#'
#' This function computes Kruskal-Wallis test p-values for each peptide comparison between treatment and control groups. It also calculates adjusted false discovery rates (FDR) using the Bonferroni method. The Kruskal-Wallis test is used to assess differences in the distributions of quantitative values across multiple groups.
#'
#' @param treatment A dataframe representing the treatment group with quantitative values for each peptide.
#' @param control A dataframe representing the control group with quantitative values for each peptide.
#'
#' @return A dataframe containing Kruskal-Wallis test p-values and corresponding Bonferroni-adjusted FDR values for each peptide comparison.
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' # Example using two dataframes for treatment and control:
#' treatment_data <- data.frame(peptide1 = c(10, 20, 30),
#' peptide2 = c(15, 25, 35))
#' control_data <- data.frame(peptide1 = c(12, 22, 32),
#' peptide2 = c(18, 28, 38))
#' result <- get_kruskal_percentile(treatment_data, control_data)
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname get_kruskal_percentile
get_kruskal_percentile <- function(treatment, control){

  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      { d <- data.frame(t = treatment[i,], c = control[i,])
      kw <- stats::kruskal.test(t ~ c, data = d) },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

    result[i] <- kw$p.value
  }
  return(data.frame(kruskal_p_val = result, kruskal_fdr = p.adjust(result, method = 'bonferroni')))
}

#' Calculate Rank Products test p-values and adjusted false discovery rates for peptide comparisons.
#'
#' This function computes p-values and adjusted false discovery rates (FDR) using the Rank Products (RankProd) statistical test for each peptide comparison between treatment and control groups. The Rank Products test is specifically designed for identifying differentially abundant features (e.g., peptides) in high-throughput data.
#'
#' @param treatment A dataframe representing the treatment group with quantitative values for each peptide.
#' @param control A dataframe representing the control group with quantitative values for each peptide.
#'
#' @return A dataframe containing Rank Products test p-values and corresponding adjusted FDR values for each peptide comparison each way.
#'
#' @import RankProd
#'
#' @examples
#' \dontrun{
#' # Example using two dataframes for treatment and control:
#' treatment_data <- data.frame(peptide1 = c(10, 20, 30),
#' peptide2 = c(15, 25, 35))
#' control_data <- data.frame(peptide1 = c(12, 22, 32),
#' peptide2 = c(18, 28, 38))
#' result <- get_rp_percentile(treatment_data, control_data)
#' }
#'
#' @keywords data
#' @family statistical analysis
#' @rdname get_rp_percentile
#' @export
get_rp_percentile <- function(treatment, control){
  nt <- dim(treatment)[2]
  nc <-  dim(control)[2]
  cl <- c(rep(1, nt), rep(0, nc) )
  d <- cbind(treatment, control)
  invisible(capture.output(r <- {RankProd::RankProducts(d, cl, logged = FALSE, plot = FALSE, na.rm = TRUE)}))
  return(
    data.frame(
      rank_prod_p2_p_val = r$pval[, 1],
      rank_prod_p1_p_val = r$pval[, 2],
      rank_prod_p2_fdr = r$pfp[, 1],
      rank_prod_p1_fdr = r$pfp[, 2]
    )
  )

}

#' Calculate Gamma Values for Peptide Data
#'
#' This function computes gamma values for each peptide based on a treatment and control group.
#' Compares a null and a full model a la sleuth for RNAseq.
#'
#' @param treatment A matrix or data frame representing the treatment group's data.
#' @param control A matrix or data frame representing the control group's data.
#'
#' @return A data frame containing gamma p-values and FDR-adjusted gamma values for each peptide.
#'
#' @details The function fits a Generalized Linear Model (GLM) using the Gamma family to compare treatment and control groups for each peptide.
#'
#' @references
#' Davison, A. C., & Hinkley, D. V. (1997). Bootstrap Methods and Their Application.
#'
#' @examples
#' \dontrun{
#' treatment <- matrix(rnorm(100), ncol = 5)
#' control <- matrix(rnorm(100), ncol = 5)
#' results <- get_gamma(treatment, control)
#' }
#' @seealso
#' \code{\link{glm}}, \code{\link{p.adjust}}
#'
get_gamma <- function(treatment, control) {

  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    df <- data.frame(
      x = c( rep("treatment", ncol(treatment) ),
             rep("control", ncol(control) )
             ),
      y = c(treatment[i,], control[i,])
    )

      full <- glm(y ~ x, data=df, family = Gamma)
      int_only <- glm(y ~ 1, data=df, family = Gamma)

      result[i] <- log_likelihood_p(int_only, full)

      #fit <- glm(y ~ x, data=df, family = Gamma)
      #result[i] <- coef(summary(fit, dispersion = est$parameters[1]))[2,4]

      ##TODO
      ## Above computes a dispersion corrected gamma glm and looks at the co-efficient. Instead compare with a null
      ## model. Intercept only model.
      ## Check out how to do likelihood ratio test between two models if model from here is better than simple intercept model then the difference is likely significant! CF limma. Sleuth.

  }

  data.frame(gamma_p_val = result,
             gamma_fdr = p.adjust(result, method = "fdr")
             )
}


get_log <- function(treatment, control) {
  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    df <- data.frame(
      x = c( rep("treatment", ncol(treatment) ),
             rep("control", ncol(control) )
      ),
      y = c(treatment[i,], control[i,])
    )

  }
}


log_likelihood_p <- function(null, full) {
  a <- logLik(null)
  b <- logLik(full)
  st <- -2 * (as.numeric(a) - as.numeric(b))

  pchisq(st, df=1, lower.tail=FALSE) #df is number of variables in full - number of variables in null ( 1 - 0 in this.)
}



#' Calculate Empirical Bayes Statistics for Treatment vs. Control Comparison
#'
#' This function computes empirical Bayes statistics for comparing treatment and control groups.
#'
#' @param treatment A matrix or data frame representing the treatment group's data.
#' @param control A matrix or data frame representing the control group's data.
#'
#' @return A data frame containing empirical Bayes p-values, FDR-adjusted p-values, and Vash statistics.
#'
#' @details The function fits a linear model using limma::lmFit and performs empirical Bayes moderation with limma::eBayes.
#' The Vashr package is used to compute additional statistics for the comparison.
#'
#' @references
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015).
#' limma powers differential expression analyses for RNA-sequencing and microarray studies.
#'
#' @examples
#' \dontrun{
#' treatment <- matrix(rnorm(100), ncol = 5)
#' control <- matrix(rnorm(100), ncol = 5)
#' results <- get_eb(treatment, control)
#' }
#'
#' @export
get_eb <- function(treatment, control) {
  m <- log(cbind(treatment, control))
  tc <- c(rep("treatment", ncol(treatment)), rep("control", ncol(control)))
  design <- model.matrix(~factor(tc))
  colnames(design) <- c("trt", "ctrl")
  fit <- limma::lmFit(m, design)
  fit <- limma::eBayes(fit)
  betahat <- fit$coefficients[, 2]
  sehat <- fit$stdev.unscaled[, 2] * fit$sigma
  fit.vash <- vashr::vash(sehat = sehat, df = fit$df.residual[1], betahat = betahat)

  data.frame(
    eb_vs_p_val = fit.vash$pvalue,
    eb_vs_fdr = fit.vash$qvalue,
    eb_p_val = fit$p.value[, 2],
    eb_fdr = p.adjust(fit$p.value[, 2], method = "fdr")
  )
}

#' Perform kmeans of a dataset using just data in selected comparisons columns, then return matrices of all clusters
#' @param l list of results, usually from `compare_many()`
#' @param cols names of comparison columns to perform the k-means with. Default is all comparisons done.
#' @param log whether to log the data default = TRUE
#' @param base base used in logging (default = 2)
#' @param sig_only return only rows with 1 or more values significant at `sig_level` of `metric`
#' @param sig_level significance level cutoff
#' @param metric the test metric used to determine significance one of:
#' `bootstrap_t_p_val`, `bootstrap_t_fdr`
#' `wilcoxon_p_val`, `wilcoxon_fdr`
#' `kruskal_p_val`,  `kruskal_fdr`
#' `rank_prod_p1_p_val`, `rank_prod_p2_p_val`, `rank_prod_p1_fdr`, `rank_prod_p2_fdr`.
#' @param k number of clusters to make
#' @param nstart nstart value for `kmeans()`
#' @param itermax  number of `kmeans()` iterations (1000)
#' @return list of matrices
#' @export
kmeans_by_selected_cols <- function(l, cols=NULL, log=TRUE, base=2, sig_only=TRUE, sig_level=0.05, metric="bootstrap_t_p_val", k=NULL, nstart=25, itermax=1000) {


  fcm <- fold_change_matrix(l, log=log, base=base, sig_only = sig_only, sig_level=sig_level, metric=metric)
  if(is.null(cols)){
    cols = colnames(fcm)
  }

  keep <- fcm[,cols]
  km <- kmeans(keep, k, nstart = nstart, iter.max=itermax)

  kl <- lapply(1:max(km$cluster), function(x){
    idx <- which(km$cluster == x)
    rows <- names(km$cluster)[idx]
    fcm[rows,]
  })
  names(kl) <- 1:max(km$cluster)
  return(kl)

}
