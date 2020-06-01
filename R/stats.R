#' get p values for contrast using normal percentile
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @param iters number of iterations to perform
#' @return dataframe with one column `norm_quantile_pval`
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
  return(data.frame(norm_quantile_pval = r ))
}

#' get p values for contrast using boostrap t test
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @param iters number of iterations to perform
#' @return dataframe with two columns `bootstrap_t_p_val` and `bootstrap_t_fdr`
#' @export
get_bootstrap_percentile <- function(treatment, control, iters){


  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      {bstrap <- MKinfer::boot.t.test(treatment[i,], control[i,], R = iters)},
      warning = function(w){ NA },
      error = function(e){
        return( NA )
      },
      finally = {})

    result[i] <- bstrap$p.value
  }
  return(data.frame(bootstrap_t_p_val = result, bootstrap_t_fdr = p.adjust(result, method = 'bonferroni')))
}

#' get p values for contrast using Wilcoxon test
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @return dataframe with two columns `wilcoxon_p_val` and `wilcoxon_fdr`
#' @export
get_wilcoxon_percentile <- function(treatment, control){

  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      {wcox <- wilcox.test(treatment[i,], control[i,]) },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

    result[i] <- wcox$p.value
  }
  return(data.frame(wilcoxon_p_val = result, wilcoxon_fdr = p.adjust(result, method = 'bonferroni')))
}

#' get p values for contrast using Kruskal-Wallis test
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @return dataframe with two columns `kruskal_p_val` and `kruskal_fdr`
#' @export
get_kruskal_percentile <- function(treatment, control){

  peptide_count <- dim(control)[1]
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      { d <- data.frame(t = treatment[i,], c = control[i,])
      kw <- kruskal.test(t ~ c, data = d) },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

    result[i] <- kw$p.value
  }
  return(data.frame(kruskal_p_val = result, kruskal_fdr = p.adjust(result, method = 'bonferroni')))
}

#' get p values for contrast using Rank Products test
#'
#' @param treatment matrix of treatment data columns
#' @param control matrix of control data columns
#' @return dataframe with four columns, two for the test each way from RankProducts `rank_prod_p1_p_val`, `rank_prod_p2_p_val` and `rank_prod_p1_fdr`, `rank_prod_p2_fdr`.
#' @export
get_rp_percentile <- function(treatment, control){
  nt <- dim(treatment)[2]
  nc <-  dim(control)[2]
  cl <- c(rep(1, nt), rep(0, nc) )
  d <- cbind(treatment, control)
  r <- RankProd::RankProducts(d, cl, logged = FALSE, plot = FALSE, na.rm = TRUE)
  return(data.frame(rank_prod_p2_p_val = r$pval[,1], rank_prod_p1_p_val = r$pval[,2], rank_prod_p2_fdr = r$pfp[,1], rank_prod_p1_fdr = r$pfp[,2]))

}
