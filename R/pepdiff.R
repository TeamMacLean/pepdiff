

# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' read data from a file
#'
#' reads data, renames columns appropriately, discards unused columns, factors and
#' reorders, discards duplicate rows
#'
#' @param file Path to the file to load - must be a csv file
#' @param treatment Column containing the treatment of the observation
#' @param bio_rep Column containing the biological replicate of the observation
#' @param tech_rep Column containing the technical replicate of the observation
#' @param quant Column containing the quantitation data
#' @param seconds Column containing timepoint of observation
#' @param gene_id Column containing the id of the gene this hit
#' @param peptide_sequence Column containing the sequence of this peptide
#' @return tibble with columns id, gene_id, peptide, treatment, seconds, bio_rep, tech_rep, quant
#' @export
import_data <- function(file,
                        treatment = "genotype",
                        bio_rep = "bio_rep",
                        tech_rep = "tech_rep",
                        quant = "total_area",
                        seconds = "seconds",
                        gene_id = "gene_id",
                        peptide = "peptide_sequence"
                        ) {

  readr::read_csv(file
                    ) %>%
    dplyr::rename(treatment = treatment,
                  bio_rep = bio_rep,
                  tech_rep = tech_rep,
                  quant = quant,
                  seconds = seconds,
                  gene_id = gene_id,
                  peptide = peptide) %>%
    dplyr::mutate(id = paste0(gene_id, ":", peptide)) %>%
    dplyr::transmute(id = as.character(id),
                     gene_id = as.character(gene_id),
                     peptide = as.character(peptide),
                     treatment = as.character(treatment),
                     seconds = as.numeric(seconds),
                     bio_rep = as.character(bio_rep),
                     tech_rep = as.character(tech_rep),
                     quant = as.numeric(quant)
                     ) %>%
    dplyr::distinct()


}


#' calculate the proportion of peptides with missing values per group in a data set.
#'
#' Group the data by treatment, seconds, bio rep and tech rep, then calculate the percent
#' of NA in each group.
#'
#' @param df dataframe with unmerged tech reps; typically from `import_data()`
#' @return grouped summary dataframe
#' @export
assess_missing <- function(df){

  dplyr::group_by( df, treatment, seconds, bio_rep, tech_rep) %>%
    dplyr::summarize(peptides = dplyr::n(), percent_missing = (sum(is.na(quant)) / dplyr::n() ) * 100 )
}

#' calculate the number of non- NA or NaNs in a vector
#'
#' @param x vector of values
#' @return number
#'
is_useable <- function(x){ !is.na(x) & !is.nan(x)}

#' calculate number of measurements of each peptide in each treatment and time
#'
#' For each peptide, works out how many biologically replicated measurements are
#' available in the different combinations of treatment and seconds
#'
#' @param df dataframe. Typically from `import_data()`
#' @return dataframe
#' @export
times_measured <- function(df){
  combine_tech_reps(df) %>%
    dplyr::group_by(gene_id, peptide, treatment, seconds) %>%
    dplyr::summarize(times_measured = sum( is_useable(mean_tr_quant))) %>%
    dplyr::arrange(desc(times_measured))
}

#' plot the count of the number of times peptides were measured.
#'
#' Calculates and plots the number of times each peptide was measured in each
#' combination of treatment and seconds and presents a summary plot
#'
#' @param df dataframe. Typically from `import_data()`
#' @return ggplot2 plot
#' @export
times_measured_plot <- function(df){
  times_measured(df) %>%
  dplyr::group_by(treatment, seconds, times_measured) %>%
  dplyr::summarize(count = n()) %>%
  ggplot2::ggplot() +
  ggplot2::aes(times_measured, treatment) +
  ggplot2::geom_tile(ggplot2::aes( fill = count)) +
  ggplot2::facet_wrap( ~ seconds) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal()
}

#' plot the representation of peptides in each group.
#'
#' Shows what proportion of the whole set of peptides is missing in each group
#' of treatment, seconds, bio rep and tech rep.
#'
#' @param df dataframe with unmerged tech reps; typically from `import_data()`
#' @return ggplot2 plot
#' @export
missing_peptides_plot <- function(df){
  assess_missing(df) %>%
    ggplot2::ggplot() +
    ggplot2::aes(bio_rep, tech_rep) +
    ggplot2::geom_tile(  ggplot2::aes(fill = percent_missing)) +
    ggplot2::facet_grid(treatment ~ seconds) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal()

}

#' combine tech replicates into one biological replicate measurement
#' @param df dataframe; typically from `import_data()`
#' @return dataframe
combine_tech_reps <- function(df){
  df %>%
    dplyr::group_by(gene_id, peptide, treatment, seconds, bio_rep) %>%
    dplyr::summarize(mean_tr_quant = mean(quant, na.rm = TRUE) )
}


#' draw density plots for data
#'
#' Plot density of quantities in data for each treatment, seconds and biological
#' replicate
#' @param df dataframe; typically from `import_data()`
#' @param log perform log transform of data
#' @param base base to use in log transform
#' @return ggplot2 plot
#' @export
plot_quant_distributions <- function(df, log = FALSE, base = 2){

  df <- combine_tech_reps(df)
  bio_rep_count <- length(unique(df$bio_rep))
  if(log){
    p <- dplyr::mutate(df, log_mean_tr_quant = log(mean_tr_quant, base = base)) %>%
      ggplot2::ggplot() +   ggplot2::aes(log_mean_tr_quant)
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::aes(mean_tr_quant)
  }
  p +
    ggplot2::geom_density(  ggplot2::aes(fill = bio_rep), alpha = I(1/bio_rep_count)) +
    ggplot2::facet_grid(treatment ~ seconds) +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::theme_minimal()

}

#' draw qqplots for data
#'
#' Plot qqplot of distribution of quantifications in data for each treatment,
#' seconds and biological replicate
#' @param df dataframe; typically from `import_data()`
#' @param log perform log transform of data
#' @param base base to use in log transform
#' @return ggplot2 plot
#' @export
norm_qqplot <- function(df, log = FALSE, base = 2){
  df <- combine_tech_reps(df)
  if(log){
    p <- dplyr::mutate(df, log_mean_tr_quant = log(mean_tr_quant, base = base)) %>%
      ggplot2::ggplot() +   ggplot2::aes(sample = log_mean_tr_quant)
  } else {
    p <-   ggplot2::ggplot(df) +
      ggplot2::aes(sample = mean_tr_quant)
  }
    p +
      ggplot2::geom_qq(  ggplot2::aes(colour = bio_rep)) +
      ggplot2::geom_qq_line(  ggplot2::aes(colour = bio_rep)) +
      ggplot2::facet_grid(treatment ~ seconds)  +
      ggplot2::scale_color_viridis_d() +
      ggplot2::theme_minimal()
}
##TODO boxcox normalisation

##TODO sample PCA

#' convert dataframe to matrix
#'
#' @param df dataframe, typically from `import_data()`
#' @return list with members `row_info` - gene ID and peptide sequence and `data`
#' a matrix version of the data in df
#'
matrix_data <- function(df){
  df <- combine_tech_reps(df)
  row_info <- dplyr::group_by(df, gene_id, peptide, treatment, seconds, bio_rep) %>%
    dplyr::summarize(col_count = dplyr::n() )

  dm <- df %>%
    tidyr::pivot_wider(names_from = c(treatment, seconds, bio_rep), values_from = mean_tr_quant) %>%
    as.matrix()

  row_info <- dm[,c("gene_id", "peptide")]
  col_info <- colnames(dm[,3:length(colnames(dm))])
  dm <-  matrix(as.numeric(dm[,3:ncol(dm)]), nrow=nrow(dm) )
  colnames(dm) <- col_info
  return(list( row_info = row_info, data = dm ))
}

#' extract data columns for specified contrast
#'
#' @param l list of data and row_info
#' @param treatment name of treatment to use
#' @param t_seconds value of seconds to use with treatment
#' @param control name of control to use
#' @param c_seconds value of seconds to use with control
#' @return list with two members, `treatment` - matrix of treatment data; `control` -
#' matrix of control data
#'
select_columns_for_contrast <- function(l, treatment = NA,
                                        t_seconds = NA,
                                        control = NA,
                                        c_seconds = NA){
  td <- paste(treatment, t_seconds, sep = "_")
  cd <- paste(control, c_seconds, sep = "_")
  t_ind <- which(stringr::str_detect(colnames(l$data), td))
  c_ind <- which(stringr::str_detect(colnames(l$data), cd))

  return(list(
    treatment = l$data[,t_ind],
    control = l$data[,c_ind]
  ))

}

#' calculate mean fold change for peptide
#'
#' for each peptide calculate the mean quantification over all experiments
#' then get the natural scale fold change.
#'
#' @param l list with members `treatment` and `control`; each of data
#' @return vector of mean fold changes `mean(treatment) / mean(control)`
mean_fold_change <- function(l){

  rowMeans(l$treatment, na.rm = TRUE) / rowMeans(l$control, na.rm = TRUE)

}

#' calculate differential

#' @export
calc_de <- function(df, iters = 1000,
                     treatment = NA,
                     t_seconds= NA,
                     control = NA,
                     c_seconds = NA){
  all_d <- matrix_data(df)
  selected_cols <- select_columns_for_contrast(all_d, treatment = treatment,
                                   t_seconds = t_seconds,
                                   control = control,
                                   c_seconds = c_seconds)

  mean_fold_change <- mean_fold_change(selected_cols)

  ##get per peptide mean then mean of those for whole data set mean, if needed
  global_mean <- mean(apply(all_d$data, MARGIN  = 1, mean, na.rm = TRUE), na.rm = TRUE)

  ## get per peptide sd, then mean of that for whole data set sd, if needed
  global_sd <- mean(apply(all_d$data, MARGIN  = 1, sd, na.rm = TRUE), na.rm = TRUE)

  ## get mean, sd for each peptide in control treatment, sample r times from that normal distribution,
  ## if mean of test is outside confidence interval then sig.
  ## mark percentile where mean of test occurs as p
  ## method 1
    p_vals_1 <- get_percentile_iterative(selected_cols, global_mean, global_sd, iters)

  ## get p_vals replacing with lowest observed value
    lowest_vals <- apply(all_d$data, MARGIN = 1, min, na.rm = TRUE)
     r <- get_percentile_lowest_observed_value_iterative(selected_cols, lowest_vals, iters)
     p_vals_2 <- r$p
     fc_2 <- r$fc
     br <- get_bootstrap_percentile(selected_cols, lowest_vals, iters)
     p_vals_3 <- br$p
     fc_3 <- br$fc
     fdr_3 <- br$fdr
     wr <- get_wilcoxon_percentile(selected_cols, lowest_vals)
     p_vals_4 <- wr$p
     fc_4 <- wr$fc
     fdr_4 <- wr$fdr
     kw <- get_kruskal_percentile(selected_cols, lowest_vals)
     p_vals_5 <- kw$p
     fc_5 <- kw$fc
     fdr_5 <- kw$fdr
     rp <- get_rp_percentile(selected_cols, lowest_vals)
     p_vals_61 <- rp$p1
     fc_6 <- rp$fc
     fdr_61 <- rp$fdr_1
     p_vals_62 <- rp$p2
     fdr_62 <- rp$fdr_2
    # return(p_vals_2)
  control_obs <- apply(selected_cols$control, MARGIN = 1, function(x){ sum(! is.na(x))})
  treatment_obs <- apply(selected_cols$treatment, MARGIN = 1, function(x){ sum(! is.na(x))})

  col_order <- c('gene_id', 'peptide', colnames(selected_cols$control), colnames(selected_cols$treatment), 'control_observations', 'treatment_observations',
                 'fc_global_mean', 'p_val_iter_global_mean',
                 'fc_lowest_observed', 'p_val_iter_lowest_obs',
                 'p_val_bootstrap','fdr_bootstrap',
                 'p_val_wilcox','fdr_wilcox',
                 'p_val_kruskal','fdr_kruskal',
                 'p_val_rp_t_gt_c','fdr_rp_t_gt_c',
                 'p_val_rp_c_gt_t','fdr_rp_c_gt_t'
                 )
  a <- tibble::tibble(
    gene_id = all_d$row_info[,'gene_id'],
    peptide = all_d$row_info[,'peptide'],
    control_observations = control_obs,
    treatment_observations = treatment_obs,
    fc_global_mean = mean_fold_change,
    p_val_iter_global_mean = p_vals_1,
    fc_lowest_observed = fc_2,
    p_val_iter_lowest_obs = p_vals_2,
#    fc_3 = fc_3,
    p_val_bootstrap = p_vals_3,
    fdr_bootstrap = fdr_3,
#    fc_4 = fc_4,
    p_val_wilcox = p_vals_4,
    fdr_wilcox = fdr_4,
#    fc_5 = fc_5,
    p_val_kruskal = p_vals_5,
    fdr_kruskal = fdr_5,
#    fc_6 = fc_6,
    p_val_rp_t_gt_c = p_vals_61,
    fdr_rp_t_gt_c = fdr_61,
    p_val_rp_c_gt_t = p_vals_62,
    fdr_rp_c_gt_t = fdr_62
  ) %>%
    dplyr::bind_cols( as.data.frame(selected_cols$control), as.data.frame(selected_cols$treatment)) %>%
    dplyr::select(col_order)

  return(a)
}

#' first developed method for getting the p estimate
#'
#' @export
get_percentile_iterative <- function(d, global_mean, global_sd, iters = 1000){

  ## vectors for control samples
  obs_count <- apply(d$control, MARGIN = 1, function(x) {sum(! is.na(x) )} )
  no_sd <- which(obs_count <= 1)
  no_mean <- which(obs_count == 0)

  peptide_means <- apply(d$control, MARGIN = 1, mean, na.rm = TRUE)
  peptide_means[no_mean] <- global_mean

  peptide_sds <- apply(d$control, MARGIN = 1, sd, na.rm = TRUE)
  peptide_sds[no_sd] <- global_sd

  ## control samples
  nobs <- dim(d$control)[1]
  m <- matrix(NA, nrow = nobs, ncol = iters)
  for (i in 1:iters) {
    m[,i] <-  rnorm(nobs, peptide_means, peptide_sds)
  }

  #vector of treatment means
  t_means <- apply(d$treatment, MARGIN = 1, mean, na.rm = TRUE)
  t_means[is.na(t_means)] <- global_mean
  find_p <- function(x, perc){ecdf(x)(perc)}
  r <- rep(NA, nobs)
  for (i in 1:nobs ){
    r[i] <- find_p(m[i,], t_means[i])
  }
  return(r)

}

#' @export
plot_result <- function(df){
  #cols <- setdiff(colnames(df), c("gene_id", "peptide", "fc", "p_val", "control_observations", "treatment_observations"))
  #tidyr::pivot_longer(df, vars(cols), names_to = "observations") %>%

  long_results(df) %>%
  dplyr::mutate(fc_lowest_observed, observations = paste(control_observations, "-", treatment_observations)) %>%
  ggplot2::ggplot() +
    ggplot2::aes(log(fc_lowest_observed, base = 2), value) +
    ggplot2::geom_point(shape = 21, ggplot2::aes(colour = observations)) +
    ggplot2::theme_minimal() +
    #ggplot2::scale_x_continuous(trans = "log2") +
    #ggplot2::coord_trans(x = "log2") +
    ggplot2::scale_color_viridis_d() +
    ggplot2::facet_grid(observations ~ test)

}

get_percentile_lowest_observed_value_iterative <- function(d, lowest_vals, iters = 1000){
  ## vectors for control samples
  replace_vals <- function(x, lowest_vals){
    to_replace <- which(is.na(x))
    x[to_replace] <- lowest_vals[to_replace]
    return(x)
  }

  control_vals <- apply(d$control, MARGIN = 2, replace_vals, lowest_vals)
  treatment_vals <- apply(d$treatment, MARGIN = 2, replace_vals, lowest_vals)

  new_fc = rowSums(treatment_vals) / rowSums(control_vals)

  peptide_means <- apply(control_vals, MARGIN = 1, mean, na.rm = TRUE)
  peptide_sds <- apply(control_vals, MARGIN = 1, sd, na.rm = TRUE)

  ## control samples
  nobs <- dim(control_vals)[1]
  m <- matrix(NA, nrow = nobs, ncol = iters)
  for (i in 1:iters) {
    m[,i] <-  rnorm(nobs, peptide_means, peptide_sds)
  }

  #vector of treatment means
  t_means <- apply(treatment_vals, MARGIN = 1, mean, na.rm = TRUE)
  # return(list(
  #   control_vals = control_vals,
  #   treatment_vals = treatment_vals,
  #   peptide_means = peptide_means,
  #   peptide_sds = peptide_sds,
  #   m = m,
  #   t_means = t_means
  # ))
  find_p <- function(x, perc){ecdf(x)(perc)}
  r <- rep(NA, nobs)
  for (i in 1:nobs ){
    r[i] <- find_p(m[i,], t_means[i])
  }
  return(list(p = r, fc = new_fc))
}

#' @export
get_bootstrap_percentile <- function(d, lowest_vals, iters){

  ## vectors for control samples
  replace_vals <- function(x, lowest_vals){
    to_replace <- which(is.na(x))
    x[to_replace] <- lowest_vals[to_replace]
    return(x)
  }

  control_vals <- apply(d$control, MARGIN = 2, replace_vals, lowest_vals)
  treatment_vals <- apply(d$treatment, MARGIN = 2, replace_vals, lowest_vals)
  new_fc = rowSums(treatment_vals) / rowSums(control_vals)

  peptide_count <- length(lowest_vals)
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      {bstrap <- MKinfer::boot.t.test(treatment_vals[i,], control_vals[i,], R = iters)},
      warning = function(w){ NA },
      error = function(e){
        return( NA )
      },
      finally = {})

    result[i] <- bstrap$p.value
  }
  return(list(p = result, fc = new_fc, fdr = p.adjust(result, method = 'bonferroni')))
}

#' @export
get_wilcoxon_percentile <- function(d, lowest_vals){

  ## vectors for control samples
  replace_vals <- function(x, lowest_vals){
    to_replace <- which(is.na(x))
    x[to_replace] <- lowest_vals[to_replace]
    return(x)
  }

  control_vals <- apply(d$control, MARGIN = 2, replace_vals, lowest_vals)
  treatment_vals <- apply(d$treatment, MARGIN = 2, replace_vals, lowest_vals)
  new_fc = rowSums(treatment_vals) / rowSums(control_vals)

  peptide_count <- length(lowest_vals)
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      {wcox <- wilcox.test(treatment_vals[i,], control_vals[i,]) },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

    result[i] <- wcox$p.value
  }
  return(list(p = result, fc = new_fc, fdr = p.adjust(result, method = 'bonferroni')))
}

#' @export
get_kruskal_percentile <- function(d, lowest_vals){
  ## vectors for control samples
  replace_vals <- function(x, lowest_vals){
    to_replace <- which(is.na(x))
    x[to_replace] <- lowest_vals[to_replace]
    return(x)
  }

  control_vals <- apply(d$control, MARGIN = 2, replace_vals, lowest_vals)
  treatment_vals <- apply(d$treatment, MARGIN = 2, replace_vals, lowest_vals)
  new_fc = rowSums(treatment_vals) / rowSums(control_vals)

  peptide_count <- length(lowest_vals)
  result <- rep(NA, peptide_count)
  for (i in 1:peptide_count){
    tryCatch(
      { d <- data.frame(t = treatment_vals[i,], c = control_vals[i,])
        kw <- kruskal.test(t ~ c, data = d) },
      warning = function(w){ list(p.value = NA) },
      error = function(e){
        return( list(p.value = NA) )
      },
      finally = {})

    result[i] <- kw$p.value
  }
  return(list(p = result, fc = new_fc, fdr = p.adjust(result, method = 'bonferroni')))
}

#' @export
get_rp_percentile <- function(d, lowest_vals){
  replace_vals <- function(x, lowest_vals){
    to_replace <- which(is.na(x))
    x[to_replace] <- lowest_vals[to_replace]
    return(x)
  }

  control_vals <- apply(d$control, MARGIN = 2, replace_vals, lowest_vals)
  treatment_vals <- apply(d$treatment, MARGIN = 2, replace_vals, lowest_vals)
  new_fc = rowSums(treatment_vals) / rowSums(control_vals)



   nt <- dim(treatment_vals)[2]
   nc <-  dim(control_vals)[2]
   cl <- c(rep(1, nt), rep(0, nc) )
  # colnames(treatment_vals) <- paste( rep("t", nt), "_", 1:nt )
  # colnames(control_vals) <- paste(rep("c", nc), "_", 1:nc )
  # d <- dplyr::bind_cols(treatment_vals, control_vals)
  d <- cbind(treatment_vals, control_vals)
  r <- RankProd::RankProducts(d, cl, logged = FALSE, plot = FALSE, na.rm = TRUE)
  return(list(p1 = r$pval[,1], p2 = r$pval[,2], fc = new_fc, fdr_1 = r$pfp[,1], fdr_2 = r$pfp[,2]))
  return(list(p1 = new_fc, p2 = new_fc, fc = new_fc, fdr_1 = new_fc, fdr_2 = new_fc))
}

#' @export
long_results <- function(r){
  r %>%
    dplyr::select(gene_id, peptide, control_observations, treatment_observations, fc_global_mean, fc_lowest_observed, dplyr::starts_with("p_val" ) ) %>%
    tidyr::pivot_longer(dplyr::starts_with("p_val"), names_to = 'test')
}

#' @export
p_value_hist <- function(r){
  r %>% long_results() %>%
    ggplot2::ggplot() +
    ggplot2::aes(value) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::facet_grid( ~ test, scales = 'free') +
    ggplot2::theme_minimal()
}

#' @export
p_value_corrs <- function(r){
  ##not useful - some tests are two side, other tests are one side
  r %>% dplyr::select(dplyr::starts_with("p_val")) %>%
    corrr::correlate() %>%
    corrr::rplot()
}

#' @export
compare_calls <- function(r){

 sets <- list(
   kruskal = which(r$p_val_kruskal <= 0.05),
   wilcox = which(r$p_val_wilcox <= 0.05),
   bootstrap = which(r$p_val_bootstrap <= 0.05),
   rp = union(which(r$p_val_rp_c_gt_t <= 0.05), which(r$p_val_rp_t_gt_c <= 0.05)),
   norm_iter_global_mean = which(r$p_val_iter_global_mean >= 0.975 | r$p_val_iter_global_mean <= 0.025),
   norm_iter_lowest_obs = which(r$p_val_iter_lowest_obs >= 0.975 | r$p_val_iter_lowest_obs <= 0.025)
 )
  UpSetR::upset(UpSetR::fromList(sets), order.by = "freq")
}

#' @export
plot_fc <- function(df, log = FALSE){
  p <- ggplot2::ggplot(df)
  if (log) {
    p <- p + ggplot2::aes(log(fc_lowest_observed, base = 2))
  } else {
    p <- p + ggplot2::aes( fc_lowest_observed )
  }
    p +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::theme_minimal()
}
#' @export
fc_qqplot <- function(df, log = FALSE){

  p <- ggplot2::ggplot(df)
  if (log){
    p <- p + ggplot2::aes(sample = log(fc_lowest_observed, base =2))
  } else {
    p <- p + ggplot2::aes(sample = fc_lowest_observed)
  }
  p +  ggplot2::geom_qq() +
    ggplot2::geom_qq_line(  ) +
    ggplot2::theme_minimal()

}

