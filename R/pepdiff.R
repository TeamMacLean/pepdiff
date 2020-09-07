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
#' @param peptide Column containing the sequence of this peptide
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
    dplyr::mutate(id = paste0(.data$gene_id, ":", .data$peptide)) %>%
    dplyr::transmute(id = as.character(.data$id),
                     gene_id = as.character(.data$gene_id),
                     peptide = as.character(.data$peptide),
                     treatment = as.character(.data$treatment),
                     seconds = as.numeric(.data$seconds),
                     bio_rep = as.character(.data$bio_rep),
                     tech_rep = as.character(.data$tech_rep),
                     quant = as.numeric(.data$quant)
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
#' @importFrom rlang .data
assess_missing <- function(df){

  dplyr::group_by( df, .data$treatment, .data$seconds, .data$bio_rep, .data$tech_rep) %>%
    dplyr::summarize(peptides = dplyr::n(), percent_missing = (sum(is.na(.data$quant)) / dplyr::n() ) * 100 )
}

#' calculate number of measurements of each peptide in each treatment and time
#'
#' For each peptide, works out how many biologically replicated measurements are
#' available in the different combinations of treatment and seconds
#'
#' @param df dataframe. Typically from `import_data()`
#' @return dataframe
#' @export
#' @importFrom rlang .data
times_measured <- function(df){
  combine_tech_reps(df) %>%
    dplyr::group_by(.data$gene_id, .data$peptide, .data$treatment, .data$seconds) %>%
    dplyr::summarize(times_measured = sum( is_useable(.data$mean_tr_quant))) %>%
    dplyr::arrange(dplyr::desc(.data$times_measured))
}


#' compare the peptide quantities in two experiments
#'
#' For each peptide this function carries out the selected tests to determine peptides
#' that are differentially abundant in the two experiments specified. Before tests are
#' performed the technical replicates are merged and lowest observed value replacement
#' for each missing biological replicate is done. All comparisins are performed as `treatment / control`
#'
#' @param df dataframe. Typically from `import_data()`
#' @param iters number of iterations to perform for iterative tests
#' @param treatment name of the experimental treatment to use as the `treatment` condition
#' @param control name of the experimental treatment to use as the `control` condition
#' @param t_seconds time point of the `treatment` condition to use
#' @param c_seconds time point of the `control` condition to use
#' @param metrics character vector of tests to use, one or more of: `norm_quantile`, `bootstrap_t`, `wilcoxon`, `kruskal-wallis`, `rank_product`
#' @return dataframe with original and replaced quantification values, natural fold change, biological replicates and p-value / fdr for each
#' @export
compare <- function(df,
                    iters = 1000,
                    treatment = NA,
                    t_seconds = NA,
                    control = NA,
                    c_seconds = NA,
                    metrics = NA){
  d <- matrix_data(df)
  selected_cols <- select_columns_for_contrast(
    d,
    treatment = treatment,
    t_seconds = t_seconds,
    control = control,
    c_seconds = c_seconds
  )

  lowest_vals <- min_peptide_values(d)

  result <- list()

  control_reps <- apply(selected_cols$control, MARGIN = 1, function(x){ sum(! is.na(x))})
  treatment_reps <- apply(selected_cols$treatment, MARGIN = 1, function(x){ sum(! is.na(x))})
  treatment <- apply(selected_cols$treatment, MARGIN = 2, replace_vals, lowest_vals)
  control <- apply(selected_cols$control, MARGIN = 2, replace_vals, lowest_vals)
  fc <- mean_fold_change(treatment, control)

  treatment_mean_count <- rowMeans(treatment, na.rm = TRUE)
  control_mean_count <- rowMeans(control, na.rm = TRUE)

  result$info <- as.data.frame(d$row_info)
  result$treatment <- as.data.frame(treatment)
  result$control <- as.data.frame(control)
  result$fold_change <- data.frame(fold_change = fc)
  result$treatment_mean_count <- data.frame(treatment_mean_count = treatment_mean_count)
  result$control_mean_count <- data.frame(control_mean_count = control_mean_count)

  result$unreplaced_treatment <- as.data.frame(selected_cols$treatment)
  names(result$unreplaced_treatment) <- paste0("unreplaced_", names(result$unreplaced_treatment))
  result$unreplaced_control <- as.data.frame(selected_cols$control)
  names(result$unreplaced_control) <- paste0("unreplaced_", names(result$unreplaced_control))

  result$treatment_replicates <- data.frame(treatment_replicates = treatment_reps)
  result$control_replicates <- data.frame(control_replicates = control_reps)

  if ("norm_quantile" %in% metrics){
    result$norm_quantile <- get_percentile_lowest_observed_value_iterative(treatment, control, iters)
  }
  if ("bootstrap_t" %in% metrics){
    result$bootstrap_t <- get_bootstrap_percentile(treatment, control, iters)
  }
  if ("wilcoxon" %in% metrics){
    result$wilcoxon <- get_wilcoxon_percentile(treatment, control)
  }
  if ("kruskal-wallis" %in% metrics){
    result$kw <- get_kruskal_percentile(treatment, control)
  }
  if ("rank_product" %in% metrics){
    result$rp <- get_rp_percentile(treatment, control)
  }
  return(tibble::as_tibble(dplyr::bind_cols(result)))
}




#' compare many combinations of treatment and control
#'
#' for each combination of treatment and control condition, runs the `compare()` function
#' and collates the results
#'
#' @param df dataframe. Typically from `import_data()`
#' @param comparison path to file or dataframe of comparisons with columns treatment, t_seconds, control, c_seconds
#' @param iters number of iterations to perform for iterative tests
#' @param metrics character vector of tests to use, one or more of: `norm_quantile`, `bootstrap_t`, `wilcoxon`, `kruskal-wallis`, `rank_product`
#' @export
compare_many <- function(df, comparison, iters = 1000, metrics = NA) {
  if (!is.data.frame(comparison)) {
    comparison <- readr::read_csv(comparison )
  }

  comparison <- dplyr::transmute_all(comparison, .funs = as.character )

  get_names <- function(r){
    a <- paste(r['treatment'], r['t_seconds'], sep = "_")
    b <- paste(r['control'], r['c_seconds'], sep = "_")
    paste(a,b, sep = "-")
  }

  do_compare <- function(r, df, iters, metrics) {
    treatment = r['treatment']
    t_seconds = r['t_seconds']
    control = r['control']
    c_seconds = r['c_seconds']
    compare(df, iters = iters, metrics = metrics,
            treatment = treatment,
            t_seconds = t_seconds,
            control = control,
            c_seconds = c_seconds)

  }

  result <- apply(comparison, MARGIN = 1, do_compare, df, iters, metrics)
  names(result) <- apply(comparison, MARGIN = 1, get_names)
  return(result)
}










