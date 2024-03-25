

#' Import and preprocess peptide quantification data from a file or data frame.
#'
#' This function reads a data file (e.g., CSV) or accepts a data frame as input and preprocesses the data by renaming columns, creating a unique identifier, and ensuring proper data types. It is designed to handle peptide quantification data and standardize the column names for treatment, biological replicate, technical replicate, peptide, gene ID, and more.
#'
#' @param file A path to a data file (e.g., CSV) or a data frame containing peptide quantification data.
#' @param treatment The column name representing treatment information in the data. Default is "genotype".
#' @param bio_rep The column name representing biological replicate information. Default is "bio_rep".
#' @param tech_rep The column name representing technical replicate information. Default is "tech_rep".
#' @param quant The column name representing quantification values. Default is "total_area".
#' @param seconds The column name representing time points (in seconds). Default is "seconds".
#' @param gene_id The column name representing gene IDs. Default is "gene_id".
#' @param peptide The column name representing peptide sequences. Default is "peptide_sequence".
#' @param quant The column name representing the peptide quantification. Default is "total_area"
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example using a CSV file:
#' data_file <- "peptide_quantification_data.csv"
#' imported_data <- import_data(data_file)
#'
#' # Example using a data frame:
#' data_frame <- data.frame(genotype = c("A", "B", "A"),
#'                          bio_rep = c(1, 2, 1),
#'                          tech_rep = c(1, 1, 2),
#'                          total_area = c(100, 150, 90),
#'                          seconds = c(0, 30, 60),
#'                          gene_id = c("Gene1", "Gene2", "Gene1"),
#'                          quant = runif(3),
#'                          peptide_sequence = c("AAA", "BBB", "CCC"))
#  imported_data <- import_data(data_frame)
#' }
#'
#' @return A preprocessed data frame with standardized column names and data types.
import_data <- function(file,
                        treatment = "genotype",
                        bio_rep = "bio_rep",
                        tech_rep = "tech_rep",
                        quant = "total_area",
                        seconds = "seconds",
                        gene_id = "gene_id",
                        peptide = "peptide_sequence"
                        ) {
    if (is.data.frame(file)){
      csv <- file
    } else {
      file_info <- file.info(file)
      if ( !is.na(file_info$size) ){
        csv <- readr::read_csv(file)
      }
    }


    csv %>%
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
count_peptides_measured <- function(df){
  combine_tech_reps(df) %>%
    dplyr::group_by(.data$gene_id, .data$peptide, .data$treatment, .data$seconds) %>%
    dplyr::summarize(peptide_measurements = sum( is_useable(.data$mean_tr_quant))) %>%
    dplyr::arrange(dplyr::desc(.data$peptide_measurements))
}


#' Compare peptide quantities between experimental conditions using selected statistical tests.
#'
#' This function compares peptide quantities between experimental conditions and applies selected statistical tests to identify differentially abundant peptides. The function calculates natural fold changes, performs various tests, and provides information about replicates and statistical power.
#'
#' @param df A dataframe containing peptide quantification data, typically obtained from `import_data()`.
#' @param iters The number of iterations to perform for iterative tests. Default is 1000.
#' @param treatment The name of the experimental treatment condition to compare.
#' @param t_seconds The time point of the treatment condition to compare.
#' @param control The name of the experimental control condition to compare.
#' @param c_seconds The time point of the control condition to compare.
#' @param tests A character vector of tests to use, one or more of: 'norm_quantile', 'bootstrap_t', 'wilcoxon', 'kruskal-wallis', 'rank_product', 'gamma', 'eb'
#' @param log whether or not to log the data
#' @param base base for logs
#' @return A data frame with information about original and replaced quantification values, natural fold changes, replicates, p-values, false discovery rate (FDR), and statistical power for each peptide.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example using imported data and performing a bootstrap t-test:
#' data_file <- "peptide_quantification_data.csv"
#' imported_data <- import_data(data_file)
#' comparison_result <- compare(imported_data, treatment = "A",
#' control = "B", tests = "bootstrap_t")
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname compare
compare <- function(df,
                    iters = 1000,
                    treatment = NA,
                    t_seconds = NA,
                    control = NA,
                    c_seconds = NA,
                    tests = c("bootstrap_t"),
                    log=FALSE,
                    base=2){
  d <- matrix_data(df, log=log, base=base)

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
  fc <- mean_fold_change(treatment, control, log)

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

  if (all(colnames(result$treatment) %in% colnames(result$control))  ) {
    warning("Duplicated names for test and control found. Presuming self vs self comparison - do you need this\nAppending _t and _c to treatment and control names.")
    colnames(result$treatment) <- paste0(colnames(result$treatment), "_t")
    colnames(result$control) <- paste0(colnames(result$control), "_c")
    colnames(result$unreplaced_treatment) <- paste0(colnames(result$unreplaced_treatment), "_t")
    colnames(result$unreplaced_treatment) <- paste0(colnames(result$unreplaced_treatment), "_c")
  }



  if ("norm_quantile" %in% tests){
    result$norm_quantile <- get_percentile_lowest_observed_value_iterative(treatment, control, iters)
  }
  if ("bootstrap_t" %in% tests){
    result$bootstrap_t <- get_bootstrap_percentile(treatment, control, iters)
  }
  if ("wilcoxon" %in% tests){
    result$wilcoxon <- get_wilcoxon_percentile(treatment, control)
  }
  if ("kruskal-wallis" %in% tests){
    result$kw <- get_kruskal_percentile(treatment, control)
  }
  if ("rank_product" %in% tests){
    result$rp <- get_rp_percentile(treatment, control)
  }

  if ("gamma" %in% tests){
    result$gamma <- get_gamma(treatment, control)
  }

  if ("eb" %in% tests){
    result$eb <- get_eb(treatment, control)
  }


  powers <- get_power(treatment, control)
  result$power <- powers$power
  result$d <- powers$d
  result$min_reps <- powers$min_reps

  return(tibble::as_tibble(dplyr::bind_cols(result, .name_repair="minimal")))
  #return(result)
}


#' Compare peptide quantities for multiple experimental conditions using selected statistical tests.
#'
#' This function compares peptide quantities between multiple experimental conditions specified in a comparison dataframe. It applies selected statistical tests to identify differentially abundant peptides for each pairwise condition comparison.
#'
#' @param df A dataframe containing peptide quantification data, typically obtained from `import_data()`.
#' @param comparison A comparison dataframe specifying the experimental conditions to compare. Each row defines a pairwise condition comparison, including treatment, control, and time points. If not already a data frame, you can provide a path to a CSV file, and it will be read into a data frame.
#' @param iters The number of iterations to perform for iterative tests. Default is 1000.
#' @param tests #' @param tests character vector of tests to use, one or more of: `norm_quantile`, `bootstrap_t`, `wilcoxon`, `kruskal-wallis`, `rank_product`, `gamma`, `eb`
#' @param log should the data be logged prior to analysis
#' @param base base to use for logging (2)
#'
#' @return A list of dataframes, where each dataframe represents the comparison results for a pairwise condition comparison. The names of the list elements are constructed based on the conditions being compared.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example using imported data and a comparison dataframe:
#' data_file <- "peptide_quantification_data.csv"
#' comparison_data <- "comparison_conditions.csv"
#' imported_data <- import_data(data_file)
#' comparison_result <- compare_many(imported_data, comparison_data,
#' tests = c("bootstrap_t", "wilcoxon"))
#' }
#' @keywords data
#' @family statistical analysis
#' @rdname compare_many
compare_many <- function(df, comparison, iters = 1000, tests = c("bootstrap_t"), log=FALSE, base=2) {
  if (!is.data.frame(comparison)) {
    comparison <- readr::read_csv(comparison )
  }

  comparison <- dplyr::transmute_all(comparison, .funs = as.character )

  get_names <- function(r){
    a <- paste(r['treatment'], r['t_seconds'], sep = "_")
    b <- paste(r['control'], r['c_seconds'], sep = "_")
    paste(a,b, sep = "-")
  }

  do_compare <- function(r, df, iters, tests) {
    treatment = r['treatment']
    t_seconds = r['t_seconds']
    control = r['control']
    c_seconds = r['c_seconds']
    compare(df, iters = iters, tests = tests,
            treatment = treatment,
            t_seconds = t_seconds,
            control = control,
            c_seconds = c_seconds,
            log=log, base=base)

  }

  result <- apply(comparison, MARGIN = 1, do_compare, df, iters, tests)
  names(result) <- apply(comparison, MARGIN = 1, get_names)
  return(result)
}


#' summarize the health of the experiment as a whole
#'
#' returns statistics on the probability of detecting effects on at least half
#' the outcomes and completeness of sampling
#' over threshold b
#'
#' @param r the results object, from `compare()` or `compare_many()`
#' @param b the statistical power at which to evaluate
#' @return i
#' @export
health <- function(r, b=0.8) {

 r <- dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::select(.data$comparison, .data$power, .data$d, .data$min_reps,.data$fold_change,dplyr::starts_with("unreplaced"))

 cc <- dplyr::select(r, dplyr::starts_with("unreplaced")) %>%
   complete.cases() %>%
   sum()

 po <- 1/2 - min(r$power, na.rm=TRUE)
 co <- cc / nrow(r)

 c(
   power_health = po,
   completeness_health = co
 )

}

#' Count Peptides Needing a give Number of Replicates
#'
#' This function takes results object and counts the number of peptides
#' that need fewer than a specified minimum number of replicates for the given power.
#'
#' @param r results usually from `compare_many()`  containing information about peptide comparisons.
#' @param m The minimum number of replicates required for a peptide at the given power.
#' @param b The power threshold for categorizing comparisons ("Sufficient Power" or "Under Power").
#'
#' @return A data frame with counts of peptides for each comparison and power category.
#'
#' @export
#'
count_peptides_needing_m_reps <- function(r, m=10, b=0.8) {

  dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::filter(min_reps < m) %>%
    dplyr::mutate(power = dplyr::if_else(power >=b, "Sufficient Power", "Under Power")) %>%
    dplyr::group_by(comparison, power) %>%
    dplyr::tally() %>%
    knitr::kable()
}




