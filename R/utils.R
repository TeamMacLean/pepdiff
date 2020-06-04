#' calculate the number of non- NA or NaNs in a vector
#'
#' @param x vector of values
#' @return vector of TRUE/FALSE indicating positions of useable stuff
#'
is_useable <- function(x){ !is.na(x) & !is.nan(x)}

#' combine tech replicates into one biological replicate measurement
#'
#' @param df dataframe; typically from `import_data()`
#' @return dataframe
#' @importFrom rlang .data
combine_tech_reps <- function(df){
  df %>%
    dplyr::group_by(.data$gene_id, .data$peptide, .data$treatment, .data$seconds, .data$bio_rep) %>%
    dplyr::summarize(mean_tr_quant = mean(.data$quant, na.rm = TRUE) )
}

#' convert dataframe to matrix
#'
#' @param df dataframe, typically from `import_data()`
#' @return list with members `row_info` - gene ID and peptide sequence and `data`
#' a matrix version of the data in df
#' @importFrom rlang .data
matrix_data <- function(df){
  df <- combine_tech_reps(df)
  row_info <- dplyr::group_by(df, .data$gene_id, .data$peptide, .data$treatment, .data$seconds, .data$bio_rep) %>%
    dplyr::summarize(col_count = dplyr::n() )

  dm <- df %>%
    tidyr::pivot_wider(names_from = c(.data$treatment, .data$seconds, .data$bio_rep), values_from = .data$mean_tr_quant) %>%
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
#' @param control vector of `control` data
#' @param treatment vector of `treatment` data
#' @return vector of mean fold changes `mean(treatment) / mean(control)`
mean_fold_change <- function(treatment, control){

  rowMeans(treatment, na.rm = TRUE) / rowMeans(control, na.rm = TRUE)

}


#' find minimal values for a peptide across replicates
#'
#' @param d matrix of peptide quantifications. rows are peptides, columns are
#' experimental conditions
#' @return vector of minimume values of length nrow(d)
#'
min_peptide_values <- function(d){
  apply(d$data, MARGIN = 1, min, na.rm = TRUE)
}

#' replace missing values in a vector with others provided
#'
#' @param x vector of values with missing values
#' @param lowest_vals vector of replacements
#' @return vector
replace_vals <- function(x, lowest_vals){
  to_replace <- which(is.na(x))
  x[to_replace] <- lowest_vals[to_replace]
  return(x)
}

#' convert wide format results table containing p-value estimates to long format
#'
#' tidies up the wide results table from `compare()` to a long format, dropping the fdr columns and quantitie columns
#'
#' @param r results dataframe typically from `compare()`
#' @return dataframe in long format missing some columns from r
#' @export
#' @importFrom rlang .data
long_results <- function(r){
  r %>%
    dplyr::select(.data$gene_id, .data$peptide, .data$treatment_replicates, .data$control_replicates, .data$fold_change, dplyr::ends_with("p_val" ) ) %>%
    tidyr::pivot_longer(dplyr::ends_with("p_val"), names_to = 'test')
}

