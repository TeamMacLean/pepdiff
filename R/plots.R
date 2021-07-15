#' plot the count of the number of times peptides were measured.
#'
#' Calculates and plots the number of times each peptide was measured in each
#' combination of treatment and seconds and presents a summary plot
#'
#' @param df dataframe. Typically from `import_data()`
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
times_measured_plot <- function(df){
  times_measured(df) %>%
    dplyr::group_by(.data$treatment, .data$seconds, .data$times_measured) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    ggplot2::ggplot() +
    ggplot2::aes(.data$times_measured, .data$treatment) +
    ggplot2::geom_tile(ggplot2::aes( fill = .data$count)) +
    ggplot2::geom_text(ggplot2::aes(label= .data$count)) +
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
#' @importFrom rlang .data
missing_peptides_plot <- function(df){
  assess_missing(df) %>%
    ggplot2::ggplot() +
    ggplot2::aes(.data$bio_rep, .data$tech_rep) +
    ggplot2::geom_tile(  ggplot2::aes(fill = .data$percent_missing)) +
    ggplot2::facet_grid(treatment ~ seconds) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal()

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
#' @importFrom rlang .data
plot_quant_distributions <- function(df, log = FALSE, base = 2){

  df <- combine_tech_reps(df)
  bio_rep_count <- length(unique(df$bio_rep))
  p <- NULL
  if(log){
    p <- dplyr::mutate(df, log_mean_tr_quant = log(.data$mean_tr_quant, base = base)) %>%
      ggplot2::ggplot() +
      ggplot2::aes(x = .data$log_mean_tr_quant) +
      ggplot2::geom_density(  ggplot2::aes(fill = .data$bio_rep), alpha = I(1/bio_rep_count))
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::aes(x = .data$mean_tr_quant) +
      ggplot2::geom_density(ggplot2::aes(fill = .data$bio_rep),alpha = I(1/bio_rep_count))
  }
  p +

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
#' @importFrom rlang .data
norm_qqplot <- function(df, log = FALSE, base = 2){
  #return(NULL)
  df <- combine_tech_reps(df)
  if(log){
    p <- dplyr::mutate(df, log_mean_tr_quant = log(.data$mean_tr_quant, base = base)) %>%
      ggplot2::ggplot() +   ggplot2::aes(sample = .data$log_mean_tr_quant)
  } else {
    p <-   ggplot2::ggplot(df) +
      ggplot2::aes(sample = .data$mean_tr_quant)
  }
  p +
    ggplot2::geom_qq(  ggplot2::aes(colour = .data$bio_rep)) +
    ggplot2::geom_qq_line(  ggplot2::aes(colour = .data$bio_rep)) +
    ggplot2::facet_grid(treatment ~ seconds)  +
    ggplot2::scale_color_viridis_d() +
    ggplot2::theme_minimal()
}

#' plot the p-values against fold change for the tests used in `compare()`
#'
#' plots fold change against p-value in all tests used, also splits data on the number
#' of biological replicates were available before value replacement for each peptide
#'
#' @param df result dataframe typically from `compare()`
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_result <- function(df){

  long_results(df) %>%
    dplyr::mutate(replicates = paste(.data$treatment_replicates, "-", .data$control_replicates)) %>%
    ggplot2::ggplot() +
    ggplot2::aes(log(.data$fold_change, base = 2), .data$value) +
    ggplot2::geom_point(shape = 21, ggplot2::aes(colour = .data$replicates)) +
    ggplot2::theme_minimal() +
    #ggplot2::scale_x_continuous(trans = "log2") +
    #ggplot2::coord_trans(x = "log2") +
    ggplot2::scale_color_viridis_d() +
    ggplot2::facet_grid(replicates ~ test)

}

#' plot histograms of p-values for each test used
#'
#' @param r list of result dataframes typically from `compare()`
#' @return ggplot2 plot
#'
#' @export
#' @importFrom rlang .data
p_value_hist <- function(r){
  r %>% long_results() %>%
    ggplot2::ggplot() +
    ggplot2::aes(.data$value) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::facet_grid( ~ test, scales = 'free') +
    ggplot2::theme_minimal()
}

#' plot histogram of fold change distribution for a comparison
#'
#' @param df result dataframe, typically from `compare()`
#' @param log log the fold change values
#' @param base base of the log, if used
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_fc <- function(df, log = FALSE, base = 2){
  p <- ggplot2::ggplot(df)
  if (log) {
    p <- p + ggplot2::aes(log(.data$fold_change, base = base))
  } else {
    p <- p + ggplot2::aes( .data$fold_change )
  }
  p +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::theme_minimal()
}

#' plot qqplot of fold changes from a comparison
#' @param df result dataframe, typically from `compare()`
#' @param log log the fold change values
#' @param base base of the log, if used
#' @export
#' @importFrom rlang .data
fc_qqplot <- function(df, log = FALSE, base = 2 ){

  p <- ggplot2::ggplot(df)
  if (log){
    p <- p + ggplot2::aes(sample = log(.data$fold_change, base = base))
  } else {
    p <- p + ggplot2::aes(sample = .data$fold_change)
  }
  p +  ggplot2::geom_qq() +
    ggplot2::geom_qq_line(  ) +
    ggplot2::theme_minimal()

}


#' compare sets of significant peptides called by the used data
#'
#' produces an UpSet plot showing intersections and set-size of the different
#' sets of significant peptides called by the methods used in the provided result
#' dataframe
#'
#' @param r result dataframe typically from `compare()`
#' @param sig significance cut-off to select peptides
#' @return UpSet plot
#'
#' @export
compare_calls <- function(r, sig = 0.05){
  upper_sig <- 1 - round(sig / 2, 4)
  lower_sig <- round(sig / 2, 4)
  sets <- list()


  if ("norm_quantile_pval" %in% names(r)){
    sets$norm_iter = which(r$norm_quantile_p_val>= upper_sig | r$norm_quantile_p_val <= lower_sig)
  }
  if ("bootstrap_t_p_val" %in% names(r)){
    sets$bootstrap_t = which(r$bootstrap_t_p_val <= sig)
  }
  if ("wilcoxon_p_val" %in% names(r)){
    sets$wilcoxon = which(r$wilcoxon_p_val <= sig)
  }
  if ("kruskal_p_val" %in% names(r)) {
    sets$kruskal = which(r$kruskal_p_val <= sig)
  }
  if ("rank_prod_p1_p_val" %in% names(r)){
    sets$rank_prod = union(which(r$rank_prod_p1_p_val <= sig), which(r$rank_prod_p2_p_val <= sig))
  }

  #return(sets)
  UpSetR::upset(UpSetR::fromList(sets), order.by = "freq")
}

#' @importFrom rlang .data
drop_columns <- function(df, sig, metric, log, base, rows_to_keep = NULL){

  small <- df %>%
    dplyr::select(.data$gene_id, .data$peptide, .data$fold_change, dplyr::ends_with('p_val'))

  if('rank_prod_p1_p_val' %in% names(small) | 'rank_prod_p2_p_val' %in% names(small)){
    small <- small %>%
      dplyr::rowwise() %>%
      dplyr::mutate(rank_product_p_val = min(c(.data$rank_prod_p1_p_val, .data$rank_prod_p2_p_val))) %>%
      dplyr::select(-.data$rank_prod_p1_p_val, -.data$rank_prod_p2_p_val)
  }

  if (is.null(rows_to_keep)){
    small <- small %>% dplyr::select(.data$gene_id, .data$peptide, .data$fold_change, dplyr::starts_with(metric)) %>%
      dplyr::filter_at(dplyr::vars(dplyr::ends_with('p_val')), dplyr::any_vars((.data$. <= sig))) %>%
      dplyr::select(.data$gene_id, .data$peptide, .data$fold_change) %>%
      dplyr::ungroup()
      return(small)
  }
  else {
    small <- small %>% dplyr::select(.data$gene_id, .data$peptide, .data$fold_change, dplyr::starts_with(metric)) %>%
      dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " ")) %>%
      dplyr::filter(.data$gene_peptide %in% rows_to_keep) %>%
      dplyr::select(.data$gene_id, .data$peptide, .data$fold_change) %>%
      dplyr::ungroup()
    return(small)
  }
}

#' makes heatmap from all experiments, pass a single metric and sig value
#' reduces dataframes and makeslong list, makes a heatmap
#' @param l list of results data frames, typically from `compare_many()`
#' @param sig significance level at which to call peptides
#' @param metric one metric to use for the source of p-values
#' @param log log the fold change values to show in the heatmap
#' @param base base to use for log transformation
#' @param col_order specify a column order for the plot, default is names(l)
#' @param row_kms k for kmeans splitting of heatmap along rows
#' @param col_kms k for kmeans splitting of heatmap on columns
#' @param pal cbrewer palette to use "RdBu"
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_heatmap <- function(l, sig = 0.05, metric = NA, log = FALSE, base = 2, col_order = NULL, rotate_x_labels = TRUE, all_points = FALSE, only_sig_points = TRUE, pal="RdBu", row_kms = NULL, col_kms = NULL, scale_min = -2, scale_max = 2) {

  if (is.null(col_order)) {
    col_order <- names(l)
  }

  filtered <- NA
  if (all_points){
    rows_to_keep <- dplyr::bind_rows(l) %>%
      dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " ))
    rows_to_keep <- unique(rows_to_keep$gene_peptide)
    filtered <- lapply(l, drop_columns, sig, metric, log, base, rows_to_keep)
  }
  else {
    filtered <- lapply(l, drop_columns, sig, metric, log, base)
  }

  if (! only_sig_points){
    rows_to_keep <- dplyr::bind_rows(filtered) %>%
      dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " ))
    rows_to_keep <- unique(rows_to_keep$gene_peptide)
    filtered <- lapply(l, drop_columns, sig, metric, log, base, rows_to_keep)

  }

  x <- dplyr::bind_rows(filtered, .id = "comparison")

  row_title <- "Comparison"
  col_title <- "Gene Peptide"
  name <- "Fold Change"

    if (log){
    x$fold_change <- log(x$fold_change, base = base)
    name <- paste("Log", base, "Fold Change")
  }

  p <- x %>%
       dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " )) %>%
    tidybulk::impute_missing_abundance(~1, .sample=comparison, .transcript=gene_peptide, .abundance=fold_change) %>%
    tidyHeatmap::heatmap(gene_peptide, comparison, fold_change,
                         column_order = col_order,
                         row_km = row_kms,
                         column_km = col_kms,
                         column_title = col_title,
                         row_title = row_title,
                         name = name,
                       palette_value =  circlize::colorRamp2(
                           seq(scale_min, scale_max, length.out = 11),
                           rev(RColorBrewer::brewer.pal(11, pal))
                         )
                         )

  return(p)
}


#' plots a pca on the treatment, seconds, bio-rep
#'
#' Performs and draws a PCA plot with four panels, PCA with sample names
#' coloured by treatment, seconds and biorep and a scree plot of the PCA dimensions
#'
#' @param df dataframe, typically from `import_data()`
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_pca <- function(df) {
  #lowest_vals <- min_peptide_values(d)
   d <- matrix_data(df)
  lowest_vals <- min_peptide_values(d)
  d$data <- apply(d$data, MARGIN = 2, replace_vals, lowest_vals)
  pca <- stats::prcomp(d$data, scale = TRUE, center = TRUE)
  sample_names <- colnames(d$data)
  pca_df <- as.data.frame(pca$rotation)
  pca_df$sample_key <- sample_names
  pca_df <- tidyr::separate(pca_df, .data$sample_key, into = c("treatment", "seconds", "bio_rep"), sep = "_", remove= FALSE ) %>%
    dplyr::mutate(sample = paste0(.data$treatment, "_", .data$seconds))

  a <-  ggplot2::ggplot(pca_df) +
    ggplot2::aes(x = .data$PC1,y = .data$PC2,label=.data$sample, color = .data$bio_rep ) +
    ggplot2::geom_text(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_viridis_d()


  b <- ggplot2::ggplot(pca_df) +
    ggplot2::aes(x = .data$PC1,y = .data$PC2,label=.data$sample, color = .data$treatment ) +
    ggplot2::geom_text(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_viridis_d()

  c <-  ggplot2::ggplot(pca_df) +
    ggplot2::aes(x = .data$PC1,y = .data$PC2,label=.data$sample, color = .data$seconds ) +
    ggplot2::geom_text(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_viridis_d()

  d <- factoextra::fviz_screeplot(pca)
  cowplot::plot_grid(a,b,c,d, nrow = 2, ncol=2)
}

#' K-means cluster the data on the samples
#'
#'Performs and draws a K-means cluster on the samples. Estimates number of clusters
#'as the product of the number of treatments and seconds. So tries to group the bio reps together
#'
#' @param df dataframe, typically from `import_data()`
#' @param nstart nstart points for `kmeans()` function
#' @param iter.max max iterations to perform for `kmeans()` function
#' @return ggplot2 plot
#' @export
plot_kmeans <- function(df, nstart = 25, iter.max = 1000){
  n <- length(unique(df$treatment)) * length(unique(df$seconds))
  if (n < 1) {
    n <- 1
  }
  d <- matrix_data(df)
  lowest_vals <- min_peptide_values(d)
  d$data <- apply(d$data, MARGIN = 2, replace_vals, lowest_vals)
  pca <- stats::prcomp(d$data, scale = TRUE, center = TRUE)
  main_components <- pca$rotation[,1:3]
  km <- stats::kmeans(main_components, n, nstart=nstart, iter.max = iter.max)
  factoextra::fviz_cluster(km, data = main_components, palette = viridis::viridis(n),
               ggtheme = ggplot2::theme_minimal(),
               main = paste0("k-Means Sample: ", n, " Cluster(s)"))
}


#' volcano plot the data
#'
#' draws a plot of peptide count against log fc at either protein or peptide level for samples
#' @param l list of results data frames, typically from `compare_many()`
#' @param metric single metric to use for volcano plot
#' @return ggplot2 plot
#' @export
#'
volcano_plot <- function(l, log = FALSE, base = 2, by="peptide", sig = 0.05, metric = NA, option="E", direction=-1  ) {
  xlabtxt <- paste0("Log ",base," Fold Change")
  ylabtxt <- paste0("-Log ",base, " P")
  dplyr::bind_rows(l, .id = "comparison")  %>%
    dplyr::select(comparison, gene_id, peptide, treatment_mean_count, control_mean_count, fold_change, dplyr::starts_with(metric)) %>%
    dplyr::rename(p_val = dplyr::ends_with('p_val') ) %>%
    dplyr::mutate(
                  gene_peptide = paste(gene_id, peptide, sep = " "),
                  log_fc = log(fold_change, base = base),
                  log_p = log(p_val, base=base) * -1,
                  signal = log((treatment_mean_count + control_mean_count), base) * -1,
                  change = dplyr::if_else(log_fc > 0 & p_val <= sig, "Up",
                            dplyr::if_else(log_fc < 0 & p_val <= sig, "Down", "None")),
                  change = forcats::fct_relevel(change, rev)
                  ) %>%
    dplyr::select(comparison, gene_peptide, log_fc, log_p, change) %>%
    ggplot2::ggplot() +
    ggplot2::aes(log_fc, log_p) +
    ggplot2::geom_point(ggplot2::aes(colour = change)) +
    ggplot2::facet_wrap( ~ comparison) +
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_viridis_d(option=option, direction=direction) +
    ggplot2::xlab(xlabtxt) +
    ggplot2::ylab(ylabtxt) +
    ggplot2::labs(colour = "Change")




}

#' converts a results object to a matrix as if for direct use in external heatmap functions
#'
list2mat <- function(r) {
  x <- dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " )) %>%
    dplyr::select(-.data$gene_id, -.data$peptide) %>%
    tidyr::pivot_wider(.data$gene_peptide, names_from = .data$comparison, values_from = .data$fold_change )

  rownames <- x$gene_peptide
  x$gene_peptide <- NULL
  x <- as.matrix(x)
  x[is.na(x)] <- 1
  rownames(x) <- rownames
  return(x)
}

#' plots a Figure of Merit curve to help estimate the number of clusters in the results
#' @param r the results object from compare_many()
#' @export
estimate_result_clusters <- function(r) {
  results_mat <- list2mat(r)
  factoextra::fviz_nbclust(results_mat, kmeans, method="wss")
}

#' @export
plot_ridges <- function(r, log = TRUE, base=2, col_order=NULL, all_points = FALSE, sig = 0.05) {

  filtered <- NA
  if (all_points){
    rows_to_keep <- dplyr::bind_rows(r) %>%
      dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " ))
    rows_to_keep <- unique(rows_to_keep$gene_peptide)
    filtered <- lapply(l, drop_columns, sig, metric, log, base, rows_to_keep)
  }
  else {
    filtered <- lapply(l, drop_columns, sig, metric, log, base)
  }


  x <- dplyr::bind_rows(filtered, .id = "comparison") %>%
    dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " ))

  if (log){
    x$fold_change <- log(x$fold_change, base = base)
  }
  if (is.null(col_order)) {
    col_order <- names(r)
  }




  ggplot2::ggplot(x, ggplot2::aes(comparison, gene_id, height=fold_change, group=gene_id, min_height=-4)) +
             ggridges::geom_ridgeline()

}
