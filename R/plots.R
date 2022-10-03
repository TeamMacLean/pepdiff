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

#' makes heatmap from all experiments, filter on a single metric and sig value
#'
#' reduces dataframes and makes long list, makes a basic heatmap.
#' Use `fold_change_matrix()` to extract data in a heatmappable format
#'
#' @param l list of results, usually from `compare_many()`
#' @param log whether to log the data
#' @param base base used in logging (default = 2)
#' @param sig_only return only rows with 1 or more values significant at `sig_level` of `metric`
#' @param sig_level significance level cutoff
#' @param metric the test metric used to determine significance one of:
#' `bootstrap_t_p_val`, `bootstrap_t_fdr`
#' `wilcoxon_p_val`, `wilcoxon_fdr`
#' `kruskal_p_val`,  `kruskal_fdr`
#' `rank_prod_p1_p_val`, `rank_prod_p2_p_val`, `rank_prod_p1_fdr`, `rank_prod_p2_fdr`.
#' @param col_order specify a column order for the plot, default is names(l)
#' @param pal cbrewer palette to use "RdBu", needs minimum 11 colours
#' @param lgd_x value to pass to ComplexHeatmap::draw for x position of legend in 'in'
#' @param lgd_y value to pass to ComplexHeatmap::draw for y position of legend in 'in'
#' @param padding vector of padding values to pass to ComplexHeatmap::draw for padding of heatmap sections
#' @param lgd_x x offset of legend placement in `in` units
#' @param lgd_y y offset of legend placement in `in` units
#' @return NULL
#' @export
#' @importFrom rlang .data
plot_heatmap <- function(l, sig_level = 0.05, metric = "bootstrap_t_fdr", log = TRUE,
                         base = 2, col_order = NULL, sig_only = TRUE, pal="RdBu",
                         lgd_x = 1.7, lgd_y = 1, padding=c(0,0,0,3)) {

  if (is.null(col_order)) {
    col_order <- names(l)
  }

  leg_title <- "Fold Change"
  if (log) leg_title <- paste("Log", base, "Fold Change")

  fcm <- fold_change_matrix(l, log=log, base=base, sig_only = sig_only, sig_level=sig_level, metric=metric)
  ul <- max(abs(fcm))
  ll <- ul * -1

  ht <- ComplexHeatmap::Heatmap(fcm, column_order = col_order,
                          col = circlize::colorRamp2(
                            seq(ll,ul, length.out = 11),
                            rev(RColorBrewer::brewer.pal(11, pal)),
                          ),
                          show_heatmap_legend = FALSE,
                          column_names_gp = grid::gpar(fontsize=24, fontface="bold"),
                          row_names_gp = grid::gpar(fontsize=15, fontface="bold")
  )
  lgd <- ComplexHeatmap::Legend(direction = "horizontal",
                                col_fun = circlize::colorRamp2(
                                  seq(ll, ul, length.out = 11),
                                  rev(RColorBrewer::brewer.pal(11, "RdBu"))
                                ),
                                legend_width = grid::unit(3, "in"),
                                title = leg_title)

  ComplexHeatmap::draw(ht, padding= grid::unit(padding, "in"))
  ComplexHeatmap::draw(lgd, x = grid::unit(as.numeric(lgd_x), "in"), y = grid::unit(as.numeric(lgd_y), "in"))

}

#' returns a matrix of fold change values
#'
#' Computes the fold change relative to the control sample and returns a
#' matrix with comparisons in columns and peptides in rows. Use this if you want data for a
#' customised heatmap
#'
#' @param l list of results, usually from `compare_many()`
#' @param log whether to log the data
#' @param base base used in logging (default = 2)
#' @param sig_only return only rows with 1 or more values significant at `sig_level` of `metric`
#' @param sig_level significance level cutoff
#' @param metric the test metric used to determine significance one of:
#' `bootstrap_t_p_val`, `bootstrap_t_fdr`
#' `wilcoxon_p_val`, `wilcoxon_fdr`
#' `kruskal_p_val`,  `kruskal_fdr`
#' `rank_prod_p1_p_val`, `rank_prod_p2_p_val`, `rank_prod_p1_fdr`, `rank_prod_p2_fdr`.
#' @return matrix
fold_change_matrix <- function(l, log=TRUE, base=2, sig_only=FALSE, sig_level=0.05, metric="bootstrap_t_fdr") {

  if (!metric %in% metrics() ){
    stop("unknown metric requested.")
  }

  fcm <- list2mat(l, column="fold_change")
  if (log) fcm <- log(fcm, base=base)
  if (!sig_only) return(fcm)


  sigr <- get_sig_rows(l, metric=metric, sig_level=sig_level)
  return( fcm[sigr,] )

}
#' reports metrics available for significance values
metrics <- function() {
  return( c("bootstrap_t_p_val", "bootstrap_t_fdr",
            "wilcoxon_p_val", "wilcoxon_fdr",
            "kruskal_p_val",  "kruskal_fdr",
            "rank_prod_p1_p_val", "rank_prod_p2_p_val", "rank_prod_p1_fdr", "rank_prod_p2_fdr")
          )
}
#' works out if a peptide has at least one significant value across the experiment
#' Composes a matrix of the `metric` significance values with peptides in rows, experiments
#' in columns and works out if each peptide row has a value below the stated cut off
#'
#' #' returns a logical vector of length equal to row number of matrix
#'
#'
#' @param l list of results, usually from `compare_many()`
#' @param sig_level significance level cutoff
#' @param metric the test metric used to determine significance one of:
get_sig_rows <- function(l, metric="bootstrap_t_pval", sig_level=0.05){

  sig_m <- list2mat(l, column=metric)
  ## work out if any rows (peptides) have any sig values <= 0.05 )
  sig_rows <- apply(sig_m, 1, function(x) {any(x <= sig_level )})

  ## deal with NAs from uncomplete bootstraps
  sig_rows[is.na(sig_rows)] <- FALSE
  return(sig_rows)
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
#' @param log log the data
#' @param base base for logging
#' @param sig_level significance cutoff for colour
#' @param metric metric to use for significance
#' @param option viridis colour scheme key to use
#' @param direction viridis colour scheme direction (1/-1)
#' @return ggplot2 plot
#' @export
#'
#' @importFrom rlang .data
volcano_plot <- function(l, log = FALSE, base = 2, sig_level = 0.05, metric = "bootstrap_t_p_val", option="E", direction=-1  ) {
  xlabtxt <- paste0("Log ",base," Fold Change")
  ylabtxt <- paste0("-Log ",base, " P")
  dplyr::bind_rows(l, .id = "comparison")  %>%
    dplyr::select(.data$comparison, .data$gene_id, .data$peptide, .data$treatment_mean_count, .data$control_mean_count, .data$fold_change, dplyr::starts_with(metric)) %>%
    dplyr::rename(p_val = dplyr::ends_with('p_val') ) %>%
    dplyr::mutate(
                  gene_peptide = paste(.data$gene_id, .data$peptide, sep = " "),
                  log_fc = log(.data$fold_change, base = base),
                  log_p = log(.data$p_val, base=base) * -1,
                  signal = log((.data$treatment_mean_count + .data$control_mean_count), base) * -1,
                  change = dplyr::if_else(.data$log_fc > 0 & .data$p_val <= sig_level, "Up",
                            dplyr::if_else(.data$log_fc < 0 & .data$p_val <= sig_level, "Down", "None")),
                  change = forcats::fct_relevel(.data$change, rev)
                  ) %>%
    dplyr::select(.data$comparison, .data$gene_peptide, .data$log_fc, .data$log_p, .data$change) %>%
    ggplot2::ggplot() +
    ggplot2::aes(.data$log_fc, .data$log_p) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$change)) +
    ggplot2::facet_wrap( ~ .data$comparison) +
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_viridis_d(option=option, direction=direction) +
    ggplot2::xlab(xlabtxt) +
    ggplot2::ylab(ylabtxt) +
    ggplot2::labs(colour = "Change")




}

#' converts a results object to a matrix as if for direct use in external heatmap functions
#' @param r results object, usually from `compare_many()`
#' @param column column from results data to put into matrix, default = "fold_change"
list2mat <- function(r,column="fold_change") {
  x <- dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::mutate(gene_peptide = paste(.data$gene_id, .data$peptide, sep = " " )) %>%
    dplyr::select(-.data$gene_id, -.data$peptide) %>%
    tidyr::pivot_wider(.data$gene_peptide, names_from = .data$comparison, values_from = .data[[column]] ) %>%
    dplyr::arrange()

  rownames <- x$gene_peptide
  x$gene_peptide <- NULL
  x <- as.matrix(x)
  x[is.na(x)] <- 1
  rownames(x) <- rownames
  return(x)
}

#' plots a Figure of Merit curve to help estimate the number of clusters in the results
#' @param r the results object from `compare_many()`
#' @export
estimate_result_clusters <- function(r) {
  results_mat <- list2mat(r)
  factoextra::fviz_nbclust(results_mat, stats::kmeans, method="wss")
}


