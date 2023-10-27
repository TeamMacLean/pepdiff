plot_peptides_measured <- function(i){
  dplyr::group_by(i, treatment, seconds, bio_rep) %>%
    dplyr::summarise(peptide_count = sum(is_useable(.data$quant))) %>%
    ggplot2::ggplot() +
    ggplot2::aes(bio_rep, peptide_count) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::facet_grid(treatment ~ seconds, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs( x= "Biological Replicate", y ="Number of Peptides with Measurements",
                   title = "Peptides observed at each time point in each treatment for each biological replicate"

    )
}


#' plot the representation of peptides in each group.
#'
#' Shows what proportion of the whole set of peptides has a NA value in each group
#' of treatment, seconds, bio rep and tech rep.
#'
#' @param i dataframe with unmerged tech reps; typically from `import_data()`
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_missing_peptides <- function(i){
  assess_missing(i) %>%
    ggplot2::ggplot() +
    ggplot2::aes(.data$bio_rep, .data$tech_rep) +
    ggplot2::geom_tile(  ggplot2::aes(fill = .data$percent_missing)) +
    ggplot2::facet_grid(treatment ~ seconds) +
    ggplot2::scale_fill_viridis_c(guide=ggplot2::guide_legend(title = "% NA")) +
    ggplot2::theme_minimal() +
    ggplot2::labs( x="Biolgical Replicate", y="Technical Replicate",
                   title = "Percent of peptides with a missing quantification value (NA) in each sample"
                   )
}

#' draw density plots for data
#'
#' Plot distribution of quantity for each treatment, seconds and biological
#' replicate. For biological replicate mean of technical replicates is taken
#' @param df dataframe; typically from `import_data()`
#' @param log perform log transform of data
#' @param base base to use in log transform
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_quant_distributions <- function(i, log = FALSE, base = 2){

  i <- combine_tech_reps(i)
  bio_rep_count <- length(unique(i$bio_rep))
  # p <- NULL
  # if(log){
  #   p <- dplyr::mutate(df, log_mean_tr_quant = log(.data$mean_tr_quant, base = base)) %>%
  #     ggplot2::ggplot() +
  #     ggplot2::aes(x = .data$log_mean_tr_quant) +
  #     ggplot2::geom_density(  ggplot2::aes(fill = .data$bio_rep), alpha = I(1/bio_rep_count))
  # } else {
  #   p <- ggplot2::ggplot(df) +
  #     ggplot2::aes(x = .data$mean_tr_quant) +
  #     ggplot2::geom_density(ggplot2::aes(fill = .data$bio_rep),alpha = I(1/bio_rep_count))
  # }
  # p +
  #
  #   ggplot2::facet_grid(treatment ~ seconds) +
  #   ggplot2::scale_fill_viridis_d() +
  #   ggplot2::theme_minimal() +
  #   ggplot2::labs(x=)

  dplyr::mutate(i, mean_tr_quant = dplyr::if_else(log, log(.data$mean_tr_quant, base = base), .data$mean_tr_quant)) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$mean_tr_quant) +
    ggplot2::geom_density(  ggplot2::aes(fill = .data$bio_rep), alpha = I(1/bio_rep_count)) +
    ggplot2::facet_grid(treatment ~ seconds) +
    ggplot2::scale_fill_viridis_d(guide = ggplot2::guide_legend(title="Biological\nReplicate")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x= "Quantity", y="density", title = "Distribution of quantity for biological replicate, time point and treatment.")

}

#' draw qqplots for data
#'
#' Plot qqplot of distribution of quantifications in data for each treatment,
#' seconds and biological replicate
#' @param i dataframe; typically from `import_data()`
#' @param log perform log transform of data
#' @param base base to use in log transform
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_norm_qq <- function(i, log = FALSE, base = 2){
  #return(NULL)
  #df <- combine_tech_reps(df)
  # if(log){
  #   p <- dplyr::mutate(df, log_mean_tr_quant = log(.data$mean_tr_quant, base = base)) %>%
  #     ggplot2::ggplot() +   ggplot2::aes(sample = .data$log_mean_tr_quant)
  # } else {
  #   p <-   ggplot2::ggplot(df) +
  #     ggplot2::aes(sample = .data$mean_tr_quant)
  # }
  combine_tech_reps(i) %>%
  dplyr::mutate(mean_tr_quant = dplyr::if_else(log, log(.data$mean_tr_quant, base = base), .data$mean_tr_quant)) %>%
    ggplot2::ggplot() +
    ggplot2::aes(sample = .data$mean_tr_quant) +
    ggplot2::geom_qq(  ggplot2::aes(colour = .data$bio_rep)) +
    ggplot2::geom_qq_line(  ggplot2::aes(colour = .data$bio_rep)) +
    ggplot2::facet_grid(treatment ~ seconds)  +
    ggplot2::scale_color_viridis_d(guide = ggplot2::guide_legend(title = "Biological\nReplicate")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y="Sampled Values", x="Theoretical Normal Values", title="Normal QQ plot of quantity in biological replicates")
}


#' plot histograms of p-values for each test used
#'
#' @param r dataframe or list of result dataframes typically from `compare()` or `compare_many()`
#' @return ggplot2 plot
#'
#' @export
#' @importFrom rlang .data
plot_p_value_dist <- function(r){

  dplyr::bind_rows(r, .id = "comparison") %>%

    dplyr::select(.data$comparison, .data$gene_id, .data$peptide, .data$treatment_replicates, .data$control_replicates, .data$fold_change, dplyr::ends_with("p_val" ) ) %>%
    tidyr::pivot_longer(dplyr::ends_with("p_val"), names_to = 'test') %>%
      ggplot2::ggplot() +
      ggplot2::aes(.data$value) +
      ggplot2::geom_histogram(bins = 20, fill="steelblue") +
      ggplot2::facet_grid( comparison ~ test, scales = 'free') +
      ggplot2::theme_minimal() +
      ggplot2::labs(x="p-value", y="Frequency")

}

#' plot histogram of fold change distribution for a comparison
#'
#' @param r result dataframe or list, typically from `compare()` or `compare_many()`
#' @param log log the fold change values
#' @param base base of the log, if used
#' @return ggplot2 plot
#' @export
#' @importFrom rlang .data
plot_fc <- function(r, log = TRUE, base = 2){
  #p <- ggplot2::ggplot(l)
  #if (log) {
  #  p <- p + ggplot2::aes(log(.data$fold_change, base = base))
  #} else {
  #  p <- p + ggplot2::aes( .data$fold_change )
  #}

  xtitle <- "Fold Change"
  if (log){ xtitle <- "Log Fold Change"}
  dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::select(.data$comparison, .data$fold_change ) %>%
    dplyr::mutate(do_log = log, fold_change = dplyr::if_else(do_log, log(.data$fold_change, base = base), .data$fold_change)) %>%
     ggplot2::ggplot() +
     ggplot2::aes(fold_change) +
     ggplot2::geom_histogram(bins = 25, fill="steelblue") +
     ggplot2::theme_minimal() +
     ggplot2::facet_wrap( ~ comparison) +
     ggplot2::labs(x=xtitle, y="Frequency")
}

#' plot qqplot of fold changes from a comparison
#' @param r result dataframe or list, typically from `compare()` or `compare_many()`
#' @param log log the fold change values
#' @param base base of the log, if used
#' @export
#' @importFrom rlang .data
plot_fc_qq <- function(r, log = TRUE, base = 2 ){

  # p <- ggplot2::ggplot(df)
  # if (log){
  #   p <- p + ggplot2::aes(sample = log(.data$fold_change, base = base))
  # } else {
  #   p <- p + ggplot2::aes(sample = .data$fold_change)
  # }
  ytitle <- "Sample Values"
  if (log){ytitle <- "Log Sample Values"}
  dplyr::bind_rows(r, .id = "comparison") %>%
    dplyr::select(.data$comparison, .data$fold_change ) %>%
    dplyr::mutate(do_log = log, fold_change = dplyr::if_else(do_log, log(.data$fold_change, base = base), .data$fold_change)) %>%
    ggplot2::ggplot() +
    ggplot2::aes(sample = fold_change) +
    ggplot2::geom_qq(colour = "steelblue", alpha=0.5) +
    ggplot2::geom_qq_line(  ) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~ comparison) +
    ggplot2::labs(y=ytitle, x="Theoretical Normal Values", title="Normal QQ plot of Fold Changes")

}


#' compare sets of significant peptides called by the used data
#'
#' produces an UpSet plot showing intersections and set-size of the different
#' sets of significant peptides called by the methods used in the provided result
#' dataframe
#'
#' @param r result dataframe or list typically from `compare()` or `compare_many()`
#' @param sig significance cut-off to select peptides
#' @return UpSet plot
#'
#' @export
plot_compared_calls <- function(r, sig = 0.05){
  upper_sig <- 1 - round(sig / 2, 4)
  lower_sig <- round(sig / 2, 4)


  mk <- function(df) {
    sets <- list()

    if ("norm_quantile_pval" %in% names(df)) {
      sets$norm_iter = which(df$norm_quantile_p_val >= upper_sig |
                               df$norm_quantile_p_val <= lower_sig)
    }
    if ("bootstrap_t_p_val" %in% names(df)) {
      sets$bootstrap_t = which(df$bootstrap_t_p_val <= sig)
    }
    if ("wilcoxon_p_val" %in% names(df)) {
      sets$wilcoxon = which(df$wilcoxon_p_val <= sig)
    }
    if ("kruskal_p_val" %in% names(df)) {
      sets$kruskal = which(df$kruskal_p_val <= sig)
    }
    if ("rank_prod_p1_p_val" %in% names(df)) {
      sets$rank_prod = union(which(df$rank_prod_p1_p_val <= sig),
                             which(df$rank_prod_p2_p_val <= sig))
    }

    #return(sets)
    if (length(sets) > 1) {
      UpSetR::upset(UpSetR::fromList(sets), order.by = "freq")
    } else {
      warning("Can't do a comparison with only one test performed")
    }
  }

  if (is.data.frame(r)){
    mk(r)
  }

  else{
    l <- lapply(r, mk)
    l <- lapply(l, ggplotify::as.grob)
    cowplot::plot_grid(plotlist=l, labels = names(l))

  }

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
#' @param r list of results, usually from `compare_many()`
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
#' @param col_fontsize size of treatment annotation labels
#' @param col_rowsize size of peptide annotation labels
#' @return NULL
#' @export
#' @importFrom rlang .data
plot_heatmap <- function(r, sig_level = 0.05, metric = "bootstrap_t_fdr", log = TRUE,
                         base = 2, col_order = NULL, sig_only = TRUE, pal="RdBu",
                         lgd_x = 1.7, lgd_y = 1, padding=c(0,0,0,3), col_fontsize=24, row_fontsize=15) {

  if (is.null(col_order)) {
    col_order <- names(r)
  }

  leg_title <- "Fold Change"
  if (log) leg_title <- paste("Log", base, "Fold Change")

  fcm <- fold_change_matrix(r, log=log, base=base, sig_only = sig_only, sig_level=sig_level, metric=metric)
  ul <- max(abs(fcm))
  ll <- ul * -1

  ht <- ComplexHeatmap::Heatmap(fcm, column_order = col_order,
                          col = circlize::colorRamp2(
                            seq(ll,ul, length.out = 11),
                            rev(RColorBrewer::brewer.pal(11, pal)),
                          ),
                          show_heatmap_legend = FALSE,
                          column_names_gp = grid::gpar(fontsize=col_fontsize, fontface="bold"),
                          row_names_gp = grid::gpar(fontsize=row_fontsize, fontface="bold")
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
plot_sample_kmeans <- function(df, nstart = 25, iter.max = 1000){
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
#' @param l list of results data frames, typically from `compare_many()` or single data frame from `compare()`
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
plot_volcano <- function(l, log = FALSE, base = 2, sig_level = 0.05, metric = "bootstrap_t_p_val", option="E", direction=-1  ) {
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

#' plots a Figure of Merit curve to help estimate the number of clusters in the results for a given metric and a
#' given significance level. Clustering is done based on Fold
#' @param r the results object from `compare()` or `compare_many()`
#' @param log log the fold changes
#' @param base base for logging
#' @param sig_level significance cutoff for colour
#' @param metric metric to use for significance
#' @export
plot_cluster_estimate <- function(r, log=TRUE, base=2, sig_only = FALSE, sig_level = 0.05, metric = "bootstrap_t_fdr") {
  results_mat <- fold_change_matrix(r, log=log, base=base, sig_only=sig_only,sig_level=sig_level, metric = metric)
  factoextra::fviz_nbclust(results_mat, stats::kmeans, method="wss")
}


#' plots the kmeans clusters
#' @param kl the kmeans clusters from `kmeans_by_selected_cols()`
#' @param col_order vector of column names in the order you want them to be drawn
#' @param logged indicates whether the fold changes are logged or not
#' @param base if logged the base that was used
#' @param pal pallete to use from RColorBrewer
#' @param nrow number of rows to plot clusters in
#' @param legx x position of legend in range (-0.5 .. 0.5)
#' @param legy y position of legend in range(-0.5 .. 0.5)
#' @param labh horizontal justification of subplot labels
#' @param labv vertical justification of subplot labels
#' @export
plot_kmeans_cluster_hmap <- function(kl, col_order=NULL, logged=TRUE, base=NULL, pal="RdBu", nrow=2, col_fontsize=6, row_fontsize=6, legx=0, legy=-0.45, labh=-1.1, labv=1.1) {


  ul <- max(unlist(lapply(kl,FUN=max)))
  ll <- ul*-1

  make_hmap <- function(mat, col_order, ll, ul, pal) {
    if (is.null(col_order)) {
      col_order <- colnames(mat)
    }

    hm <- ComplexHeatmap::Heatmap(mat,column_order = col_order,
                                      col = circlize::colorRamp2(
                                        seq(ll,ul, length.out = 11),
                                        rev(RColorBrewer::brewer.pal(11, pal)),
                                      ),
                                      show_heatmap_legend = FALSE,
                                  column_names_gp = grid::gpar(fontsize=col_fontsize ),
                                  row_names_gp = grid::gpar(fontsize=row_fontsize))

    grid::grid.grabExpr(
      ComplexHeatmap::draw(hm)
      )

    }

  hms <- lapply(kl, function(mat){ make_hmap(mat, col_order, ll, ul, pal)} )


  titletex <- "Fold Change"
  if (logged){
    titletex <- paste0("Log", base, " Fold Change")
  }

  lgd <- grid::grid.grabExpr(

    ComplexHeatmap::draw(ComplexHeatmap::Legend(direction = "horizontal",
                                col_fun = circlize::colorRamp2(
                                  seq(ll, ul, length.out = 11),
                                  rev(RColorBrewer::brewer.pal(11, pal))
                                ),
                                title = titletex)
    )
  )
  cowplot::plot_grid(plotlist=hms, nrow=nrow,labels = names(hms), hjust=labh, vjust=labv) +
    cowplot::draw_grob(lgd, x=legx, y = legy)

}

plot_kmeans_cluster_profile <- function(kl){

}

