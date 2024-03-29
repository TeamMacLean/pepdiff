% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{get_sig_rows}
\alias{get_sig_rows}
\title{works out if a peptide has at least one significant value across the experiment
Composes a matrix of the `metric` significance values with peptides in rows, experiments
in columns and works out if each peptide row has a value below the stated cut off}
\usage{
get_sig_rows(l, metric = "bootstrap_t_pval", sig_level = 0.05)
}
\arguments{
\item{l}{list of results, usually from `compare_many()`}

\item{metric}{the test metric used to determine significance one of:}

\item{sig_level}{significance level cutoff}
}
\description{
#' returns a logical vector of length equal to row number of matrix
}
