% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{volcano_plot}
\alias{volcano_plot}
\title{volcano plot the data}
\usage{
volcano_plot(
  l,
  log = FALSE,
  base = 2,
  by = "peptide",
  sig = 0.05,
  metric = NA,
  option = "E",
  direction = -1
)
}
\arguments{
\item{l}{list of results data frames, typically from `compare_many()`}

\item{metric}{single metric to use for volcano plot}
}
\value{
ggplot2 plot
}
\description{
draws a plot of peptide count against log fc at either protein or peptide level for samples
}
