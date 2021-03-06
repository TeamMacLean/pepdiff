% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{norm_qqplot}
\alias{norm_qqplot}
\title{draw qqplots for data}
\usage{
norm_qqplot(df, log = FALSE, base = 2)
}
\arguments{
\item{df}{dataframe; typically from `import_data()`}

\item{log}{perform log transform of data}

\item{base}{base to use in log transform}
}
\value{
ggplot2 plot
}
\description{
Plot qqplot of distribution of quantifications in data for each treatment,
seconds and biological replicate
}
