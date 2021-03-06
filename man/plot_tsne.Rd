% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tsne.R
\name{plot_tsne}
\alias{plot_tsne}
\title{Plot a tSNE colored by feature.}
\usage{
plot_tsne(data, x = "tSNE_1", y = "tSNE_2", feature,
  legend_name = "", c_cols = pretty_palette, d_cols = NULL,
  pt_size = 0.25, scale_limits = NULL)
}
\arguments{
\item{data}{input data}

\item{x}{x variable}

\item{y}{y variable}

\item{feature}{feature to color by}

\item{legend_name}{legend name to display, defaults to no name}

\item{c_cols}{character vector of colors to built color gradient
for continuous values. defaults to \code{\link[clustifyR:pretty_palette]{clustifyR::pretty_palette}}}

\item{d_cols}{character vector of colors for discrete values.
defaults to RColorBrewer paired palette}

\item{pt_size}{point size}

\item{scale_limits}{defaults to min = 0, max = max(data$x),
otherwise a two-element numeric vector indicating min and max to plot}
}
\description{
Plot a tSNE colored by feature.
}
