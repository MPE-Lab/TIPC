% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trend_plot_hexLen.R
\name{trend_plot_hexLen}
\alias{trend_plot_hexLen}
\title{Trend plot of grid size effect}
\usage{
trend_plot_hexLen(
  root_dir = NULL,
  hex_len_range = NULL,
  hex_len_stepsize = NULL,
  output_dir = NULL,
  pdf_width = 18,
  pdf_height = 12
)
}
\arguments{
\item{root_dir}{A directory path pointing to the parent directory of the
ouptput folders containing TIPC metric counts generated by
\code{\link[TIPC]{count_TIPC_cat}} at various hexagonal lengths.}

\item{hex_len_range}{A vector of 2 integer number indicating the range of
hexagonal length. If unspecified (NULL), all found TIPC result folders of
different hexagonal lengths will be processed.}

\item{hex_len_stepsize}{An integer indicating the step size of hexagonal
length. If unspecified (NULL), all found TIPC result folders of
different hexagonal lengths will be processed.}

\item{output_dir}{A directory path for saving output plot; if NULL, plots will
be saved in root_dir}

\item{pdf_width, pdf_height}{the width and height of the graphics region in inches.
The default values are 18 and 12; inherited from \code{\link[grDevices]{pdf}}}
}
\description{
Plotting the changes in TIPC metric count as the grid size i.e. hexagon
length increases; the 5 shifting directional data will
  be collapsed specific hex_len.
}
