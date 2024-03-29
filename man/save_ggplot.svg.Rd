% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{save_ggplot.svg}
\alias{save_ggplot.svg}
\title{save_ggplot.svg}
\usage{
save_ggplot.svg(
  plot,
  folder = NULL,
  name = "plot",
  wide = TRUE,
  quality = "4K",
  dpi = 80,
  custom = NULL
)
}
\arguments{
\item{plot}{'ggplot2' object. Object to plot and save.}

\item{folder}{Character. Folder path as character type.}

\item{name}{Character. File name.}

\item{wide}{Logical. If TRUE, widescreen format (16:9) is used, in other case (4:3) format.}

\item{quality}{Character. One of: "HD", "FHD", "2K", "4K", "8K"}

\item{dpi}{Numeric. Dpi value for the image.}

\item{custom}{Numeric vector. Custom size of the image. Numeric vector of width and height.}
}
\value{
Generate as many plot images as list objects in the specific folder or working directory.
}
\description{
Allows to save 'ggplot2' objects in .svg format based on an specific resolution.
}
\examples{
\dontrun{
library(ggplot2)
data(iris)
g <- ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species))
g <- g + geom_point(size = 4)
save_ggplot.svg(g)
}
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
