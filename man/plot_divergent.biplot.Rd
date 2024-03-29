% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{plot_divergent.biplot}
\alias{plot_divergent.biplot}
\title{plot_divergent.biplot}
\usage{
plot_divergent.biplot(
  X,
  Y,
  NAMEVAR1,
  NAMEVAR2,
  BREAKTIME,
  x.text = "N. of Samples"
)
}
\arguments{
\item{X}{Numeric matrix or data.frame. Explanatory variables with "NAMEVAR1" and "NAMEVAR2"
variables. "NAMEVAR1" must be a factor variable.}

\item{Y}{Numeric matrix or data.frame. Response variables. Object must have two columns named as
"time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
event observations.}

\item{NAMEVAR1}{Character. Factor variable name (must be located in colnames(X) and have to have
two levels).}

\item{NAMEVAR2}{Character. Numerical variable name (must be located in colnames(X)).}

\item{BREAKTIME}{Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points.
If BREAKTIME = NULL, "n.breaks" is used (default: NULL).}

\item{x.text}{Character. Title for X axis.}
}
\value{
A 'ggplot2' two side bar plot. X axis represent the number of samples per each NAMEVAR1
factor levels and Y axis, the X NAMEVAR2 numerical variables categorize in groups of breaks.
}
\description{
Generates a divergent biplot visualizing the distribution of a qualitative variable
against a quantitative variable, further categorized by an event matrix.
}
\details{
The function \code{plot_divergent.biplot} is designed to offer a comprehensive visualization
of the relationship between a qualitative and a quantitative variable, while also taking into
account an associated event matrix. The qualitative variable, denoted by "NAMEVAR1", is expected
to be a factor with two levels, and the quantitative variable, "NAMEVAR2", is numerically
represented. The event matrix, "Y", consists of two columns: "time" and "event". The "event"
column indicates whether an observation is censored or an event, represented by binary values
(0/1 or FALSE/TRUE).

The function processes the input data to categorize the quantitative variable into groups based
on the specified "BREAKTIME" parameter. Each group represents a range of values for the quantitative
variable. The resulting plot displays the number of samples for each level of the qualitative
variable on the X-axis, while the Y-axis represents the categorized groups of the quantitative
variable. The bars in the plot are further colored based on the event type, providing a clear
distinction between censored and event observations.
}
\examples{
\dontrun{
  NAMEVAR1 = "sex"
  NAMEVAR2 = "age"
  plot_divergent.biplot(X, Y, NAMEVAR1, NAMEVAR2, BREAKTIME = 5, x.text = "N. of Patients")
}
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
