% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{plot_observation.eventHistogram}
\alias{plot_observation.eventHistogram}
\title{plot_observation.eventHistogram}
\usage{
plot_observation.eventHistogram(
  observation,
  model,
  time = NULL,
  type = "lp",
  size = 3,
  color = "red"
)
}
\arguments{
\item{observation}{Numeric matrix or data.frame. New explanatory variables (raw data) for one
observation. Qualitative variables must be transform into binary variables.}

\item{model}{Coxmos model.}

\item{time}{Numeric. Time point where the AUC will be evaluated (default: NULL).}

\item{type}{Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").}

\item{size}{Numeric. Point size (default: 3).}

\item{color}{String. R Color.}
}
\value{
A ggplot object representing a histogram of the predicted event values based on the
provided Coxmos model.
}
\description{
Generates a histogram plot for observation event data based on a given Coxmos model. The
function visualizes the distribution of predicted values and highlights the prediction for a
specific observation.
}
\details{
The \code{plot_observation.eventHistogram} function is designed to provide a visual representation
of the distribution of predicted event values based on a Coxmos model. The function takes in observation
data, a specified time point, and a Coxmos model to compute the prediction. The resulting histogram
plot displays the distribution of these predictions, with a specific emphasis on the prediction
for the provided observation data. The prediction is represented as a point on the histogram, allowing
for easy comparison between the specific observation's prediction and the overall distribution of
predictions. The type of prediction ("lp", "risk", "expected", or "survival") can be specified,
offering flexibility in the kind of insights one wishes to derive from the visualization. The
appearance of the point representing the observation's prediction can be customized using the \code{size}
and \code{color} parameters.
}
\examples{
data("X_proteomic")
data("Y_proteomic")
set.seed(123)
index_train <- caret::createDataPartition(Y_proteomic$event, p = .5, list = FALSE, times = 1)
X_train <- X_proteomic[index_train,1:50]
Y_train <- Y_proteomic[index_train,]
X_test <- X_proteomic[-index_train,1:50]
Y_test <- Y_proteomic[-index_train,]
coxEN.model <- coxEN(X_train, Y_train, x.center = TRUE, x.scale = TRUE)
observation = X_test[1,,drop=FALSE]
plot_observation.eventHistogram(observation = observation, model = coxEN.model, time = NULL)
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
