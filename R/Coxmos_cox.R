#### ### ##
# METHODS #
#### ### ##

#' cox
#' @description
#' The `cox` function conducts a Cox proportional hazards regression analysis, a type of survival
#' analysis. It is designed to handle right-censored data and is built upon the `coxph` function from
#' the `survival` package. The function returns an object of class "Coxmos" with the attribute model
#' labeled as "cox".
#'
#' @details
#' The Cox proportional hazards regression model is a linear model that describes the relationship
#' between the hazard rate and one or more predictor variables. The function provided here offers
#' several preprocessing steps to ensure the quality and robustness of the model.
#'
#' The function allows for the centering and scaling of predictor variables, which can be essential
#' for the stability and interpretability of the model. It also provides options to remove variables
#' with near-zero or zero variance, which can be problematic in regression analyses. Such variables
#' offer little to no information and can lead to overfitting.
#'
#' Another notable feature is the ability to remove non-significant predictors from the final model
#' through a backward selection process. This ensures that only variables that contribute significantly
#' to the model are retained.
#'
#' The function also checks for the minimum number of events per variable (EPV) to ensure the
#' robustness of the model. If the specified EPV is not met, the function can either halt the
#' computation or proceed based on user preference.
#'
#' It's important to note that while this function is tailored for standard Cox regression, it might
#' not be suitable for high-dimensional data. In such cases, users are advised to consider alternative
#' methods like `coxEN()` or PLS-based Cox methods.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance
#' variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will
#' be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance
#' filtering (default: NULL).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant
#' variables/components in final cox model will be removed until all variables are significant by
#' forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold
#' (default: 0.05).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param FORCE Logical. In case the MIN_EPV is not meet, it allows to compute the model (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "Coxmos" and model "cox". The class contains the following elements:
#'
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(data)}: normalized Y matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix
#'  }
#' \code{survival_model}: List of survival model information
#' \itemize{
#'  \item \code{fit}: coxph object.
#'  \item \code{AIC}: AIC of cox model.
#'  \item \code{BIC}: BIC of cox model.
#'  \item \code{lp}: linear predictors for train data.
#'  \item \code{coef}: Coefficients for cox model.
#'  \item \code{YChapeau}: Y Chapeau residuals.
#'  \item \code{Yresidus}: Y residuals.
#' }
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{nsv}: Variables removed by remove_non_significant if any.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{nz_coeffvar}: Variables removed by coefficient variation near zero.
#'
#' \code{removed_variables_correlation}: Variables removed by being high correlated with other
#' variables.
#'
#' \code{class}: Model class.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @references
#' \insertRef{Cox_1972}{Coxmos}
#' \insertRef{Concato_1995}{Coxmos}
#' \insertRef{survival_package}{Coxmos}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cox(X, Y)
#' cox(X, Y, x.center = TRUE, x.scale = TRUE)
#' }

cox <- function (X, Y,
                 x.center = TRUE, x.scale = FALSE,
                 remove_near_zero_variance = TRUE, remove_zero_variance = FALSE, toKeep.zv = NULL,
                 remove_non_significant = FALSE, alpha = 0.05,
                 MIN_EPV = 5, FORCE = FALSE, returnData = TRUE, verbose = FALSE){

  t1 <- Sys.time()
  FREQ_CUT <- 95/5
  y.center = y.scale = FALSE

  #### Check values classes and ranges
  params_with_limits <- list("alpha" = alpha)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("MIN_EPV" = MIN_EPV)
  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose)
  check_class(logical_params, class = "logical")

  #### Check rownames
  lst_check <- checkXY.rownames(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Check colnames in X for Illegal Chars (affect cox formulas)
  old_colnames <- colnames(X)
  colnames(X) <- transformIllegalChars(old_colnames)
  if(all(old_colnames %in% colnames(X))){
    #colnames changed
    FLAG_COLNAMES = TRUE
  }else{
    FLAG_COLNAMES = FALSE
  }

  #### REQUIREMENTS
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### ZERO VARIANCE - ALWAYS
  if(remove_near_zero_variance || remove_zero_variance){
    lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                            remove_near_zero_variance = remove_near_zero_variance,
                                            remove_zero_variance = remove_zero_variance,
                                            toKeep.zv = toKeep.zv,
                                            freqCut = FREQ_CUT)
    X <- lst_dnz$X
    variablesDeleted <- lst_dnz$variablesDeleted
  }else{
    variablesDeleted <- NULL
  }

  #### COEF VARIATION
  lst_dnzc <- deleteNearZeroCoefficientOfVariation(X = X)
  X <- lst_dnzc$X
  variablesDeleted_cvar <- lst_dnzc$variablesDeleted

  #### MAX PREDICTORS
  if(!check.maxPredictors.cox(X, Y, MIN_EPV, FORCE)){
    return(NA)
  }

  #### SCALING
  lst_scale <- XY.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  #### INITIALISING VARIABLES
  d <- as.data.frame(cbind(Xh, Yh)) #data
  best_cox <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = d,
                      ties = "efron",
                      singular.ok = TRUE,
                      robust = TRUE,
                      nocenter = rep(1, ncol(Xh)),
                      model = TRUE, x = TRUE)
    },
    # Specifying error message
    error = function(e){
      message(paste0("COX: ", e))
      # invisible(gc())
      return(NA)
    }
  )

  cox_trycatch <- function(X){
    tryCatch(
      # Specifying expression
      expr = {
        survival::coxph(formula = survival::Surv(time,event) ~ .,
                        data = X,
                        ties = "efron",
                        singular.ok = TRUE,
                        robust = TRUE,
                        nocenter = rep(1, ncol(X)-2),
                        model = TRUE, x = TRUE)
      },
      # Specifying error message
      error = function(e){
        message(paste0("COX: ", e))
        # invisible(gc())
        return(NA)
      }
    )
  }

  if(all(is.na(best_cox))){
    # Probably by "data contains an infinite predictor"
    message("It is possible that your data has NAs. Cox algorithm cannot manage NAs. Try to clean the data.")
    # lst_cox_uni <- purrr::map(colnames(Xh), ~cox_trycatch(as.data.frame(cbind(Xh[,.], Yh))))
    # names(lst_cox_uni) = colnames(Xh)
    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")

    func_call <- match.call()

    # invisible(gc())
    return(cox_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                          Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                          survival_model = NULL,
                          call = func_call,
                          X_input = if(returnData) X_original else NA,
                          Y_input = if(returnData) Y_original else NA,
                          nsv = NULL,
                          nzv = variablesDeleted,
                          removed_variables_correlation = NULL,
                          class = pkg.env$cox,
                          time = time)))
  }

  # RETURN a MODEL with ALL significant Variables from complete, deleting one by one
  removed_variables <- NULL
  removed_variables_cor <- NULL
  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  if(all(c("time", "event") %in% colnames(d))){
    lst_model <- removeNAorINFcoxmodel(model = best_cox, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAorINFcoxmodel(model = best_cox, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
  }

  best_cox <- lst_model$model
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(d))){
      lst_rnsc <- removeNonSignificativeCox(cox = best_cox, alpha = alpha, cox_input = d, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = best_cox, alpha = alpha, cox_input = cbind(d, Yh), time.value = NULL, event.value = NULL)
    }

    best_cox <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  if(all(c("time", "event") %in% colnames(d))){
    lst_model <- removeNAorINFcoxmodel(model = best_cox, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAorINFcoxmodel(model = best_cox, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
  }
  best_cox <- lst_model$model
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  if(isa(best_cox,"coxph")){
    survival_model <- getInfoCoxModel(best_cox)
  }else{
    survival_model <- NULL
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(cox_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                        Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                        survival_model = survival_model,
                        call = func_call,
                        X_input = if(returnData) X_original else NA,
                        Y_input = if(returnData) Y_original else NA,
                        nsv = removed_variables,
                        removed_variables_correlation = removed_variables_cor,
                        nzv = variablesDeleted,
                        nz_coeffvar = variablesDeleted_cvar,
                        class = pkg.env$cox,
                        time = time)))
}


### ## ##
# CLASS #
### ## ##

cox_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$cox)
  return(model)
}
