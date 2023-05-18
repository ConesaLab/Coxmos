#### ### ##
# METHODS #
#### ### ##

#' sPLSDA-COX
#' @description This function performs a sparse partial least squares discriminant analysis Cox (sPLSDA-COX) by dynamic variable selection methodology.
#' The function returns a HDcox model with the attribute model as "sPLSDA-COX".
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 10).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of variables to use (default: 10).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of variables to use (default: 1000).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables. If only two cut points are selected, minimum and maximum size are used (default: 5)
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to continue evaluating higher values in the multiple tested parameters. If it is not reached for next 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops (default: 0.01).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param EVAL_METHOD Character. If EVAL_METHOD = "AUC", AUC metric will be use to compute the best number of variables. In other case, c-index metrix will be used (default: "AUC").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "sPLS-DACOX-Dynamic". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: sPLS weights
#'  \item \code{(W.star)}: sPLS W* vector
#'  \item \code{(loadings)}: sPLS loadings
#'  \item \code{(scores)}: sPLS scores/variates
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#'
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix'
#'  }
#'
#' \code{survival_model}: List of survival model information.
#' \itemize{
#'  \item \code{fit}: coxph object.
#'  \item \code{AIC}: AIC of cox model.
#'  \item \code{BIC}: BIC of cox model.
#'  \item \code{lp}: linear predictors for train data.
#'  \item \code{coef}: Coefficients for cox model.
#'  \item \code{YChapeau}: Y Chapeau residuals.
#'  \item \code{Yresidus}: Y residuals.
#' }
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{n.varX}: Number of Variables selected in each PLS component.
#'
#' \code{var_by_component}: Variables selected in each PLS component.
#'
#' \code{plot_accuracyPerVariable}: If NULL vector is selected, return a plot for understanding the number of variable selection.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{alpha}: alpha value selected
#'
#' \code{removed_variables_cox}: Variables removed by sparse penalty.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{class}: Model class.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' splsdacox_dynamic(X, Y)
#' splsdacox_dynamic(X, Y, n.comp = 3, vector = NULL, x.center = TRUE, x.scale = TRUE)
#' }

splsdacox_dynamic <- function (X, Y,
                               n.comp = 4, vector = NULL,
                               MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                               MIN_AUC_INCREASE = 0.01,
                               x.center = TRUE, x.scale = FALSE,
                               remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                               remove_non_significant = F, alpha = 0.05, tol = 1e-10,
                               EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                               times = NULL, max_time_points = 15,
                               MIN_EPV = 5, returnData = T, verbose = F){

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### Check values classes and ranges
  params_with_limits <- list("alpha" = alpha, "MIN_AUC_INCREASE" = MIN_AUC_INCREASE)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("n.comp" = n.comp, "MIN_NVAR" = MIN_NVAR, "MAX_NVAR" = MAX_NVAR, "n.cut_points" = n.cut_points,
                  "max_time_points" = max_time_points,
                  "MIN_EPV" = MIN_EPV, "tol" = tol, "max.iter" = max.iter)
  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose)
  check_class(logical_params, class = "logical")

  character_params <- list("EVAL_METHOD" = EVAL_METHOD, "pred.method" = pred.method)
  check_class(character_params, class = "character")

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
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = FREQ_CUT)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  #### SCALING
  lst_scale <- XY.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  #### MAX PREDICTORS
  n.comp <- check.maxPredictors(X, Y, MIN_EPV, n.comp)

  #### ### ### ### ### ### ### ### ### ###
  # DIVIDE Y VENCERAS - BEST VECTOR SIZE #
  #### ### ### ### ### ### ### ### ### ###
  plotVAR <- NULL
  DR_coxph <- NULL

  if(is.null(vector)){
    lst_BV <- getBestVector(Xh, DR_coxph, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                            cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F,
                            mode = "splsda", times = times, max_time_points = max_time_points, verbose = verbose)
    keepX <- lst_BV$best.keepX
    plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
  }else{
    if(is.numeric(vector)){
      keepX <- vector
      if(length(keepX)>1){
        message("keepX must be a number, not a vector. Maximum value will be selected for compute the sPLS model.")
        keepX <- max(keepX)
      }

      if(keepX>ncol(X)){
        message("keepX must be a lesser than the number of columns in X. The value will be updated to that one.")
        keepX <- ncol(X)
      }
    }else{
      message("Vector does not has the proper structure. Optimizing best n.variables by using your vector as start vector.")
      lst_BV <- getBestVector(Xh, DR_coxph, Yh, n.comp, max.iter, vector = NULL, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                             EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F, mode = "splsda", times = times, max_time_points = max_time_points, verbose = verbose)
      keepX <- lst_BV$best.keepX
      plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
    }
  }

  #### ### ### ### ### ### ### ### ### ### ### ###
  ### ##             PLSDA-COX             ###  ##
  #### ### ### ### ### ### ### ### ### ### ### ###

  splsda <- mixOmics::splsda(Xh, Yh[,"event"], scale=F, ncomp = n.comp, keepX = rep(keepX, n.comp), max.iter = max.iter, near.zero.var = T)

  #last model includes all of them
  tt_splsDR = data.matrix(splsda$variates$X)
  ww_splsDR = data.matrix(splsda$loadings$X)
  pp_splsDR = data.matrix(splsda$mat.c)

  colnames(tt_splsDR) <- paste0("comp_", 1:n.comp)

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ### #
  ##  ##              PLS-COX            ##  ##
  #### ### ### ### ### ### ### ### ### ### ### #
  d <- as.data.frame(cbind(tt_splsDR, Yh))
  cox_model <- NULL
  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = d,
                      ties = "efron",
                      singular.ok = T,
                      robust = T,
                      nocenter = rep(1, ncol(d)),
                      model=T, x = T)
    },
    # Specifying error message
    error = function(e){
      message(paste0("splsdacox_dynamic: ", e))
      # invisible(gc())
      return(NA)
    }
  )

  # RETURN a MODEL with ALL significant Variables from complete, deleting one by one
  removed_variables <- NULL
  removed_variables_cor <- NULL
  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  if(all(c("time", "event") %in% colnames(d))){
    lst_model <- removeNAcoxmodel(model = cox_model$fit, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAcoxmodel(model = cox_model$fit, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
  }
  cox_model$fit <- lst_model$model
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(d))){
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$fit, alpha = alpha, cox_input = d, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$fit, alpha = alpha, cox_input = cbind(d, Yh), time.value = NULL, event.value = NULL)
    }

    cox_model$fit <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  survival_model = NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  #get W.star
  W <- ww_splsDR
  P <- pp_splsDR

  if(is.null(P) | is.null(W)){
    message(paste0(pkg.env$splsdacox_dynamic, " model cannot be computed because P or W vectors are NULL. Returning NA."))
    # invisible(gc())
    return(NA)
  }

  #W.star
  #sometimes solve(t(P) %*% W)
  #system is computationally singular: reciprocal condition number = 6.24697e-18
  PW <- tryCatch(expr = {solve(t(P) %*% W, tol = tol)},
                 error = function(e){
                   if(verbose){
                     message(e$message)
                   }
                   NA
                 })

  if(all(is.na(PW))){
    message(paste0(pkg.env$splsdacox_dynamic, " model cannot be computed due to solve(t(P) %*% W). Multicollineality could be present in your data. Optional (not recommended): Reduce 'tol' parameter to fix it. Returning NA."))
    # invisible(gc())
    return(NA)
  }

  # What happen when you cannot compute W.star but you have P and W?
  W.star <- W %*% PW

  Ts <- tt_splsDR

  func_call <- match.call()

  rownames(Ts) <- rownames(X)
  #rownames(P) <- rownames(W) <-  rownames(W.star) <- colnames(X) #as some variables cannot be selected, that name does not work

  colnames(Ts) <- colnames(P) <- colnames(W) <-  colnames(W.star) <- paste0("comp_", 1:n.comp)

  # variable per component
  var_by_component = list()
  for(cn in colnames(W)){
    var_by_component[[cn]] <- rownames(W)[W[,cn]!=0]
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(splsdacox_dynamic_class(list(X = list("data" = if(returnData) X_norm else NA, "weightings" = W, "W.star" = W.star, "loadings" = P, "scores" = Ts, "x.mean" = xmeans, "x.sd" = xsds),
                                      Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                                      survival_model = survival_model,
                                      n.comp = n.comp, #number of components
                                      n.varX = keepX,
                                      var_by_component = var_by_component,
                                      plot_accuracyPerVariable = plotVAR,
                                      call = func_call,
                                      X_input = if(returnData) X_original else NA,
                                      Y_input = if(returnData) Y_original else NA,
                                      alpha = alpha,
                                      removed_variables_cox = removed_variables,
                                      nzv = variablesDeleted,
                                      class = pkg.env$splsdacox_dynamic,
                                      time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation splsdacox_dynamic
#' @description plsdacox_dynamic cross validation model
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of variables to use (default: 10).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of variables to use (default: 1000).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables. If only two cut points are selected, minimum and maximum size are used (default: 5)
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to continue evaluating higher values in the multiple tested parameters. If it is not reached for next 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops (default: 0.01).
#' @param EVAL_METHOD Character. If EVAL_METHOD = "AUC", AUC metric will be use to compute the best number of variables. In other case, c-index metrix will be used (default: "AUC").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_c.index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test observations. Once all have their linear predictors, the evaluation is perform across all the observations together (default: FALSE).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "HDcox" and model "cv.sPLS-DACOX-Dynamic".
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha information).
#'
#' \code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
#' \code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.
#'
#' \code{opt.comp}: Optimal component selected by the best_model.
#' \code{opt.nvar}: Optimal number of variables selected by the best_model.
#'
#' \code{plot_AIC}: AIC plot by each hyper-parameter.
#' \code{plot_c_index}: C-Index plot by each hyper-parameter.
#' \code{plot_BRIER}: Brier Score plot by each hyper-parameter.
#' \code{plot_AUC}: AUC plot by each hyper-parameter.
#'
#' \code{class}: Cross-Validated model class.
#'
#' \code{lst_train_indexes}: List (of lists) of indexes for the observations used in each run/fold for train the models.
#' \code{lst_test_indexes}: List (of lists) of indexes for the observations used in each run/fold for test the models.
#'
#' \code{time}: time consumed for running the cross-validated function.
#' @export
#'
#' @examples
#' \dontrun{
#' cv.splsdacox_dynamic_model <- cv.splsdacox_dynamic(X, Y, max.ncomp = 10, vector = NULL,
#' x.center = TRUE, x.scale = TRUE)
#' splsdacox_model <- splsdacox(X, Y, n.comp = cv.splsdacox_dynamic_model$opt.comp,
#' vector = cv.splsdacox_dynamic_model$opt.nvar, x.center = TRUE, x.scale = TRUE)
#' }

cv.splsdacox_dynamic <- function(X, Y,
                        max.ncomp = 10, vector = NULL,
                        n_run = 5, k_folds = 10,
                        x.center = TRUE, x.scale = FALSE,
                        remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                        remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                        MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                        MIN_AUC_INCREASE = 0.01,
                        EVAL_METHOD = "AUC",
                        w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                        MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                        pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                        max.iter = 500,
                        MIN_EPV = 5, return_models = F, returnData = F, tol = 1e-10,
                        PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### ### ###
  # WARNINGS #
  #### ### ###

  #### Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### Check values classes and ranges
  params_with_limits <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                 "w_AIC" = w_AIC, "w_c.index" = w_c.index, "w_AUC" = w_AUC, "w_BRIER" = w_BRIER)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("max.ncomp" = max.ncomp, "MIN_NVAR" = MIN_NVAR, "MAX_NVAR" = MAX_NVAR, "n.cut_points" = n.cut_points,
                  "n_run" = n_run, "k_folds" = k_folds, "max_time_points" = max_time_points,
                  "MIN_COMP_TO_CHECK" = MIN_COMP_TO_CHECK, "MIN_EPV" = MIN_EPV, "seed" = seed, "tol" = tol)
  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_variance_at_fold_level" = remove_variance_at_fold_level,
                      "remove_non_significant_models" = remove_non_significant_models,
                      "remove_non_significant" = remove_non_significant,
                      "return_models" = return_models,"returnData" = returnData, "verbose" = verbose, "PARALLEL" = PARALLEL)
  check_class(logical_params, class = "logical")

  character_params <- list("EVAL_METHOD" = EVAL_METHOD, "pred.attr" = pred.attr, "pred.method" = pred.method)
  check_class(character_params, class = "character")

  #### Check cv-folds
  lst_checkFR <- checkFoldRuns(Y, n_run, k_folds, fast_mode)
  n_run <- lst_checkFR$n_run
  fast_mode <- lst_checkFR$fast_mode

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars(X)

  #### REQUIREMENTS
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  check.cv.weights(c(w_AIC, w_c.index, w_BRIER, w_AUC))
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### MAX PREDICTORS
  max.ncomp <- check.ncomp(X, max.ncomp)
  max.ncomp <- check.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)
  if(MIN_COMP_TO_CHECK >= max.ncomp){
    MIN_COMP_TO_CHECK = max.ncomp-1
  }

  #### REQUIREMENTS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
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

  #### #
  # CV #
  #### #
  lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds, seed = seed) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test
  k_folds <- lst_data$k_folds

  lst_train_indexes <- lst_data$lst_train_index
  message(print(lst_train_indexes)) #!!!
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###

  total_models <- 1 * k_folds * n_run

  comp_model_lst <- get_HDCOX_models2.0(method = pkg.env$splsdacox_dynamic,
                                        lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                        max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL, max.variables = NULL, vector = vector,
                                        n_run = n_run, k_folds = k_folds,
                                        MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE, EVAL_METHOD = EVAL_METHOD,
                                        n.cut_points = n.cut_points,
                                        x.center = x.center, x.scale = x.scale,
                                        y.center = y.center, y.scale = y.scale,
                                        remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                        alpha = alpha, MIN_EPV = MIN_EPV,
                                        remove_non_significant = remove_non_significant, tol = tol, max.iter = max.iter,
                                        returnData = returnData, total_models = total_models,
                                        PARALLEL = PARALLEL, verbose = verbose)

  if(all(is.null(comp_model_lst))){
    message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdacox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdacox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst, alpha = alpha,
                                                    max.ncomp = max.ncomp, eta.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models, verbose = verbose)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdacox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdacox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # EVALUATING BRIER SCORE #
  #### ### ### ### ### ### #
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_comp_flag <- F

  if(TRUE){ #compute always BRIER SCORE
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_BRIER(comp_model_lst = comp_model_lst,
                                       fast_mode = fast_mode,
                                       lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                       df_results_evals = df_results_evals, times = times,
                                       pred.method = pred.method, pred.attr = pred.attr,
                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       w_BRIER = w_BRIER, method.train = pkg.env$splsdacox_dynamic, PARALLEL = F, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    #total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp)#inside get_COX_evaluation_AUC

    #times should be the same for all folds
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, method.train = pkg.env$splsdacox_dynamic, PARALLEL = F, verbose = verbose)

    if(is.null(df_results_evals_comp)){
      df_results_evals_comp <- lst_df$df_results_evals_comp
    }else{
      df_results_evals_comp$AUC <- lst_df$df_results_evals_comp$AUC
    }

    if(is.null(df_results_evals_run)){
      df_results_evals_run <- lst_df$df_results_evals_run
    }else{
      df_results_evals_run$AUC <- lst_df$df_results_evals_run$AUC
    }

    if(is.null(df_results_evals_fold)){
      df_results_evals_fold <- lst_df$df_results_evals_fold
    }else{
      df_results_evals_fold$AUC <- lst_df$df_results_evals_fold$AUC
    }

    optimal_comp_index <- lst_df$optimal_comp_index
    optimal_comp_flag <- lst_df$optimal_comp_flag
  }

  #### ### ### #
  # BEST MODEL #
  #### ### ### #

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_BRIER, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[optimal_comp_index,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  ### ### #
  # PLOTS #
  ### ### #
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_BRIER = w_BRIER, max.ncomp = max.ncomp,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER", x.text = "Component")

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp

  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_BRIER <- lst_EVAL_PLOTS$ggp_BRIER
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  #### ### #
  # RETURN #
  #### ### #
  best_model_info$n.var <- as.numeric(as.character(best_model_info$n.var)) #just in case be a factor

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  if(return_models){
    return(cv.splsdacox_dynamic_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class= pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.splsdacox_dynamic_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsdacox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }

}

### ## ##
# CLASS #
### ## ##

splsdacox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$splsdacox_dynamic)
  return(model)
}

cv.splsdacox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.splsdacox_dynamic)
  return(model)
}
