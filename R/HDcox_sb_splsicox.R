#### ### ##
# METHODS #
#### ### ##

#' SB.sPLS-ICOX
#' @description This function performs a single-block sparse partial least squares individual Cox (SB.sPLS-ICOX).
#' The function returns a HDcox model with the attribute model as "SB.sPLS-ICOX".
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 10).
#' @param spv_penalty Numeric. Penalty for variable selection for the individual cox models. Variables with a lower P-Value than "spv_penalty" in the individual cox analysis will be keep for the sPLS-ICOX approach (default: 1).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "sb.splsicox". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: PLS weights
#'  \item \code{(weightings_norm)}: PLS normalize weights
#'  \item \code{(W.star)}: PLS W* vector
#'  \item \code{(scores)}: PLS scores/variates
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
#'  \item \code{(dr.mean)}: mean values for deviance residuals Y matrix
#'  \item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix'
#'  }
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
#' \code{list_pls_models}: List of sPLS-ICOX models computed for each block.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
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
#' sb.splsicox(X, Y)
#' sb.splsicox(X, Y, n.comp = 3, spv_penalty = 0.5, x.center = TRUE, x.scale = TRUE)
#' }

sb.splsicox <- function (X, Y,
                        n.comp = 4, spv_penalty = 1,
                        x.center = TRUE, x.scale = FALSE,
                        y.center = FALSE, y.scale = FALSE,
                        remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                        remove_non_significant = F, alpha = 0.05, tol = 1e-10,
                        MIN_EPV = 5, returnData = T, verbose = F){

  t1 <- Sys.time()

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### Check values classes and ranges
  lst_01 <- list("alpha" = alpha, "eta" = spv_penalty)
  check_min0_max1_variables(lst_01)

  lst_num <- list("n.comp" = n.comp,
                  "MIN_EPV" = MIN_EPV, "tol" = tol)
  check_class(lst_num, class = "numeric")

  lst_logical <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                      "y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose)
  check_class(lst_logical, class = "logical")

  #### REQUIREMENTS
  lst_check <- checkXY.mb.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  checkY.colnames(Y)

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
                                            remove_near_zero_variance = remove_near_zero_variance,
                                            remove_zero_variance = remove_zero_variance,
                                            toKeep.zv = toKeep.zv,
                                            freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  #### SCALING
  lst_scale <- XY.mb.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  #### MAX PREDICTORS
  n.comp <- check.mb.maxPredictors(X, Y, MIN_EPV, n.comp, verbose = verbose)

  # CREATE INDIVIDUAL MODELS
  lst_sb.pls <- list()
  for(b in names(Xh)){
    lst_sb.pls[[b]] <- splsicox(X = Xh[[b]], Y = Yh, n.comp = n.comp, spv_penalty = spv_penalty,
                               x.scale = F, x.center = F, y.scale = F, y.center = F,
                               remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL, #zero_var already checked
                               remove_non_significant = remove_non_significant, alpha = alpha, tol = tol,
                               returnData = F, verbose = verbose)
  }

  # CHECK ALL MODELS SAME COMPONENTS
  aux_ncomp <- purrr::map(lst_sb.pls, ~.$n.comp)

  # CREATE COMBINE MODEL
  data <- NULL
  cn.merge <- NULL
  for(b in names(Xh)){
    if(!is.null(lst_sb.pls[[b]]$survival_model)){
      data <- cbind(data, lst_sb.pls[[b]]$X$scores)
      cn.merge <- c(cn.merge, paste0(colnames(lst_sb.pls[[b]]$X$scores), "_", b))
    }else{
      next
    }
  }

  #colnames(data) <- apply(expand.grid(colnames(lst_sb.pls[[1]]$X$scores), names(Xh)), 1, paste, collapse="_")
  colnames(data) <- cn.merge
  cox_model <- cox(X = data, Y = Yh,
                   x.center = F, x.scale = F,
                   y.center = F, y.scale = F,
                   remove_near_zero_variance = F, remove_zero_variance = F,
                   remove_non_significant = remove_non_significant, FORCE = T)

  #### ### #
  # RETURN #
  #### ### #
  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(sb.splsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "x.mean" = xmeans, "x.sd" = xsds),
                                Y = list("data" = lst_sb.pls[[1]]$Y$data, "y.mean" = ymeans, "y.sd" = ysds),
                                survival_model = cox_model$survival_model,
                                list_pls_models = lst_sb.pls,
                                n.comp = n.comp, #number of components used, but could be lesser than expected because not computed models
                                spv_penalty = spv_penalty,
                                call = func_call,
                                X_input = if(returnData) X_original else NA,
                                Y_input = if(returnData) Y_original else NA,
                                nzv = variablesDeleted,
                                class = pkg.env$sb.splsicox,
                                time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation cv.sb.splsicox
#' @description cv.sb.splsicox cross validation model
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param spv_penalty.list Numeric vector. Penalty for variable selection for the individual cox models. Variables with a lower P-Value than "spv_penalty" in the individual cox analysis will be keep for the sPLS-ICOX approach (default: seq(0.1,1,0.1)).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_c.index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to continue evaluating higher values in the multiple tested parameters. If it is not reached for next 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops (default: 0.01).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test observations. Once all have their linear predictors, the evaluation is perform across all the observations together (default: FALSE).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "HDcox" and model "cv.SB.sPLS-ICOX".
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha information).
#'
#' \code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
#' \code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.
#'
#' \code{opt.comp}: Optimal component selected by the best_model.
#' \code{opt.spv_penalty}: Optimal penalty value selected by the best_model.
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
#' cv.sb.splsicox_model <- cv.sb.splsicox(X, Y, max.ncomp = 10,
#' spv_penalty.list = seq(0.1,1,0.1), x.center = TRUE, x.scale = TRUE)
#' sb.splsicox_model <- sb.splsicox(X, Y, n.comp = cv.sb.splsicox_model$opt.comp,
#' spv_penalty = cv.sb.splsicox_model$opt.spv_penalty, x.center = TRUE, x.scale = TRUE)
#' }


cv.sb.splsicox <- function(X, Y,
                          max.ncomp = 10, spv_penalty.list = seq(0.1,1,0.1),
                          n_run = 5, k_folds = 10,
                          x.center = TRUE, x.scale = FALSE,
                          y.center = FALSE, y.scale = FALSE,
                          remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                          remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                          w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                          MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                          pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                          MIN_EPV = 5, return_models = F, returnData = F, tol = 1e-10,
                          PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  #### ### ###
  # WARNINGS #
  #### ### ###

  #### Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### Check values classes and ranges
  lst_01 <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                 "w_AIC" = w_AIC, "w_c.index" = w_c.index, "w_AUC" = w_AUC, "w_BRIER" = w_BRIER)
  check_min0_max1_variables(lst_01)

  lst_num <- list("max.ncomp" = max.ncomp, "spv_penalty.list" = spv_penalty.list,
                  "n_run" = n_run, "k_folds" = k_folds, "max_time_points" = max_time_points,
                  "MIN_COMP_TO_CHECK" = MIN_COMP_TO_CHECK, "MIN_EPV" = MIN_EPV, "seed" = seed, "tol" = tol)
  check_class(lst_num, class = "numeric")

  lst_logical <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                      "y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_variance_at_fold_level" = remove_variance_at_fold_level,
                      "remove_non_significant_models" = remove_non_significant_models,
                      "remove_non_significant" = remove_non_significant,
                      "return_models" = return_models,"returnData" = returnData, "verbose" = verbose, "PARALLEL" = PARALLEL)
  check_class(lst_logical, class = "logical")

  lst_character <- list("pred.attr" = pred.attr, "pred.method" = pred.method)
  check_class(lst_character, class = "character")

  #### Check cv-folds
  n_run <- checkFoldRuns(Y, n_run, k_folds)

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars.mb(X)

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_BRIER, w_AUC))
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### ZERO VARIANCE - ALWAYS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
                                               remove_near_zero_variance = remove_near_zero_variance,
                                               remove_zero_variance = remove_zero_variance,
                                               toKeep.zv = toKeep.zv,
                                               freqCut = 95/5)
    X <- lst_dnz$X
    variablesDeleted <- lst_dnz$variablesDeleted
  }else{
    variablesDeleted <- NULL
  }

  #### MAX PREDICTORS
  max.ncomp <- check.mb.ncomp(X, max.ncomp)
  max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)
  if(MIN_COMP_TO_CHECK >= max.ncomp){
    MIN_COMP_TO_CHECK = max.ncomp-1
  }

  #### #
  # CV #
  #### #
  set.seed(seed)
  lst_data <- splitData_Iterations_Folds.mb(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test
  k_folds <- lst_data$k_folds

  lst_train_indexes <- lst_data$lst_train_index
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  #total_models <- 1 * k_folds * n_run
  total_models <- max.ncomp * k_folds * n_run * length(spv_penalty.list)

  lst_model <- get_HDCOX_models2.0(method = pkg.env$sb.splsicox,
                                   lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                   max.ncomp = max.ncomp, eta.list = spv_penalty.list, EN.alpha.list = NULL,
                                   n_run = n_run, k_folds = k_folds,
                                   remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                   remove_non_significant = remove_non_significant, alpha = alpha, MIN_EPV = MIN_EPV, tol = tol,
                                   x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
                                   total_models = total_models, PARALLEL = PARALLEL, verbose = verbose)

  # lst_model <- get_HDCOX_models(method = pkg.env$sb.splsicox,
  #                               lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
  #                               max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL,
  #                               n_run = n_run, k_folds = k_folds,
  #                               x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
  #                               total_models = total_models)

  if(all(is.na(unlist(lst_model)))){
    message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.sb.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = lst_model, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.sb.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class= pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = lst_model, alpha = alpha,
                                                    max.ncomp = max.ncomp, eta.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models, verbose = verbose)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.sb.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = lst_model, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.sb.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class= pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
  optimal_eta_index <- NULL
  optimal_eta <- NULL

  if(TRUE){ #compute always BRIER SCORE
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_BRIER_sPLS(comp_model_lst = lst_model,
                                            lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                            df_results_evals = df_results_evals, times = times,
                                            pred.method = pred.method, pred.attr = pred.attr,
                                            max.ncomp = max.ncomp, eta.list = spv_penalty.list, n_run = n_run, k_folds = k_folds,
                                            MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                            w_BRIER = w_BRIER, method.train = pkg.env$sb.splsicox, PARALLEL = F, verbose = verbose)

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

    lst_df <- get_COX_evaluation_AUC_sPLS(comp_model_lst = lst_model,
                                          lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                          df_results_evals = df_results_evals, times = times,
                                          fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                          max.ncomp = max.ncomp, eta.list = spv_penalty.list, n_run = n_run, k_folds = k_folds,
                                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          w_AUC = w_AUC, method.train = pkg.env$sb.splsicox, PARALLEL = F, verbose = verbose)

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
    optimal_eta <- lst_df$optimal_eta
    optimal_eta_index <- lst_df$optimal_eta_index
  }

  #### ### ### #
  # BEST MODEL #
  #### ### ### #

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_BRIER, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  #### ###
  # PLOT #
  #### ###
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_BRIER = w_BRIER, max.ncomp = max.ncomp, eta.list = NULL,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER", x.text = "Component")

  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_BRIER <- lst_EVAL_PLOTS$ggp_BRIER
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp

  #### ### #
  # RETURN #
  #### ### #

  message(paste0("Best model obtained.\n"))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.sb.splsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = lst_model, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.spv_penalty = best_model_info$eta, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.sb.splsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.spv_penalty = best_model_info$eta, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.sb.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }
}

#' Cross validation fast.cv.sb.splsicox
#' @description fast.cv.sb.splsicox cross validation model
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param spv_penalty.list Numeric vector. Penalty for variable selection for the individual cox models. Variables with a lower P-Value than "spv_penalty" in the individual cox analysis will be keep for the sPLS-ICOX approach (default: seq(0.1,1,0.1)).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_c.index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to continue evaluating higher values in the multiple tested parameters. If it is not reached for next 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops (default: 0.01).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test observations. Once all have their linear predictors, the evaluation is perform across all the observations together (default: FALSE).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "HDcox" and model "sb.splsicox". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: PLS weights
#'  \item \code{(weightings_norm)}: PLS normalize weights
#'  \item \code{(W.star)}: PLS W* vector
#'  \item \code{(scores)}: PLS scores/variates
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
#'  \item \code{(dr.mean)}: mean values for deviance residuals Y matrix
#'  \item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix'
#'  }
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
#' \code{list_pls_models}: List of sPLS-ICOX models computed for each block.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
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
#' sb.splsicox_model <- fast.cv.sb.splsicox(X, Y, max.ncomp = 10,
#' spv_penalty.list = seq(0.1,1,0.1), x.center = TRUE, x.scale = TRUE)
#' }

fast.cv.sb.splsicox <- function(X, Y,
                               max.ncomp = 10, spv_penalty.list = seq(0.1,1,0.1),
                               n_run = 5, k_folds = 10,
                               x.center = TRUE, x.scale = FALSE,
                               y.center = FALSE, y.scale = FALSE,
                               remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                               remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                               w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                               MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                               pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                               MIN_EPV = 5, returnData = T, return_models = F, tol = 1e-10,
                               PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  #### Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### Check values classes and ranges
  lst_01 <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                 "w_AIC" = w_AIC, "w_c.index" = w_c.index, "w_AUC" = w_AUC, "w_BRIER" = w_BRIER)
  check_min0_max1_variables(lst_01)

  lst_num <- list("max.ncomp" = max.ncomp, "spv_penalty.list" = spv_penalty.list,
                  "n_run" = n_run, "k_folds" = k_folds, "max_time_points" = max_time_points,
                  "MIN_COMP_TO_CHECK" = MIN_COMP_TO_CHECK, "MIN_EPV" = MIN_EPV, "seed" = seed, "tol" = tol)
  check_class(lst_num, class = "numeric")

  lst_logical <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                      "y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_variance_at_fold_level" = remove_variance_at_fold_level,
                      "remove_non_significant_models" = remove_non_significant_models,
                      "remove_non_significant" = remove_non_significant,
                      "return_models" = return_models,"returnData" = returnData, "verbose" = verbose, "PARALLEL" = PARALLEL)
  check_class(lst_logical, class = "logical")

  lst_character <- list("pred.attr" = pred.attr, "pred.method" = pred.method)
  check_class(lst_character, class = "character")

  #### Check cv-folds
  n_run <- checkFoldRuns(Y, n_run, k_folds)

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars.mb(X)

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_BRIER, w_AUC))
  max.ncomp <- check.mb.ncomp(X, max.ncomp)

  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### REQUIREMENTS
  lst_check <- checkXY.mb.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  checkY.colnames(Y)

  #### SCALE
  if(length(x.center)==1){
    x.center <- rep(x.center, length(names(X)))
    names(x.center) <- names(X)
  }
  if(length(x.scale)==1){
    x.scale <- rep(x.scale, length(names(X)))
    names(x.scale) <- names(X)
  }

  #### ZERO VARIANCE - ALWAYS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
                                               remove_near_zero_variance = remove_near_zero_variance,
                                               remove_zero_variance = remove_zero_variance,
                                               toKeep.zv = toKeep.zv,
                                               freqCut = 95/5)
    X <- lst_dnz$X
    variablesDeleted <- lst_dnz$variablesDeleted
  }else{
    variablesDeleted <- NULL
  }

  #### SCALING
  lst_scale <- XY.mb.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  #### MAX PREDICTORS
  max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)
  if(MIN_COMP_TO_CHECK >= max.ncomp){
    MIN_COMP_TO_CHECK = max.ncomp-1
  }

  # CREATE INDIVIDUAL MODELS
  lst_sb.pls <- list()
  for(b in names(Xh)){

    message(paste0("\nRunning cross validation ", pkg.env$sb.splsicox, " for block: ", b, "\n"))

    cv.splsdrcox_res <- cv.splsicox(X = Xh[[b]], Y = Yh,
                                   max.ncomp = max.ncomp,
                                   n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                   w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                   MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                   x.scale = x.scale[[b]], x.center = x.center[[b]], y.scale = y.scale, y.center = y.center,
                                   remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                   remove_non_significant = remove_non_significant,
                                   fast_mode = fast_mode, return_models = return_models, tol = tol,
                                   MIN_EPV = MIN_EPV, verbose = verbose,
                                   pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

    lst_sb.pls[[b]] <- splsicox(X = Xh[[b]],
                               Y = Yh,
                               n.comp = cv.splsdrcox_res$opt.comp,
                               remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                               remove_non_significant = remove_non_significant, alpha = alpha, tol = tol,
                               returnData = F,
                               x.center = x.center[[b]], x.scale = x.scale[[b]],
                               y.scale = y.scale, y.center = y.center,
                               MIN_EPV = MIN_EPV, verbose = verbose)
  }

  # CHECK ALL MODELS SAME COMPONENTS
  aux_ncomp <- purrr::map(lst_sb.pls, ~.$n.comp)

  # CREATE COMBINE MODEL
  data <- NULL
  cn.merge <- NULL
  for(b in names(Xh)){
    if(!is.null(lst_sb.pls[[b]]$survival_model)){
      data <- cbind(data, lst_sb.pls[[b]]$X$scores)
      cn.merge <- c(cn.merge, paste0(colnames(lst_sb.pls[[b]]$X$scores), "_", b))
    }else{
      next
    }
  }

  #colnames(data) <- apply(expand.grid(colnames(lst_sb.pls[[1]]$X$scores), names(Xh)), 1, paste, collapse="_")
  colnames(data) <- cn.merge
  cox_model <- cox(X = data, Y = Yh,
                   x.center = F, x.scale = F,
                   y.center = F, y.scale = F,
                   remove_non_significant = remove_non_significant,
                   FORCE = T)

  # RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  # already performed in cox() function
  # removed_variables <- NULL
  # if(remove_non_significant){
  #   if(all(c("time", "event") %in% colnames(data))){
  #     lst_rnsc <- removeNonSignificativeCox(cox = cox_model$survival_model$fit, alpha = alpha, cox_input = data, time.value = NULL, event.value = NULL)
  #   }else{
  #     lst_rnsc <- removeNonSignificativeCox(cox = cox_model$survival_model$fit, alpha = alpha, cox_input = cbind(data, Yh), time.value = NULL, event.value = NULL)
  #   }
  #
  #   cox_model$survival_model$fit <- lst_rnsc$cox
  #   removed_variables <- lst_rnsc$removed_variables
  # }

  if(remove_non_significant){
    removed_variables <- cox_model$nsv
  }else{
    removed_variables <- NULL
  }

  #### ### #
  # RETURN #
  #### ### #
  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(sb.splsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "x.mean" = xmeans, "x.sd" = xsds),
                                Y = list("data" = lst_sb.pls[[1]]$Y$data, "y.mean" = ymeans, "y.sd" = ysds),
                                survival_model = cox_model$survival_model,
                                list_pls_models = lst_sb.pls,
                                n.comp = aux_ncomp, #number of components used, but could be lesser than expected because not computed models
                                call = func_call,
                                X_input = if(returnData) X_original else NA,
                                Y_input = if(returnData) Y_original else NA,
                                alpha = alpha,
                                removed_variables_cox = removed_variables,
                                class = pkg.env$sb.splsicox,
                                time = time)))
}

### ## ##
# CLASS #
### ## ##

sb.splsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$sb.splsicox)
  return(model)
}

cv.sb.splsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.sb.splsicox)
  return(model)
}
