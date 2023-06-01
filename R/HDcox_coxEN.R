#### ### ##
# METHODS #
#### ### ##

#' coxEN
#' @description This function performs a cox elastic net model (based on glmnet R package).
#' The function returns a HDcox model with the attribute model as "coxEN".
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param EN.alpha Numeric. Elastic net mixing parameter. If EN.alpha = 1 is the lasso penalty, and EN.alpha = 0 the ridge penalty (default: 0.5).
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If MIN_EPV is not meet, the value will be change automatically (default: 20).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables/components in final cox model will be removed until all variables are significant by forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "coxEN". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
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
#' \code{opt.lambda}: Optimal lambda computed by the model with maximum % Var from glmnet function.
#'
#' \code{EN.alpha}: EN.alpha selected
#'
#' \code{n.var}: Number of variables selected
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{convergence_issue}: If any convergence issue has been found.
#'
#' \code{alpha}: alpha value selected
#'
#' \code{selected_variables_cox}: Variables selected to enter the cox model.
#'
#' \code{removed_variables_cox}: Variables removed by EN penalty.
#'
#' \code{removed_variables_correlation}: Variables removed by being high correlated with other variables.
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
#' coxEN(X, Y)
#' coxEN(X, Y, EN.alpha = 0.75, x.center = TRUE, x.scale = TRUE)
#' }

coxEN <- function(X, Y,
                  EN.alpha = 0.5, max.variables = 15,
                  x.center = TRUE, x.scale = FALSE,
                  remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                  remove_non_significant = F, alpha = 0.05,
                  MIN_EPV = 5, returnData = T, verbose = F){

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### Check values classes and ranges
  params_with_limits <- list("alpha" = alpha, "EN.alpha" = EN.alpha)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("max.variables" = max.variables,
                  "MIN_EPV" = MIN_EPV)
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
  max.variables <- check.ncomp(X, max.variables)
  max.variables <- check.maxPredictors(X, Y, MIN_EPV, max.variables, verbose = verbose)

  #### REMOVING event in time 0 patients
  if(any(Yh[,"time"]==0)){
    if(verbose){
      message("Some patients get the event at time 0. Those patients will be removed for the analysis")
    }
    pat_names <- rownames(Yh)[Yh[,"time"]==0]
    Yh <- Yh[!rownames(Yh) %in% pat_names,]
    Xh <- Xh[rownames(Xh) %in% rownames(Yh),]
  }

  #### INITIALISING VARIABLES
  best_cox <- NULL
  best_lambda <- NULL
  selected_variables <- NULL
  problem <- FALSE #convergence issues

  #A small detail in the Cox model: if death times are tied with censored times, we assume the censored times occurred
  #just before the death times in computing the Breslow approximation; if users prefer the usual convention of after,
  #they can add a small number to all censoring times to achieve this effect.
  for(t in unique(Yh[,"time"])){
    set_pat <- Yh[Yh[,"time"]==t,,drop=F]
    if(nrow(set_pat)>1 & length(unique(set_pat[,"event"]))>1){
      names_censored <- names(which(set_pat[,"event"]==0))
      Yh[names_censored,"time"] <- Yh[names_censored,"time"] + 0.0001
    }
  }

  # BEST lambda
  # I cannot add the limit for maximum number of variables bc it fails
  # cvfit <- cv.glmnet(x = Xh, y = survival::Surv(time = Yh[,"time"], event = Yh[,"event"]),
  #                    family = "cox", type.measure = "C",
  #                    alpha = EN.alpha, dfmax = max.variables,
  #                    standardize = F, nlambda=200)

  #I cannot add the limit for maximum number of variables bc it fails
  #pmax = max.variables
  EN_cox <- tryCatch(
    # Specifying expression
    # pmax - coeffitients to be non-zero
    expr = {
      glmnet::glmnet(x = Xh, y = survival::Surv(time = Yh[,"time"], event = Yh[,"event"]),
                     family = "cox", alpha = EN.alpha, standardize = T, nlambda = 300, pmax = max.variables)
    },
    # Specifying error message
    error = function(e){
      message(paste0("coxEN: ", e))
      # invisible(gc())
      return(NA)
    },
    warning = function(e){
      if(verbose){
        message("Model probably has a convergence issue...\n")
      }
      suppressWarnings(
        # dfmax - variables in the model
        res <- glmnet::glmnet(x = Xh, y = survival::Surv(time = Yh[,"time"], event = Yh[,"event"]),
                              family = "cox", alpha = EN.alpha, standardize = F, nlambda = 300, dfmax = max.variables)
      )
      list(res = res, problem = T)
    }
  )

  #only if problems
  if(all(isa(EN_cox, "list"))){
    problem = EN_cox$problem
    EN_cox = EN_cox$res
  }

  if(all(class(EN_cox) %in% c("coxnet", "glmnet"))){
    best_lambda <- EN_cox$lambda[which.max(EN_cox$dev.ratio)]

    if(best_lambda!=Inf){
      coef.matrix <- as.matrix(coef(EN_cox, s = best_lambda))
      selected_variables <- rownames(coef.matrix)[which(coef.matrix != 0)]

      d <- as.data.frame(cbind(Xh[,selected_variables,drop=F], Yh)) #data

      best_cox <- tryCatch(
        # Specifying expression
        expr = {
          survival::coxph(formula = survival::Surv(time,event) ~ .,
                          data = d,
                          ties = "efron",
                          singular.ok = T,
                          robust = T,
                          nocenter = rep(1, ncol(Xh)),
                          model=T, x = T)
        },
        # Specifying error message
        error = function(e){
          message(paste0("COX: ", e))
          # invisible(gc())
          return(NA)
        }
      )
    }else{
      best_cox = NA
      best_lambda = NA
      selected_variables = NA

      func_call <- match.call()

      t2 <- Sys.time()
      time <- difftime(t2,t1,units = "mins")

      survival_model <- NULL

      return(coxEN_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                              Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                              survival_model = survival_model,
                              EN.alpha = EN.alpha,
                              n.var = max.variables,
                              #alpha = alpha,
                              call = func_call,
                              X_input = if(returnData) X_original else NA,
                              Y_input = if(returnData) Y_original else NA,
                              nzv = variablesDeleted,
                              selected_variables = selected_variables,
                              removed_variables = removed_variables,
                              removed_variables_correlation = removed_variables_cor,
                              opt.lambda = best_lambda,
                              convergence_issue = problem,
                              class = pkg.env$coxEN,
                              time = time)))
    }

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

  if(all(is.na(best_cox)) || (problem & all(best_cox$linear.predictors==0))){
    best_cox = NA
    best_lambda = NA
    selected_variables = NA

    func_call <- match.call()

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")

    survival_model <- NULL

    # invisible(gc())
    return(coxEN_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                            Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                            survival_model = survival_model,
                            EN.alpha = EN.alpha,
                            n.var = max.variables,
                            #alpha = alpha,
                            call = func_call,
                            X_input = if(returnData) X_original else NA,
                            Y_input = if(returnData) Y_original else NA,
                            nzv = variablesDeleted,
                            selected_variables = selected_variables,
                            removed_variables = removed_variables,
                            removed_variables_correlation = removed_variables_cor,
                            opt.lambda = best_lambda,
                            convergence_issue = problem,
                            class = pkg.env$coxEN,
                            time = time)))
  }

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(d))){
      lst_rnsc <- removeNonSignificativeCox(cox = best_cox, alpha = alpha, cox_input = d, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = best_cox, alpha = alpha, cox_input = cbind(d, Yh), time.value = NULL, event.value = NULL)
    }

    best_cox <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
    selected_variables <- selected_variables[!selected_variables %in% lst_rnsc$removed_variables]
  }

  if(class(best_cox)[[1]] %in% "coxph.null"){
    survival_model <- NULL
  }else if(isa(best_cox,"coxph")){
    survival_model <- getInfoCoxModel(best_cox)
  }else{
    survival_model <- NULL
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(coxEN_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                          Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                          survival_model = survival_model,
                          opt.lambda = best_lambda,
                          EN.alpha = EN.alpha,
                          n.var = max.variables,
                          call = func_call,
                          X_input = if(returnData) X_original else NA,
                          Y_input = if(returnData) Y_original else NA,
                          convergence_issue = problem,
                          alpha = alpha,
                          selected_variables_cox = selected_variables,
                          removed_variables_cox = removed_variables,
                          nzv = variablesDeleted,
                          class = pkg.env$coxEN,
                          time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' coxEN Cross-Validation
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param EN.alpha.list Numeric vector. Elastic net mixing parameter values to test in cross validation. EN.alpha = 1 is the lasso penalty, and EN.alpha = 0 the ridge penalty (default: seq(0,1,0.1)).
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If MIN_EPV is not meet, the value will be change automatically (default: 20).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
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
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "HDcox" and model "cv.coxEN". The class contains the following elements:
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha information).
#'
#' \code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
#' \code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.
#'
#' \code{opt.EN.alpha}: Optimal EN.alpha value selected by the best_model.
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
#' cv.coxEN_model <- cv.coxEN(X, Y, EN.alpha.list = c(0.1,0.5,0.1), x.center = TRUE, x.scale = TRUE)
#' coxEN_model <- coxEN(X, Y, EN.alpha = cv.coxEN_model$opt.EN.alpha, x.center = TRUE, x.scale = TRUE)
#' }

cv.coxEN <- function(X, Y,
                     EN.alpha.list = seq(0,1,0.1),
                     max.variables = 15,
                     n_run = 5, k_folds = 10,
                     x.center = TRUE, x.scale = FALSE,
                     remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                     remove_non_significant = F, alpha = 0.05,
                     w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                     MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                     pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                     MIN_EPV = 5, return_models = F, returnData = F,
                     PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()
  FREQ_CUT <- 95/5
  y.center = y.scale = FALSE

  #### ### ###
  # WARNINGS #
  #### ### ###

  #### Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### Check values classes and ranges
  params_with_limits <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                        "w_AIC" = w_AIC, "w_c.index" = w_c.index, "w_AUC" = w_AUC, "w_BRIER" = w_BRIER)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("EN.alpha.list" = EN.alpha.list, "max.variables" = max.variables,
                  "n_run" = n_run, "k_folds" = k_folds, "max_time_points" = max_time_points,
                  "MIN_COMP_TO_CHECK" = MIN_COMP_TO_CHECK, "MIN_EPV" = MIN_EPV, "seed" = seed)
  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_variance_at_fold_level" = remove_variance_at_fold_level,
                      "remove_non_significant" = remove_non_significant,
                      "return_models" = return_models,"returnData" = returnData, "verbose" = verbose, "PARALLEL" = PARALLEL)
  check_class(logical_params, class = "logical")

  character_params <- list("pred.attr" = pred.attr, "pred.method" = pred.method)
  check_class(character_params, class = "character")

  #### Check cv-folds
  lst_checkFR <- checkFoldRuns(Y, n_run, k_folds, fast_mode)
  n_run <- lst_checkFR$n_run
  fast_mode <- lst_checkFR$fast_mode

  #### Check rownames
  lst_check <- checkXY.rownames(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars(X)

  #### REQUIREMENTS
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  check.cv.weights(c(w_AIC, w_c.index, w_BRIER, w_AUC))

  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### MAX PREDICTORS
  max.variables <- check.ncomp(X, max.variables)
  max.variables <- check.maxPredictors(X, Y, MIN_EPV, max.variables, verbose = verbose)

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
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  total_models <- k_folds * n_run * length(EN.alpha.list)
  comp_model_lst <- get_HDCOX_models2.0(method = pkg.env$coxEN,
                                        lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                        max.ncomp = NULL, eta.list = NULL, EN.alpha.list = EN.alpha.list, max.variables = max.variables, vector = NULL,
                                        n_run = n_run, k_folds = k_folds,
                                        MIN_NVAR = NULL, MAX_NVAR = NULL, MIN_AUC_INCREASE = NULL, EVAL_METHOD = NULL,
                                        n.cut_points = NULL,
                                        x.center = x.center, x.scale = x.scale,
                                        y.center = y.center, y.scale = y.scale,
                                        remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                        alpha = alpha, MIN_EPV = MIN_EPV,
                                        remove_non_significant = remove_non_significant, tol = NULL, max.iter = NULL,
                                        returnData = returnData, total_models = total_models,
                                        PARALLEL = PARALLEL, verbose = verbose)

  count_problems = 0
  count_null = 0
  count = 0
  for(l in names(comp_model_lst)){
    for(r in names(comp_model_lst[[l]])){
      for(f in names(comp_model_lst[[l]][[r]])){
        if(is.null(comp_model_lst[[l]][[r]][[f]]$convergence_issue)){
          count_null = count_null + 1
        }else if(comp_model_lst[[l]][[r]][[f]]$convergence_issue){
          count_problems = count_problems + 1
        }else{
          count = count + 1
        }
      }
    }
  }

  if(count_problems>0){
    message(paste0("There are a total of ", count_problems, " out of ", total_models, " models that could present convergence issues.\n"))
  }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst, alpha = alpha,
                                                    max.ncomp = EN.alpha.list, eta.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = F, verbose = verbose) #deletion at variable level for EN

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.coxEN_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.EN.alpha = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.coxEN, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.coxEN_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.EN.alpha = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.coxEN, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
                                       max.ncomp = EN.alpha.list, n_run = n_run, k_folds = k_folds,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       w_BRIER = w_BRIER, method.train = pkg.env$coxEN, PARALLEL = F, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    #total_models <- ifelse(!fast_mode, n_run * length(EN.alpha.list), k_folds * n_run * length(EN.alpha.list))#inside get_COX_evaluation_AUC

    #times should be the same for all folds
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = EN.alpha.list, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, method.train = pkg.env$coxEN, PARALLEL = F, verbose = verbose)

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
  df_results_evals_comp <- cv.getScoreFromWeight(lst_cox_mean = df_results_evals_comp, w_AIC, w_c.index, w_BRIER, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER")

  flag_no_models = F
  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==EN.alpha.list[[optimal_comp_index]],, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    if(all(is.nan(df_results_evals_comp[,"score"]))){
      #message("No models computed. All of them have problems in convergency. Probably due to a high number of variables.")
      best_model_info <- df_results_evals_comp[1,,drop=F]
      flag_no_models = T
    }else{
      best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
      best_model_info <- as.data.frame(best_model_info)
    }
  }

  #### ###
  # PLOT #
  #### ###
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_BRIER = w_BRIER, max.ncomp = EN.alpha.list, eta.list = NULL,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER", x.text = "EN.alpha Penalization")

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp
  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_BRIER <- lst_EVAL_PLOTS$ggp_BRIER
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  #### ### ### ### ### ### #
  # CHANGE 1s COLUMN_NAME  #
  #### ### ### ### ### ### #
  colnames(best_model_info)[1] <- "EN.alpha"
  colnames(df_results_evals_fold)[1] <- "EN.alpha"
  colnames(df_results_evals_run)[1] <- "EN.alpha"
  colnames(df_results_evals_comp)[1] <- "EN.alpha"

  ### ## ###
  # RETURN #
  ### ## ###

  if(!flag_no_models){
    message(paste0("Best model obtained."))
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  if(return_models){
    return(cv.coxEN_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.EN.alpha = best_model_info$EN.alpha, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.coxEN, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.coxEN_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.EN.alpha = best_model_info$EN.alpha, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.coxEN, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }
}

### ## ##
# CLASS #
### ## ##

coxEN_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$coxEN)
  return(model)
}

cv.coxEN_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$cv.coxEN)
  return(model)
}
