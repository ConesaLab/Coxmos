#### ### ##
# METHODS #
#### ### ##

#' sPLS-DRCOX
#' @description Performs a sPLS-DRCOX_dynamic model (mixing sPLS-DRCOX and mixOmics custom variable selection).
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
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
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
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "sPLS-DRCOX-Dynamic". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: PLS weights
#'  \item \code{(weightings_norm)}: PLS normalize weights
#'  \item \code{(W.star)}: PLS W* vector
#'  \item \code{(scores)}: PLS scores/variates
#'  \item \code{(E)}: error matrices
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
#' \code{eta}: Penalty value selected.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{n.varX}: Number of Variables selected in each PLS component.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{beta_matrix}: PLS beta matrix
#'
#' \code{R2}: PLS R2
#'
#' \code{SCR}: PLS SCR
#'
#' \code{SCT}: PLS SCT
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export

splsdrcox_dynamic <- function (X, Y,
                                n.comp = 4, vector = NULL,
                                MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                MIN_AUC_INCREASE = 0.01,
                                x.center = TRUE, x.scale = FALSE,
                                y.center = FALSE, y.scale = FALSE,
                                remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                                remove_non_significant = F, alpha = 0.05, tol = 1e-15,
                                EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                                times = NULL, max_time_points = 15,
                                MIN_EPV = 5, returnData = T, PARALLEL = F, verbose = F){

  t1 <- Sys.time()

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### REQUIREMENTS
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  checkY.colnames(Y)

  # if(k_folds.mixOmics<2){
  #   message("MixOmics k_folds must be 2 as minimum.")
  #   k_folds.mixOmics = 2
  # }

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = 95/5)
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

  E <- list()
  R2 <- list()
  SCR <- list()
  SCT <- list()

  XXNA <- is.na(Xh) #T is NA
  YNA <- is.na(Y) #T is NA

  #### ### ### ### ### ### ### ### ### ### ### ###
  ### ###             sPLS-COX             ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###

  #2. Surv function - NULL model
  coxDR <- survival::coxph(survival::Surv(time = time, event = event, type = "right") ~ 1, as.data.frame(Xh))

  #3. Residuals - Default is deviance because eval type="deviance"
  DR_coxph <- residuals(coxDR, type = "deviance") #"martingale", "deviance", "score", "schoenfeld", "dfbeta"', "dfbetas", "scaledsch" and "partial"

  #### ### ### ### ### ### ### ### ### ### ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  #### ### ### ### ### ### ### ### ### ### ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###

  #4. sPLS Algorithm
  n_obs <- nrow(Xh)
  n_var <- ncol(Xh)
  n_dr <- ncol(DR_coxph)

  if(is.null(n_dr))
    n_dr=1

  #Norm Y
  mu <- mean(DR_coxph) #equivalent because Y it is not normalized
  DR_coxph <- scale(DR_coxph, center = mu, scale = FALSE) #center DR to DR / patients
  DR_coxph_ori <- DR_coxph

  #Norm X - Do not center again
  # x_center_by_nobs <- colMeans(Xh, na.rm = T)
  # Xh <- scale(Xh, center = x_center_by_nobs, scale = FALSE)
  flag = T
  cv.spls <- NA

  #### ### ### ### ### ### ### ### ### ###
  # DIVIDE Y VENCERAS - BEST VECTOR SIZE #
  #### ### ### ### ### ### ### ### ### ###

  if(is.null(vector)){
    lst_BV <- getBestVector(Xh, DR_coxph, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                           EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F, mode = "spls", times = times, max_time_points = max_time_points, verbose = verbose)
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
                                 EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F, mode = "spls", times = times, max_time_points = max_time_points, verbose = verbose)
        keepX <- lst_BV$best.keepX
        plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
    }
  }

  #### ### ### ### ### ### ### ### ### ### ###
  ### ##             sPLS              ###  ##
  #### ### ### ### ### ### ### ### ### ### ###

  spls <- mixOmics::spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = rep(keepX, n.comp), scale = F)

  # PREDICTION

  # both functions work fine (predict and predict.mixOmixs.pls)
  # predplsfit <- predict.mixOmixs.pls(spls, newdata=Xh[,rownames(spls$loadings$X),drop=F])
  # predplsfit <- predict(spls, newdata=Xh[,rownames(spls$loadings$X),drop=F])

  # sometimes solve(t(P) %*% W) in predict can cause an error
  # system is computationally singular: reciprocal condition number = 6.24697e-18
  predplsfit <- tryCatch(expr = {predict(spls, newdata=Xh[,rownames(spls$loadings$X),drop=F])},
                 error = function(e){
                   if(verbose){
                     message(e$message)
                   }
                   NA
                 })

  if(!all(is.na(predplsfit))){
    # R2 calculation
    for(h in 1:n.comp){
      E[[h]] <- DR_coxph_ori - predplsfit$predict[,,h]
      SCR[[h]] = sum(apply(E[[h]],2,function(x) sum(x**2)))
      SCT[[h]] = sum(apply(as.matrix(DR_coxph_ori),2,function(x) sum(x**2))) #equivalent sum((DR_coxph_ori - mean(DR_coxph_ori))**2)
      R2[[h]] = 1 - (SCR[[h]]/SCT[[h]]) #deviance residuals explanation
    }
  }else{
    E <- NULL
    SCR <- NULL
    SCT <- NULL
    R2 <- NULL
  }

  #last model includes all of them
  tt_splsDR = spls$variates$X
  ww_splsDR = spls$loadings$X
  rr_splsDR = spls$loadings.star
  pp_splsDR = spls$mat.c

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ## ### ### ### ### #
  ### ##              PLS-COX            ### ##
  #### ### ### ### ### ### ## ### ### ### ### #
  n.comp_used <- ncol(tt_splsDR) #can be lesser than expected because we have lesser variables to select because penalization
  n.varX_used <- keepX

  d <- as.data.frame(tt_splsDR[,,drop=F])
  rownames(d) <- rownames(X)
  colnames(d) <- paste0("comp_", 1:n.comp_used)
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
      message(paste0("splsdrcox_dynamic: ",conditionMessage(e)))
      invisible(gc())
      return(NA)
    }
  )

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  removed_variables <- NULL
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(d))){
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$fit, alpha = alpha, cox_input = d, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$fit, alpha = alpha, cox_input = cbind(d, Yh), time.value = NULL, event.value = NULL)
    }

    cox_model$fit <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  survival_model <- NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  #get W.star
  W <- ww_splsDR
  P <- pp_splsDR

  if(is.null(P) | is.null(W)){
    message(paste0(pkg.env$splsdrcox_dynamic," model cannot be computed because P or W vectors are NULL. Returning NA."))
    invisible(gc())
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
    message(paste0(pkg.env$splsdrcox_dynamic," model cannot be computed due to solve(t(P) %*% W). Reduce 'tol' parameter to fix it. Returning NA."))
    invisible(gc())
    return(NA)
  }

  W.star <- W %*% PW
  Ts <- tt_splsDR

  func_call <- match.call()

  rownames(Ts) <- rownames(X)
  rownames(P) <- rownames(W) <-  rownames(W.star) <- rownames(ww_splsDR)

  colnames(Ts) <- colnames(P) <- colnames(W) <-  colnames(W.star)<- paste0("comp_", 1:n.comp_used)

  # variable per component
  var_by_component = list()
  for(cn in colnames(W)){
    var_by_component[[cn]] <- rownames(W)[W[,cn]!=0]
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(splsdrcox_dynamic_class(list(X = list("data" = if(returnData) X_norm else NA, "loadings" = P, "weightings" = W, "W.star" = W.star, "scores" = Ts, "E" = E, "x.mean" = xmeans, "x.sd" = xsds),
                                      Y = list("deviance_residuals" = if(returnData) DR_coxph_ori else NA, "dr.mean" = NULL, "dr.sd" = NULL, #deviance_residuals object already centered
                                                "data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                                      survival_model = survival_model,
                                      n.comp = n.comp_used, #number of components
                                      n.varX = n.varX_used,
                                      var_by_component = var_by_component, #variables selected for each component
                                      plot_accuracyPerVariable = plotVAR,
                                      call = func_call,
                                      X_input = if(returnData) X_original else NA,
                                      Y_input = if(returnData) Y_original else NA,
                                      R2 = R2,
                                      SCR = SCR,
                                      SCT = SCT,
                                      alpha = alpha,
                                      removed_variables_cox = removed_variables,
                                      nzv = variablesDeleted,
                                      class = pkg.env$splsdrcox_dynamic,
                                      time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation sPLS-DRCOX
#' @description sPLS-DRCOX cross validation model
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL).
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
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation. A non-significant model is a model with at least one component/variable with a P-Value higher than the alpha cutoff.
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
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "HDcox" and model "cv.sPLS-DRCOX-Dynamic".
#' @export

cv.splsdrcox_dynamic <- function (X, Y,
                                  max.ncomp = 10, n_run = 5, k_folds = 10,
                                  vector = NULL,
                                  x.center = TRUE, x.scale = FALSE,
                                  y.center = FALSE, y.scale = FALSE,
                                  remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                                  remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                                  MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                  MIN_AUC_INCREASE = 0.01,
                                  EVAL_METHOD = "AUC",
                                  w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                                  MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                  pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                                  MIN_EPV = 5, return_models = F, returnData = F, tol = 1e-15,
                                  PARALLEL = F, verbose = F, seed = 123){
  t1 <- Sys.time()

  #### ### ###
  # WARNINGS #
  #### ### ###

  #Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #Illegal chars in colnames
  X <- checkColnamesIllegalChars(X)

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_BRIER, w_AUC))
  max.ncomp <- check.ncomp(X, max.ncomp)

  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### MAX PREDICTORS
  max.ncomp <- check.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)

  #### REQUIREMENTS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                            remove_near_zero_variance = remove_near_zero_variance,
                                            remove_zero_variance = remove_zero_variance,
                                            toKeep.zv = toKeep.zv,
                                            freqCut = 95/5)
    X <- lst_dnz$X
    variablesDeleted <- lst_dnz$variablesDeleted
  }else{
    variablesDeleted <- NULL
  }

  #### #
  # CV #
  #### #
  set.seed(seed)
  lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test

  lst_train_indexes <- lst_data$lst_train_index
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  total_models <- 1 * k_folds * n_run

  comp_model_lst  <- get_HDCOX_models2.0(method = pkg.env$splsdrcox_dynamic,
                                         lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                         max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL, n_run = n_run, k_folds = k_folds,
                                         vector = vector,
                                         MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                         MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                         EVAL_METHOD = EVAL_METHOD,
                                         x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
                                         remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                         remove_non_significant = remove_non_significant, tol = tol,
                                         total_models = total_models, PARALLEL = PARALLEL, verbose = verbose)

  if(all(is.null(comp_model_lst))){
    message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdrcox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdrcox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
      return(cv.splsdrcox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdrcox_dynamic_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
    lst_df <- get_COX_evaluation_BRIER(comp_model_lst = comp_model_lst,
                                       lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                       df_results_evals = df_results_evals, times = times,
                                       pred.method = pred.method, pred.attr = pred.attr,
                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       w_BRIER = w_BRIER, method.train = pkg.env$splsdrcox_dynamic, PARALLEL = F, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    #total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp) #inside get_COX_evaluation_AUC

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
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, method.train = pkg.env$splsdrcox_dynamic, PARALLEL = F, verbose = verbose)

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
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  #### ###
  # PLOT #
  #### ###
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

  invisible(gc())
  if(return_models){
    return(cv.splsdrcox_dynamic_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.splsdrcox_dynamic_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_model_info$n.var, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsdrcox_dynamic, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }
}

### ## ##
# CLASS #
### ## ##

splsdrcox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$splsdrcox_dynamic)
  return(model)
}

cv.splsdrcox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.splsdrcox_dynamic)
  return(model)
}
