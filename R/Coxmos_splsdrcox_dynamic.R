#### ### ##
# METHODS #
#### ### ##

#' sPLS-DRCOX Dynamic
#' @description
#' The sPLS-DRCOX Dynamic function conducts a sparse partial least squares deviance residual Cox
#' regression analysis using a dynamic variable selection approach. This method is particularly
#' useful for high-dimensional survival data where variable selection is crucial. The function returns
#' a model of class "Coxmos" with the attribute model specified as "sPLS-DRCOX".
#'
#' @details
#' The function employs a sparse partial least squares (sPLS) approach combined with deviance
#' residuals from a Cox model to handle survival data. The dynamic variable selection methodology
#' ensures that only the most relevant predictors are included in the model, enhancing interpretability
#' and potentially improving predictive performance.
#'
#' The input matrices X and Y represent the explanatory and response variables, respectively. It is
#' essential to note that qualitative variables in X should be transformed into binary format. The
#' response matrix Y should have two columns named "time" and "event", where the "event" column can
#' take values 0/1 or FALSE/TRUE, representing censored and event observations.
#'
#' Several parameters allow users to fine-tune the model. For instance, n.comp determines the number
#' of latent components for the PLS model, and vector aids in computing the optimal number of variables.
#' Parameters like MIN_NVAR and MAX_NVAR define the range for computing cut points to select the best
#' number of variables. The function also provides options for data preprocessing, such as centering
#' and scaling the X matrix and removing variables with near-zero or zero variance.
#'
#' The evaluation metric for determining the best number of variables can be specified using the
#' EVAL_METHOD parameter. The function supports various evaluation algorithms for assessing model
#' performance, as indicated by the pred.method parameter.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 10).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 10).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use (default: 1000).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
#' as n.cut_points^n.blocks models will be computed as minimum (default: 5).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01).
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
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param EVAL_METHOD Character. If EVAL_METHOD = "AUC", AUC metric will be use to compute the best
#' number of variables. In other case, c-index metric will be used (default: "AUC").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "Coxmos" and model "sPLS-DRCOX-Dynamic". The class contains the following
#' elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: PLS weights
#'  \item \code{(W.star)}: PLS W* vector
#'  \item \code{(loadings)}: sPLS loadings
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
#' \code{n.comp}: Number of components selected.
#'
#' \code{n.varX}: Number of Variables selected in each PLS component.
#'
#' \code{var_by_component}: Variables selected in each PLS component.
#'
#' \code{plot_accuracyPerVariable}: If NULL vector is selected, return a plot for understanding the
#' number of variable selection.
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
#' \code{alpha}: alpha value selected
#'
#' \code{nsv}: Variables removed by cox alpha cutoff.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{nz_coeffvar}: Variables removed by coefficient variation near zero.
#'
#' \code{class}: Model class.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @references
#' \insertRef{Bastien_2008}{Coxmos}
#' \insertRef{Bastien_2015}{Coxmos}
#' \insertRef{MixOmics}{Coxmos}
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' X <- X_proteomic[,1:50]
#' Y <- Y_proteomic
#' splsdrcox_dynamic(X, Y, n.comp = 3, vector = NULL, x.center = TRUE, x.scale = TRUE)

splsdrcox_dynamic <- function (X, Y,
                               n.comp = 4, vector = NULL,
                               MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                               MIN_AUC_INCREASE = 0.01,
                               x.center = TRUE, x.scale = FALSE,
                               remove_near_zero_variance = TRUE, remove_zero_variance = TRUE,
                               toKeep.zv = NULL,
                               remove_non_significant = FALSE, alpha = 0.05,
                               EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                               times = NULL, max_time_points = 15,
                               MIN_EPV = 5, returnData = TRUE, verbose = FALSE){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

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

  #### Check rownames
  lst_check <- checkXY.rownames(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### REQUIREMENTS
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  # if(k_folds.mixOmics<2){
  #   message("MixOmics k_folds must be 2 as minimum.")
  #   k_folds.mixOmics = 2
  # }

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

  #### COEF VARIATION
  lst_dnzc <- deleteNearZeroCoefficientOfVariation(X = X)
  X <- lst_dnzc$X
  variablesDeleted_cvar <- lst_dnzc$variablesDeleted

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

  XXNA <- is.na(Xh) #TRUE is NA
  YNA <- is.na(Y) #TRUE is NA

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

  if(is.null(n_dr)){
    n_dr=1
  }

  #Norm Y
  mu <- mean(DR_coxph) #equivalent because Y it is not normalized
  DR_coxph <- scale(DR_coxph, center = mu, scale = FALSE) #center DR to DR / patients
  DR_coxph_ori <- DR_coxph

  #Norm X - Do not center again
  # x_center_by_nobs <- colMeans(Xh, na.rm = TRUE)
  # Xh <- scale(Xh, center = x_center_by_nobs, scale = FALSE)
  flag = TRUE
  cv.spls <- NA

  #### ### ### ### ### ### ### ### ### ###
  # DIVIDE Y VENCERAS - BEST VECTOR SIZE #
  #### ### ### ### ### ### ### ### ### ###
  plotVAR <- NULL
  if(is.null(vector)){
    lst_BV <- getBestVector(Xh, DR_coxph, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                           EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = FALSE, mode = "spls", times = times, max_time_points = max_time_points, verbose = verbose)
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
                                 EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = FALSE, mode = "spls", times = times, max_time_points = max_time_points, verbose = verbose)
        keepX <- lst_BV$best.keepX
        plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
    }
  }

  #### ### ### ### ### ### ### ### ### ### ###
  ### ##             sPLS              ###  ##
  #### ### ### ### ### ### ### ### ### ### ###

  spls <- mixOmics::spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = rep(keepX, n.comp), scale = FALSE)

  # PREDICTION

  # both functions work fine (predict and predict_mixOmixs.pls)
  # predplsfit <- predict_mixOmixs.pls(spls, newdata=Xh[,rownames(spls$loadings$X),drop = FALSE])
  # predplsfit <- predict(spls, newdata=Xh[,rownames(spls$loadings$X),drop = FALSE])

  # sometimes solve(t(P) %*% W) in predict can cause an error
  # system is computationally singular: reciprocal condition number = 6.24697e-18
  predplsfit <- tryCatch(expr = {predict(spls, newdata=Xh[,rownames(spls$loadings$X),drop = FALSE])}, #mixomics
                         error = function(e){
                           if(verbose){
                             message("Predicting values using a pseudo-inverse matrix...\n")
                           }
                           # Estimation matrix W, P and C
                           predict <- NULL
                           Pmat = crossprod(spls$X, spls$variates$X)
                           Cmat = crossprod(spls$Y, spls$variates$X)
                           Wmat = spls$loadings

                           PW <- list()
                           for(i in 1:n.comp){
                             PW[[i]] <- tryCatch(expr = {MASS::ginv(t(Pmat[,1:i]) %*% Wmat$X[,1:i])},
                                                 error = function(e){
                                                   if(verbose){
                                                     message(e$message)
                                                   }
                                                   NA
                                                 })
                           }


                           Ypred = lapply(1:n.comp, function(x){Xh %*% Wmat$X[, 1:x] %*% PW[[x]] %*% t(Cmat)[1:x, ]})
                           Ypred = sapply(Ypred, function(x){x}, simplify = "array")
                           predict = array(Ypred, c(nrow(spls$X), ncol(spls$Y), n.comp)) # in case one observation and only one Y, we need array() to keep it an array with a third dimension being ncomp

                           predplsfit <- list()
                           predplsfit$predict <- predict
                           predplsfit
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

  d <- as.data.frame(tt_splsDR[,,drop = FALSE])
  rownames(d) <- rownames(X)
  colnames(d) <- paste0("comp_", 1:n.comp_used)
  cox_model <- NULL

  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = d,
                      ties = "efron",
                      singular.ok = TRUE,
                      robust = TRUE,
                      nocenter = rep(1, ncol(d)),
                      model = TRUE, x = TRUE)
    },
    # Specifying error message
    error = function(e){
      message(paste0("splsdrcox_dynamic: ",conditionMessage(e)))
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
    lst_model <- removeNAorINFcoxmodel(model = cox_model$fit, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAorINFcoxmodel(model = cox_model$fit, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
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

  survival_model <- NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  #get W.star
  W <- ww_splsDR
  P <- pp_splsDR

  if(is.null(P) | is.null(W)){
    message(paste0(pkg.env$splsdrcox_dynamic," model cannot be computed because P or W vectors are NULL. Returning NA."))
    # invisible(gc())
    return(NA)
  }

  #W.star
  #sometimes solve(t(P) %*% W)
  #system is computationally singular: reciprocal condition number = 6.24697e-18
  # PW <- tryCatch(expr = {solve(t(P) %*% W, tol = tol)},
  #                error = function(e){
  #                  if(verbose){
  #                    message(e$message)
  #                  }
  #                  NA
  #                })
  PW <- tryCatch(expr = {MASS::ginv(t(P) %*% W)},
                 error = function(e){
                   if(verbose){
                     message(e$message)
                   }
                   NA
                 })

  if(all(is.na(PW))){
    message(paste0(pkg.env$splsdrcox_dynamic," model cannot be computed due to ginv(t(P) %*% W). Multicollineality could be present in your data. Returning NA."))
    # invisible(gc())
    return(NA)
  }

  W.star <- W %*% PW
  Ts <- tt_splsDR

  rownames(Ts) <- rownames(X)
  rownames(P) <- rownames(W) <-  rownames(W.star) <- rownames(ww_splsDR)

  colnames(Ts) <- colnames(P) <- colnames(W) <-  colnames(W.star)<- paste0("comp_", 1:n.comp_used)

  # variable per component
  var_by_component = list()
  for(cn in colnames(W)){
    var_by_component[[cn]] <- rownames(W)[W[,cn]!=0]
  }

  #if we filter some components
  if(n.comp_used != length(names(cox_model$fit$coefficients))){
    if(verbose){
      message(paste0("Updating vectors. Final model select ", length(names(cox_model$fit$coefficients))," components instead of ", n.comp,"."))
    }
    #update all values
    which_to_keep <- which(colnames(W) %in% names(cox_model$fit$coefficients))

    W <- W[,names(cox_model$fit$coefficients),drop = FALSE]
    W.star = W.star[,names(cox_model$fit$coefficients),drop = FALSE]
    P = P[,names(cox_model$fit$coefficients),drop = FALSE]
    Ts = Ts[,names(cox_model$fit$coefficients),drop = FALSE]

    E = E[which_to_keep]
    n.comp = ncol(max(which_to_keep))
  }

  func_call <- match.call()

  if(!returnData){
    survival_model <- removeInfoSurvivalModel(survival_model)
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(splsdrcox_dynamic_class(list(X = list("data" = if(returnData) X_norm else NA,
                                               "weightings" = W, #used for computed number of variables, bc mixomics do not put 0 in loadings
                                               "W.star" = W.star,
                                               "loadings" = P,
                                               "scores" = Ts,
                                               "E" = if(returnData) E else NA,
                                               "x.mean" = xmeans, "x.sd" = xsds),
                                      Y = list("deviance_residuals" = if(returnData) DR_coxph_ori else NA,
                                               "dr.mean" = NULL, "dr.sd" = NULL, #deviance_residuals object already centered
                                               "data" = Yh,
                                               "y.mean" = ymeans, "y.sd" = ysds),
                                      survival_model = survival_model,
                                      n.comp = n.comp_used, #number of components
                                      n.varX = n.varX_used,
                                      var_by_component = var_by_component, #variables selected for each component
                                      plot_accuracyPerVariable = plotVAR,
                                      call = if(returnData) func_call else NA,
                                      X_input = if(returnData) X_original else NA,
                                      Y_input = if(returnData) Y_original else NA,
                                      R2 = R2,
                                      SCR = SCR,
                                      SCT = SCT,
                                      alpha = alpha,
                                      nsv = removed_variables,
                                      nzv = variablesDeleted,
                                      nz_coeffvar = variablesDeleted_cvar,
                                      class = pkg.env$splsdrcox_dynamic,
                                      time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation sPLS-DRCOX
#' @description The function cv.splsdrcox_dynamic conducts a cross-validation for the sPLS-DRCOX model,
#' which is a specialized model tailored for survival analysis. The function aims to optimize the model's
#' performance by determining the best number of PLS components and variables through cross-validation.
#'
#' @details
#' The cv.splsdrcox_dynamic function is designed to perform cross-validation for the sPLS-DRCOX model,
#' a specialized model for survival analysis. The function's primary objective is to identify the
#' optimal number of PLS components and variables that yield the best model performance.
#'
#' The function accepts both numeric matrices and data frames for explanatory (X) and response (Y)
#' variables. It is essential to ensure that qualitative variables in X are transformed into binary
#' format. The response variable Y should have two columns: "time" and "event". The event column
#' should contain binary values, where 0/1 or FALSE/TRUE represent censored and event observations,
#' respectively.
#'
#' The cross-validation process is controlled by several parameters, including the maximum number of
#' PLS components (max.ncomp), the number of runs (n_run), and the number of folds (k_folds). The
#' function also provides options for data preprocessing, such as centering and scaling of the X matrix,
#' and removal of variables with near-zero or zero variance.
#'
#' Significance testing is incorporated into the model evaluation process. Users can specify the alpha
#' threshold (alpha) for determining significance. Non-significant models or variables can be optionally
#' removed from the evaluation based on user-defined criteria.
#'
#' The function also offers flexibility in model evaluation metrics. Users can choose between different
#' metrics such as AUC, AIC, C-Index, and Brier Score. The importance of each metric in the evaluation
#' can be controlled using weights (w_AIC, w_c.index, w_AUC, w_BRIER).
#'
#' For computational efficiency, the function provides an option to run the cross-validation in parallel
#' (PARALLEL). Additionally, verbose logging can be enabled to display extra messages during the execution.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event
#' observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation
#' (default: 8).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL).
#' @param n_run Numeric. Number of runs for cross validation (default: 3).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance
#' variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will
#' be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance
#' filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero
#' variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE,
#' non-significant models are removed before computing the evaluation. A non-significant model is a
#' model with at least one component/variable with a P-Value higher than the alpha cutoff.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant
#' variables/components in final cox model will be removed until all variables are significant by
#' forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 10).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use (default: 1000).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
#' as n.cut_points^n.blocks models will be computed as minimum (default: 5).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops (default: 0.01).
#' @param EVAL_METHOD Character. If EVAL_METHOD = "AUC", AUC metric will be use to compute the best
#' number of variables. In other case, c-index metric will be used (default: "AUC").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_c.index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is
#' reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the
#' 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC
#' improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the
#' evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following:
#' "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated simultaneously.
#' If fast_mode = FALSE, for each run, all linear predictors are computed for test observations. Once
#' all have their linear predictors, the evaluation is perform across all the observations together
#' (default: FALSE).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your
#' total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "Coxmos" and model "cv.sPLS-DRCOX-Dynamic".
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha
#' information).
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
#' \code{lst_train_indexes}: List (of lists) of indexes for the observations used in each run/fold
#' for train the models.
#' \code{lst_test_indexes}: List (of lists) of indexes for the observations used in each run/fold
#' for test the models.
#'
#' \code{time}: time consumed for running the cross-validated function.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_proteomic$event, p = .5, list = FALSE, times = 1)
#' X_train <- X_proteomic[index_train,1:20]
#' Y_train <- Y_proteomic[index_train,]
#' cv.splsdrcox_dynamic_model <- cv.splsdrcox_dynamic(X_train, Y_train, max.ncomp = 1, vector = NULL,
#' n_run = 1, k_folds = 2, x.center = TRUE, x.scale = TRUE)

cv.splsdrcox_dynamic <- function (X, Y,
                                  max.ncomp = 8, vector = NULL,
                                  n_run = 3, k_folds = 10,
                                  x.center = TRUE, x.scale = FALSE,
                                  remove_near_zero_variance = TRUE, remove_zero_variance = TRUE,
                                  toKeep.zv = NULL, remove_variance_at_fold_level = FALSE,
                                  remove_non_significant_models = FALSE, remove_non_significant = FALSE,
                                  alpha = 0.05,
                                  MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                  MIN_AUC_INCREASE = 0.01,
                                  EVAL_METHOD = "AUC",
                                  w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL,
                                  max_time_points = 15,
                                  MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                  pred.attr = "mean", pred.method = "cenROC", fast_mode = FALSE,
                                  max.iter = 200,
                                  MIN_EPV = 5, return_models = FALSE, returnData = FALSE,
                                  PARALLEL = FALSE, verbose = FALSE, seed = 123){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

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

  #### Check rownames
  lst_check <- checkXY.rownames(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars(X)

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :penalty]"

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

  #### COEF VARIATION
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnzc <- deleteNearZeroCoefficientOfVariation(X = X)
    X <- lst_dnzc$X
    variablesDeleted_cvar <- lst_dnzc$variablesDeleted
  }else{
    variablesDeleted_cvar <- NULL
  }

  #### #
  # CV #
  #### #
  # lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds, seed = seed) #FOR TEST
  # lst_X_train <- lst_data$lst_X_train
  # lst_Y_train <- lst_data$lst_Y_train
  # lst_X_test <- lst_data$lst_X_test
  # lst_Y_test <- lst_data$lst_Y_test
  # k_folds <- lst_data$k_folds
  #
  # lst_train_indexes <- lst_data$lst_train_index
  # lst_test_indexes <- lst_data$lst_test_index

  lst_data <- splitData_Iterations_Folds_indexes(Y, n_run = n_run, k_folds = k_folds, seed = seed) #FOR TEST

  lst_train_indexes <- lst_data$lst_train_index
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  total_models <- 1 * k_folds * n_run

  comp_model_lst  <- get_Coxmos_models2.0(method = pkg.env$splsdrcox_dynamic,
                                         X_train = X, Y_train = Y,
                                         lst_X_train = lst_train_indexes, lst_Y_train = lst_train_indexes,
                                         max.ncomp = max.ncomp, penalty.list = NULL, EN.alpha.list = NULL, max.variables = NULL, vector = vector,
                                         n_run = n_run, k_folds = k_folds,
                                         MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE, EVAL_METHOD = EVAL_METHOD,
                                         n.cut_points = n.cut_points,
                                         x.center = x.center, x.scale = x.scale,
                                         y.center = y.center, y.scale = y.scale,
                                         remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = FALSE, toKeep.zv = NULL,
                                         alpha = alpha, MIN_EPV = MIN_EPV,
                                         remove_non_significant = remove_non_significant, tol = tol,
                                         max.iter = max.iter, times = times, pred.method = pred.method, max_time_points = max_time_points,
                                         returnData = returnData, total_models = total_models,
                                         PARALLEL = PARALLEL, verbose = verbose)

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
                                                    max.ncomp = max.ncomp, penalty.list = NULL, n_run = n_run, k_folds = k_folds,
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
  optimal_comp_flag <- FALSE
  optimal_eta_index <- NULL
  optimal_eta <- NULL

  if(TRUE){ #compute always BRIER SCORE
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL = FALSE
    lst_df <- get_COX_evaluation_BRIER(comp_model_lst = comp_model_lst,
                                       fast_mode = fast_mode,
                                       X_test = X, Y_test = Y,
                                       lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                       df_results_evals = df_results_evals, times = times,
                                       pred.method = pred.method, pred.attr = pred.attr,
                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       w_BRIER = w_BRIER, method.train = pkg.env$splsdrcox_dynamic, PARALLEL = FALSE, verbose = verbose)

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

    #As we are measuring just one evaluator and one method - PARALLEL = FALSE
    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     X_test = X, Y_test = Y,
                                     lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, method.train = pkg.env$splsdrcox_dynamic, PARALLEL = FALSE, verbose = verbose)

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
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index,, drop = FALSE][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = TRUE)),, drop = FALSE][1,]
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

  # invisible(gc())
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
