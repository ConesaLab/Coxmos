#### ### ##
# METHODS #
#### ### ##

#' PLS-ICOX
#' @description Performs a PLS-ICOX model (based on plsRcox R package idea).
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param n.comp Numeric. Number of principal components to compute in the PLS model.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W)
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param returnData Logical. Return original and normalized X and Y matrices.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "PLS-ICOX". The class contains the following elements:
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
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export

plsicox <- function (X, Y,
                     n.comp = 2,
                     x.center = TRUE, x.scale = FALSE,
                     y.center = FALSE, y.scale = FALSE,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                     remove_non_significant = F, alpha = 0.05, tol = 1e-15,
                     MIN_EPV = 5, returnData = T, verbose = F){

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

  #### INITIALISING VARIABLES
  Ts <- NULL
  W <- NULL
  W_norm <- NULL
  P <- NULL
  E <- list(Xh)

  XXNA <- is.na(Xh) #T is NA
  YNA <- is.na(Yh) #T is NA

  #### ### ### ### ### ### ### ### ### ### ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  #### ### ### ### ### ### ### ### ### ### ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###

  #Update NAs by 0s
  Xh[XXNA] <- 0

  stopped = F
  for(h in 1:n.comp){

    # Break if residuals is not the solution because in some cases, it is not a problem to compute wh vector. !!!!

    # #residuals^2
    # residuals <- sqrt(colSums(Xh^2, na.rm=TRUE))
    #
    # #break iteration
    # if(any(residuals<tol)){
    #   bad_var <- paste0(names(residuals)[residuals < tol], collapse = ", ")
    #   message(paste0(paste0("Individual COX model cannot be computed for variables (", bad_var, "). Stopped at component ", h, ".")))
    #   stopped = T
    #   break
    # }

    #### ### ### ### ### ### ### ### ### ### ### #
    #                                            #
    #     Weight computation for each model      #
    #                                            #
    #### ### ### ### ### ### ### ### ### ### ### #

    #### ### ### ##
    ##  PLS-COX  ##
    #### ### ### ##

    #2. wh <- individual cox regression vector
    Xh[XXNA] <- NA
    Xh[XXNA] <- 0

    #Sometimes, fit can not be compute by NA at cox calculus, we cannot avoid printing the NAs matrix... !!!!
    wh <- tryCatch(
      # Specifying expression
      expr = {
        as.matrix(apply(Xh, 2, function(x){
          eps = 1e-14
          control <- survival::coxph.control(eps = eps, toler.chol = .Machine$double.eps^0.90,
                                             iter.max = 220, toler.inf = sqrt(eps), outer.max = 100, timefix = TRUE)
          fit <- survival::coxph(survival::Surv(time = time, event = event, type = "right") ~ ., as.data.frame(cbind(Ts,x)), control = control)
          fit$coefficients["x"] #cause variable of study is called 'x' and we extract new coefficient taking into account components already computed
        }))
      },
      # Specifying error message
      error = function(e){
        message(paste0("plsicox: ", e))
        invisible(gc())
        return(NA)
      }
    )

    if(all(is.na(wh))){
      message(paste0("Individual COX models cannot be computed for each variable. Stopped at component ", h, "."))
      stopped = T
      break
    }

    if(any(is.na(wh))){
      message(paste0(paste0("Individual COX model cannot be computed for variables (", paste0(rownames(wh)[is.na(wh)], collapse = ", ") ,")."), " Stopped at component ", h, "."))
      stopped = T
      break
    }

    # if(h!=1){
    #   wh <- wh[h,]
    # }
    wh <- wh[,1]

    #3. wh <- wh / ||wh||
    wh_norm <- wh/as.vector(sqrt(sum(wh^2))) #as.vector(sqrt(t(wh) %*% wh))

    #4. t = Xh wh / wh'wh
    #4. t = Xh wh_norm (solo si wh ya normalizado)
    #normalization for NAs
    Xh[XXNA] <- 0
    th <- (Xh %*% wh_norm)/((!XXNA) %*% wh_norm^2)

    #th <- t(lm(t(Xh)~0 + wh_norm)$coefficients)/((!XXNA)%*%(wh_norm^2))

    #deberia ser...
    #th <- t(lm(t(Xh)~0+wh)$coefficients)
    #th <- th/as.vector(sqrt(sum(th^2)))

    #5. p
    #ph <- t(Xh) %*% th/as.vector(t(th) %*% th)
    #normalization for NAs
    ph <- t(t(th) %*% Xh) / (t((!XXNA)) %*% th^2)

    # temppp <- rep(0,res$nc)
    # for (jj in 1:(res$nc)) {
    #   temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod((!XXNA)[,jj],temptt^2))
    # }

    #6. Residuals
    #res$residXX <- XXwotNA-temptt%*%temppp #residuals to the new matrix to get next components
    Xh <- Xh - th %*% t(ph)

    Ts <- cbind(Ts, th)
    colnames(Ts) <- paste0("comp_", 1:h)
    P <- cbind(P, ph)
    colnames(P) <- paste0("comp_", 1:h)
    W <- cbind(W, wh)
    colnames(W) <- paste0("comp_", 1:h)
    W_norm <- cbind(W_norm, wh_norm)
    colnames(W_norm) <- paste0("comp_", 1:h)
    E[[h]] <- Xh
  }

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ##
  ###              PLS-COX                 ##
  #### ### ### ### ### ### ### ### ### ### ##

  if(stopped & h==1){ #if it is the first component, no cox model has been computed

    func_call <- match.call()
    invisible(gc())

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")

    return(plsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "loadings" = NULL, "weightings" = NULL, "weightings_norm" = NULL, "W.star" = NULL, "scores" = NULL, "E" = NULL, "x.mean" = xmeans, "x.sd" = xsds),
                              Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                              beta_matrix = NULL, #NEED TO BE COMPUTED
                              survival_model = NULL,
                              n.comp = n.comp,
                              call = func_call,
                              X_input = if(returnData) X_original else NA,
                              Y_input = if(returnData) Y_original else NA,
                              nzv = variablesDeleted,
                              class = pkg.env$plsicox,
                              time = time)))
  }

  cox_model = NULL

  d <- as.data.frame(cbind(as.matrix(Ts)))
  h <- ncol(Ts)
  aux <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = d[,1:h,drop=F],
                      ties = "efron",
                      singular.ok = T,
                      robust = T,
                      nocenter = rep(1, ncol(d[,1:h,drop=F])),
                      model=T)
    },
    # Specifying error message
    error = function(e){
      message(e$message)
      invisible(gc())
      return(NA)
    }
  )

  # keep at least one component
  while(all(is.na(aux)) & h>1){
    h <- h-1
    aux <- tryCatch(
      # Specifying expression
      expr = {
        survival::coxph(formula = survival::Surv(time,event) ~ .,
                        data = d[,1:h,drop=F],
                        ties = "efron",
                        singular.ok = T,
                        robust = T,
                        nocenter = rep(1, ncol(d[,1:h,drop=F])),
                        model=T)
      },
      # Specifying error message
      error = function(e){
        if(verbose){
          message(e$message)
        }
        invisible(gc())
        return(NA)
      }
    )
  }

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  removed_variables <- NULL
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(d))){
      lst_rnsc <- removeNonSignificativeCox(cox = aux, alpha = alpha, cox_input = d, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = aux, alpha = alpha, cox_input = cbind(d, Yh), time.value = NULL, event.value = NULL)
    }

    aux <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  cox_model <- NULL
  cox_model$fit <- aux

  #we cannot compute all components
  if(h != n.comp & !all(is.na(cox_model$fit))){
    if(verbose){
      message(paste0("Model cannot be computed for all components. Final model select ", h," components instead of ", n.comp,"."))
    }
    #update all values
    W <- W[,1:h,drop=F]
    W_norm = W_norm[,1:h,drop=F]
    P = P[,1:h,drop=F]
    Ts = Ts[,1:h,drop=F]
    E = E[1:h]
    n.comp = ncol(Ts)
  }

  survival_model = NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  if(is.null(P) | is.null(W)){
    message(paste0(pkg.env$plsicox, " model cannot be computed because P or W vectors are NULL. Returning NA."))
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
    message(paste0(pkg.env$plsicox, " model cannot be computed due to solve(t(P) %*% W). Reduce 'tol' parameter to fix it. Returning NA."))
    invisible(gc())
    return(NA)
  }

  # What happen when you cannot compute W.star but you have P and W?
  W.star <- W %*% PW

  rownames(Ts) <- rownames(X)
  rownames(P) <- rownames(W_norm) <- rownames(W) <-  rownames(W.star) <- colnames(Xh)

  if(stopped){
    colnames(Ts) <- colnames(P) <- colnames(W_norm) <- colnames(W) <-  colnames(W.star) <- paste0("comp_", 1:h)
  }else{
    colnames(Ts) <- colnames(P) <- colnames(W_norm) <- colnames(W) <-  colnames(W.star) <- paste0("comp_", 1:n.comp)
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(plsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "loadings" = P, "weightings" = W, "weightings_norm" = W_norm, "W.star" = W.star, "scores" = Ts, "E" = E, "x.mean" = xmeans, "x.sd" = xsds),
                            Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                            survival_model = survival_model,
                            n.comp = n.comp,
                            call = func_call,
                            X_input = if(returnData) X_original else NA,
                            Y_input = if(returnData) Y_original else NA,
                            alpha = alpha,
                            removed_variables_cox = removed_variables,
                            nzv = variablesDeleted,
                            class = pkg.env$plsicox,
                            time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation PLS-ICOX
#' @description PLS-ICOX cross validation model
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation.
#' @param n_run Number. Number of runs for cross validation.
#' @param k_folds Number. Number of folds for cross validation.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation. A non-significant model is a model with at least one component/variable with a P-Value higher than the alpha cutoff.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC evaluator. All three weights must sum 1 (default: 0).
#' @param w_c.index Numeric. Weight for C-Index evaluator. All three weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All three weights must sum 1 (default: 1).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 15 points will be selected equally distributed.
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different EN.alpha.list to continue evaluating. If not reached for the next MIN_COMP_TO_CHECK penalties and the minimum MIN_AUC is reach, the evaluation stop.
#' @param MIN_AUC Numeric. Minimum AUC desire.
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties to check whether the AUC improves.
#' @param pred.attr Character. Method for average the AUC. Must be one of the following: "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC method for evaluation. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC")
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test observations. Once all have their linear predictors, the evaluation is perform across all the observations together (default: FALSE).
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param return_models Logical. Return all models computed in cross validation.
#' @param tol Numeric. Tolerance for solving: solve(t(P) %*% W)
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for perform the runs/folds divisions.
#'
#' @return Instance of class "HDcox" and model "cv.PLS-ICOX".
#' @export

cv.plsicox <- function (X, Y,
                        max.ncomp = 10,
                        n_run = 10, k_folds = 10,
                        x.center = TRUE, x.scale = FALSE,
                        y.center = FALSE, y.scale = FALSE,
                        remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                        remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                        w_AIC = 0, w_c.index = 0, w_AUC = 1, times = NULL,
                        MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                        pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                        MIN_EPV = 5, return_models = F, tol = 1e-15,
                        PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  #### ### ###
  # WARNINGS #
  #### ### ###

  #Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_AUC))
  max.ncomp <- check.ncomp(X, max.ncomp)
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### MAX PREDICTORS
  max.ncomp <- check.maxPredictors(X, Y, MIN_EPV, max.ncomp)

  #### REQUIREMENTS
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  #### #
  # CV #
  #### #

  set.seed(seed)
  lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###

  total_models <- 1 * k_folds * n_run #with greatest component we have all of them

  comp_model_lst <- get_HDCOX_models2.0(method = pkg.env$plsicox,
                                       lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                       max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL, n_run = n_run, k_folds = k_folds,
                                       x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
                                       remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL,
                                       remove_non_significant = remove_non_significant,
                                       total_models = total_models, tol = tol, PARALLEL = PARALLEL, verbose = verbose)

  if(all(is.na(unlist(comp_model_lst)))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.plsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.plsicox, time = time)))
    }else{
      return(cv.plsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.plsicox, time = time)))
    }
  }

  # comp_model_lst <- get_HDCOX_models(method = pkg.env$plsicox,
  #                                    lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
  #                                    max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL, n_run = n_run, k_folds = k_folds,
  #                                    x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
  #                                    remove_near_zero_variance = remove_near_zero_variance,
  #                                    total_models = total_models)

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
      return(cv.plsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.plsicox, time = time)))
    }else{
      return(cv.plsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.plsicox, time = time)))
    }
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_comp_flag <- NULL

  if(w_AUC!=0){
    #total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp)#inside get_COX_evaluation_AUC

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, #total_models = total_models,
                                     method.train = pkg.env$plsicox, PARALLEL = F)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
    optimal_comp_index <- lst_df$optimal_comp_index
    optimal_comp_flag <- lst_df$optimal_comp_flag
  }else{
    df_results_evals_fold <- df_results_evals
  }

  #### ### ### #
  # BEST MODEL #
  #### ### ### #
  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_AUC, colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[optimal_comp_index,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  #### ###
  # PLOTS #
  #### ###
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, max.ncomp = max.ncomp,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", x.text = "Component")
  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp
  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  #### ### #
  # RETURN #
  #### ### #

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.plsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, class = pkg.env$cv.plsicox, time = time)))
  }else{
    return(cv.plsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, class = pkg.env$cv.plsicox, time = time)))
  }

}

### ## ##
# CLASS #
### ## ##

plsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$plsicox)
  return(model)
}

cv.plsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.plsicox)
  return(model)
}
