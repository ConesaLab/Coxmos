#### ### ##
# METHODS #
#### ### ##

#' sPLS-ICOX
#' @description This function performs a sparse partial least squares individual Cox (sPLS-ICOX) (based on plsRcox R package).
#' The function returns a HDcox model with the attribute model as "sPLS-ICOX".
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 10).
#' @param spv_penalty Numeric. Penalty for variable selection for the individual cox models. Variables with a lower P-Value than "spv_penalty" in the individual cox analysis will be keep for the sPLS-ICOX approach (default: 1).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
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
#' @return Instance of class "HDcox" and model "sPLS-ICOX". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: sPLS weights
#'  \item \code{(weightings_norm)}: sPLS normalize weights
#'  \item \code{(W.star)}: sPLS W* vector
#'  \item \code{(loadings)}: sPLS loadings
#'  \item \code{(scores)}: sPLS scores/variates
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
#' \code{var_by_component}: Variables selected by each component.
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
#' splsicox(X, Y)
#' splsicox(X, Y, n.comp = 3, spv_penalty = 0.5, x.center = TRUE, x.scale = TRUE)
#' }

splsicox <- function(X, Y,
                     n.comp = 4, spv_penalty = 1,
                     x.center = TRUE, x.scale = FALSE,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                     remove_non_significant = F, alpha = 0.05, tol = 1e-10,
                     MIN_EPV = 5, returnData = T, verbose = F){

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### Check values classes and ranges
  params_with_limits <- list("alpha" = alpha, "eta" = spv_penalty)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("n.comp" = n.comp,
                  "MIN_EPV" = MIN_EPV, "tol" = tol)
  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                      "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                      "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose)
  check_class(logical_params, class = "logical")

  #### REQUIREMENTS
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

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
  if(length(XXNA)>0){
    Xh[XXNA] <- 0
  }

  var_by_component <- list()
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
    if(length(XXNA)>0){
      Xh[XXNA] <- 0
    }

    #Sometimes, fit can not be compute by NA at cox calculus, we cannot avoid printing the NAs matrix... !!!!
    # returning wh[,1] coefficients and wh[,2] p-values
    wh <- getIndividualCox(data = cbind(Xh, Yh), time_var = "time", event_var = "event", score_data = Ts)
    # wh2 <- tryCatch(
    #   # Specifying expression
    #   expr = {
    #     as.data.frame(t(apply(Xh, 2, function(x){
    #       eps = 1e-14
    #       control <- survival::coxph.control(eps = eps, toler.chol = .Machine$double.eps^0.90,
    #                                          iter.max = 220, toler.inf = sqrt(eps), outer.max = 100, timefix = TRUE)
    #       fit <- survival::coxph(survival::Surv(time = time,
    #                                             event = event,
    #                                             type = "right") ~ ., as.data.frame(cbind(Ts,x)),
    #                              control = control,
    #                              singular.ok = T)
    #
    #       if(length(getPvalFromCox(fit))==1){
    #         aux <- c(fit$coefficients["x"], getPvalFromCox(fit))
    #       }else{
    #         aux <- c(fit$coefficients["x"], getPvalFromCox(fit)["x"]) #cause variable of study is called 'x' and we extract new coefficient taking into account components already computed
    #       }
    #       aux
    #     })))
    #   },
    #   # Specifying error message
    #   error = function(e){
    #     message(paste0("splsicox: ", e))
    #     invisible(gc())
    #     return(NA)
    #     #if error we could return beta=0 (no significant) instead a NA!!!
    #   }
    # )

    # if(any(is.na(wh))){
    #   message(paste0(paste0("Individual COX model cannot be computed for variables (", paste0(rownames(wh)[is.na(wh[,1])], collapse = ", ") ,").")))
    #
    #   #wh <- wh[-which(is.na(wh[,1])),]
    #   #replace for beta of 0, and p-value of 1
    #   wh[which(is.na(wh[,1])),] <- c(rep(0, length(rownames(wh)[is.na(wh[,1])])), rep(1, length(rownames(wh)[is.na(wh[,1])])))
    #
    #   #stopped = T
    #   #break
    # }

    if(all(is.na(wh))){
      message(paste0("Stopping at component ", h-1, ": The weight vector could not be computed.."))
      h = h-1
      stopped = T
      break
    }

    ### ## ## ##
    #filter variables by p-val cutoff
    ### ## ## ##
    index2zero <- which(wh[,2,drop=T]>spv_penalty)
    index2keep <- which(wh[,2,drop=T]<=spv_penalty)

    if(length(index2keep)==0){
      if(verbose){
        message(paste0("Stopping at component ", h-1, ": No significant variables found in component ", h,"."))
      }
      h = h-1
      break
    }

    wh[index2zero,1] <- 0
    var_by_component[[h]] <- rownames(wh)[index2keep]

    nm <- rownames(wh)
    nm_keep <- var_by_component[[h]]
    wh <- wh[,1,drop=T]
    names(wh) <- nm

    sub_wh <- wh[nm_keep]

    #3. wh <- wh / ||wh||
    wh_norm <- data.frame(wh)
    wh_norm[,1] <- 0
    sub_wh_norm <- sub_wh/as.vector(sqrt(sum(sub_wh^2))) #as.vector(sqrt(t(wh) %*% wh))

    wh_norm[nm_keep,] <- sub_wh_norm

    #4. t = Xh wh / wh'wh
    #4. t = Xh wh_norm (solo si wh ya normalizado)
    #normalization for NAs
    if(length(XXNA)>0){
      Xh[XXNA[,nm_keep]] <- 0
    }

    if(length(XXNA)>0){
      # th <- (Xh[,nm_keep,drop=F] %*% sub_wh_norm)/((!XXNA[,nm_keep,drop=F]) %*% sub_wh_norm^2) # do not remember why
      th <- (Xh[,nm_keep,drop=F] %*% sub_wh_norm)
    }else{
      th <- (Xh[,nm_keep,drop=F] %*% sub_wh_norm)
    }

    #th <- t(lm(t(Xh)~0 + sub_wh_norm)$coefficients)/((!XXNA)%*%(sub_wh_norm^2))

    #deberia ser...
    #th <- t(lm(t(Xh)~0+wh)$coefficients)
    #th <- th/as.vector(sqrt(sum(th^2)))

    #5. p
    #ph <- t(Xh) %*% th/as.vector(t(th) %*% th)
    #normalization for NAs

    ph <- data.frame(wh)
    ph[,1] <- 0

    if(length(XXNA)>0){
      # sub_ph <- t(t(th) %*% Xh[,nm_keep,drop=F]) / (t((!XXNA[,nm_keep,drop=F])) %*% th^2) # do not remember why
      sub_ph <- t(t(th) %*% Xh[,nm_keep,drop=F]) / (as.vector(t(th) %*% th))
    }else{
      sub_ph <- t(t(th) %*% Xh[,nm_keep,drop=F]) / (as.vector(t(th) %*% th))
    }

    ph[nm_keep,] <- sub_ph
    # temppp <- rep(0,res$nc)
    # for (jj in 1:(res$nc)) {
    #   temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod((!XXNA)[,jj],temptt^2))
    # }

    #6. Residuals
    #res$residXX <- XXwotNA-temptt%*%temppp #residuals to the new matrix to get next components
    Xh_sub <- Xh[,nm_keep,drop=F] - th %*% t(sub_ph)
    Xh[,nm_keep] <- Xh_sub

    Ts <- cbind(Ts, th)
    P <- cbind(P, as.matrix(ph))
    W <- cbind(W, wh)
    W_norm <- cbind(W_norm, as.matrix(wh_norm))
    E[[h]] <- Xh

  }

  # Problems computing firts component
  if(h==0){ #no significant individual cox model at first component
    func_call <- match.call()
    invisible(gc())

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")

    return(splsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "weightings" = NULL, "weightings_norm" = NULL, "W.star" = NULL, "loadings" = NULL, "scores" = NULL, "E" = NULL, "x.mean" = xmeans, "x.sd" = xsds),
                               Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                               beta_matrix = NULL, #NEED TO BE COMPUTED
                               survival_model = NULL,
                               n.comp = h,
                               spv_penalty = spv_penalty,
                               var_by_component = var_by_component, #variables selected for each component
                               call = func_call,
                               X_input = if(returnData) X_original else NA,
                               Y_input = if(returnData) Y_original else NA,
                               nzv = variablesDeleted,
                               class = pkg.env$splsicox,
                               time = time)))
  }

  ### DELETE VARIABLES ALL 0s (in all components)
  rm_var <- names(which(rowSums(W)==0))
  W <- W[!rownames(W) %in% rm_var,,drop=F]
  W_norm <- W_norm[!rownames(W_norm) %in% rm_var,,drop=F]
  P <- P[!rownames(P) %in% rm_var,,drop=F]

  colnames(Ts) <- paste0("comp_", 1:h)
  colnames(P) <- paste0("comp_", 1:h)
  colnames(W) <- paste0("comp_", 1:h)
  colnames(W_norm) <- paste0("comp_", 1:h)
  names(var_by_component) <- paste0("comp_", 1:h)

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ##
  ###              PLS-COX                 ##
  #### ### ### ### ### ### ### ### ### ### ##

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
                      model=T, x = T)
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
                        model=T, x = T)
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

  # RETURN a MODEL with ALL significant Variables from complete, deleting one by one
  removed_variables <- NULL
  removed_variables_cor <- NULL
  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  if(all(c("time", "event") %in% colnames(d))){
    lst_model <- removeNAcoxmodel(model = aux, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAcoxmodel(model = aux, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
  }
  aux <- lst_model$model
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
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
    message(paste0(pkg.env$splsicox, " model cannot be computed because P or W vectors are NULL. Returning NA."))
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
    message(paste0(pkg.env$splsicox, " model cannot be computed due to solve(t(P) %*% W). Reduce 'tol' parameter to fix it. Returning NA."))
    invisible(gc())
    return(NA)
  }

  # What happen when you cannot compute W.star but you have P and W?
  W.star <- W %*% PW

  rownames(Ts) <- rownames(X)
  #rownames(P) <- rownames(W_norm) <- rownames(W) <-  rownames(W.star) <- colnames(Xh)

  if(stopped){
    colnames(Ts) <- colnames(P) <- colnames(W_norm) <- colnames(W) <-  colnames(W.star) <- paste0("comp_", 1:h)
  }else{
    colnames(Ts) <- colnames(P) <- colnames(W_norm) <- colnames(W) <-  colnames(W.star) <- paste0("comp_", 1:n.comp)
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(splsicox_class(list(X = list("data" = if(returnData) X_norm else NA, "weightings" = W, "weightings_norm" = W_norm, "W.star" = W.star, "loadings" = P, "scores" = Ts, "E" = E, "x.mean" = xmeans, "x.sd" = xsds),
                            Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                            survival_model = survival_model,
                            n.comp = h,
                            spv_penalty = spv_penalty,
                            var_by_component = var_by_component, #variables selected for each component
                            call = func_call,
                            X_input = if(returnData) X_original else NA,
                            Y_input = if(returnData) Y_original else NA,
                            alpha = alpha,
                            removed_variables_cox = removed_variables,
                            nzv = variablesDeleted,
                            class = pkg.env$splsicox,
                            time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' sPLS-ICOX Cross-Validation
#' @description This function performs cross-validated sparse partial least squares Cox (sPLS-ICOX).
#' The function returns the optimal number of components and the optimal sparsity penalty value based on cross-validation.
#' The performance could be based on multiple metrics as Area Under the Curve (AUC), Brier Score or C-Index.
#' Furthermore, the user could establish more than one metric simultaneously.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param spv_penalty.list Numeric vector. Penalty for variable selection for the individual cox models. Variables with a lower P-Value than "spv_penalty" in the individual cox analysis will be keep for the sPLS-ICOX approach (default: seq(0.1,1,0.1)).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation. A non-significant model is a model with at least one component/variable with a P-Value higher than the alpha cutoff.
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
#' @return Instance of class "HDcox" and model "cv.sPLS-ICOX".
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
#' cv.splsicox_model <- cv.splsicox(X, Y, max.ncomp = 10, spv_penalty.list = c(0.1,0.5,0.8),
#' x.center = TRUE, x.scale = TRUE)
#' splsicox_model <- splsicox(X, Y, n.comp = cv.splsicox_model$opt.comp,
#' spv_penalty = cv.splsicox_model$opt.spv_penalty, x.center = TRUE, x.scale = TRUE)
#' }

cv.splsicox <- function (X, Y,
                        max.ncomp = 10, spv_penalty.list = seq(0.1,1,0.1),
                        n_run = 5, k_folds = 10,
                        x.center = TRUE, x.scale = FALSE,
                        remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                        remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                        w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                        MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                        pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
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

  numeric_params <- list("max.ncomp" = max.ncomp, "spv_penalty.list" = spv_penalty.list,
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

  character_params <- list("pred.attr" = pred.attr, "pred.method" = pred.method)
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
  max.ncomp <- check.maxPredictors(X, Y, MIN_EPV, max.ncomp)
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

  set.seed(seed)
  lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
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

  total_models <- 1 * k_folds * n_run * length(spv_penalty.list) #with greatest component we have all of them
  t1 <- Sys.time()
  #we use spv_penalty.list as penalty parameter (sPLS-DRCOX and sPLS-ICOX)
  lst_model <- get_HDCOX_models2.0(method = pkg.env$splsicox,
                                   lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                   max.ncomp = max.ncomp, eta.list = spv_penalty.list, EN.alpha.list = NULL, max.variables = NULL, vector = NULL,
                                   n_run = n_run, k_folds = k_folds,
                                   MIN_NVAR = NULL, MAX_NVAR = NULL, MIN_AUC_INCREASE = NULL, EVAL_METHOD = NULL,
                                   n.cut_points = NULL,
                                   x.center = x.center, x.scale = x.scale,
                                   y.center = y.center, y.scale = y.scale,
                                   remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = F, toKeep.zv = NULL,
                                   alpha = alpha, MIN_EPV = MIN_EPV,
                                   remove_non_significant = remove_non_significant, tol = tol, max.iter = NULL,
                                   returnData = returnData, total_models = total_models,
                                   PARALLEL = PARALLEL, verbose = verbose)

  comp_model_lst = lst_model$comp_model_lst
  info = lst_model$info

  t2 <- Sys.time()
  t2-t1
  if(all(is.null(comp_model_lst))){
    message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, nzv = variablesDeleted, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, nzv = variablesDeleted, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run * length(spv_penalty.list)
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst, alpha = alpha,
                                                    max.ncomp = max.ncomp, eta.list = spv_penalty.list, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models, verbose = verbose)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsicox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.spv_penalty = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
    lst_df <- get_COX_evaluation_BRIER_sPLS(comp_model_lst = comp_model_lst,
                                            fast_mode = fast_mode,
                                            lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                            df_results_evals = df_results_evals, times = times,
                                            pred.method = pred.method, pred.attr = pred.attr,
                                            max.ncomp = max.ncomp, eta.list = spv_penalty.list, n_run = n_run, k_folds = k_folds,
                                            MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                            w_BRIER = w_BRIER, method.train = pkg.env$splsicox, PARALLEL = F, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    total_models <- ifelse(!fast_mode, n_run * max.ncomp * length(spv_penalty.list), k_folds * n_run * max.ncomp * length(spv_penalty.list))
    #total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp)#inside get_COX_evaluation_AUC

    #times should be the same for all folds
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC_sPLS(comp_model_lst = comp_model_lst,
                                          lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                          df_results_evals = df_results_evals, times = times,
                                          fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                          max.ncomp = max.ncomp, eta.list = spv_penalty.list, n_run = n_run, k_folds = k_folds,
                                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          w_AUC = w_AUC, method.train = pkg.env$splsicox, PARALLEL = F, verbose = verbose)

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
  df_results_evals_comp <- cv.getScoreFromWeight(lst_cox_mean = df_results_evals_comp, w_AIC, w_c.index, w_BRIER, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER")

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
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_BRIER = w_BRIER, max.ncomp = max.ncomp, eta.list = spv_penalty.list,
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

  df_results_evals$eta <- as.numeric(as.character(df_results_evals$eta))
  df_results_evals_run$eta <- as.numeric(as.character(df_results_evals_run$eta))
  df_results_evals_comp$eta <- as.numeric(as.character(df_results_evals_comp$eta))

  message(paste0("Best model obtained."))

  #### ### ### ### ##
  # Change eta name #
  #### ### ### ### ##
  colnames(best_model_info)[which(colnames(best_model_info)=="eta")] <- "spv_penalty"
  colnames(df_results_evals)[which(colnames(df_results_evals)=="eta")] <- "spv_penalty"
  colnames(df_results_evals_run)[which(colnames(df_results_evals_run)=="eta")] <- "spv_penalty"
  colnames(df_results_evals_comp)[which(colnames(df_results_evals_comp)=="eta")] <- "spv_penalty"
  ggp_AUC <- ggp_AUC + guides(color=guide_legend(title="spv_penalty"))
  ggp_BRIER <- ggp_BRIER + guides(color=guide_legend(title="spv_penalty"))
  ggp_c_index <- ggp_c_index + guides(color=guide_legend(title="spv_penalty"))
  ggp_AIC <- ggp_AIC + guides(color=guide_legend(title="spv_penalty"))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.splsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.spv_penalty = best_model_info$spv_penalty, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.splsicox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.spv_penalty = best_model_info$spv_penalty, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsicox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }

}

### ## ##
# CLASS #
### ## ##

splsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$splsicox)
  return(model)
}

cv.splsicox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.splsicox)
  return(model)
}
