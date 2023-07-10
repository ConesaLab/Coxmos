#### ### ##
# METHODS #
#### ### ##

#' sPLS-DRCOX
#' @description This function performs a sparse partial least squares deviance residual Cox (sPLS-DRCOX) (based on plsRcox R package).
#' The function returns a Coxmos model with the attribute model as "sPLS-DRCOX".
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 10).
#' @param eta Numeric (0-1). Penalty for sPLS. If eta = 0 no penalty is applied and 1 maximum penalty (no variables are selected). Equal or greater than 1 cannot be selected (default: 0.5).
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
#' @return Instance of class "Coxmos" and model "sPLS-DRCOX". The class contains the following elements:
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
#'  \item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
#'  \item \code{(dr.mean)}: mean values for deviance residuals Y matrix
#'  \item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: sPLS weights
#'  \item \code{(loadings)}: sPLS loadings
#'  \item \code{(scores)}: sPLS scores/variates
#'  \item \code{(ratio)}: r value for the sPLS model (used to perform predictions)
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
#' \code{eta}: Penalty value selected.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{var_by_component}: Variables selected in each PLS component.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{B.hat}: sPLS beta matrix
#'
#' \code{R2}: sPLS R2
#'
#' \code{SCR}: sPLS SCR
#'
#' \code{SCT}: sPLS SCT
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
#' @export
#'
#' @examples
#' \dontrun{
#' splsdrcox(X, Y)
#' splsdrcox(X, Y, n.comp = 3, eta = 0.25, x.center = TRUE, x.scale = TRUE)
#' }

splsdrcox <- function (X, Y,
                      n.comp = 4, eta = 0.5,
                      x.center = TRUE, x.scale = FALSE,
                      remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                      remove_non_significant = F, alpha = 0.05,
                      MIN_EPV = 5, returnData = T, verbose = F){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### Check values classes and ranges
  params_with_limits <- list("eta" = eta)
  check_min0_less1_variables(params_with_limits)

  params_with_limits <- list("alpha" = alpha)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("n.comp" = n.comp,
                  "MIN_EPV" = MIN_EPV, "tol" = tol)
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

  XXNA <- is.na(Xh) #T is NA
  YNA <- is.na(Y) #T is NA

  #### MAX PREDICTORS
  n.comp <- check.maxPredictors(X, Y, MIN_EPV, n.comp)

  #### ### ### ### ### ### ### ### ### ### ### ###
  ### ###             sPLS-COX             ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###

  #2. Surv function - NULL model
  coxDR <- survival::coxph(survival::Surv(time = time, event = event, type = "right") ~ 1, as.data.frame(Xh))

  #3. Residuals - Default is deviance because eval type="deviance"
  DR_coxph <- residuals(coxDR, type = "deviance") #"martingale", "deviance", "score", "schoenfeld", "dfbeta"', "dfbetas", "scaledsch" and "partial"

  #### ### ### ### ### ### #### ### ### ### ### #
  #### ### ### ### ### ### #### ### ### ### ### #
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  #### ### ### ### ### ### #### ### ### ### ### #
  #### ### ### ### ### ### #### ### ### ### ### #

  #4. sPLS Algorithm
  n_var <- ncol(Xh)
  n_dr <- ncol(DR_coxph)

  if(is.null(n_dr)){
    n_dr=1
  }

  #CENTER DEVIANCE RESIUDALS
  mu <- mean(DR_coxph) #equivalent because Y it is not normalized
  DR_coxph <- scale(DR_coxph, center = mu, scale = FALSE) #center DR to DR / patients
  DR_coxph_ori <- DR_coxph

  #### INITIALISING VARIABLES
  beta_pls <- matrix(0, n_var, n_dr)
  beta_matrix <- list()
  var_by_component <- list()
  var_by_component_nzv <- list()
  plsfit_list <- list()
  predplsfit_list <- list()
  Z_list <- list()
  ww_list <- list()

  last.pls <- NULL

  E <- list()
  R2 <- list()
  SCR <- list()
  SCT <- list()

  for(h in 1:n.comp) {
    #4.1
    #t(Xh) * Deviance Residuals
    #por cada variable, multiplicamos su valor por el residuo de cada paciente (aunque este sea el total de todo)
    #Ej. Sex = 1.1724346  1.1724346  1.1724346  1.1724346 -0.8528423  1.1724346
    #Residual = -0.1754595 -0.004732449 0.8676313 -0.757973 -0.3859317 -0.004732449
    #Res = 0.2408943 (Valor de relacion entre los residuos y las variables)
    #valores cercanos a cero es porque no hay residuos

    Xh[XXNA] <- 0
    #for matrix multiplications, NA as 0 is the same as not using that patient
    Z <- t(Xh) %*% DR_coxph #this Z is not the same as the paper
    Z_list[[h]] <- Z

    Xh[XXNA] <- NA

    #4.2 Get selected variables using eta: what <- spls.dv(Z, eta, kappa, eps, maxstep)
    # spls DEVIANCE RESIDUALS
    Z_median <- median(abs(Z))
    Z_norm <- Z/Z_median #normalizar respecto la mediana

    # ww2 <- matrix(0, n_var, 1)
    # rownames(ww2) <- colnames(Xh)
    # lambda = 0.2
    # Z2_l <- (abs(Z2) - lambda/2)
    # ww2[Z2_l >= 0] <- Z2_l[Z2_l>=0] * sign(Z2)[Z2_l>=0]

    ww <- matrix(0, n_var, 1)
    rownames(ww) <- colnames(Xh)
    if(eta < 1) {#threshold para penalizacion?
      Z_mod <- abs(Z_norm) - eta * max(abs(Z_norm)) #keep those variables greater than eta * max value

      if(sum(Z_mod >= 0)==1){ #only one variable does not allow to compute a PLS model (take the another one)
        Z_mod_neg <- Z_mod[which(Z_mod < 0),,drop=F]
        cn_extra <- rownames(Z_mod_neg)[which.max(Z_mod_neg)]
        Z_mod[cn_extra,] <- 0.01
      }

      ww[Z_mod >= 0] <- Z_mod[Z_mod >= 0] * (sign(Z_mod))[Z_mod >= 0] #keeping the sign
    }else{
      stop_quietly("eta should be a value between [0, 1), default is 0.5")
    }

    ww_list[[h]] <- ww

    #4.3 Get variables greater than 0 and variables already selected (beta !=0)
    if(h==1){
      A <- rownames(ww[which(ww != 0),,drop=F])
    }else{
      A <- unique(c(rownames(ww[which(ww != 0),,drop=F]), names(beta_matrix[,h-1])[which(beta_matrix[, h-1]!=0)]))
    }

    Xa <- Xh[,A,drop = FALSE]

    #4.4 Run standard PLS with new components and the deviance residuals - Filter near zero variables
    nZ <- caret::nearZeroVar(Xa, saveMetrics = T) #to check if we have to do some changes in the data
    td <- rownames(nZ[nZ$nzv==T,])

    #Do not delete
    # if(any(mustKeep %in% td)){
    #   td <- td[-which(td %in% mustKeep)]
    # }

    lstDeleted <- td
    if(length(lstDeleted)>0 & ncol(Xa)>2){
      Xa <- Xa[,!colnames(Xa) %in% lstDeleted, drop=F]
    }
    A_nzv <- colnames(Xa)

    #### ### ### #
    # PREDICTION #
    #### ### ### #

    #But always using the complete residuals
    # plsfit3 <- ropls::opls(x = Xa, y = DR_coxph_ori, predI = min(h, ncol(Xa)), scaleC = "none")
    # plsfit2 <- mixOmics::pls(X = Xa, Y = DR_coxph_ori, ncomp = min(h, ncol(Xa)), scale = F)
    plsfit <- tryCatch(
      # Specifying expression
      expr = {
        pls2(X = Xa, Y = DR_coxph_ori, n.comp = min(h, ncol(Xa)),
             x.center = F, x.scale = F, y.center = F, y.scale = F,
             it = 100, tol.W.star = tol, verbose = verbose)
      },
      # Specifying error message
      error = function(e){
        message(paste0("splsdrcox: ", e$message, "\n"))
        NA
      }
    )

    if(all(is.na(plsfit))){
      break
    }

    #plsfit$X$weightings_norm;plsfit2$loadings$X;plsfit3@weightMN

    XaNA <- is.na(Xa) #T is NA
    Xa[XaNA] <- 0
    predplsfit <- Xa[,rownames(plsfit$X$loadings),drop=F] %*% plsfit$B[,,drop=F]
    Xa[XaNA] <- NA

    #### ### ### ### ###
    # UPDATING RESULTS #
    #### ### ### ### ###
    predplsfit_list[[h]] <- predplsfit

    beta_pls <- matrix(0, n_var, n_dr)
    rownames(beta_pls) <- colnames(Xh)
    #beta_pls[A_nzv,] <- matrix(data = plsfit$B[,min(h, ncol(Xa)),drop=F], nrow = length(A_nzv), ncol = n_dr) #n.comp from new pls and not all
    beta_pls[A_nzv,] <- matrix(data = plsfit$B[,,drop=F], nrow = length(A_nzv), ncol = n_dr) #n.comp from new pls and not all
    beta_matrix <- cbind(beta_matrix, beta_pls) #res

    var_by_component[[h]] <- A
    var_by_component_nzv[[h]] <- A_nzv

    #### ### ### ### ##
    # UPDATING VALUES #
    #### ### ### ### ##
    #DR_coxph <- DR_coxph_ori - predplsfit[,min(h, ncol(Xa)),drop=F] #for manual prediction
    DR_coxph <- DR_coxph_ori - predplsfit[,,drop=F] #for manual prediction

    #R2 calculation
    #E[[h]] = DR_coxph_ori - predplsfit$predict[,,plsfit$n.comp] #same formula, but adding components
    E[[h]] = DR_coxph #same formula, but adding components

    SCR[[h]] = sum(apply(E[[h]],2,function(x) sum(x**2)))
    SCT[[h]] = sum(apply(as.matrix(DR_coxph_ori),2,function(x) sum(x**2))) #equivalent sum((DR_coxph_ori - mean(DR_coxph_ori))**2)
    R2[[h]] = 1 - (SCR[[h]]/SCT[[h]]) #deviance residuals explanation

    last.pls <- plsfit
  }

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ### ##
  ### ###              PLS-COX            ### ###
  #### ### ### ### ### ### ### ### ### ### ### ##

  #n.comp_used <- ncol(tt_splsDR) #can be lesser than expected because we have lesser variables to select because penalization
  n.comp_used <- ncol(last.pls$X$scores)

  #d <- as.data.frame(last.pls$X$scores[,,drop=F])
  d <- as.data.frame(last.pls$X$scores[,,drop=F])
  rownames(d) <- rownames(X)
  colnames(d) <- paste0("comp_", 1:n.comp_used)

  aux <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = d[,1:n.comp_used,drop=F],
                      ties = "efron",
                      singular.ok = T,
                      robust = T,
                      nocenter = rep(1, ncol(d[,1:n.comp_used,drop=F])),
                      model=T, x = T)
    },
    # Specifying error message
    error = function(e){
      if(verbose){
        message(e)
      }
      # invisible(gc())
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
          message(e)
        }
        # invisible(gc())
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
    lst_model <- removeNAorINFcoxmodel(model = aux, data = d, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAorINFcoxmodel(model = aux, data = cbind(d, Yh), time.value = NULL, event.value = NULL)
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
  names(var_by_component_nzv) <- paste0("comp_", 1:n.comp_used)

  #we cannot compute all components
  if(h != n.comp & !all(is.na(cox_model$fit))){
    if(verbose){
      message(paste0("Model cannot be computed for all components. Final model select ", h," components instead of ", n.comp,"."))
    }
    #update all values
    last.pls$X$weightings <- last.pls$X$weightings[,1:h,drop=F]
    last.pls$X$W.star = last.pls$X$W.star[,1:h,drop=F]
    last.pls$X$loadings = last.pls$X$loadings[,1:h,drop=F]
    last.pls$X$scores = last.pls$X$scores[,1:h,drop=F]
    last.pls$Y$weightings = last.pls$Y$weightings[,1:h,drop=F]
    last.pls$Y$loadings = last.pls$Y$loadings[,1:h,drop=F]
    last.pls$Y$scores = last.pls$Y$scores[,1:h,drop=F]
    last.pls$Y$ratio = last.pls$Y$ratio[,1:h,drop=F]
    var_by_component_nzv = var_by_component_nzv[1:h] #variables selected for each component
    names(var_by_component_nzv) <- paste0("comp_", 1:h)
    last.pls$B = last.pls$B[,,drop=F] #only final coefficients

    n.comp_used <- h
  }

  #or if we filter some components
  if(h != length(names(cox_model$fit$coefficients))){
    if(verbose){
      message(paste0("Updating vectors. Final model select ", length(names(cox_model$fit$coefficients))," components instead of ", n.comp,"."))
    }
    #update all values
    which_to_keep <- which(colnames(last.pls$X$weightings) %in% names(cox_model$fit$coefficients))

    last.pls$X$weightings <- last.pls$X$weightings[,names(cox_model$fit$coefficients),drop=F]
    last.pls$X$W.star = last.pls$X$W.star[,names(cox_model$fit$coefficients),drop=F]
    last.pls$X$loadings = last.pls$X$loadings[,names(cox_model$fit$coefficients),drop=F]
    last.pls$X$scores = last.pls$X$scores[,names(cox_model$fit$coefficients),drop=F]
    last.pls$Y$weightings = last.pls$Y$weightings[,names(cox_model$fit$coefficients),drop=F]
    last.pls$Y$loadings = last.pls$Y$loadings[,names(cox_model$fit$coefficients),drop=F]
    last.pls$Y$scores = last.pls$Y$scores[,names(cox_model$fit$coefficients),drop=F]
    last.pls$Y$ratio = last.pls$Y$ratio[,names(cox_model$fit$coefficients),drop=F]
    var_by_component_nzv = var_by_component_nzv[which_to_keep] #variables selected for each component
    names(var_by_component_nzv) <- paste0("comp_", which_to_keep)
    last.pls$B = last.pls$B[,,drop=F] #only final coefficients

    n.comp_used <- ncol(max(which_to_keep))
  }

  survival_model <- NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  func_call <- match.call()

  if(!returnData){
    survival_model <- removeInfoSurvivalModel(survival_model)
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(splsdrcox_class(list(X = list("data" = if(returnData) X_norm else NA,
                                      "weightings" = if(returnData) last.pls$X$weightings else NA,
                                      "W.star" = last.pls$X$W.star,
                                      "loadings" = last.pls$X$loadings,
                                      "scores" = last.pls$X$scores,
                                      "E" = if(returnData) E else NA,
                                      "x.mean" = xmeans, "x.sd" = xsds),
                             Y = list("deviance_residuals" = if(returnData) DR_coxph_ori else NA,
                                      "dr.mean" = mu,
                                      "dr.sd" = NULL, #deviance_residuals object already centered
                                      "data" = Yh,
                                      "weightings" = if(returnData) last.pls$Y$weightings else NA,
                                      "loadings" = if(returnData) last.pls$Y$loadings else NA,
                                      "scores" = if(returnData) last.pls$Y$scores else NA,
                                      "ratio" = if(returnData) last.pls$Y$ratio else NA,
                                      "y.mean" = ymeans, "y.sd" = ysds),
                             survival_model = survival_model,
                             eta = eta,
                             n.comp = n.comp_used, #number of components
                             var_by_component = var_by_component_nzv, #variables selected for each component
                             call = if(returnData) func_call else NA,
                             X_input = if(returnData) X_original else NA,
                             Y_input = if(returnData) Y_original else NA,
                             B.hat = last.pls$B,
                             R2 = R2,
                             SCR = SCR,
                             SCT = SCT,
                             alpha = alpha,
                             nsv = removed_variables,
                             nzv = variablesDeleted,
                             nz_coeffvar = variablesDeleted_cvar,
                             class = pkg.env$splsdrcox,
                             time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' sPLS-DRCOX Cross-Validation
#' @description sPLS-DRCOX cross validation model
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation (default: 10).
#' @param eta.list Numeric vector. Vector of penalty values. Penalty for sPLS. If eta = 0 no penalty is applied and 1 maximum penalty (no variables are selected). Equal or greater than 1 cannot be selected (default: seq(0.1,0.9,0.1)).
#' @param n_run Numeric. Number of runs for cross validation (default: 5).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation. A non-significant model is a model with at least one component/variable with a P-Value higher than the alpha cutoff. @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
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
#' @return Instance of class "Coxmos" and model "cv.sPLS-DRCOX".
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha information).
#'
#' \code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
#' \code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.
#'
#' \code{opt.comp}: Optimal component selected by the best_model.
#' \code{opt.eta}: Optimal eta/penalty selected by the best_model.
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
#' cv.splsdrcox_model <- cv.splsdrcox(X, Y, max.ncomp = 10,
#' eta.list = seq(0.1,0.9,0.1), x.center = TRUE, x.scale = TRUE)
#' splsdrcox_model <- splsdrcox(X, Y, n.comp = cv.splsdrcox_model$opt.comp,
#' eta = cv.splsdrcox_model$opt.eta, x.center = TRUE, x.scale = TRUE)
#' }

cv.splsdrcox <- function (X, Y,
                         max.ncomp = 10, eta.list = seq(0.1,0.9,0.1),
                         n_run = 5, k_folds = 10,
                         x.center = TRUE, x.scale = FALSE,
                         remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL, remove_variance_at_fold_level = F,
                         remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                         w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
                         MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                         pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                         MIN_EPV = 5, return_models = F, returnData = F,
                         PARALLEL = F, verbose = F, seed = 123){
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
  params_with_limits <- list("eta.list" = eta.list)
  check_min0_less1_variables(params_with_limits)

  params_with_limits <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                 "w_AIC" = w_AIC, "w_c.index" = w_c.index, "w_AUC" = w_AUC, "w_BRIER" = w_BRIER)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("max.ncomp" = max.ncomp,
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

  #### FIX possible SEQ() problems
  eta.list <- as.character(eta.list)
  eta.list <- as.numeric(eta.list)

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
  #total_models <- 1 * k_folds * n_run * length(eta.list)
  total_models <- max.ncomp * k_folds * n_run * length(eta.list)

  lst_model <- get_HDCOX_models2.0(method = pkg.env$splsdrcox,
                                   X_train = X, Y_train = Y,
                                   lst_X_train = lst_train_indexes, lst_Y_train = lst_train_indexes,
                                   max.ncomp = max.ncomp, eta.list = eta.list, EN.alpha.list = NULL, max.variables = NULL, vector = NULL,
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

  if(all(is.null(comp_model_lst))){
    message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run * length(eta.list)
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst, alpha = alpha,
                                                    max.ncomp = max.ncomp, eta.list = eta.list, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models, verbose = verbose)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AIC = NULL, plot_c_index = NULL, plot_BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
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
                                            X_test = X, Y_test = Y,
                                            lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                            df_results_evals = df_results_evals, times = times,
                                            pred.method = pred.method, pred.attr = pred.attr,
                                            max.ncomp = max.ncomp, eta.list = eta.list, n_run = n_run, k_folds = k_folds,
                                            MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                            w_BRIER = w_BRIER, method.train = pkg.env$splsdrcox, PARALLEL = F, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    total_models <- ifelse(!fast_mode, n_run * max.ncomp * length(eta.list), k_folds * n_run * max.ncomp * length(eta.list))

    #times should be the same for all folds
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC_sPLS(comp_model_lst = comp_model_lst,
                                          X_test = X, Y_test = Y,
                                          lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                          df_results_evals = df_results_evals, times = times,
                                          fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                          max.ncomp = max.ncomp, eta.list = eta.list, n_run = n_run, k_folds = k_folds,
                                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          w_AUC = w_AUC, method.train = pkg.env$splsdrcox, PARALLEL = F, verbose = verbose)

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
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index & df_results_evals_comp[,"eta"]==optimal_eta,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  #### ###
  # PLOT #
  #### ###
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_BRIER = w_BRIER, max.ncomp = max.ncomp, eta.list = eta.list,
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

  df_results_evals$eta <- as.numeric(as.character(df_results_evals$eta))
  df_results_evals_run$eta <- as.numeric(as.character(df_results_evals_run$eta))
  df_results_evals_comp$eta <- as.numeric(as.character(df_results_evals_comp$eta))

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  if(return_models){
    return(cv.splsdrcox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.eta = best_model_info$eta, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class= pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.splsdrcox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.eta = best_model_info$eta, plot_AIC = ggp_AIC, plot_c_index = ggp_c_index, plot_BRIER = ggp_BRIER, plot_AUC = ggp_AUC, class = pkg.env$cv.splsdrcox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }
}

#### ### ### #
# PREDICTION #
#### ### ### #

predict.mixOmixs.pls <- function(object, newdata){


  if (missing(newdata))
    stop_quietly("No new data available.")

  if(!class(object) %in% c("mixo_pls", "mixo_spls"))
    stop_quietly("Object must be a 'mixo_pls' object.")

  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
  n.comp = object$ncomp

  x.center <- attr(X,"scaled:center")
  x.scale <- attr(X,"scaled:scale")

  if(is.null(x.center)){
    x.center <- F
  }
  if(is.null(x.scale)){
    x.scale <- F
  }

  newdata = scale(newdata, center = x.center, scale = x.scale)

  if(ncol(newdata) != p){
    stop_quietly(paste0("'newdata' must be a numeric matrix with ncol = ", p, " or a vector of length = ", p, "."))
  }

  X_ww = object$loadings$X #loading x
  Y_ww = object$loadings$Y #loading y

  X_sco = object$variates$X[,drop=F] #variate x
  Y_sco = object$variates$Y[,drop=F] #variate Y

  pp = object$mat.c #matrix of coefficients from the regression of X / residual matrices X on the X-variates, to be used internally by predict.

  newdata = as.matrix(newdata)
  I=diag(1, nrow=n.comp)
  B.hat = array(0, dim = c(p, q, n.comp)) #beta predictor for new data
  Y.hat = array(0, dim = c(nrow(newdata), q, n.comp))

  #Ww <- X_ww %*% solve(t(pp) %*% X_ww)
  Ww <- X_ww %*% MASS::ginv(t(pp) %*% X_ww)
  B <- Ww %*% t(Y_ww)

  Q <- crossprod(Y, X_sco)
  P <- crossprod(X, X_sco)
  #Ww.mod <- X_ww %*% solve(t(P) %*% X_ww)
  Ww.mod <- X_ww %*% MASS::ginv(t(P) %*% X_ww)

  Ypred = lapply(1 : n.comp, function(x){newdata %*% Ww.mod[,1:x] %*% t(Q)[1:x,]})

  aux = NULL
  for(c in 1:length(Ypred)){
    aux = cbind(aux,Ypred[[c]])
  }

  Ypred = aux

  y.center <- attr(Y,"scaled:center")
  y.scale <- attr(Y,"scaled:scale")

  if(is.null(y.center)){
    y.center <- 0
  }
  if(is.null(y.scale)){
    y.scale <- 1
  }

  predict = Ypred * y.scale + y.center

  return(list(predict = predict, B.hat = B))
}

### ## ##
# CLASS #
### ## ##

splsdrcox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$splsdrcox)
  return(model)
}

cv.splsdrcox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.splsdrcox)
  return(model)
}

### ###
# PLS #
### ###

# NA VALUES AS 0, then NA again
pls2 <- function(X, Y, n.comp, x.center = T, x.scale = F, y.center = T, y.scale = F, it = 500, tol = 1e-20, tol.W.star = 1e-20, verbose = F){

  if(n.comp >= nrow(X)) {
    n.comp <- qr(X)$rank-1
  }

  if(is.data.frame(X)) X <- as.matrix(X)
  if(is.data.frame(Y)) Y <- as.matrix(Y)

  xmeans = NULL
  xsds = NULL
  ymeans = NULL
  ysds = NULL

  # Aready centered?
  if(is.null(attr(X, "scaled:center")) | !is.null(attr(X, "scaled:scale"))){
    # Centering and/or scaling
    if(x.center | x.scale){
      Xh <- scale(X, center = x.center, scale = x.scale)
      if(x.center) xmeans <- attr(Xh, "scaled:center")
      if(x.scale) xsds <- attr(Xh, "scaled:scale")
    }else{
      Xh <- X
    }
  }else{
    Xh <- X
    xmeans <- attr(Xh, "scaled:center")
    xsds <- attr(Xh, "scaled:scale")
  }

  if(is.null(attr(Y, "scaled:scale")) |  !is.null(attr(Y, "scaled:scale"))){
    if(y.center | y.scale){
      Yh <- scale(Y, center = y.center, scale = y.scale)
      if(y.center) ymeans <- attr(Yh, "scaled:center")
      if(y.scale) ysds <- attr(Yh, "scaled:scale")
    }else{
      Yh <- Y
    }
  }else{
    Yh <- Y
    ymeans <- attr(Yh, "scaled:center")
    ysds <- attr(Yh, "scaled:scale")
  }

  Yh_ori <- Yh

  Ts <- NULL
  W <- NULL
  Q <- NULL
  U <- NULL
  P <- NULL
  R <- NULL
  C <- NULL
  W <- NULL
  Wnorm <- NULL
  Cnorm <- NULL

  XXNA <- is.na(Xh)
  Xh[XXNA] = 0

  for (h in 1:n.comp) { #for component h to a
    perm <- 0
    nr <- 0
    #1. Select the first u vector with maximum variance in Y
    uh = Yh[,max(which(apply(Yh,2,var) == max(apply(Yh,2,var))))] #select only 1 column with max var in Y
    ende <- FALSE
    while (!ende & nr <= it) { #firs iteration
      nr <- nr + 1 #iteration counts

      #2. wh = uh'Xh / uh'uh
      wh_num <- t(uh) %*% Xh #if any NA, it will be NA
      wh_den <- as.vector(t(uh) %*% uh)
      wh <- t(wh_num / wh_den)

      #3. whn = w' / ||w'||
      whn <- wh / as.vector(sqrt(sum(wh^2, na.rm = T)))

      #4. t = Xh wh / wh'wh
      #4. t = Xh whn
      th <- Xh %*% whn

      #5. c = t'Yh / t't
      #5. q = t'Yh / t't
      ch_num <- t(th) %*% Yh
      ch_den <- as.vector(t(th) %*% th)
      ch <- t(ch_num/ch_den)

      chn <- ch / as.vector(sqrt(sum(ch^2, na.rm = T)))

      #6. uhnew = Yhc / c'c
      #6. uhnew = Yhq / q'q
      uhnew_num <- as.numeric(Yh %*% ch)
      uhnew_den <- as.numeric(t(ch) %*% ch)
      uhnew <- uhnew_num / uhnew_den

      deltau <- uhnew - uh
      unorm <- sqrt(sum(deltau^2))

      if (unorm < tol) {
        ende <- TRUE
      }
      uh <- uhnew
    }

    #7. ph <- t'X / t't
    ph <- t(t(th) %*% Xh / as.vector(t(th) %*% th))
    qh <- t(t(uh) %*% Yh / as.vector(t(uh) %*% uh))

    rh <- t(uh) %*% th/as.vector(t(th) %*% th)

    Xh <- Xh - th %*% t(ph)
    Yh <- Yh - th %*% t(ch) * as.vector(rh)

    Ts <- cbind(Ts, th)
    Q <- cbind(Q, qh)
    U <- cbind(U, uh)
    P <- cbind(P, ph)
    R <- cbind(R, rh)
    C <- cbind(C, ch)
    W <- cbind(W, wh)
    Wnorm <- cbind(Wnorm, whn)
    Cnorm <- cbind(Cnorm, chn)

    if(is.null(P) | is.null(W)){
      message(paste0(pkg.env$splsdrcox, " model cannot be computed because P or W vectors are NULL. Returning NA."))
      # invisible(gc())
      return(NA)
    }

    #system is computationally singular: reciprocal condition number = 6.24697e-18
    # PW <- tryCatch(expr = {solve(t(P) %*% W, tol = tol.W.star)},
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
      message(paste0(pkg.env$splsdrcox, " model cannot be computed due to ginv(t(P) %*% W). Multicollineality could be present in your data. Returning NA."))
      # invisible(gc())
      return(NA)
    }

    # What happen when you cannot compute W.star but you have P and W?
    W.star <- W %*% PW
  }

  Xh[XXNA] <- NA

  # B <- lapply(1:n.comp, function(x){W.star[,1:x,drop=F] %*% t(C[,1:x,drop=F]) %*% R[,1:x,drop=F]})
  # B <- lapply(1:n.comp, function(x){W.star[,1:x,drop=F] %*% t(C[,1:x,drop=F])})

  B <- W.star %*% t(C)
  Brc <- W.star %*% apply(C, 1, function(x){R * x}) #apply gets the transpose matrix
  Bq <- W.star %*% t(Q)

  # aux <- NULL
  # for(c in 1:length(B)){
  #   aux <- cbind(aux, B[[c]][,c,drop=F])
  # }
  # B <- aux

  colnames(W) <- paste0("comp_",1:n.comp)
  colnames(Wnorm) <- paste0("comp_",1:n.comp)
  colnames(W.star) <- paste0("comp_",1:n.comp)
  colnames(P) <- paste0("comp_",1:n.comp)
  colnames(Ts) <- paste0("comp_",1:n.comp)
  colnames(C) <- paste0("comp_",1:n.comp)
  colnames(Cnorm) <- paste0("comp_",1:n.comp)
  colnames(Q) <- paste0("comp_",1:n.comp)
  colnames(U) <- paste0("comp_",1:n.comp)
  colnames(R) <- paste0("comp_",1:n.comp)
  colnames(B) <- colnames(Brc) <- colnames(Bq) <- "coefficient"

  # colnames(B) <- paste0("comp_",1:n.comp)
  # colnames(Brc) <- paste0("comp_",1:n.comp)
  # colnames(Bq) <- paste0("comp_",1:n.comp)

  func_call <- match.call()

  list(X = list("data" = Xh,
                "weightings" = W,
                "weightings_norm" = Wnorm,
                "W.star" = W.star,
                "loadings" = P,
                "scores" = Ts,
                "x.mean" = xmeans,
                "x.sd" = xsds),
       Y = list("data" = Yh_ori,
                "weightings" = C,
                "weightings_norm" = Cnorm,
                "loadings" = Q,
                "scores" = U,
                "ratio" = R,
                "y.mean" = ymeans,
                "y.sd" = ysds),
       B = B,
       Bq = Bq,
       Brc = Brc,
       n.comp = n.comp, #number of components
       call = func_call,
       X_input = X,
       Y_input = Y)
}
