#### ### ##
# METHODS #
#### ### ##

#' sPLS-DRCOX
#' @description Performs a sPLS-DRCOX model (based on plsRcox R package idea).
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param n.comp Numeric. Number of principal components to compute in the PLS model.
#' @param eta Numeric [0-1). Penalty for sPLS. Mean 0 no penalty and 1 maximum penalty. Greater than 1 cannot be selected.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param returnData Logical. Return original and normalized X and Y matrices.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "sPLS-DRCOX". The class contains the following elements:
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
#' \code{var_by_component}: Variables selected in each PLS component.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{B.hat}: PLS beta matrix
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

splsdrcox <- function (X, Y,
                      n.comp = 2, eta = 0.5,
                      x.center = TRUE, x.scale = FALSE,
                      y.center = FALSE, y.scale = FALSE,
                      remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
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

  XXNA <- is.na(Xh) #T is NA
  YNA <- is.na(Y) #T is NA

  ####MAX PREDICTORS
  n.comp <- check.maxPredictors(X, Y, MIN_EPV, n.comp)

  ##############################################
  ######             sPLS-COX             ######
  ##############################################

  #2. Surv function - NULL model
  coxDR <- survival::coxph(survival::Surv(time = time, event = event, type = "right") ~ 1, as.data.frame(Xh))

  #3. Residuals - Default is deviance because eval type="deviance"
  DR_coxph <- residuals(coxDR, type = "deviance") #"martingale", "deviance", "score", "schoenfeld", "dfbeta"', "dfbetas", "scaledsch" and "partial"

  ################################################
  ################################################
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  ################################################
  ################################################

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

    ##############
    # PREDICTION #
    ##############

    #But always using the complete residuals
    # plsfit3 <- ropls::opls(x = Xa, y = DR_coxph_ori, predI = min(h, ncol(Xa)), scaleC = "none")
    # plsfit2 <- mixOmics::pls(X = Xa, Y = DR_coxph_ori, ncomp = min(h, ncol(Xa)), scale = F)
    plsfit <- tryCatch(
      # Specifying expression
      expr = {
        pls2(X = Xa, Y = DR_coxph_ori, n.comp = min(h, ncol(Xa)),
             x.center = F, x.scale = F, y.center = F, y.scale = F,
             it = 10, tol = 1e-10)
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

    ####################
    # UPDATING RESULTS #
    ####################
    predplsfit_list[[h]] <- predplsfit

    beta_pls <- matrix(0, n_var, n_dr)
    rownames(beta_pls) <- colnames(Xh)
    #beta_pls[A_nzv,] <- matrix(data = plsfit$B[,min(h, ncol(Xa)),drop=F], nrow = length(A_nzv), ncol = n_dr) #n.comp from new pls and not all
    beta_pls[A_nzv,] <- matrix(data = plsfit$B[,,drop=F], nrow = length(A_nzv), ncol = n_dr) #n.comp from new pls and not all
    beta_matrix <- cbind(beta_matrix, beta_pls) #res

    var_by_component[[h]] <- A
    var_by_component_nzv[[h]] <- A_nzv

    ###################
    # UPDATING VALUES #
    ###################
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

  ##############################################
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  ##############################################

  ##############################################
  ######              PLS-COX            ######
  ##############################################
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
                      model=T)
    },
    # Specifying error message
    error = function(e){
      if(verbose){
        message(e)
      }
      invisible(gc())
      return(NA)
    }
  )

  while(all(is.na(aux)) & h>0){
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
          message(e)
        }
        invisible(gc())
        return(NA)
      }
    )
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

  survival_model <- NULL
  if(!length(cox_model$fit) == 1){
    survival_model <- getInfoCoxModel(cox_model$fit)
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(splsdrcox_class(list(X = list("data" = if(returnData) X_norm else NA,
                                      "weightings" = last.pls$X$weightings,
                                      "W.star" = last.pls$X$W.star,
                                      "loadings" = last.pls$X$loadings,
                                      "scores" = last.pls$X$scores,
                                      "E" = E,
                                      "x.mean" = xmeans,
                                      "x.sd" = xsds),
                             Y = list("deviance_residuals" = if(returnData) DR_coxph_ori else NA,
                                      "dr.mean" = mu,
                                      "dr.sd" = NULL, #deviance_residuals object already centered
                                      "data" = Yh,
                                      "weightings" = last.pls$Y$weightings,
                                      "loadings" = last.pls$Y$loadings,
                                      "scores" = last.pls$Y$scores,
                                      "ratio" = last.pls$Y$ratio,
                                      "y.mean" = ymeans,
                                      "y.sd" = ysds),
                             survival_model = survival_model,
                             eta = eta,
                             n.comp = n.comp_used, #number of components
                             var_by_component = var_by_component_nzv, #variables selected for each component
                             call = func_call,
                             X_input = if(returnData) X_original else NA,
                             Y_input = if(returnData) Y_original else NA,
                             B.hat = last.pls$B,
                             R2 = R2,
                             SCR = SCR,
                             SCT = SCT,
                             nzv = variablesDeleted,
                             time = time)))
}

splsdrcox.modelPerComponent <- function (X, Y,
                                        n.comp = 4, eta = 0.5,
                                        x.center = TRUE, x.scale = FALSE,
                                        y.center = FALSE, y.scale = FALSE,
                                        remove_near_zero_variance = T, toKeep.zv = NULL,
                                        returnData = T, remove_zero_variance = F, verbose = F){

  #### REQUIREMENTS
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  checkY.colnames(Y)

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

  XXNA <- is.na(Xh) #T is NA
  YNA <- is.na(Y) #T is NA

  ##############################################
  ######             sPLS-COX             ######
  ##############################################

  #2. Surv function - NULL model
  coxDR <- survival::coxph(survival::Surv(time = time, event = event, type = "right") ~ 1, as.data.frame(Xh))

  #3. Residuals - Default is deviance because eval type="deviance"
  DR_coxph <- residuals(coxDR, type = "deviance") #"martingale", "deviance", "score", "schoenfeld", "dfbeta"', "dfbetas", "scaledsch" and "partial"

  ################################################
  ################################################
  ##                                            ##
  ##  Beginning of the loop for the components  ##
  ##                                            ##
  ################################################
  ################################################

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
  lst_HDcox_spls <- list()

  beta_pls <- matrix(0, n_var, n_dr)
  beta_matrix <- list()
  var_by_component <- list()
  var_by_component_nzv <- list()
  plsfit_list <- list()
  predplsfit_list <- list()
  Z_list <- list()
  ww_list <- list()

  lst_pls <- NULL
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

    ##############
    # PREDICTION #
    ##############

    #But always using the complete residuals
    #plsfit2 <- mixOmics::pls(Xa, DR_coxph_ori, ncomp = min(h, ncol(Xa)), scale = F)
    plsfit <- pls2(X = Xa, Y = DR_coxph_ori, n.comp = min(h, ncol(Xa)),
                   x.center = F, x.scale = F, y.center = F, y.scale = F,
                   it = 1, tol = 1e-10)
    plsfit_list[[h]] <- plsfit

    XaNA <- is.na(Xa) #T is NA
    Xa[XaNA] <- 0
    predplsfit <- Xa[,rownames(plsfit$X$loadings),drop=F] %*% plsfit$B[,,drop=F]
    Xa[XaNA] <- NA

    ####################
    # UPDATING RESULTS #
    ####################
    predplsfit_list[[h]] <- predplsfit

    beta_pls <- matrix(0, n_var, n_dr)
    rownames(beta_pls) <- colnames(Xh)
    #beta_pls[A_nzv,] <- matrix(data = plsfit$B[,min(h, ncol(Xa)),drop=F], nrow = length(A_nzv), ncol = n_dr) #n.comp from new pls and not all
    beta_pls[A_nzv,] <- matrix(data = plsfit$B[,,drop=F], nrow = length(A_nzv), ncol = n_dr) #just return the best prediction (all components, so you do not need to select the last one)
    beta_matrix <- cbind(beta_matrix, beta_pls) #res

    var_by_component[[h]] <- A
    var_by_component_nzv[[h]] <- A_nzv

    ###################
    # UPDATING VALUES #
    ###################
    #DR_coxph <- DR_coxph_ori - predplsfit[,min(h, ncol(Xa)),drop=F] #for manual prediction
    DR_coxph <- DR_coxph_ori - predplsfit[,,drop=F] #just return the best prediction (all components, so you do not need to select the last one)

    #R2 calculation
    #E[[h]] = DR_coxph_ori - predplsfit$predict[,,plsfit$n.comp] #same formula, but adding components
    E[[h]] = DR_coxph #same formula, but adding components

    SCR[[h]] = sum(apply(E[[h]],2,function(x) sum(x**2)))
    SCT[[h]] = sum(apply(as.matrix(DR_coxph_ori),2,function(x) sum(x**2))) #equivalent sum((DR_coxph_ori - mean(DR_coxph_ori))**2)
    R2[[h]] = 1 - (SCR[[h]]/SCT[[h]]) #deviance residuals explanation

    last.pls <- plsfit

    ##############################################
    #                                            #
    #      Computation of the coefficients       #
    #      of the model with kk components       #
    #                                            #
    ##############################################

    ##############################################
    ######              PLS-COX            ######
    ##############################################

    n.comp_used <- ncol(last.pls$X$scores)

    d <- as.data.frame(last.pls$X$scores[,,drop=F])
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
                        model=T)
      },
      # Specifying error message
      error = function(e){
        message(e)
        invisible(gc())
        return(NA)
      }
    )

    survival_model <- NULL
    if(!length(cox_model$fit) == 1){
      survival_model <- getInfoCoxModel(cox_model$fit)
    }

    func_call <- match.call()

    lst_HDcox_spls[[h]] <- splsdrcox_class(list(X = list("data" = if(returnData) X_norm else NA,
                                                        "weightings" = last.pls$X$weightings,
                                                        "W.star" = last.pls$X$W.star,
                                                        "loadings" = last.pls$X$loadings,
                                                        "scores" = last.pls$X$scores,
                                                        "E" = E,
                                                        "x.mean" = xmeans,
                                                        "x.sd" = xsds),
                                               Y = list("deviance_residuals" = if(returnData) DR_coxph_ori else NA,
                                                        "dr.mean" = mu,
                                                        "dr.sd" = NULL, #deviance_residuals object already centered
                                                        "data" = Yh,
                                                        "weightings" = last.pls$Y$weightings,
                                                        "loadings" = last.pls$Y$loadings,
                                                        "scores" = last.pls$Y$scores,
                                                        "ratio" = last.pls$Y$ratio,
                                                        "y.mean" = ymeans,
                                                        "y.sd" = ysds),
                                               survival_model = survival_model,
                                               eta = eta,
                                               n.comp = n.comp_used, #number of components
                                               var_by_component = var_by_component_nzv, #variables selected for each component
                                               call = func_call,
                                               X_input = if(returnData) X_original else NA,
                                               Y_input = if(returnData) Y_original else NA,
                                               B.hat = last.pls$B,
                                               R2 = R2,
                                               SCR = SCR,
                                               SCT = SCT))

  }

  names(lst_HDcox_spls) <- paste0("comp_",1:length(lst_HDcox_spls)) #max comp

  invisible(gc())
  return(lst_HDcox_spls) #list of spls_models \ one per component 1:X

}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation sPLS-DRCOX
#' @description sPLS-DRCOX cross validation model
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation.
#' @param eta.list Numeric vector. Vector of penalty values.
#' @param n_run Number. Number of runs for cross validation.
#' @param k_folds Number. Number of folds for cross validation.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation.#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
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
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for perform the runs/folds divisions.
#'
#' @return Instance of class "HDcox" and model "cv.sPLS-DRCOX".
#' @export

cv.splsdrcox <- function (X, Y,
                         max.ncomp = 10, eta.list = seq(0.1,0.9,0.1),
                         n_run = 10, k_folds = 10,
                         x.center = TRUE, x.scale = FALSE,
                         y.center = FALSE, y.scale = FALSE,
                         remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                         remove_non_significant_models = F, alpha = 0.05,
                         w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
                         MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                         pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                         MIN_EPV = 5, return_models = F,
                         PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  ############
  # WARNINGS #
  ############

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

  ####MAX PREDICTORS
  max.ncomp <- check.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)

  #### REQUIREMENTS
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  ######
  # CV #
  ######
  set.seed(seed)
  lst_data <- splitData_Iterations_Folds(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test

  ################
  # TRAIN MODELS #
  ################
  #total_models <- 1 * k_folds * n_run * length(eta.list)
  total_models <- max.ncomp * k_folds * n_run * length(eta.list)

  lst_model <- get_HDCOX_models2.0(method = pkg.env$splsdrcox,
                                lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                max.ncomp = max.ncomp, eta.list = eta.list, EN.alpha.list = NULL,
                                n_run = n_run, k_folds = k_folds,
                                x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
                                remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL,
                                total_models = total_models, PARALLEL = PARALLEL, verbose = verbose)

  # lst_model <- get_HDCOX_models(method = pkg.env$splsdrcox,
  #                               lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
  #                               max.ncomp = max.ncomp, eta.list = eta.list, EN.alpha.list = NULL,
  #                               n_run = n_run, k_folds = k_folds,
  #                               x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
  #                               total_models = total_models)

  comp_model_lst = lst_model$comp_model_lst
  info = lst_model$info

  ##########################
  # BEST MODEL FOR CV DATA #
  ##########################
  total_models <- max.ncomp * k_folds * n_run * length(eta.list)
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst,
                                                    max.ncomp = max.ncomp, eta.list = eta.list, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, time = time)))
    }else{
      return(cv.splsdrcox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = NULL, opt.comp = NULL, opt.eta = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, time = time)))
    }
  }

  ##################
  # EVALUATING AUC #
  ##################
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_eta_index <- NULL
  optimal_eta <- NULL
  optimal_comp_flag <- NULL

  if(w_AUC!=0){
    total_models <- ifelse(!fast_mode, n_run * max.ncomp * length(eta.list), k_folds * n_run * max.ncomp * length(eta.list))
    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC_sPLS(comp_model_lst = comp_model_lst,
                                          lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                          df_results_evals = df_results_evals, times = times,
                                          fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                          max.ncomp = max.ncomp, eta.list = eta.list, n_run = n_run, k_folds = k_folds,
                                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          w_AUC = w_AUC, total_models = total_models, method.train = "spls", PARALLEL = F)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
    optimal_comp_index <- lst_df$optimal_comp_index
    optimal_comp_flag <- lst_df$optimal_comp_flag
    optimal_eta <- lst_df$optimal_eta
    optimal_eta_index <- lst_df$optimal_eta_index
  }else{
    df_results_evals_fold <- df_results_evals
  }

  ##############
  # BEST MODEL #
  ##############

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index & df_results_evals_comp[,"eta"]==optimal_eta,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  ########
  # PLOT #
  ########
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, max.ncomp = max.ncomp, eta.list = eta.list,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", x.text = "Component")

  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp

  ##########
  # RETURN #
  ##########

  df_results_evals$eta <- as.numeric(as.character(df_results_evals$eta))
  df_results_evals_run$eta <- as.numeric(as.character(df_results_evals_run$eta))
  df_results_evals_comp$eta <- as.numeric(as.character(df_results_evals_comp$eta))

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.splsdrcox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.eta = best_model_info$eta, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, time = time)))
  }else{
    return(cv.splsdrcox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.eta = best_model_info$eta, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, time = time)))
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

  Ww <- X_ww %*% solve(t(pp) %*% X_ww)
  B <- Ww %*% t(Y_ww)

  Q <- crossprod(Y, X_sco)
  P <- crossprod(X, X_sco)
  Ww.mod <- X_ww %*% solve(t(P) %*% X_ww)

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
pls2 <- function(X, Y, n.comp, x.center = T, x.scale = F, y.center = T, y.scale = F, it = 500, tol = 1e-20){

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

    W.star <- W %*% solve(t(P) %*% W, tol = 1e-20)
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
