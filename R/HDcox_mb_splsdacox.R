#### ### ##
# METHODS #
#### ### ##

#' mb.splsdacox
#' @description Performs a mb.splsdacox model.
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param n.comp Numeric. Number of principal components to compute in the PLS model.
#' @param vector Numeric vector. Used for computing best number of variables. If NULL, an automatic detection is perform.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param MIN_NVAR Numeric If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param MAX_NVAR Numeric If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param n.cut_points Numeric. Number of start cut points for look the optimal number of variable. 2 cut points mean start with the minimum and maximum. 3 start with minimum, maximum and middle point...(default: 3)
#' @param MIN_AUC_INCREASE Numeric If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param EVAL_METHOD Numeric. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param pred.method Character. AUC method for evaluation. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC")
#' @param max.iter Maximum number of iterations for PLS convergence.
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param returnData Logical. Return original and normalized X and Y matrices.
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "MB.sPLS-DACOX". The class contains the following elements:
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
#' \code{mb.model}: List of splsdacox_mixOmics models computed for each block.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{n.varX}: Number of variables selected for each block.
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
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export

mb.splsdacox <- function (X, Y,
                          n.comp = 4, vector = NULL,
                          x.center = TRUE, x.scale = FALSE,
                          y.center = FALSE, y.scale = FALSE,
                          remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                          remove_non_significant = T,
                          MIN_NVAR = 10, MAX_NVAR = 10000, n.cut_points = 5,
                          MIN_AUC_INCREASE = 0.01,
                          EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                          MIN_EPV = 5, returnData = T, PARALLEL = F, verbose = F){

  t1 <- Sys.time()

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

  E <- list()
  R2 <- list()
  SCR <- list()
  SCT <- list()

  XXNA <- purrr::map(Xh, ~is.na(.)) #T is NA
  YNA <- is.na(Y) #T is NA

  #### ### ### ### ### ### ### ### ### ### ### ##
  ### ##         MB.sPLSDA-COX             ### ##
  #### ### ### ### ### ### ### ### ### ### ### ##

  # set up a full design where every block is connected
  design = matrix(1, ncol = length(Xh), nrow = length(Xh),
                  dimnames = list(c(names(Xh)), c(names(Xh))))
  diag(design) =  0

  #### ### ### ### ### ### ### ### ### ###
  # DIVIDE Y VENCERAS - BEST VECTOR SIZE #
  #### ### ### ### ### ### ### ### ### ###

  DR_coxph = NULL #not used in plsda

  if(is.null(vector)){
    keepX <- getBestVectorMB(Xh, DR_coxph, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                             EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F, mode = "splsda", verbose = verbose)
  }else{
    if(isa(vector, "list")){
      keepX <- vector
      #if list, but not n.comp length... and just one value in each block
      if(!all(unlist(purrr::map(keepX, ~length(.)==n.comp))) & all(unlist(purrr::map(keepX, ~length(.)==1)))){
        keepX <- purrr::map(keepX, ~rep(., n.comp))
      }else if(!all(unlist(purrr::map(keepX, ~length(.)==1)))){
        #more than one value... just take the first one
        keepX <- purrr::map(keepX, ~rep(.[[1]], n.comp))
      }
    }else{
      #vector is the same length of blocks in X (so each value correspond to each vector)
      if(length(vector)==length(X)){
        keepX <- list()
        for(e in 1:length(vector)){
          keepX[[e]] <- rep(vector[[e]], n.comp)
        }
        names(keepX) <- names(X)
      }else{
        message("Vector does not has the proper structure. Optimizing best n. variables by using your vector as start vector.")
        keepX <- getBestVectorMB(Xh, DR_coxph, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                                 EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = F, mode = "splsda", verbose = verbose)
      }
    }
  }

  mb.splsda <- mixOmics::block.splsda(Xh, Yh[,"event"], scale=F, ncomp = n.comp, keepX = keepX, max.iter = max.iter, near.zero.var = F, all.outputs = T)

  #PREDICTION
  predplsfit <- predict(mb.splsda, newdata=Xh)

  for(block in names(predplsfit$predict)){
    E[[block]] <- list()
    SCR[[block]] <- list()
    SCT[[block]] <- list()
    for(h in 1:n.comp){
      E[[block]][[h]] <- Yh[,"event"] - predplsfit$predict[[block]][,,h]

      SCR[[block]][[h]] = sum(apply(E[[block]][[h]],2,function(x) sum(x**2)))
      SCT[[block]][[h]] = sum(apply(as.matrix(Yh[,"event"]),2,function(x) sum(x**2))) #equivalent sum((Yh[,"event"] - mean(Yh[,"event"]))**2)

      R2[[block]][[h]] = 1 - (SCR[[block]][[h]]/SCT[[block]][[h]]) #deviance residuals explanation
    }
    R2[[block]] = mb.splsda$prop_expl_var[[block]]
  }

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ### ###
  ### ###         MB:sPLSDA-COX            ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###
  n.comp_used <- ncol(mb.splsda$variates$Y) #can be lesser than expected because we have lesser variables to select because penalization

  n.varX_used <- list()
  for(i in names(Xh)){
    aux <- list()
    for(j in 1:n.comp){
      aux[[j]] <- rownames(mb.splsda$loadings[[i]][which(mb.splsda$loadings[[i]][,j]!=0),j,drop=F])
    }
    names(aux) <- colnames(mb.splsda$loadings[[i]])
    n.varX_used[[i]] <- aux
  }

  data <- as.data.frame(mb.splsda$variates[[1]][,,drop=F])
  for(b in names(Xh)[2:length(Xh)]){
    data <- cbind(data, as.data.frame(mb.splsda$variates[[b]][,,drop=F]))
  }

  update_colnames <- paste0("comp_", 1:ncol(mb.splsda$variates[[1]]))
  colnames(data) <- apply(expand.grid(update_colnames, names(Xh)), 1, paste, collapse="_")
  cox_model <- cox(X = data, Y = Yh, x.center = F, x.scale = F, y.center = F, y.scale = F, remove_non_significant = remove_non_significant, FORCE = T)

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  removed_variables <- NULL
  if(remove_non_significant){
    lst_rnsc <- removeNonSignificativeCox(cox = cox_model$fit, alpha = alpha, cox_input = cbind(data, Yh))

    cox_model$fit <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  survival_model <- cox_model$survival_model

  #get W.star
  Tmat <- Pmat <- Cmat <- Wmat <- W.star <- B.hat <- list()
  for(i in 1:length(Xh)){
    #select just features != 0 (selected features)
    names <- purrr::map(1:n.comp_used, ~rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,.,drop=F]!=0)])
    all_names <- unique(unlist(names))

    aux_Pmat = matrix(data = 0, nrow = ncol(Xh[[i]]), ncol = n.comp_used)
    rownames(aux_Pmat) <- colnames(Xh[[i]])
    colnames(aux_Pmat) <- colnames(mb.splsda$loadings[[i]])

    for(c in 1:n.comp_used){
      names <- rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,c,drop=F]!=0)]
      aux <- crossprod(Xh[[i]][,names,drop=F], mb.splsda$variates[[i]][,c])
      aux_Pmat[names,c] = aux
    }

    Pmat[[i]] = aux_Pmat
    Cmat[[i]] = crossprod(Yh[,"event"], mb.splsda$variates[[i]])
    Wmat[[i]] = mb.splsda$loadings[[i]]
    Tmat[[i]] = mb.splsda$variates[[i]]

    colnames(Wmat[[i]]) <- paste0("comp_", 1:ncol(Wmat[[i]]))
    colnames(Pmat[[i]]) <- paste0("comp_", 1:ncol(Pmat[[i]]))
    colnames(Tmat[[i]]) <- paste0("comp_", 1:ncol(Tmat[[i]]))

    # W.star[[i]] <- lapply(1:n.comp, function(x){Wmat[[i]][,1:x,drop=F] %*% solve(t(Pmat[[i]][,1:x,drop=F]) %*% Wmat[[i]][, 1:x,drop=F])})
    # B.hat[[i]] <- lapply(1:n.comp, function(x){W.star[[i]][[x]][,1:x,drop=F] %*% t(Cmat[[i]][,1:x,drop=F])})

    aux_W.star = matrix(data = 0, nrow = ncol(Xh[[i]]), ncol = n.comp_used)
    rownames(aux_W.star) <- colnames(Xh[[i]])
    colnames(aux_W.star) <- colnames(mb.splsda$loadings[[i]])

    for(c in 1:n.comp_used){
      names <- rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,c,drop=F]!=0)]
      aux <- Wmat[[i]][names,c,drop=F] %*% solve(t(Pmat[[i]][names,c,drop=F]) %*% Wmat[[i]][names,c,drop=F])
      aux_W.star[names,c] = aux
    }

    W.star[[i]] <- aux_W.star
    B.hat[[i]] <- W.star[[i]] %*% t(Cmat[[i]][,1:n.comp_used,drop=F])

    colnames(W.star[[i]]) <- paste0("comp_", 1:ncol(W.star[[i]]))
  }

  # #get W.star
  # Tmat <- Pmat <- Cmat <- Wmat <- W.star <- B.hat <- list()
  # for(i in 1:length(Xh)){
  #   Pmat[[i]] = crossprod(Xh[[i]], mb.splsda$variates[[i]])
  #   Cmat[[i]] = crossprod(Yh[,"event"], mb.splsda$variates[[i]])
  #   Wmat[[i]] = mb.splsda$loadings[[i]]
  #   Tmat[[i]] = mb.splsda$variates[[i]]
  #
  #   colnames(Wmat[[i]]) <- paste0("comp_", 1:ncol(Wmat[[i]]))
  #   colnames(Pmat[[i]]) <- paste0("comp_", 1:ncol(Pmat[[i]]))
  #   colnames(Tmat[[i]]) <- paste0("comp_", 1:ncol(Tmat[[i]]))
  #
  #   # W.star[[i]] <- lapply(1:n.comp, function(x){Wmat[[i]][,1:x,drop=F] %*% solve(t(Pmat[[i]][,1:x,drop=F]) %*% Wmat[[i]][, 1:x,drop=F])})
  #   # B.hat[[i]] <- lapply(1:n.comp, function(x){W.star[[i]][[x]][,1:x,drop=F] %*% t(Cmat[[i]][,1:x,drop=F])})
  #
  #   W.star[[i]] <- Wmat[[i]][,1:n.comp_used,drop=F] %*% solve(t(Pmat[[i]][,1:n.comp_used,drop=F]) %*% Wmat[[i]][, 1:n.comp_used,drop=F])
  #   B.hat[[i]] <- W.star[[i]] %*% t(Cmat[[i]][,1:n.comp_used,drop=F])
  # }

  names(Pmat) <- names(Xh)
  names(Cmat) <- names(Xh)
  names(Wmat) <- names(Xh)
  names(Tmat) <- names(Xh)
  names(W.star) <- names(Xh)
  names(B.hat) <- names(Xh)

  #MIX Omics, a la hora de generar los nuevos scores para nuevas X (o las mismas de entrenamiento),
  #a parte de realizar la multiplicacion X*W.STAR, realiza luego una normalización de los scores en base a la norma de la propia X usada,
  #de esa manera, en el multiblock de SPLS los resultados no coinciden con los de la funcion predict de MIXOMICS. La siguiente linea es
  #la que se ejecuta una vez realizado el calculo de los nuevos SCORES.

  # head(predplsfit$variates$genes)
  # head(mb.splsda$X$genes %*% W.star[[1]][[n.comp]])
  # head(mb.splsda$X$genes %*% W.star[[1]][[n.comp]])
  # head(mb.splsda$X$genes %*% Wmat[[1]] %*% solve(t(Pmat[[1]]) %*% Wmat[[1]]))
  # #
  # Pmat[[1]] = crossprod(Xh$genes, tt_mbsplsDR[[1]])
  # Wmat[[1]] = mb.splsda$AVE$AVE_inner
  # head(mb.splsda$X$genes %*% Wmat[[1]] %*% solve(t(Pmat[[1]]) %*% Wmat[[1]]))
  #
  # new_t <- mb.splsda$X$genes %*% W.star$genes[[n.comp_used]]
  # new_t2 <- matrix(data = sapply(1:ncol(new_t),
  #                                function(x) {new_t[, x] * apply(mb.splsda$variates$genes, 2,
  #                                                                function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(Xh$genes), ncol = ncol(new_t))
  # head(new_t2)

  # Si lo aplicamos a SPLS normal, también falla el cáclulo de la W*. Puede ser que sea debido a que los cálculos de
  # los loadings de X se estén realizando con la normalización de la C y por tanto la corrección de la norma soluciona el problema.
  # Sin embargo, hubiera sido más sencillo trabajar directamente con una metodología correcta. En mi caso, si utilizo mixomics, debo usar
  # su función siempre para predecir los scores de las nuevas X y NO LO ESTOY HACIENDO!

  #get W.star
  W <- Wmat
  P <- Pmat
  W.star <- W.star
  B.hat <- B.hat
  Ts <- Tmat

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(mb.splsdacox_class(list(X = list("data" = if(returnData) X_norm else NA, "loadings" = P, "weightings" = W, "W.star" = W.star, "scores" = Ts, "E" = E, "x.mean" = xmeans, "x.sd" = xsds),
                                 Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                                 survival_model = survival_model,
                                 mb.model = mb.splsda,
                                 n.comp = n.comp_used, #number of components
                                 n.varX = n.varX_used,
                                 call = func_call,
                                 X_input = if(returnData) X_original else NA,
                                 Y_input = if(returnData) Y_original else NA,
                                 B.hat = B.hat,
                                 R2 = R2,
                                 SCR = SCR,
                                 SCT = SCT,
                                 alpha = alpha,
                                 removed_variables_cox = removed_variables,
                                 nzv = variablesDeleted,
                                 class = pkg.env$mb.splsdacox,
                                 time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation cv.mb.splsdacox
#' @description cv.mb.splsdacox cross validation model
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation.
#' @param vector Numeric vector. Used for computing best number of variables. If NULL, an automatic detection is perform.
#' @param n_run Number. Number of runs for cross validation.
#' @param k_folds Number. Number of folds for cross validation.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE, non-significant models are removed before computing the evaluation. A non-significant model is a model with at least one component/variable with a P-Value higher than the alpha cutoff. @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param MIN_NVAR Numeric If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param MAX_NVAR Numeric If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param n.cut_points Numeric. Number of start cut points for look the optimal number of variable. 2 cut points mean start with the minimum and maximum. 3 start with minimum, maximum and middle point...(default: 3)
#' @param EVAL_METHOD Numeric. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
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
#' @return Instance of class "HDcox" and model "cv.MB.sPLS-DACOX".
#' @export

cv.mb.splsdacox <- function(X, Y,
                            max.ncomp = 10, vector = NULL,
                            n_run = 5, k_folds = 10,
                            x.center = TRUE, x.scale = FALSE,
                            y.center = FALSE, y.scale = FALSE,
                            remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                            remove_non_significant_models = F, remove_non_significant = F, alpha = 0.05,
                            MIN_NVAR = 10, MAX_NVAR = 10000, n.cut_points = 5, EVAL_METHOD = "cenROC",
                            w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
                            MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                            pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                            MIN_EPV = 5, return_models = F,
                            PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  #### ### ###
  # WARNINGS #
  #### ### ###

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_AUC))
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
                                            remove_near_zero_variance = remove_near_zero_variance,
                                            remove_zero_variance = remove_zero_variance,
                                            toKeep.zv = toKeep.zv,
                                            freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  max.ncomp <- check.mb.ncomp(X, max.ncomp)

  #### MAX PREDICTORS
  max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)

  #### #
  # CV #
  #### #
  set.seed(seed)
  lst_data <- splitData_Iterations_Folds.mb(X, Y, n_run = n_run, k_folds = k_folds) #FOR TEST
  lst_X_train <- lst_data$lst_X_train
  lst_Y_train <- lst_data$lst_Y_train
  lst_X_test <- lst_data$lst_X_test
  lst_Y_test <- lst_data$lst_Y_test

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  total_models <- 1 * k_folds * n_run
  #total_models <- max.ncomp * k_folds * n_run

  lst_model <- get_HDCOX_models2.0(method = pkg.env$mb.splsdacox, vector = vector,
                                lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
                                max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL,
                                n_run = n_run, k_folds = k_folds, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
                                remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL,
                                alpha = alpha, MIN_EPV = MIN_EPV,
                                remove_non_significant = remove_non_significant,
                                total_models = total_models, PARALLEL = PARALLEL, verbose = verbose)

  # lst_model <- get_HDCOX_models(method = pkg.env$mb.splsdacox, vector = vector,
  #                               lst_X_train = lst_X_train, lst_Y_train = lst_Y_train,
  #                               max.ncomp = max.ncomp, eta.list = NULL, EN.alpha.list = NULL,
  #                               n_run = n_run, k_folds = k_folds,
  #                               x.center = x.center, x.scale = x.scale, y.center = y.center, y.scale = y.scale,
  #                               total_models = total_models)

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = lst_model,
                                                    max.ncomp = max.ncomp, eta.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = lst_model, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.mb.splsdacox, time = time)))
    }else{
      return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, class = pkg.env$cv.mb.splsdacox, time = time)))
    }
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_eta_index <- NULL
  optimal_eta <- NULL
  optimal_comp_flag <- NULL

  if(w_AUC!=0){
    total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp)

    lst_df <- get_COX_evaluation_AUC(comp_model_lst = lst_model,
                                     lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, total_models = total_models, method.train = pkg.env$mb.splsdacox, PARALLEL = F)

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

  #### ### ### #
  # BEST MODEL #
  #### ### ### #

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index,, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  best_n_var <- list()
  aux_n_var <- as.numeric(strsplit(as.character(best_model_info$n.var), "_")[[1]])
  for(e in 1:length(aux_n_var)){
    best_n_var[[e]] <- aux_n_var[[e]]
  }
  names(best_n_var) <- names(X)

  #### ###
  # PLOT #
  #### ###
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, max.ncomp = max.ncomp, eta.list = NULL,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", x.text = "Component")

  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp

  #### ### #
  # RETURN #
  #### ### #

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.mb.splsdacox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = lst_model, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_n_var, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, class = pkg.env$cv.mb.splsdacox, time = time)))
  }else{
    return(cv.mb.splsdacox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_n_var, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, class = pkg.env$cv.mb.splsdacox, time = time)))
  }
}

### ## ##
# CLASS #
### ## ##

mb.splsdacox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$mb.splsdacox)
  return(model)
}

cv.mb.splsdacox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.mb.splsdacox)
  return(model)
}
