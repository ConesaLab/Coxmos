checkColnamesIllegalChars.mb <- function(X){

  for(block in names(X)){
    new_cn_X <- deleteIllegalChars(colnames(X[[block]]))

    if(length(unique(new_cn_X)) == length(unique(colnames(X[[block]])))){
      colnames(X[[block]]) <- new_cn_X
    }else{
      stop(paste0("When deleting illegal chars, some colnames in X (", block,") get the same name. Update manually the colnames to avoid the next chars: ", paste0(pkg.env$IllegalChars, collapse = " ")))
    }
  }

  return(X)

}

#' getEPV.mb
#' @description Provides a quantitative assessment of the dataset by computing the Events per Variable
#' (EPV) metric for multi-block data, which gauges the proportionality between observed events and the
#' number of explanatory variables.
#'
#' @details In the realm of survival analysis, the balance between observed events and explanatory
#' variables is paramount. The `getEPV` function serves as a tool for researchers to ascertain this
#' balance, which can be pivotal in determining the robustness and interpretability of subsequent
#' statistical models. By evaluating the ratio of events in the `Y` matrix to the variables in the `X`
#' matrix, the function yields the EPV metric. It is of utmost importance that the `Y` matrix encompasses
#' two distinct columns, namely "time" and "event". The latter, "event", should strictly encapsulate
#' binary values, delineating censored (either 0 or FALSE) and event (either 1 or TRUE) observations.
#' To ensure the integrity of the data and the precision of the computation, the function is equipped
#' with an error mechanism that activates if the "event" column remains undetected.
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform
#' into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event
#' observations.
#'
#' @return Return the EPV value for a specific X (explanatory variables) and Y (time and censored variables) data.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' data("Y_multiomic")
#' X <- X_multiomic
#' Y <- Y_multiomic
#' getEPV.mb(X,Y)

getEPV.mb <- function(X,Y){
  EPV_lst <- list()
  for(b in names(X)){
    if("event" %in% colnames(Y)){
      EPV_lst[[b]] <- getEPV(X[[b]], Y)
    }else{
      stop("Column event has not been detected in Y matrix.")
    }
  }

  return(EPV_lst)
}

#' deleteZeroOrNearZeroVariance.mb
#' @description Provides a robust mechanism to filter out variables from a dataset that exhibit zero
#' or near-zero variance, thereby enhancing the quality and interpretability of subsequent statistical
#' analyses.
#'
#' @details The `deleteZeroOrNearZeroVariance` function is an indispensable tool in the preprocessing
#' phase of statistical modeling. In many datasets, especially high-dimensional ones, certain variables
#' might exhibit zero or near-zero variance. Such variables can be problematic as they offer limited
#' information variance and can potentially distort the results of statistical models, leading to
#' issues like overfitting. By leveraging the `caret::nearZeroVar()` function, this tool offers a
#' rigorous method to identify and exclude these variables. Users are afforded flexibility in their
#' choices, with options to remove only zero variance variables, near-zero variance variables, or
#' both. The function also provides the capability to set a frequency cutoff, `freqCut`, which
#' determines the threshold for near-zero variance based on the ratio of the most frequent value to
#' the second most frequent value. For scenarios where certain variables are deemed essential and
#' should not be removed regardless of their variance, the `toKeep.zv` parameter allows users to specify
#' a list of such variables.
#' @param X List of numeric matrices or data.frame. Explanatory variables. Qualitative variables must
#' be transform into binary variables.
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance
#' variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will
#' be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance
#' filtering (default: NULL).
#' @param freqCut Numeric. Cutoff for the ratio of the most common value to the second most common
#' value (default: 95/5).
#'
#' @return A list of two objects.
#' \code{X}: A list with as many blocks as X input, but with the variables filtered.
#' \code{variablesDeleted}: A list with as many blocks as X input, with the name of the variables
#' that have been removed.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' X <- X_multiomic
#' filter <- deleteZeroOrNearZeroVariance.mb(X, remove_near_zero_variance = TRUE)

deleteZeroOrNearZeroVariance.mb <- function(X, remove_near_zero_variance = FALSE,
                                            remove_zero_variance = TRUE,
                                            toKeep.zv = NULL, freqCut = 95/5){

  auxX <- X

  variablesDeleted <- NULL
  if(remove_near_zero_variance){
    lst.zv <- purrr::map(auxX, ~deleteZeroVarianceVariables(data = ., info = TRUE, mustKeep = toKeep.zv, freqCut = freqCut))
    variablesDeleted <- purrr::map(lst.zv, ~.$variablesDeleted[,1])
    if(any(unlist(lapply(variablesDeleted, is.null)))){ #if any not null
      for(n in names(variablesDeleted)){
        if(is.null(variablesDeleted[[n]])){
          next
        }else{
          auxX[[n]] <- auxX[[n]][,!colnames(auxX[[n]]) %in% variablesDeleted[[n]]]
        }
      }
    }
  }else if(remove_zero_variance){
    lst.zv <- purrr::map(auxX, ~deleteZeroVarianceVariables(data = ., info = TRUE, mustKeep = toKeep.zv, onlyZero = TRUE))
    variablesDeleted <- purrr::map(lst.zv, ~.$variablesDeleted[,1])
    if(any(unlist(lapply(variablesDeleted, is.null)))){ #if any not null
      for(n in names(variablesDeleted)){
        if(is.null(variablesDeleted[[n]])){
          next
        }else{
          auxX[[n]] <- auxX[[n]][,!colnames(auxX[[n]]) %in% variablesDeleted[[n]]]
        }
      }
    }
  }

  return(list(X = auxX, variablesDeleted = variablesDeleted))

}

#' deleteNearZeroCoefficientOfVariation.mb
#' @description Filters out variables from a dataset that exhibit a coefficient of variation below a
#' specified threshold, ensuring the retention of variables with meaningful variability.
#'
#' @details The `deleteNearZeroCoefficientOfVariation` function is a pivotal tool in data preprocessing,
#' especially when dealing with high-dimensional datasets. The coefficient of variation (CoV) is a
#' normalized measure of data dispersion, calculated as the ratio of the standard deviation to the mean.
#' In many scientific investigations, variables with a low CoV might be considered as offering limited
#' discriminative information, potentially leading to noise in subsequent statistical analyses. By
#' setting a threshold through the `LIMIT` parameter, this function provides a systematic approach to
#' identify and exclude variables that do not meet the desired variability criteria. The underlying
#' rationale is that variables with a CoV below the set threshold might not contribute significantly
#' to the variability of the dataset and could be redundant or even detrimental for certain analyses.
#' The function returns a modified dataset, a list of deleted variables, and the computed coefficients
#' of variation for each variable. This comprehensive output ensures that researchers are well-informed
#' about the preprocessing steps and can make subsequent analytical decisions with confidence.
#'
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must
#' be transform into binary variables.
#' @param LIMIT Numeric. Cutoff for minimum variation. If coefficient is lesser than the limit, the
#' variables are removed because not vary enough (default: 0.1).
#'
#' @return A list of three objects.
#' \code{X}: A list with as many blocks as X input, but with the variables filtered.
#' \code{variablesDeleted}: A list with as many blocks as X input, with the name of the variables that have been removed.
#' \code{coeff_variation}: A list with as many blocks as X input, with the coefficient of variation per variable.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' X <- X_multiomic
#' filter <- deleteNearZeroCoefficientOfVariation.mb(X, LIMIT = 0.1)
#'
deleteNearZeroCoefficientOfVariation.mb <- function(X, LIMIT = 0.1){
  variablesDeleted <- list()
  coef_list <- list()
  newX <- X
  for(b in names(newX)){
    cvar <- apply(newX[[b]], 2, function(x){sd(x)/abs(mean(x))})
    variablesDeleted[[b]] <- names(cvar)[which(cvar <= LIMIT)]
    if(length(variablesDeleted[[b]])>0){
      newX[[b]] <- newX[[b]][,!colnames(newX[[b]]) %in% variablesDeleted[[b]]]
    }else{
      variablesDeleted[[b]] <- 0
    }
    coef_list[[b]] <- cvar
  }

  variablesDeleted <- purrr::map(variablesDeleted, ~{.==0;NULL})

  return(list("X" = newX, "variablesDeleted" = variablesDeleted, "coeff_variation" = coef_list))
}

checkXY.rownames.mb <- function(X, Y, verbose = TRUE){
  # Check if X and Y are matrices
  if (!isa(X, "list")){
    if(verbose){
      stop("X data is not a list\n")
    }
    X <- data.matrix(X)
  }else{
    for(b in names(X)){
      if(!nrow(X[[b]])==nrow(Y)){
        stop(paste0("X[",b,"]", " and Y have different number of observations."))
      }

      if(is.null(rownames(X[[b]])) & is.null(rownames(Y))){
        if(verbose){
          message(paste0("Rownames of X[",b,"] and Y are NULL. Named from 1 to ", nrow(X[[b]]), "."))
        }
        rownames(Y) <- rownames(X[[b]]) <- 1:nrow(X[[b]])
      }else if(is.null(rownames(X[[b]]))){
        rownames(X[[b]]) <- rownames(Y)
      }else{
        rownames(Y) <- rownames(X[[b]])
      }
    }
  }

  return(list(X = X, Y = Y))
}

checkXY.mb.class <- function(X, Y, verbose = TRUE){
  # Check if X and Y are matrices
  if (!isa(X, "list")){
    if(verbose){
      stop("X data is not a list\n")
    }
    X <- data.matrix(X)
  }else{
    if(all(unlist(lapply(X, class)) %in% "data.frame")){
      X <- lapply(X, data.matrix)
    }else if(all(unlist(lapply(X, function(x){class(x)[[1]]})) %in% c("data.frame", "matrix"))){
      for(block in names(X)){
        if(is.data.frame(X[[block]])[[1]]){
          X[[block]] <- data.matrix(X[[block]])
        }
      }
    }else if(!all(unlist(lapply(X, class)) %in% c("matrix", "array"))){
      stop("X elements are not a matrix or a data.frame\n")
    }
  }

  if (!is.matrix(Y)) {
    if(is.data.frame(Y)){
      if(verbose){
        message("Y data is not a matrix, applying data.matrix\n")
      }
      Y <- data.matrix(Y)
    }else{
      stop("Y data is not a matrix or a data.frame")
    }
  }

  for(i in names(X)){
    if(any(rowSums(is.na(X[[i]]))>0)){
      stop(paste0("Block ", i, " have observations with NAs. MB:PLS models cannot manage NAs yet."))
    }
  }

  return(list(X = X, Y = Y))
}

check.mb.ncomp <- function(X, max.ncomp, verbose = FALSE){

  if(length(max.ncomp)>1){
    stop("max.ncomp must be a number. Not a vector.")
  }

  d <- min(unlist(lapply(X, ncol)))
  if(max.ncomp > d){
    if(verbose){
      message(paste0("Number of components used by all blocks cannot be grater than the number of variables. Updated to ", ncol(X), "."))
    }
    max.ncomp <- d
  }
  return(max.ncomp)
}

check.mb.maxPredictors <- function(X, Y, MIN_EPV, max.variables, verbose = FALSE){
  max_n_predictors <- purrr::map(X, ~getMaxNPredictors(n.var = ncol(.), Y, MIN_EPV))
  max_n_predictors <- max(unlist(max_n_predictors))
  max.variables.res <- max.variables

  if(max_n_predictors==0){
    stop_quietly(paste0("Less than ", MIN_EPV, " events. Stop program."))
  }

  message(paste0("As we are working with a multiblock approach with ", length(X),
                 " blocks, a maximum of ", round(max_n_predictors/length(X)), " components could be use."))

  max_n_predictors = round(max_n_predictors/length(X))

  if(max.variables > max_n_predictors){
    if(verbose){
      message(paste0("As we have a low number of events, we should work with a maximum of ", max_n_predictors,
                     " variables/components in final cox model."))
    }
    max.variables.res = max_n_predictors
  }

  return(max.variables.res)
}

## Divide the data into 4 lists, two train and two test
## Each list has n_run lists of k_folds train/test subsets.
## n_run = 5 and k_folds = 10 -> a total of 50 models
## X must be a list of matrix
## Y must be a data.frame with "time" and "event" columns

splitData_Iterations_Folds.mb <- function(X, Y, n_run, k_folds, seed = 123){

  set.seed(seed)

  if(!is.numeric(n_run) & n_run > 0){
    stop_quietly("Parameter 'n_run' must be a numeric greater or equal than 1.")
  }
  if(!is.numeric(k_folds) & k_folds > 1){
    stop_quietly("Parameter 'k_folds' must be a numeric greater or equal than 2.") #At least two folds for train/test
  }

  if(k_folds > nrow(Y)){
    warning(paste0("Parameter 'k_folds' cannot be greater than the number of observations (changed to ", nrow(Y), ").\n")) #Folds as observation as maximum
    k_folds <- nrow(Y)
  }

  if(!any(c("event", "status") %in% colnames(Y))){
    stop_quietly("Y data.frame must contain the colname 'event' or 'status'.")
  }else if("status" %in% colnames(Y)){
    colnames(Y)[colnames(Y)=="status"] <- "event"
  }

  lst_X_train <- list()
  lst_Y_train <- list()

  lst_X_test <- list()
  lst_Y_test <- list()

  lst_obs_index_train <- list()
  lst_obs_index_test <- list()

  for(i in 1:n_run){
    testIndex <- caret::createFolds(Y[,"event"],
                                     k = k_folds,
                                     list = TRUE)

    #for each fold, take the others as train
    lst_X_data_train_aux <- purrr::map(X, ~lapply(testIndex, function(ind, dat) dat[-ind,], dat = .))
    lst_Y_data_train <- lapply(testIndex, function(ind, dat) dat[-ind,], dat = Y)

    #for each fold, take just one as test
    lst_X_data_test_aux <- purrr::map(X, ~lapply(testIndex, function(ind, dat) dat[ind,], dat = .))
    lst_Y_data_test <- lapply(testIndex, function(ind, dat) dat[ind,], dat = Y)

    #### CHANGE ORDER
    lst_X_data_train <- list()
    lst_X_data_test <- list()
    for(f in names(lst_X_data_train_aux[[1]])){
      lst_X_data_train[[f]] <- list()
      lst_X_data_test[[f]] <- list()
      for(b in names(X)){
        lst_X_data_train[[f]][[b]] <- lst_X_data_train_aux[[b]][[f]]
        lst_X_data_test[[f]][[b]] <- lst_X_data_test_aux[[b]][[f]]
      }
    }

    lst_X_train[[i]] <- lst_X_data_train
    lst_Y_train[[i]] <- lst_Y_data_train

    lst_X_test[[i]] <- lst_X_data_test
    lst_Y_test[[i]] <- lst_Y_data_test

    lst_obs_index_test[[i]] <- testIndex
    aux <- 1:nrow(Y)
    lst_obs_index_train[[i]] <- lapply(testIndex, function(ind, aux) aux[-ind], aux = aux)
  }

  names(lst_X_train) <- paste0("run", 1:n_run)
  names(lst_Y_train) <- paste0("run", 1:n_run)
  names(lst_X_test) <- paste0("run", 1:n_run)
  names(lst_Y_test) <- paste0("run", 1:n_run)

  names(lst_obs_index_train) <- paste0("run", 1:n_run)
  names(lst_obs_index_test) <- paste0("run", 1:n_run)

  return(list(lst_X_train = lst_X_train, lst_Y_train = lst_Y_train, lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, lst_train_index = lst_obs_index_train, lst_test_index = lst_obs_index_test, k_folds = k_folds))
}

XY.mb.scale <- function(X, Y, x.center, x.scale, y.center, y.scale){
  xmeans <- NULL
  xsds <- NULL
  ymeans <- NULL
  ysds <- NULL

  #CENTERING X
  if(length(x.center)>1){
    if(length(x.center)==length(X) & all(names(X) %in% names(x.center))){
      Xh <- purrr::map2(.x = X, .y = names(X), ~scale(., center = x.center[[.y]], scale = FALSE))
    }else{
      stop("Names in x.center do not match with names in X.")
    }
  }else{
    Xh <- purrr::map(X, ~scale(., center = x.center, scale = FALSE))
  }

  #SCALING X
  if(length(x.scale)>1){
    if(length(x.scale)==length(X) & all(names(X) %in% names(x.scale))){
      Xh <- purrr::map2(.x = Xh, .y = names(X), ~scale(., center = FALSE, scale = x.scale[[.y]]))
    }else{
      stop("Names in x.scale do not match with names in X.")
    }
  }else{
    Xh <- purrr::map(Xh, ~scale(., center = FALSE, scale = x.scale))
  }

  xmeans <- purrr::map(Xh, ~attr(., "scaled:center"))
  xsds <- purrr::map(Xh, ~attr(., "scaled:scale"))

  #SCALING Y
  if(y.center | y.scale){
    Yh_time <- scale(Y[,"time", drop = FALSE], center = y.center, scale = y.scale)
    ymeans <- attr(Y, "scaled:center")
    ysds <- attr(Y, "scaled:scale")
    Yh <- Y
    Yh[,"time"] <- Yh_time
  }else{
    Yh <- Y
  }

  return(list(Xh = Xh, xmeans = xmeans, xsds = xsds, Yh = Yh, ymeans = ymeans, ysds = ysds))
}

Xh_XXNA <- function(Xh, XXNA, value = 0){
  for(omic in names(Xh)){
    Xh[[omic]][XXNA[[omic]]] <- value
  }
  return(Xh)
}

getCIndex_AUC_CoxModel_block.spls <- function(Xh, DR_coxph_ori, Yh, n.comp, keepX, scale = FALSE,
                                              near.zero.var = FALSE, EVAL_EVALUATOR = "cenROC",
                                              max.iter = 100, verbose = FALSE, times = NULL,
                                              max_time_points = 15){
  model <- mixOmics::block.spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = keepX, scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDR = model$variates[names(Xh)]

  d <- as.data.frame(tt_mbsplsDR[[1]][,,drop = FALSE])
  for(b in names(Xh)[2:length(Xh)]){
    d <- cbind(d, as.data.frame(tt_mbsplsDR[[b]][,,drop = FALSE]))
  }

  colnames(d) <- apply(expand.grid(colnames(tt_mbsplsDR[[1]]), names(Xh)), 1, paste, collapse="_")

  cox_model <- NULL
  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = cbind(d, Yh),
                      ties = "efron",
                      singular.ok = TRUE,
                      robust = TRUE,
                      nocenter = rep(1, ncol(d)),
                      model = TRUE, x = TRUE)
    },
    # Specifying error message
    error = function(e){
      message(paste0("mb_splsdrcox_dynamic: ",conditionMessage(e)))
      return(NA)
    }
  )

  lp <- getLinealPredictors(cox = cox_model$fit, data = d)

  if(is.null(times)){
    times <- getTimesVector(Yh, max_time_points)
  }

  #C-index and AUC
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Yh, times = times, bestModel = NULL, eval = "mean", method = EVAL_EVALUATOR, PARALLEL = FALSE, verbose = verbose)

  #BRIER
  #lst_BRIER_values <- survAUC_BRIER_LP(lp = lp$fit, Y = Yh, lp_new = lp$fit, Y_test = Yh)
  lst_BRIER_values <- SURVCOMP_BRIER_LP(lp_train = lp$fit, Y_train = Yh, lp_test = lp$fit, Y_test = Yh)

  return(list("c_index" = cox_model$fit$concordance["concordance"], "AUC" = lst_AUC_values$AUC, "BRIER" = lst_BRIER_values$ierror))
}

getCIndex_AUC_CoxModel_block.splsda <- function(Xh, Yh, n.comp, keepX, scale = FALSE,
                                                near.zero.var = FALSE,
                                                EVAL_EVALUATOR = "cenROC", max.iter = 100,
                                                verbose = verbose, times = NULL,
                                                max_time_points = 15){
  model <- mixOmics::block.splsda(X = Xh, Y = Yh[,"event"], ncomp = n.comp, keepX = keepX, scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDR = model$variates[names(Xh)]

  d <- as.data.frame(tt_mbsplsDR[[1]][,,drop = FALSE])
  for(b in names(Xh)[2:length(Xh)]){
    d <- cbind(d, as.data.frame(tt_mbsplsDR[[b]][,,drop = FALSE]))
  }

  colnames(d) <- apply(expand.grid(colnames(tt_mbsplsDR[[1]]), names(Xh)), 1, paste, collapse="_")

  cox_model <- NULL
  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = cbind(d, Yh),
                      ties = "efron",
                      singular.ok = TRUE,
                      robust = TRUE,
                      nocenter = rep(1, ncol(d)),
                      model = TRUE, x = TRUE)
    },
    # Specifying error message
    error = function(e){
      message(paste0("mb_splsda_dynamic: ",conditionMessage(e)))
      return(NA)
    }
  )

  if(is.null(times)){
    times <- getTimesVector(Yh, max_time_points)
  }

  lp <- getLinealPredictors(cox = cox_model$fit, data = d)

  #C-index and AUC
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Yh, times = times, bestModel = NULL, eval = "mean", method = EVAL_EVALUATOR, PARALLEL = FALSE, verbose = verbose)

  #BRIER
  #lst_BRIER_values <- survAUC_BRIER_LP(lp = lp$fit, Y = Yh, lp_new = lp$fit, Y_test = Yh)
  lst_BRIER_values <- SURVCOMP_BRIER_LP(lp_train = lp$fit, Y_train = Yh, lp_test = lp$fit, Y_test = Yh)

  return(list("c_index" = cox_model$fit$concordance["concordance"], "AUC" = lst_AUC_values$AUC, "BRIER" = lst_BRIER_values$ierror))
}

getVarExpModel_block.spls <- function(Xh, DR_coxph_ori, n.comp, keepX, scale = FALSE){
  model <- mixOmics::block.spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = keepX, scale = scale)
  var_exp <- data.frame(lapply(model$prop_expl_var, sum))
}

getBestVectorMB <- function(Xh, DR_coxph = NULL, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE,
                            MIN_NVAR = 10, MAX_NVAR = 10000, cut_points = 5, EVAL_METHOD = "AUC",
                            EVAL_EVALUATOR = "cenROC", PARALLEL = FALSE, mode = "spls", times = NULL,
                            max_time_points = 15, verbose = FALSE){

  if(!mode %in% c("spls", "splsda")){
    stop("Mode must be one of: 'spls' or 'splsda'")
  }

  if(!EVAL_METHOD %in% c("AUC", "BRIER", "c_index")){
    stop("Evaluation method must be one of: 'AUC', 'BRIER' or 'c_index'")
  }

  max_ncol <- purrr::map(Xh, ~ncol(.))

  if(is.null(vector)){
    #vector <- purrr::map(names(Xh), ~c(min(MIN_NVAR, max_ncol[[.]]), (max_ncol[[.]]+min(MIN_NVAR, max_ncol[[.]]))/2, min(max_ncol[[.]], MAX_NVAR)))
    vector <- purrr::map(names(Xh), ~getVectorCuts(vector = c(min(MIN_NVAR, max_ncol[[.]]):min(max_ncol[[.]], MAX_NVAR)), cut_points = cut_points, verbose = verbose))
    names(vector) <- names(Xh)
  }else{
    #check vector is a list, and each value is less than the max.variables of that block
    if(!isa(vector[1],"list")){
      aux <- list()
      for(b in names(Xh)){
        aux[[b]] <- min(vector, dim(Xh[[b]])) #to not overpass the dimension
      }
      vector <- aux
      message(paste0("The initial vector is: ", paste0(vector, collapse = ", ")))
    }else{
      if(!all(names(vector) %in% names(Xh))){
        message("Your vector not contains the block names. A start vector is created using your vector.")
        vector <- purrr::map(names(Xh), ~c(vector))
        names(vector) <- names(Xh)
        message(paste0("The initial vector is: ", paste0(vector, collapse = ", ")))
      }
    }
  }

  if(verbose){
    message(paste0("Original vector: "))
    message(paste0("Block ", names(vector), ": ", paste0(vector), "\n"))
  }

  lst_mb.spls <- list()
  count = 1
  var_exp = NULL

  #if n_col is minimum than MIN_NVAR, values could be the same, so delete duplicates
  vector <- purrr::map(vector, ~unique(.))

  l <- vector
  all_comb <- expand.grid(l)

  ## ALL POSSIBLE COMBINATIONS WITH N_VECTOR AND M_BLOCKS
  ## So we are trying all different number of variables per block
  list_KeepX <- list()
  for(i in 1:nrow(all_comb)){
    keepX = list()
    iter_name = NULL
    for(k in names(Xh)){
      keepX[[k]] = rep(all_comb[[k]][[i]], n.comp)
      iter_name = c(iter_name, keepX[[k]][[1]])
    }
    list_KeepX[[paste0(iter_name, collapse = "_")]] <- keepX
  }

  if(PARALLEL){
    n_cores <- future::availableCores() - 1

    if(.Platform$OS.type == "unix") {
      future::plan("multicore", workers = min(length(list_KeepX), n_cores))
    }else{
      future::plan("multisession", workers = min(length(list_KeepX), n_cores))
    }

    t1 <- Sys.time()
    if(mode %in% "spls"){
      lst_cox_value <- furrr::future_map(list_KeepX, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
    }else{
      lst_cox_value <- furrr::future_map(list_KeepX, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
    }
    t2 <- Sys.time()
    future::plan("sequential")
  }else{
    t1 <- Sys.time()
    if(mode %in% "spls"){
      lst_cox_value <- purrr::map(list_KeepX, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
    }else{
      lst_cox_value <- purrr::map(list_KeepX, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
    }
    t2 <- Sys.time()
  }

  df_cox_value <- NULL
  for(i in 1:length(lst_cox_value)){
    if(EVAL_METHOD=="AUC"){
      df_cox_value <- rbind(df_cox_value, lst_cox_value[[i]]$AUC)
    }else if(EVAL_METHOD=="c_index"){
      df_cox_value <- rbind(df_cox_value, lst_cox_value[[i]]$c_index)
    }else if(EVAL_METHOD=="BRIER"){
      df_cox_value <- rbind(df_cox_value, lst_cox_value[[i]]$BRIER)
    }
  }
  if(EVAL_METHOD=="BRIER"){
    df_cox_value <- 1 - df_cox_value #maximize 1-brier
  }
  rownames(df_cox_value) <- names(list_KeepX)

  index <- which.max(df_cox_value) #MAX CONCORDANCE

  #Select best keepX
  keepX <- list_KeepX[[index]]
  FLAG = TRUE
  cont = 0

  best_c_index <- df_cox_value[index]
  best_keepX <- keepX

  if(verbose){
    message(paste0("First selection: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), "Pred. Value (", EVAL_METHOD ,"): ", round(best_c_index, 4), "\n")
  }

  ori_vector <- vector
  aux_vector <- vector
  p_val <- rep(NA, length(vector))
  names(p_val) <- vector
  p_val <- df_cox_value[,1]

  while(FLAG){
    cont = cont + 1

    if(verbose){
      message(paste0("Iteration: ", cont, "\n"))
    }

    new_vector <- list()

    #before_vector - always go two sides
    for(b in names(best_keepX)){
      aux <- best_keepX[[b]][[1]]
      index <- which(aux_vector[[b]] < aux)
      if(length(index)==0){
        value = aux_vector[[b]][which(aux_vector[[b]] == aux)]
      }else{
        value = round(mean(c(aux, aux_vector[[b]][index][[length(aux_vector[[b]][index])]]))) #last smaller
      }
      new_vector[[b]] <- unique(c(new_vector[[b]], aux, value))
    }

    #next_vector - always go two sides
    for(b in names(best_keepX)){
      aux <- best_keepX[[b]][[1]]
      index <- which(aux_vector[[b]] > aux)
      if(length(index)==0){
        value = aux_vector[[b]][which(aux_vector[[b]] == aux)]
      }else{
        value = round(mean(c(aux, aux_vector[[b]][index][[1]]))) #first greater
      }
      new_vector[[b]] <- unique(c(new_vector[[b]], aux, value))
    }

    # If two blocks and the first with only one value:
    # We do not need to test the second block with the same value: clinic = 9, miRNA = 124, 105, 142 :: because 9,124 was tested
    # however, if multiple blocks, clinic = 9, miRNA = 124, 105, 142, protein = 124, 105, 142 :: we want to test 124 miRNA with 105 proteomic
    # so, we should keep it. Summary, if all block length 1 minus one, delete, in other case, keep. Or just do not test with the same values
    # that have been already tested

    lst_length <- unlist(lapply(new_vector, length))
    if(sum(lst_length==1)==1){ #if just one block with multiples values
      for(b in names(best_keepX)){
        #first value already tested if multiple
        if(length(new_vector[[b]])>1){
          new_vector[[b]] <- new_vector[[b]][-1]
        }
      }
    }

    ## if vector == new vector - means a custom vector has been used
    if(identical(new_vector,vector)){
      break
    }

    if(verbose){
      message(paste0("Testing: \n"), paste0("Block ", names(best_keepX), ": ", new_vector, "\n"))
    }

    all_comb <- expand.grid(new_vector)

    ### update all vector
    aux_vector <- purrr::map(names(new_vector), ~unique(c(aux_vector[[.]], new_vector[[.]])))
    names(aux_vector) <- names(new_vector)
    aux_vector <- purrr::map(names(new_vector), ~aux_vector[[.]][order(aux_vector[[.]])])
    names(aux_vector) <- names(new_vector)

    p_val <- c(p_val, rep(NA, prod(unlist(purrr::map(new_vector, length))))) #longitud es el el numero nuevo de posibles combinaciones posibles
    names(p_val)[is.na(p_val)] <- apply(all_comb, 1, function(x){paste0(unlist(x), collapse = "_")})

    ### OTHER KEEP_VECTOR
    list_KeepX_aux <- list()
    for(i in 1:nrow(all_comb)){
      keepX_aux = list()
      iter_name = NULL
      for(k in names(Xh)){
        keepX_aux[[k]] = rep(all_comb[[k]][[i]], n.comp)
        iter_name = c(iter_name, keepX_aux[[k]][[1]])
      }

      if(paste0(iter_name, collapse = "_") %in% names(list_KeepX_aux)){
        next #one already computed
      }

      list_KeepX_aux[[paste0(iter_name, collapse = "_")]] <- keepX_aux
    }

    #### ###
    ## remove any keepX_aux that has been already tested (when at least 2 omics with more than 2 vectors)
    if(cont==1){
      if(any(names(list_KeepX_aux) %in% rownames(df_cox_value))){
        index_tested <- which(names(list_KeepX_aux) %in% rownames(df_cox_value))
        list_KeepX_aux <- list_KeepX_aux[-index_tested]
      }
    }else{
      if(names(list_KeepX_aux) %in% rownames(df_cox_value_aux)){
        index_tested <- which(names(list_KeepX_aux) %in% rownames(df_cox_value_aux))
        list_KeepX_aux <- list_KeepX_aux[-index_tested]
      }
    }

    ## Compute next c_index
    if(PARALLEL){
      n_cores <- future::availableCores() - 1

      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(list_KeepX_aux), n_cores))
      }else{
        future::plan("multisession", workers = min(length(list_KeepX_aux), n_cores))
      }

      t1 <- Sys.time()
      if(mode %in% "spls"){
        lst_cox_value <- furrr::future_map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
      }else{
        lst_cox_value <- furrr::future_map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
      }
      t2 <- Sys.time()
      future::plan("sequential")
    }else{
      t1 <- Sys.time()
      if(mode %in% "spls"){
        lst_cox_value <- purrr::map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
      }else{
        lst_cox_value <- purrr::map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = FALSE, near.zero.var = FALSE, max.iter = max.iter, verbose = verbose, times = times, max_time_points = max_time_points), .progress = FALSE)
      }
      t2 <- Sys.time()
    }

    df_cox_value_aux <- NULL
    for(i in 1:length(lst_cox_value)){
      if(EVAL_METHOD=="AUC"){
        df_cox_value_aux <- rbind(df_cox_value_aux, lst_cox_value[[i]]$AUC)
      }else if(EVAL_METHOD=="c_index"){
        df_cox_value_aux <- rbind(df_cox_value_aux, lst_cox_value[[i]]$c_index)
      }else if(EVAL_METHOD=="BRIER"){
        df_cox_value_aux <- rbind(df_cox_value_aux, lst_cox_value[[i]]$BRIER)
      }
    }
    if(EVAL_METHOD=="BRIER"){
      df_cox_value_aux <- 1 - df_cox_value_aux #maximize 1-brier
    }
    rownames(df_cox_value_aux) <- names(list_KeepX_aux)

    #index <- which.max(rowSums(df_cox_value_aux)) #MAX VAR_MEDIA
    #index <- which.max(df_cox_value_aux[,"Y"]) #MAX Y?
    index <- which.max(df_cox_value_aux) #MAX CONCORDANCE
    best_c_index_aux <- df_cox_value_aux[index]
    p_val[rownames(df_cox_value_aux)] <- df_cox_value_aux

    if(best_c_index >= best_c_index_aux | best_c_index_aux-best_c_index <= MIN_AUC_INCREASE){
      FLAG = FALSE
      if(verbose){
        message(paste0("End: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value (", EVAL_METHOD ,"): ", round(best_c_index, 4), "\n"))
      }
    }else{
      best_c_index <- best_c_index_aux
      best_keepX <- list_KeepX_aux[[index]]
      if(verbose){
        message(paste0("New Vector: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value (", EVAL_METHOD ,"): ", round(best_c_index_aux, 4), "n"))
      }
    }
  }

  aux_df <- t(data.frame(sapply(names(p_val), strsplit, "_")))
  aux_df <- as.data.frame(aux_df)
  rownames(aux_df) <- names(p_val)
  colnames(aux_df) <- names(aux_vector)
  for(cn in colnames(aux_df)){
    aux_df[,cn] <- as.numeric(aux_df[,cn])
    aux_df <- aux_df[order(aux_df[,cn]),]
  }

  p_val <- p_val[rownames(aux_df)]

  keepX <- best_keepX
  return(list(best.keepX = keepX, p_val = p_val))
}
