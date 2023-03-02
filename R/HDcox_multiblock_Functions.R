#' deleteZeroOrNearZeroVariance.mb
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param freqCut Numeric. Cutoff for the ratio of the most common value to the second most common value (default: 95/5).
#'
#' @export

deleteZeroOrNearZeroVariance.mb <- function(X, remove_near_zero_variance = F, remove_zero_variance = T, toKeep.zv = NULL, freqCut = 95/5){

  auxX <- X

  variablesDeleted <- NULL
  if(remove_near_zero_variance){
    lst.zv <- purrr::map(auxX, ~deleteZeroVarianceVariables(data = ., info = T, mustKeep = toKeep.zv, freqCut = freqCut))
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
    lst.zv <- purrr::map(auxX, ~deleteZeroVarianceVariables(data = ., info = T, mustKeep = toKeep.zv, onlyZero = T))
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

checkXY.mb.class <- function(X, Y, verbose = T){
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

check.mb.ncomp <- function(X, max.ncomp, verbose = F){

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

check.mb.maxPredictors <- function(X, Y, MIN_EPV, max.variables, verbose = F){
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
  if(!any(c("event", "status") %in% colnames(Y))){
    stop_quietly("Y data.frame must contain the colname 'event' or 'status'.")
  }else if("status" %in% colnames(Y)){
    colnames(Y)[colnames(Y)=="status"] <- "event"
  }

  lst_X_train <- list()
  lst_Y_train <- list()

  lst_X_test <- list()
  lst_Y_test <- list()

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

  return(list(lst_X_train = lst_X_train, lst_Y_train = lst_Y_train, lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, lst_train_index = lst_obs_index_train, lst_test_index = lst_obs_index_test))
}

XY.mb.scale <- function(X, Y, x.center, x.scale, y.center, y.scale){
  xmeans <- NULL
  xsds <- NULL
  ymeans <- NULL
  ysds <- NULL

  #CENTERING X
  if(length(x.center)>1){
    if(length(x.center)==length(X) & all(names(X) %in% names(x.center))){
      Xh <- purrr::map2(.x = X, .y = names(X), ~scale(., center = x.center[[.y]], scale = F))
    }
  }else{
    Xh <- purrr::map(X, ~scale(., center = x.center, scale = F))
  }

  #SCALING X
  if(length(x.scale)>1){
    if(length(x.scale)==length(X) & all(names(X) %in% names(x.scale))){
      Xh <- purrr::map2(.x = Xh, .y = names(X), ~scale(., center = F, scale = x.scale[[.y]]))
    }
  }else{
    Xh <- purrr::map(X, ~scale(., center = F, scale = x.scale))
  }

  xmeans <- purrr::map(Xh, ~attr(., "scaled:center"))
  xsds <- purrr::map(Xh, ~attr(., "scaled:scale"))

  #SCALING Y
  if(y.center | y.scale){
    Yh_time <- scale(Y[,"time", drop=F], center = y.center, scale = y.scale)
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

getCIndex_AUC_CoxModel_block.spls <- function(Xh, DR_coxph_ori, Yh, n.comp, keepX, scale = F, near.zero.var = F, EVAL_EVALUATOR = "cenROC", max.iter = 100, verbose = F, times = NULL, max_time_points = 15){
  model <- mixOmics::block.spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = keepX, scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDR = model$variates[names(Xh)]

  d <- as.data.frame(tt_mbsplsDR[[1]][,,drop=F])
  for(b in names(Xh)[2:length(Xh)]){
    d <- cbind(d, as.data.frame(tt_mbsplsDR[[b]][,,drop=F]))
  }

  colnames(d) <- apply(expand.grid(colnames(tt_mbsplsDR[[1]]), names(Xh)), 1, paste, collapse="_")

  cox_model <- NULL
  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = cbind(d, Yh),
                      ties = "efron",
                      singular.ok = T,
                      robust = T,
                      nocenter = rep(1, ncol(d)),
                      model=T)
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

  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Yh, times = times, bestModel = NULL, eval = "mean", method = EVAL_EVALUATOR, PARALLEL = FALSE, verbose = verbose)

  return(list("c_index" = cox_model$fit$concordance["concordance"], "AUC" = lst_AUC_values$AUC))
}

getCIndex_AUC_CoxModel_block.splsda <- function(Xh, Yh, n.comp, keepX, scale = F, near.zero.var = F, EVAL_EVALUATOR = "cenROC", max.iter = 100, verbose = verbose, times = NULL, max_time_points = 15){
  model <- mixOmics::block.splsda(X = Xh, Y = Yh[,"event"], ncomp = n.comp, keepX = keepX, scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDR = model$variates[names(Xh)]

  d <- as.data.frame(tt_mbsplsDR[[1]][,,drop=F])
  for(b in names(Xh)[2:length(Xh)]){
    d <- cbind(d, as.data.frame(tt_mbsplsDR[[b]][,,drop=F]))
  }

  colnames(d) <- apply(expand.grid(colnames(tt_mbsplsDR[[1]]), names(Xh)), 1, paste, collapse="_")

  cox_model <- NULL
  cox_model$fit <- tryCatch(
    # Specifying expression
    expr = {
      survival::coxph(formula = survival::Surv(time,event) ~ .,
                      data = cbind(d, Yh),
                      ties = "efron",
                      singular.ok = T,
                      robust = T,
                      nocenter = rep(1, ncol(d)),
                      model=T)
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
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Yh, times = times, bestModel = NULL, eval = "mean", method = EVAL_EVALUATOR, PARALLEL = FALSE, verbose = verbose)

  return(list("c_index" = cox_model$fit$concordance["concordance"], "AUC" = lst_AUC_values$AUC))
}

getVarExpModel_block.spls <- function(Xh, DR_coxph_ori, n.comp, keepX, scale = F){
  model <- mixOmics::block.spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = keepX, scale = scale)
  var_exp <- data.frame(lapply(model$prop_expl_var, sum))
}

getBestVectorMB <- function(Xh, DR_coxph = NULL, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = 10, MAX_NVAR = 10000, cut_points = 5, EVAL_METHOD = "AUC", EVAL_EVALUATOR = "cenROC", PARALLEL = F, mode = "spls", times = NULL, max_time_points = 15, verbose = F){

  if(!mode %in% c("spls", "splsda")){
    stop("Mode must be one of: 'spls' or 'splsda'")
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
        aux[[b]] <- vector
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
      lst_cox_value <- furrr::future_map(list_KeepX, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
    }else{
      lst_cox_value <- furrr::future_map(list_KeepX, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
    }
    t2 <- Sys.time()
    future::plan("sequential")
  }else{
    t1 <- Sys.time()
    if(mode %in% "spls"){
      lst_cox_value <- purrr::map(list_KeepX, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
    }else{
      lst_cox_value <- purrr::map(list_KeepX, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
    }
    t2 <- Sys.time()
  }

  df_cox_value <- NULL
  for(i in 1:length(lst_cox_value)){
    if(EVAL_METHOD=="AUC"){
      df_cox_value <- rbind(df_cox_value, lst_cox_value[[i]]$AUC)
    }else{
      df_cox_value <- rbind(df_cox_value, lst_cox_value[[i]]$c_index)
    }
  }
  rownames(df_cox_value) <- names(list_KeepX)

  index <- which.max(df_cox_value) #MAX CONCORDANCE

  #Select best keepX
  keepX <- list_KeepX[[index]]
  FLAG = T
  cont = 0

  best_c_index <- df_cox_value[index]
  best_keepX <- keepX

  if(verbose){
    message(paste0("First selection: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), "Pred. Value: ", round(best_c_index, 4), "\n")
  }

  ori_vector <- vector
  aux_vector <- vector

  while(FLAG){
    cont = cont + 1

    if(verbose){
      message(paste0("Iteration: ", cont, "\n"))
    }

    new_vector <- list()

    #before_vector
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

    if(verbose){
      message(paste0("Testing: \n"), paste0("Block ", names(best_keepX), ": ", new_vector, "\n"))
    }

    all_comb <- expand.grid(new_vector)

    ### update all vector
    aux_vector <- purrr::map(names(new_vector), ~unique(c(aux_vector[[.]], new_vector[[.]])))
    names(aux_vector) <- names(new_vector)
    aux_vector <- purrr::map(names(new_vector), ~aux_vector[[.]][order(aux_vector[[.]])])
    names(aux_vector) <- names(new_vector)

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
        lst_cox_value <- furrr::future_map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
      }else{
        lst_cox_value <- furrr::future_map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
      }
      t2 <- Sys.time()
      future::plan("sequential")
    }else{
      t1 <- Sys.time()
      if(mode %in% "spls"){
        lst_cox_value <- purrr::map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
      }else{
        lst_cox_value <- purrr::map(list_KeepX_aux, ~getCIndex_AUC_CoxModel_block.splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, verbose = verbose, time = time, max_time_points = max_time_points), .progress = F)
      }
      t2 <- Sys.time()
    }

    df_cox_value_aux <- NULL
    for(i in 1:length(lst_cox_value)){
      if(EVAL_METHOD=="AUC"){
        df_cox_value_aux <- rbind(df_cox_value_aux, lst_cox_value[[i]]$AUC)
      }else{
        df_cox_value_aux <- rbind(df_cox_value_aux, lst_cox_value[[i]]$c_index)
      }
    }

    rownames(df_cox_value_aux) <- names(list_KeepX_aux)
    #index <- which.max(rowSums(df_cox_value_aux)) #MAX VAR_MEDIA
    #index <- which.max(df_cox_value_aux[,"Y"]) #MAX Y?
    index <- which.max(df_cox_value_aux) #MAX CONCORDANCE
    best_c_index_aux <- df_cox_value_aux[index]

    if(best_c_index >= best_c_index_aux | best_c_index_aux-best_c_index <= MIN_AUC_INCREASE){
      FLAG = F
      if(verbose){
        message(paste0("End: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value: ", round(best_c_index, 4), "\n"))
      }
    }else{
      best_c_index <- best_c_index_aux
      best_keepX <- list_KeepX_aux[[index]]
      if(verbose){
        message(paste0("New Vector: \n"), paste0(paste0("Block ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value: ", round(best_c_index_aux, 4), "n"))
      }
    }
  }

  keepX <- best_keepX
  return(keepX)
}
