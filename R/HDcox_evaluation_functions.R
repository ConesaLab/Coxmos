cv.getScoreFromWeight <- function(lst_cox_mean, w_AIC, w_c.index, w_AUC, colname_AIC = "AIC_mean", colname_c_index = "c_index_mean", colname_AUC = "AUC_mean"){

  if(nrow(lst_cox_mean)!=1){
    if(w_AUC!=0){
      aux <- scale(lst_cox_mean[,c(colname_AIC, colname_c_index, colname_AUC), drop=F])
      #c_index and AUC max is better
      aux[,1] <- aux[,colname_AIC]*-1 #min AIC better
    }else{ #if AUC was no evaluated, it has not to be considered
      aux <- scale(lst_cox_mean[,c(colname_AIC, colname_c_index), drop=F])
      aux[,1] <- aux[,colname_AIC]*-1 #min AIC better
    }

  }else{
    if(w_AUC!=0){
      aux <- lst_cox_mean[,c(colname_AIC, colname_c_index, colname_AUC), drop=F]
    }else{
      aux <- lst_cox_mean[,c(colname_AIC, colname_c_index), drop = F]
    }
  }

  for(cn in colnames(aux)){
    if(all(is.nan(aux[,cn]))){
      #all values the same, same models across penalties or components...
      if(cn == "AIC"){
        aux[,"AIC"] <- lst_cox_mean$AIC
      }else if(cn == "c_index"){
        aux[,"c_index"] <- lst_cox_mean$c_index
      }else if(cn == "AUC"){
        aux[,"AUC"] <- lst_cox_mean$AUC
      }
    }
  }

  if(w_AUC!=0){
    score = aux %*% c(w_AIC, w_c.index, w_AUC)
  }else{
    score = aux %*% c(w_AIC, w_c.index)
    lst_cox_mean <- removeColumn(lst_cox_mean, colname_AUC)
  }

  lst_cox_mean[,"score"] <- score

  return(lst_cox_mean)
}

removeColumn <- function(data, lst_cn){
  return(deleteColumn(data, lst_cn))
}

deleteColumn <- function(data, lst_cn){
  if(length(lst_cn)==0)
    return(data)

  for(cn in lst_cn){
    if(is.data.frame(data) & cn %in% colnames(data)){
      data[,cn] <- NULL
    }else if(is.matrix(data) & cn %in% colnames(data)){
      data <- data[,-cn]
    }
  }
  return(data)
}

splitData_Iterations_Folds <- function(X, Y, n_run, k_folds, seed = 123){

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
    set.seed(seed+i)
    trainIndex <- caret::createFolds(Y[,"event"],
                                     k = k_folds,
                                     list = TRUE)

    #for each fold, take the others as train
    lst_X_data_train <- lapply(trainIndex, function(ind, dat) dat[-ind,], dat = X)
    lst_Y_data_train <- lapply(trainIndex, function(ind, dat) dat[-ind,], dat = Y)

    #for each fold, take just one as test
    lst_X_data_test <- lapply(trainIndex, function(ind, dat) dat[ind,], dat = X)
    lst_Y_data_test <- lapply(trainIndex, function(ind, dat) dat[ind,], dat = Y)

    lst_X_train[[i]] <- lst_X_data_train
    lst_Y_train[[i]] <- lst_Y_data_train

    lst_X_test[[i]] <- lst_X_data_test
    lst_Y_test[[i]] <- lst_Y_data_test
  }

  names(lst_X_train) <- paste0("run", 1:n_run)
  names(lst_Y_train) <- paste0("run", 1:n_run)
  names(lst_X_test) <- paste0("run", 1:n_run)
  names(lst_Y_test) <- paste0("run", 1:n_run)

  return(list(lst_X_train = lst_X_train, lst_Y_train = lst_Y_train, lst_X_test = lst_X_test, lst_Y_test = lst_Y_test))
}

getAUC_from_LP_2.0 <- function(linear.predictors, Y, times, bestModel = NULL, method = "cenROC", eval = "median", PARALLEL = F, verbose = F){

  # if(!method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("Method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("Method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  if(!all(c("time", "event") %in% colnames(Y))){
    stop_quietly("Data.frame Y must contain the columns time and event for COX model.")
  }

  #### TIMES SHOULD BE SOMETHING NOW; BUT IN CASE WILL BE NULL...
  if(is.null(times)){
    max_time_points = 15
    times <- getVectorOfTime(Y, max_time_points)
  }

  if(is.null(linear.predictors)){
    return(NA)
  }else if(!isa(linear.predictors, "list")){
    aux = linear.predictors
    linear.predictors = NULL
    linear.predictors$fit = aux
  }else if(!"fit" %in% names(linear.predictors)){
    stop_quietly("fit must be a list inside linea predictors object.")
  }

  #order Y
  ids <- names(linear.predictors$fit)[names(linear.predictors$fit) %in% rownames(Y)]
  Y <- Y[ids,]

  if(method %in% c(pkg.env$AUC_nsROC, pkg.env$AUC_cenROC)){
    Y <- data.matrix(Y)
  }

  AUC = NULL #area under the curve for each timepoint
  out = NULL
  if(method == pkg.env$AUC_risksetROC){
    times.aux <- times #times[times.vector]
    if(PARALLEL){
      n_cores <- max(future::availableCores() - 1, 1)
      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(times.aux), n_cores))
      }else{
        future::plan("multisession", workers = min(length(times.aux), n_cores))
      }
      out <- furrr::future_map(times.aux, ~risksetROC_tryCatch(marker = linear.predictors$fit, Stime = Y[,"time"], status = Y[,"event"],
                                                               predict.time = ., verbose = verbose))
      future::plan("sequential")
    }else{
      out <- purrr::map(times.aux, ~risksetROC_tryCatch(linear.predictors$fit, Stime = Y[,"time"], status = Y[,"event"], predict.time = ., verbose = verbose))
    }

    lst.AUC <- getAUC_vector(output = out, method = method, eval = eval)
    AUC <- lst.AUC$AUC
    AUC.vector <- lst.AUC$AUC.vector

    if(all(is.na(AUC)) & all(is.na(AUC.vector))){
      AUC <- rep(NA, length(times.aux))
      AUC.vector <- rep(NA, length(times.aux))
    }

  }else if(method == pkg.env$AUC_survivalROC){
    #https://cran.r-project.org/web/packages/survivalROC/
    #needs at least 2 events per time and time can not be day 0
    times.vector <- timesAsumption_AUC_Eval(Y, times)
    times.aux <- times #times[times.vector]
    if(!all(times.vector==F)){
      times_run <- times.aux[which(times.vector==T)]
      if(PARALLEL){
        n_cores <- max(future::availableCores() - 1, 1)
        if(.Platform$OS.type == "unix") {
          future::plan("multicore", workers = min(length(times.aux), n_cores))
        }else{
          future::plan("multisession", workers = min(length(times.aux), n_cores))
        }
        out <- furrr::future_map(times_run, ~survivalROC_tryCatch(Stime = Y[,"time"], status = Y[,"event"], marker = linear.predictors$fit,
                                                                  predict.time = ., method = "NNE",
                                                                  span = 0.25 * nrow(Y)^(-0.20), verbose = verbose))
        future::plan("sequential")
      }else{
        out <- purrr::map(times_run, ~survivalROC_tryCatch(Stime = Y[,"time"], status = Y[,"event"], marker = linear.predictors$fit,
                                                           predict.time = ., method = "NNE",
                                                           span = 0.25 * nrow(Y)^(-0.20)))
      }

      lst.AUC <- getAUC_vector(output = out, method = method, eval = eval)

      aux_auc.v <- rep(NA, length(times))
      aux_auc.v[which(times.vector==T)] <- lst.AUC$AUC.vector

      AUC <- lst.AUC$AUC
      AUC.vector <- aux_auc.v
    }else{
      AUC <- rep(NA, length(times.vector))
      AUC.vector <- rep(NA, length(times.vector))
    }
  }else if(method == pkg.env$AUC_cenROC){
    #https://cran.r-project.org/web/packages/cenROC/
    #needs at least 2 events per time and time can not be day 0
    #Y is a matrix
    times.vector <- timesAsumption_AUC_Eval(Y, times, method = method)
    times.aux <- times #times[times.vector]
    if(!all(times.vector==F)){
      times_run <- times.aux[which(times.vector==T)]
      #METHOD
      if(PARALLEL){
        n_cores <- max(future::availableCores() - 1, 1)
        if(.Platform$OS.type == "unix") {
          future::plan("multicore", workers = min(length(times.aux), n_cores))
        }else{
          future::plan("multisession", workers = min(length(times.aux), n_cores))
        }
        out <- furrr::future_map(times_run, ~cenROC_tryCatch(Y = Y[,"time"], censor = Y[,"event"], M = linear.predictors$fit,
                                                             t = ., method = "tra", ktype  ="normal", alpha = 0.05, plot = F, verbose = verbose)) #problems in simulated data
        future::plan("sequential")
      }else{
        out <- purrr::map(times_run, ~cenROC_tryCatch(Y = Y[,"time"], censor = Y[,"event"], M = linear.predictors$fit,
                                                      t = ., method = "tra", ktype  ="normal", alpha = 0.05, plot = F, verbose = verbose)) #problems in simulated data
      }

      lst.AUC <- getAUC_vector(output = out, method = method, eval = eval)

      aux_auc.v <- rep(NA, length(times))
      aux_auc.v[which(times.vector==T)] <- lst.AUC$AUC.vector

      AUC <- lst.AUC$AUC
      AUC.vector <- aux_auc.v
    }else{
      AUC <- rep(NA, length(times.vector))
      AUC.vector <- rep(NA, length(times.vector))
    }
  }else if(method == pkg.env$AUC_nsROC){
    #https://cran.r-project.org/web/packages/nsROC/nsROC.pdf
    #needs at least 2 events per time and time can not be day 0
    #Y is a matrix
    #cumulative/dynamic approach
    times.vector <- timesAsumption_AUC_Eval(Y, times)
    times.aux <- times #times[times.vector]
    if(!all(times.vector==F)){
      times_run <- times.aux[which(times.vector==T)]
      #METHOD
      if(PARALLEL){
        n_cores <- max(future::availableCores() - 1, 1)
        if(.Platform$OS.type == "unix") {
          future::plan("multicore", workers = min(length(times.aux), n_cores))
        }else{
          future::plan("multisession", workers = min(length(times.aux), n_cores))
        }
        out <- furrr::future_map(times_run, ~nsROC_tryCatch(stime = Y[,"time"],
                                                            status = Y[,"event"],
                                                            marker = linear.predictors$fit,
                                                            predict.time = .,
                                                            method = "wKM", # Cox, KM, wKM (last with kernel-weighted)
                                                            kernel = "normal", #normal, Epanechnikov, other (only if wKM)
                                                            ci = F, #a lot of time and problems with NAs in bootstraping
                                                            boot.n = 200,
                                                            conf.level = 0.95,
                                                            seed = 123, verbose = verbose))
        future::plan("sequential")
      }else{
        out <- purrr::map(times_run, ~nsROC_tryCatch(stime = Y[,"time"],
                                                     status = Y[,"event"],
                                                     marker = linear.predictors$fit,
                                                     predict.time = .,
                                                     method = "wKM", # Cox, KM, wKM (last with kernel-weighted)
                                                     kernel = "normal", #normal, Epanechnikov, other (only if wKM)
                                                     ci = F, #a lot of time and problems with NAs in bootstraping
                                                     boot.n = 200,
                                                     conf.level = 0.95,
                                                     seed = 123, verbose = verbose))
      }

      lst.AUC <- getAUC_vector(output = out, method = method, eval = eval)

      aux_auc.v <- rep(NA, length(times))
      aux_auc.v[which(times.vector==T)] <- lst.AUC$AUC.vector

      AUC <- lst.AUC$AUC
      AUC.vector <- aux_auc.v
    }else{
      AUC <- rep(NA, length(times.vector))
      AUC.vector <- rep(NA, length(times.vector))
    }
  }else if(method == pkg.env$AUC_smoothROCtime_C){
    #https://cran.r-project.org/web/packages/smoothROCtime/
    times.vector <- timesAsumption_AUC_Eval(Y, times)
    times.aux <- times #times[times.vector]
    if(!all(times.vector == FALSE)){
      times_run <- times.aux[which(times.vector==T)]
      d <- cbind(Y[,"time"], Y[,"event"], linear.predictors$fit)
      out <- smoothROCtime_tryCatch(data = d, times = times_run, tcr = "C", meth = "1", verbose = verbose) #Cumulative/Dynamic with smooth method

      lst.AUC <- getAUC_vector(output = out, method = method, eval = eval, times = times_run, times.vector = times.vector)

      aux_auc.v <- rep(NA, length(times))
      aux_auc.v[which(times.vector==T)] <- lst.AUC$AUC.vector

      AUC <- lst.AUC$AUC
      AUC.vector <- aux_auc.v
    }else{
      AUC <- rep(NA, length(times.vector))
      AUC.vector <- rep(NA, length(times.vector))
    }

  }else if(method == pkg.env$AUC_smoothROCtime_I){
    #https://cran.r-project.org/web/packages/smoothROCtime/
    times.vector <- timesAsumption_AUC_Eval(Y, times)
    times.aux <- times #times[times.vector]
    if(!all(times.vector == FALSE)){
      times_run <- times.aux[which(times.vector==T)]
      d <- cbind(Y[,"time"], Y[,"event"], linear.predictors$fit)

      out <- smoothROCtime_tryCatch(data = d, times = times.aux, tcr = "I", meth = "1", verbose = verbose) #Cumulative/Dynamic with smooth method

      lst.AUC <- getAUC_vector(output = out, method = method, eval = eval, times = times_run, times.vector = NULL)

      aux_auc.v <- rep(NA, length(times))
      aux_auc.v[which(times.vector==T)] <- lst.AUC$AUC.vector

      AUC <- lst.AUC$AUC
      AUC.vector <- aux_auc.v
    }else{
      AUC <- rep(NA, length(times.vector))
      AUC.vector <- rep(NA, length(times.vector))
    }
  }else{
    stop_quietly("No available method selected.")
  }

  names(AUC.vector) <- times

  #!!!! This it's just a patch, I do not know why some NaN are present
  if(any(is.nan(AUC.vector))){
    AUC.vector[is.nan(AUC.vector)] <- NA
    if(eval=="mean"){
      AUC <- mean(AUC.vector, na.rm = T)
    }else{
      AUC <- mean(AUC.vector, na.rm = T)
    }
  }

  #!!!! This it's just a patch, I do not know why smoothROCtime_I can generate negative values
  if(any(AUC.vector[!is.na(AUC.vector)]<0)){
    vector_nna <- AUC.vector[!is.na(AUC.vector)]
    vector_nna[vector_nna<0] <- NA

    AUC.vector[!is.na(AUC.vector)] <- vector_nna

    if(eval=="mean"){
      AUC <- mean(AUC.vector, na.rm = T)
    }else{
      AUC <- mean(AUC.vector, na.rm = T)
    }
  }

  return(list(lp.used = linear.predictors, AUC = AUC, AUC.vector = AUC.vector, times = times.aux))
}

getAUC_vector <-function(output, method, eval, times = NULL, times.vector = NULL){

  if(all(is.na(output)) || is.null(output)){
    AUC <- NA
    AUC.vector <- NA
    return(list(AUC = AUC, AUC.vector = AUC.vector))
  }

  if(method %in% c(pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I)){
    AUC <- NULL
    AUC.vector <- NULL
    if(all(is.na(output))){
      AUC <- NA
      AUC.vector <- NA
    }else{
      subjects <- table(output$t)[[1]]
      times_used <- as.numeric(unique(output$t))

      index <- NULL
      for(i in 1:length(times_used)){
        index <- c(index, 1 + (i-1)*subjects)
      }

      AUC.vector <- rep(NA, length(times))
      #AUC.vector[times.vector] <- as.numeric(output$auc[index])
      AUC.vector <- as.numeric(output$auc[index])
      AUC <- ifelse(eval=="median", median(AUC.vector, na.rm = T), mean(AUC.vector, na.rm = T))
    }
  }else if(method %in% pkg.env$AUC_nsROC){
    AUC <- NULL
    AUC.vector <- NULL
    for(t in 1:length(output)){
      if(length(output[[t]]) > 1 && !all(is.na(output[[t]]))){
        AUC.vector <- c(AUC.vector, output[[t]]$auc)
      }else{
        AUC.vector <- c(AUC.vector, NA)
      }
    }
    if(!all(is.na(AUC.vector))){
      AUC <- ifelse(eval=="median", median(AUC.vector, na.rm = T), mean(AUC.vector, na.rm = T))
    }else{
      AUC <- NA
      AUC.vector <- NA
    }
  }else if(method %in% c(pkg.env$AUC_cenROC)){
    AUC <- NULL
    AUC.vector <- NULL
    for(t in 1:length(output)){
      if(length(output[[t]]) > 1 && !all(is.na(output[[t]]))){
        AUC.vector <- c(AUC.vector, output[[t]]$AUC$AUC)
      }else{
        AUC.vector <- c(AUC.vector, NA)
      }
    }
    if(!all(is.na(AUC.vector))){
      AUC <- ifelse(eval=="median", median(AUC.vector, na.rm = T), mean(AUC.vector, na.rm = T))
    }else{
      AUC <- NA
      AUC.vector <- NA
    }
  }else if(method %in% c(pkg.env$AUC_survivalROC)){
    AUC <- NULL
    AUC.vector <- NULL
    for(t in 1:length(output)){
      if(length(output[[t]]) > 1 && !all(is.na(output[[t]]))){
        AUC.vector <- c(AUC.vector, output[[t]]$AUC)
      }else{
        AUC.vector <- c(AUC.vector, NA)
      }
    }
    if(!all(is.na(AUC.vector))){
      AUC <- ifelse(eval=="median", median(AUC.vector, na.rm = T), mean(AUC.vector, na.rm = T))
    }else{
      AUC <- NA
      AUC.vector <- NA
    }
  }else if(method %in% pkg.env$AUC_risksetROC){
    AUC <- NULL
    AUC.vector <- NULL
    for(t in 1:length(output)){
      if(length(output[[t]]) > 1 && !all(is.na(output[[t]]))){
        AUC.vector <- c(AUC.vector, output[[t]]$AUC)
      }else{
        AUC.vector <- c(AUC.vector, NA)
      }
    }
    if(!all(is.na(AUC.vector))){
      AUC <- ifelse(eval=="median", median(AUC.vector, na.rm = T), mean(AUC.vector, na.rm = T))
    }else{
      AUC <- NA
      AUC.vector <- NA
    }
  }else{
    stop_quietly("No available method selected.")
  }

  return(list(AUC = AUC, AUC.vector = AUC.vector))
}

timesAsumption_AUC_Eval <- function(Y, times, method = NULL){
  #which times can be computed
  res <- NULL
  for(t in times){
    if(!sum(Y[(Y[,"event"]==1 | Y[,"event"]==TRUE),"time"]<=t)<2 & t != 0){
      if(method==pkg.env$AUC_cenROC){
        #cenROC needs also two events after the time
        if(sum(Y[(Y[,"event"]==1 | Y[,"event"]==TRUE),"time"]>=t)>2 & t != 0){
          res <- c(res, TRUE)
        }else{
          res <- c(res, FALSE)
        }
      }else{
        res <- c(res, TRUE)
      }
    }else{
      res <- c(res, FALSE)
    }
  }
  return(res)
}

smoothROCtime_tryCatch <- function(data, times, tcr = "C", meth = "1", verbose = F){
  invisible(utils::capture.output(out <- tryCatch(
    # Specifying expression
    expr = {
      smoothROCtime::stRoc(data = as.data.frame(data), t = times, tcr = tcr, meth = meth, bw="naive.pdf")
    },
    # Specifying error message
    error = function(e){
      #message(paste0(e$message, " for tcr ", tcr, " and method 'smoothROCtime'"))
      NA
    }
  )))

  FLAG_AGAIN = F
  if(all(is.na(out))){
    FLAG_AGAIN = T
  }else if(any(as.numeric(out$auc)>1)){
    FLAG_AGAIN = T
  }

  if(FLAG_AGAIN){
    if(verbose){
      message(paste0("Default Bandwith method 'naive.pdf' failed. Changed to: 'Hpi'."))
    }
    invisible(utils::capture.output(out <- tryCatch(
      # Specifying expression
      expr = {
        smoothROCtime::stRoc(data = as.data.frame(data), t = times, tcr = tcr, meth = meth, bw="Hpi")
      },
      # Specifying error message
      error = function(e){
        message(paste0(e$message, " for tcr ", tcr, " and method 'smoothROCtime'."))
        NA
      }
    )))
  }

  if(all(is.na(out))){
    message("Both methods failed. Returning NA.")
    out$auc = NA
  }else if(any(as.numeric(out$auc)>1)){
    message("Both methods failed. Returning NA.")
    out$auc = NA
  }

  return(out)
}

nsROC_tryCatch <- function(stime, status, marker, predict.time, method, kernel, ci, boot.n, conf.level, seed, verbose = F){
  invisible(utils::capture.output(out <- tryCatch(
    # Specifying expression
    expr = {
      nsROC::cdROC(stime = stime,
                   status = status,
                   marker = marker,
                   predict.time = predict.time,
                   method = method, # Cox, KM, wKM (last with kernel-weighted)
                   kernel = kernel, #normal, Epanechnikov, other (only if wKM)
                   ci = ci, #a lot of time and problems with NAs in bootstraping
                   boot.n = boot.n,
                   conf.level = conf.level,
                   seed = seed)
    },
    # Specifying error message
    error = function(e){
      message(paste0(e$message, " for time ", predict.time, " and method '",pkg.env$AUC_nsROC,"'"))
      NA
    }
  )))
  return(out)
}

cenROC_tryCatch <- function(Y, censor, M, t, method, ktype, alpha, plot, verbose = F){
  invisible(utils::capture.output(out <- tryCatch(
    # Specifying expression
    expr = {
      # cenROC::cenROC(Y = Y, censor = censor, M = M,
      #                t = t, method = method, ktype = ktype, alpha = alpha, plot = plot)

      #our function #extracted from GitHub
      cenROC(Y = Y, censor = censor, M = M,
             t = t, method = method, ktype = ktype, alpha = alpha, plot = plot)
    },
    # Specifying error message
    error = function(e){
      message(paste0(e$message, " for time ", t, " and method '",pkg.env$AUC_cenROC,"'"))
      NA
    },
    # warning
    warning = function(w){
      # message(paste0("Problem in time: ", t))
      # message(paste0("Message: ", w$message))
      NA
    }
  )))
  return(out)
}

survivalROC_tryCatch <- function(Stime, status, marker, predict.time, method, span, verbose = F){
  invisible(utils::capture.output(out <- tryCatch(
    # Specifying expression
    expr = {
      survivalROC::survivalROC(Stime = Stime, status = status, marker = marker,
                               predict.time = predict.time, method = method,
                               span = span)
    },
    # Specifying error message
    error = function(e){
      message(paste0(e$message, " for time ", predict.time, " and method '",pkg.env$AUC_survivalROC,"'"))
      NA
    },
    warning = function(e){
      NA
    }
  )))
  return(out)
}

risksetROC_tryCatch <- function(marker, Stime, status, predict.time, verbose = F){
  invisible(utils::capture.output(out <- tryCatch(
    # Specifying expression
    expr = {
      risksetROC::CoxWeights(marker = marker, Stime = Stime, status = status,
                             predict.time = predict.time)
    },
    # Specifying error message
    error = function(e){
      message(paste0(e$message, " for time ", predict.time, " and method '",pkg.env$AUC_risksetROC,"'"))
      NA
    }
  )))
  return(out)
}

getAUC_from_LP <- function(linear.predictors, Y, times, bestModel = NULL, method = "cenROC"){
  # method = c("risksetROC", "survivalROC")
  if(!all(c("time", "event") %in% colnames(Y))){
    stop_quietly("Data.frame Y must contain the columns time and event for COX model.")
  }

  if(is.null(linear.predictors)){
    return(NA)
  }else if(!isa(linear.predictors, "list")){
    aux = linear.predictors
    linear.predictors = NULL
    linear.predictors$fit = aux
  }else if(!"fit" %in% names(linear.predictors)){
    stop_quietly("fit must be a list inside linea predictors object.")
  }

  #order Y
  Y <- Y[names(linear.predictors$fit),]

  if(method %in% c(pkg.env$AUC_nsROC, pkg.env$AUC_cenROC)){
    Y <- data.matrix(Y)
  }

  AUC.cox = NULL #area under the curve for each timepoint

  #### ### ### ### ###
  #### NO PARALLEL ###
  #### ### ### ### ###
  t1 <- Sys.time()

  if(method %in% c(pkg.env$AUC_risksetROC, pkg.env$AUC_survivalROC, pkg.env$AUC_cenROC, pkg.env$AUC_nsROC)){
    for(t in times){
      #https://cran.r-project.org/web/packages/risksetROC/risksetROC.pdf
      if(method == pkg.env$AUC_risksetROC){
        # out <- risksetROC::CoxWeights(linear.predictors$fit, Stime = Y$time, status = Y$event,
        #                               predict.time = t)

        out <- risksetROC_tryCatch(linear.predictors$fit, Stime = Y$time, status = Y$event,
                                   predict.time = t)

        AUC <- out$AUC
      }else if(method == pkg.env$AUC_survivalROC){
        #https://cran.r-project.org/web/packages/survivalROC/
        #needs at least 2 events per time and time can not be day 0
        if(!sum(Y[Y$event==1,]$time<=t)<2 & t != 0){
          # out <- survivalROC::survivalROC(Stime = Y$time, status = Y$event, marker = linear.predictors$fit,
          #                                 predict.time = t, method = "NNE",
          #                                 span = 0.25 * nrow(Y)^(-0.20))

          out <- survivalROC_tryCatch(Stime = Y$time, status = Y$event, marker = linear.predictors$fit,
                                      predict.time = t, method = "NNE",
                                      span = 0.25 * nrow(Y)^(-0.20))

          AUC <- out$AUC
        }else{
          AUC <- NA
        }
      }else if(method == pkg.env$AUC_cenROC){
        #https://cran.r-project.org/web/packages/cenROC/
        #needs at least 2 events per time and time can not be day 0
        #Y is a matrix
        if(!sum(Y[Y[,"event"]==1,"time"]<=t)<2 & t != 0){
          # out <- tryCatch(
          #   # Specifying expression
          #   expr = {
          #     cenROC::cenROC(Y = Y[,"time"], censor = Y[,"event"], M = linear.predictors$fit,
          #                    t = t, method = "tra", ktype  ="normal", alpha = 0.05, plot = F) #problems in simulated data
          #   },
          #   # Specifying error message
          #   error = function(e){
          #     return(NULL)
          #   }
          # )
          out <- cenROC_tryCatch(Y = Y[,"time"], censor = Y[,"event"], M = linear.predictors$fit,
                                 t = t, method = "tra", ktype  ="normal", alpha = 0.05, plot = F, verbose = F)

          if(is.null(out)){
            AUC <- NA
          }else if(is.na(out)){
            AUC <- NA
          }else{
            AUC <- out$AUC$AUC
          }
        }else{
          AUC <- NA
        }
      }else if(method == pkg.env$AUC_nsROC){
        #https://cran.r-project.org/web/packages/nsROC/nsROC.pdf
        #needs at least 2 events per time and time can not be day 0
        #Y is a matrix
        #cumulative/dynamic approach
        if(!sum(Y[Y[,"event"]==1,"time"]<=t)<2 & t != 0){
          # out <- nsROC::cdROC(stime = Y[,"time"],
          #                     status = Y[,"event"],
          #                     marker = linear.predictors$fit,
          #                     predict.time = t,
          #                     method = "wKM", # Cox, KM, wKM (last with kernel-weighted)
          #                     kernel = "normal", #normal, Epanechnikov, other (only if wKM)
          #                     ci = F, #a lot of time and problems with NAs in bootstraping
          #                     boot.n = 200,
          #                     conf.level = 0.95,
          #                     seed = 123)

          out <- nsROC_tryCatch(stime = Y[,"time"],
                                status = Y[,"event"],
                                marker = linear.predictors$fit,
                                predict.time = t,
                                method = "wKM", # Cox, KM, wKM (last with kernel-weighted)
                                kernel = "normal", #normal, Epanechnikov, other (only if wKM)
                                ci = F, #a lot of time and problems with NAs in bootstraping
                                boot.n = 200,
                                conf.level = 0.95,
                                seed = 123)

          AUC <- out$auc
        }else{
          AUC <- NA
        }
      }
      AUC.cox <- c(AUC.cox, AUC)
    }

  }else if(method %in% c(pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I)){

    if(method == pkg.env$AUC_smoothROCtime_C){
      #https://cran.r-project.org/web/packages/smoothROCtime/
      new_times <- sapply(times, function(x){!sum(Y[Y$event==1,]$time<x)<2})

      if(!all(new_times == FALSE)){
        d <- cbind(Y$time, Y$event, linear.predictors$fit)
        #default naive.pdf for banwidth matrix

        # invisible(utils::capture.output(out <- tryCatch(
        #   # Specifying expression
        #   expr = {
        #     #smoothROCtime::stRoc(data = d, t = times[new_times], tcr = "C", meth = "2") #Cumulative/Dynamic with p-kernel method #too much problems
        #     smoothROCtime::stRoc(data = d, t = times[new_times], tcr = "C", meth = "1") #Cumulative/Dynamic with smooth method
        #   },
        #   # Specifying error message
        #   error = function(e){
        #     NULL
        #   }
        # )))

        out <- smoothROCtime_tryCatch(data = d, times = times[new_times], tcr = "C", meth = "1")

        if(is.null(out)){
          AUC <- NA
        }else if(is.na(out)){
          AUC <- NA
        }else{
          AUC <- NULL
          number_values <- table(out$t)[1]
          for(i in seq(1,length(out$auc), number_values)){ #get the auc for each time (remember each value is rep)
            AUC <- c(AUC, out$auc[i])
          }

          AUC.cox <- new_times
          AUC.cox[AUC.cox==T] <- as.numeric(AUC)
          AUC.cox[AUC.cox==0] <- NA
        }
      }else{
        AUC <- NA
      }

    }else if(method == pkg.env$AUC_smoothROCtime_I){
      #https://cran.r-project.org/web/packages/smoothROCtime/
      new_times <- sapply(times, function(x){!sum(Y[Y$event==1,]$time<x)<2})

      if(!all(new_times == FALSE)){
        d <- cbind(Y$time, Y$event, linear.predictors$fit)

        # invisible(utils::capture.output(out <- tryCatch(
        #   # Specifying expression
        #   expr = {
        #     smoothROCtime::stRoc(data = d, t = times[new_times], tcr = "I") #Incident/Dynamic
        #   },
        #   # Specifying error message
        #   error = function(e){
        #     NULL
        #   }
        # )))

        out <- smoothROCtime_tryCatch(data = d, times = times[new_times], tcr = "I")

        if(is.null(out)){
          AUC <- NA
        }else if(is.na(out)){
          AUC <- NA
        }else{
          AUC <- NULL
          number_values <- table(out$t)[1]
          for(i in seq(1,length(out$auc), number_values)){ #get the auc for each time (remember each value is rep)
            AUC <- c(AUC, out$auc[i])
          }

          AUC.cox <- new_times
          AUC.cox[AUC.cox==T] <- as.numeric(AUC)
          AUC.cox[AUC.cox==0] <- NA
        }
      }else{
        AUC <- NA
      }
    }

  }else{
    stop_quietly("No available method selected.")
  }

  t2 <- Sys.time()
  t2-t1

  return(list(lp.used = linear.predictors, AUC = AUC.cox))
}
