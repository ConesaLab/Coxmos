getSurvivalSubset <- function(X, Y, event.val = TRUE, EPV, p.censored, n.patients = NULL, num.subsets = 1, perc.test = 0.2, seed = 123){

  set.seed(seed)

  if(!"event" %in% colnames(Y) & !"status" %in% colnames(Y)){
    stop_quietly("'event' or 'status' column is not present in Y matrix/data.frame.")
  }else if(!"event" %in% colnames(Y)){
    colnames(Y)[colnames(Y)=="status"] <- "event"
  }

  if(!is.numeric(EPV)){
    stop_quietly("EPV must be a numeric variable.")
  }

  if(EPV < 5){
    message("Your EPV is low, considered increasing if you goal is to use the model as a predictor.")
  }

  if(!all(p.censored >= 0) | !all(p.censored <= 1)){
    stop_quietly("The p.censored parameter must be between 0 and less than 1.")
  }

  #both list
  if(length(EPV)>1 & length(p.censored)>1 & !length(EPV) == length(p.censored)){
    stop_quietly("Length of both list are not equal.")
  }

  #p.censored not list
  if(length(EPV)>1 & !length(p.censored)>1){
    p.censored <- rep(p.censored, length(EPV))
  }

  if(is.numeric(Y$event)){
    if(all(Y$event %in% c(0,1))){
      Y$event <- ifelse(Y$event==0, FALSE, TRUE)
    }
  }

  #EPV not list
  # if(!length(EPV)>1 & length(p.censored)>1){
  #   EPV <- rep(EPV, length(p.censored))
  # }

  aux_EPV = list()
  #We are going to fix the number of patients n=1000 and only modify the p.censored
  if(!is.null(n.patients)){

    #update n.patients to 0.8 train
    n.patients = floor(n.patients * (1-perc.test))

    EPV = NULL
    if(length(n.patients)>1){ #change patients each time

      if(length(p.censored)>1){
        stop_quietly("If n.patients is a list, just one value for p.censored must be provided.")
      }

      #EPV = dynamic EPV
      EPV = c(EPV, n.patients * (1-p.censored) / ncol(X))
      for(k in 1:length(EPV)){
        aux_censored <- list()
        for(j in 1:1){ #one element per EPV (just one set of censored patients meet the number of patients)
          aux_times <- list()
          for(i in 1:num.subsets){
            n.events = n.patients[[k]] * (1-p.censored[[j]])
            n.censored = n.patients[[k]]-n.events

            Y_event = Y[Y[,"event"]==event.val,]
            Y_censored = Y[Y[,"event"]!=event.val,]

            #event
            rows.event <- sample(nrow(Y_event), size = n.events)
            if(!length(rows.event)==n.events){
              rows.event <- sample(nrow(Y_event), size = n.events+1)
            }
            Y_subset <- Y_event[rows.event,]

            #censored
            rows.censored <- sample(nrow(Y_censored), size = n.censored)
            if(!length(rows.censored)==as.integer(n.censored)){
              rows.censored <- sample(nrow(Y_censored), size = n.censored+1)
            }
            Y_subset <- rbind(Y_subset, Y_censored[rows.censored,])
            #random rows
            rows <- sample(nrow(Y_subset))
            Y_subset <- Y_subset[rows,]
            #X data
            X_subset <- X[rownames(Y_subset),]

            #is Y factor an event?
            if(!is.factor(Y_subset[,"event"]) & !is.logical(Y_subset[,"event"]) & is.numeric(Y_subset[,"event"])){
              message("Y attribute is not logical, factor or numeric in 'event' column.")
            }

            #### ### ### ###
            #### ADD TEST ##
            #### ### ### ###

            #Y without train
            Y_event_aux <- Y_event[!rownames(Y_event) %in% rownames(Y_subset),]
            Y_censored_aux <- Y_censored[!rownames(Y_censored) %in% rownames(Y_subset),]

            #samples
            rows.event_test <- sample(nrow(Y_event_aux), size = n.events*perc.test) # 20% of train as test (test size)
            rows.censored_test <- sample(nrow(Y_censored_aux), size = n.censored*perc.test)

            Y_subset_test <- rbind(Y_event_aux[rows.event_test,], Y_censored_aux[rows.censored_test,])
            X_subset_test <- X[rownames(Y_subset_test),]

            #add to list
            aux_times[[i]] = list(X = X_subset, Y = Y_subset, X_test = X_subset_test, Y_test = Y_subset_test)
          }
          names(aux_times) <- paste0("set", 1:num.subsets)
          aux_censored[[j]] <- aux_times
        }
        names(aux_censored) <- paste0("censored", p.censored[[j]])
        aux_EPV[[k]] <- aux_censored
      }
    }else{
      #EPV = dynamic EPV
      EPV = n.patients * (1-p.censored) / ncol(X)
      #1 p.censored for each EPV
      for(k in 1:length(p.censored)){
        aux_censored <- list()
        for(j in 1:1){ #one element per EPV (just one set of censored patients meet the number of patients)
          aux_times <- list()
          for(i in 1:num.subsets){
            n.events = n.patients * (1-p.censored[[k]])
            n.censored = n.patients-n.events

            Y_event = Y[Y[,"event"]==event.val,]
            Y_censored = Y[Y[,"event"]!=event.val,]

            #event
            rows.event <- sample(nrow(Y_event), size = n.events)
            if(!length(rows.event)==n.events){
              rows.event <- sample(nrow(Y_event), size = n.events+1)
            }
            Y_subset <- Y_event[rows.event,]

            #censored
            rows.censored <- sample(nrow(Y_censored), size = n.censored)
            if(!length(rows.censored)==as.integer(n.censored)){
              rows.censored <- sample(nrow(Y_censored), size = n.censored+1)
            }
            Y_subset <- rbind(Y_subset, Y_censored[rows.censored,])
            #random rows
            rows <- sample(nrow(Y_subset))
            Y_subset <- Y_subset[rows,]
            #X data
            X_subset <- X[rownames(Y_subset),]

            #is Y factor an event?
            if(!is.factor(Y_subset[,"event"]) & !is.logical(Y_subset[,"event"]) & is.numeric(Y_subset[,"event"])){
              message("Y attribute is not logical, factor or numeric in 'event' column.")
            }

            #### ### ### ###
            #### ADD TEST ##
            #### ### ### ###

            #Y without train
            Y_event_aux <- Y_event[!rownames(Y_event) %in% rownames(Y_subset),]
            Y_censored_aux <- Y_censored[!rownames(Y_censored) %in% rownames(Y_subset),]

            #samples
            if(n.events*perc.test>=nrow(Y_event_aux)){
              rows.event_test <- 1:nrow(Y_event_aux) #if we require more, takes all patients
            }else{
              rows.event_test <- sample(1:nrow(Y_event_aux), size = n.events*perc.test) # 20% of train as test (test size)
            }

            if(n.censored*perc.test >= nrow(Y_censored_aux)){
              rows.censored_test <- 1:nrow(Y_censored_aux)
            }else{
              rows.censored_test <- sample(1:nrow(Y_censored_aux), size = n.censored*perc.test)
            }

            Y_subset_test <- rbind(Y_event_aux[rows.event_test,], Y_censored_aux[rows.censored_test,])
            X_subset_test <- X[rownames(Y_subset_test),]

            #add to list
            aux_times[[i]] = list(X = X_subset, Y = Y_subset, X_test = X_subset_test, Y_test = Y_subset_test)
          }
          names(aux_times) <- paste0("set", 1:num.subsets)
          aux_censored[[j]] <- aux_times
        }
        names(aux_censored) <- paste0("censored", p.censored[[k]])
        aux_EPV[[k]] <- aux_censored
      }
    }

  }else{
    for(k in 1:length(EPV)){
      bad_p.censored = NULL
      aux_censored <- list()
      for(j in 1:length(p.censored)){
        aux_times <- list()
        for(i in 1:num.subsets){
          n.events = EPV[[k]] * ncol(X)
          n.censored = floor(-1 * p.censored[[j]] * n.events) / (p.censored[[j]]-1) #round floor

          Y_event = Y[Y[,"event"]==event.val,]
          Y_censored = Y[Y[,"event"]!=event.val,]

          if(n.events > nrow(Y_event)){
            message("In (EPV ",EPV[k], ", %Cens: ", p.censored[j], " Set:", num.subsets[i],") - Number of events needed are larger than the number of event available (", n.events, ">", nrow(Y_event),")")
            bad_p.censored <- c(bad_p.censored, p.censored[j])
            next
          }else if(n.censored > nrow(Y_censored)){
            message("In (EPV ",EPV[k], ", %Cens: ", p.censored[j], " Set:", num.subsets[i],") - Number of censored needed are larger than the number of censored available (", n.censored, ">", nrow(Y_censored),")")
            bad_p.censored <- c(bad_p.censored, p.censored[j])
            next
          }

          #event
          rows.event <- sample(nrow(Y_event), size = n.events)
          Y_subset <- Y_event[rows.event,]
          #censored
          rows.censored <- sample(nrow(Y_censored), size = n.censored)
          Y_subset <- rbind(Y_subset, Y_censored[rows.censored,])
          #random rows
          rows <- sample(nrow(Y_subset))
          Y_subset <- Y_subset[rows,]
          #X data
          X_subset <- X[rownames(Y_subset),]

          #is Y factor an event?
          if(!is.factor(Y_subset[,"event"]) & !is.logical(Y_subset[,"event"]) & is.numeric(Y_subset[,"event"])){
            message("Y attribute is not logical, factor or numeric in 'event' column.")
          }

          #### ### ### ###
          #### ADD TEST ##
          #### ### ### ###

          #Y without train
          Y_event_aux <- Y_event[!rownames(Y_event) %in% rownames(Y_subset),]
          Y_censored_aux <- Y_censored[!rownames(Y_censored) %in% rownames(Y_subset),]

          #samples
          if(n.events*perc.test>=nrow(Y_event_aux)){
            rows.event_test <- 1:nrow(Y_event_aux) #if we require more, takes all patients
          }else{
            rows.event_test <- sample(1:nrow(Y_event_aux), size = n.events*perc.test) # 20% of train as test (test size)
          }

          if(n.censored*perc.test >= nrow(Y_censored_aux)){
            rows.censored_test <- 1:nrow(Y_censored_aux)
          }else{
            rows.censored_test <- sample(1:nrow(Y_censored_aux), size = n.censored*perc.test)
          }

          Y_subset_test <- rbind(Y_event_aux[rows.event_test,], Y_censored_aux[rows.censored_test,])
          X_subset_test <- X[rownames(Y_subset_test),]

          #add to list
          aux_times[[i]] = list(X = X_subset, Y = Y_subset, X_test = X_subset_test, Y_test = Y_subset_test)
        }
        if(length(aux_times)==0){
          next #more censored or events that we have
        }
        names(aux_times) <- paste0("set", 1:num.subsets)
        aux_censored[[j]] <- aux_times
      }
      if(length(aux_censored)==0){
        next #more censored or events that we have
      }
      names(aux_censored) <- paste0("censored", p.censored[!p.censored %in% bad_p.censored])
      aux_EPV[[k]] <- aux_censored
    }
  }

  names(aux_EPV) <- paste0("EPV", round(EPV, 4)) #more numbers for omic datasets, the get the same name for different EPV

  return(aux_EPV)
}

#' super.trainAllModels
#'
#' @param lst_subdata List of data (X, Y)
#' @param methods Methods to use for train
#' @param comp_calculation Compute optimal number of components by "manual" or "auto".
#' @param ncomp N.components for manual detection.
#' @param EN.alpha Penalization for manual detection,
#' @param eta Eta for manual detection,
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final cox model. Used to restrict the number of variables/components can be computed in final cox models. If the minimum is not meet, the model cannot be computed (default: 5).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#'
#'  @export

super.trainAllModels <- function(lst_subdata, methods,
                                 comp_calculation = "manual",
                                 ncomp = 5, EN.alpha = 0.5, eta = 0.5,
                                 x.center = T, x.scale = T,
                                 MIN_EPV = 0, PARALLEL = T){

  #test
  y.center = y.scale = FALSE
  lst_params <- list()

  for(e in names(lst_subdata)){ #EPV
    for(c in names(lst_subdata[[e]])){ #% censored
      for(s in names(lst_subdata[[e]][[c]])){ #sets
        name_m = paste0(e,"_",c,"_",s)

        l <- list()
        l[["EPV"]] <- e
        l[["censored"]] <- c
        l[["set"]] <- s

        lst_params[[name_m]] <- l
      }
    }
  }

  t1 <- Sys.time()
  if(PARALLEL){
    if(.Platform$OS.type == "unix") {
      future::plan("multicore", workers = min(future::availableCores()-1, length(lst_params)))
    }else{
      future::plan("multisession", workers = min(future::availableCores()-1, length(lst_params)))
    }


    lst_subdata_models <- furrr::future_map(lst_params, ~train_all_models2.5(
                                                           lst_X_train = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$X,
                                                           lst_Y_train = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$Y,
                                                           methods = methods,
                                                           comp_calculation = comp_calculation,
                                                           ncomp = ncomp, EN.alpha = EN.alpha, eta = eta,
                                                           x.center = x.center, x.scale = x.scale,
                                                           y.center = y.center, y.scale = y.scale, MIN_EPV = MIN_EPV),
                                             .options = furrr_options(seed = T))

    future::plan("sequential")
  }else{
    lst_subdata_models <- purrr::map(lst_params, ~train_all_models2.5(
                                     lst_X_train = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$X,
                                     lst_Y_train = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$Y,
                                     methods = methods,
                                     comp_calculation = comp_calculation,
                                     ncomp = ncomp, EN.alpha = EN.alpha, eta = eta,
                                     x.center = x.center, x.scale = x.scale,
                                     y.center = y.center, y.scale = y.scale, MIN_EPV = MIN_EPV))
  }

  t2 <- Sys.time()
  t_training <- t2-t1

  return(lst_subdata_models)
}

#' super.evalAllModels
#'
#' @param lst_subdata List of data
#' @param lst_subdata_models List of model per each data
#' @param lst_evaluations List of which evaluators to use: "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I", "risksetROC"
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param progress_bar Logical. If progress_bar = TRUE, progress bar is shown (default = FALSE).
#'
#' @export
super.evalAllModels <- function(lst_subdata, lst_subdata_models, lst_evaluations,
                                pred.attr = "mean", times = NULL, max_time_points = 15,
                                PARALLEL = F, progress_bar = F){

  #test
  lst_models_to_test <- list()

  for(e in names(lst_subdata)){ #EPV
    for(c in names(lst_subdata[[e]])){ #% censored
      for(s in names(lst_subdata[[e]][[c]])){ #sets
        for(eval in lst_evaluations){
          name_m = paste0(e,"_",c,"_",s)

          l <- list()
          l[["EPV"]] <- e
          l[["censored"]] <- c
          l[["set"]] <- s
          l[["evaluator"]] <- eval
          l[["model"]] <- name_m

          name = paste0(e,"_",c,"_",s,"_",eval)
          lst_models_to_test[[name]] <- l
        }
      }
    }
  }

  #PARALLEL
  t1 <- Sys.time()
  if(PARALLEL){
    if(.Platform$OS.type == "unix") {
      future::plan("multicore", workers = min(future::availableCores()-1, length(lst_models_to_test)))
    }else{
      future::plan("multisession", workers = min(future::availableCores()-1, length(lst_models_to_test)))
    }

    eval_results <- furrr::future_map(lst_models_to_test, ~eval_HDcox_models(lst_models = lst_subdata_models[[.$model]],
                                                                          X_test = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$X_test,
                                                                          Y_test = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$Y_test,
                                                                          pred.method = .$evaluator,
                                                                          pred.attr = pred.attr,
                                                                          times = times, max_time_points = max_time_points,
                                                                          PARALLEL = F,
                                                                          progress_bar = progress_bar, verbose = F))

    future::plan("sequential")

  }else{
    eval_results <- purrr::map(lst_models_to_test, ~eval_HDcox_models(lst_models = lst_subdata_models[[.$model]],
                                                                       X_test = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$X_test,
                                                                       Y_test = lst_subdata[[.$EPV]][[.$censored]][[.$set]]$Y_test,
                                                                       pred.method = .$evaluator,
                                                                       pred.attr = pred.attr,
                                                                       times = times, max_time_points = max_time_points,
                                                                       PARALLEL = F,
                                                                       progress_bar = progress_bar, verbose = F))
  }

  t2 <- Sys.time()
  t_evaluation <- t2-t1
  return(eval_results)

}

#' super.evalResults2DataFrame
#'
#' @param eval_results Eval results
#'
#' @export
super.evalResults2DataFrame <- function(eval_results){
  ### I need to transform the results into one data.frame including the EPV, %censored, set
  df.results_subdata = NULL

  for(n in names(eval_results)){
    txt <- strsplit(n, "_")[[1]]

    if(length(txt)==5){
      eval = paste0(txt[[4]], "_", txt[[5]])
    }else{
      eval = txt[[4]]
    }

    e <- gsub("EPV", "", txt[[1]])
    p <- gsub("censored", "", txt[[2]])
    s <- gsub("set", "", txt[[3]])
    eval <- eval

    row <- c(e,p,s,eval)
    rows <- NULL
    for(nr in 1:ifelse(!is.null(nrow(eval_results[[n]]$df)), nrow(eval_results[[n]]$df), 1)){
      rows <- rbind(rows, row)
    }

    if(!is.null(eval_results[[n]]$df)){
      auc <- eval_results[[n]]$df
    }else{
      auc <- matrix(rep(NA, 7), nrow = 1)
    }

    rows <- as.data.frame(cbind(rows, auc))
    rownames(rows) <- NULL

    if(!is.null(eval_results[[n]]$df)){
      rows$AUC <- as.numeric(rows$AUC)
    }

    colnames(rows) <- colnames(df.results_subdata) #in some cases, smoothROCtime_I return all NULL

    df.results_subdata <- rbind(df.results_subdata, rows)
    colnames(df.results_subdata) <- c("EPV", "p.censored", "set", "eval", "method",
                                      "training.time","evaluating.time", "AIC", "c.index", "time", "AUC")
  }

  lvl <- as.numeric(unique(df.results_subdata$set))
  lvl <- lvl[order(lvl)]
  df.results_subdata$set <- factor(df.results_subdata$set, levels = lvl)

  lvl <- as.numeric(unique(df.results_subdata$EPV))
  lvl <- lvl[order(lvl)]
  df.results_subdata$EPV <- factor(df.results_subdata$EPV, levels = lvl)

  lvl <- as.numeric(unique(df.results_subdata$p.censore))
  lvl <- lvl[order(lvl)]
  df.results_subdata$p.censored <- factor(df.results_subdata$p.censored, levels = lvl)

  df.results_subdata$eval <- factor(df.results_subdata$eval)
  df.results_subdata$method <- factor(df.results_subdata$method)
  df.results_subdata$time <- factor(df.results_subdata$time)

  df.results_subdata$training.time <- as.numeric(df.results_subdata$training.time)
  df.results_subdata$evaluating.time <- as.numeric(df.results_subdata$evaluating.time)
  df.results_subdata$AIC <- as.numeric(df.results_subdata$AIC)
  df.results_subdata$c.index <- as.numeric(df.results_subdata$c.index)
  df.results_subdata$AUC <- as.numeric(df.results_subdata$AUC)

  df.results_subdata <- droplevels.data.frame(df.results_subdata)

  return(df.results_subdata)
}

train_all_models2.5 <- function(lst_X_train, lst_Y_train,
                                methods = c("cox", "coxSW", "coxEN", "sPLS-ICOX", "sPLS-DRCOX", "sPLS-DRCOX-Dynamic", "sPLS-DACOX-Dynamic"),
                                ncomp = 5, EN.alpha = 0.5, eta = 0.5, comp_calculation = "manual",
                                n_run = 2, k_folds = 10, fast_mode = F, pred.method = "cenROC",
                                x.center = TRUE, x.scale = FALSE,
                                y.center = FALSE, y.scale = FALSE, MIN_EPV = 0,
                                vector = NULL,
                                MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                EVAL_METHOD = "cenROC"){

  X_train <- lst_X_train
  Y_train <- lst_Y_train

  if(ncol(X_train) == 0  | ncol(Y_train) == 0){
    stop("X_train or Y_train have 0 dimensions.")
  }

  #Algorithms Parameters
  x.center = x.center
  x.scale = x.scale
  y.center = y.center
  y.scale = y.scale
  max.ncomp = ncomp

  #Weights Parameters
  w_AIC = 0
  w_c.index = 0
  w_BRIER = 0
  w_AUC = 1

  #Eval stop detection
  MIN_AUC_INCREASE = 0.01
  MIN_AUC = 0.75
  MIN_COMP_TO_CHECK = 3

  #others
  return_models = F
  MIN_EPV = MIN_EPV
  pred.attr = "mean"
  seed = 123
  alpha = 0.05
  remove_non_significant = F #cox and coxEN variables
  remove_non_significant_models = F  #plscox methods
  times = NULL

  #mixomics num. of variables to study
  test.keepX <- NULL
  max.var <- ncol(X_train)
  if(max.var > 20){
    test.keepX <- c(test.keepX, seq(1, max.var, ceiling(max.var/20)))
  }else{
    test.keepX <- c(test.keepX, seq(1, max.var, 1))
  }

  if(!comp_calculation %in% c("auto", "manual")){
    stop_quietly(paste0("Parameter 'comp_calculation' must be one of 'auto' or 'manual' options."))
  }

  auto=F
  if(comp_calculation == "auto"){
    auto=T
  }

  #### ### ### #
  ### MODELS ###
  #### ### ### #

  lst_res <- list()
  for(m in methods){
    lst_res[[m]] <- list()
  }

  #COX
  if(pkg.env$cox %in% methods){

    res_cox <- cox(X = X_train, Y = Y_train,
                   x.center = x.center, x.scale = x.scale,
                   #y.center = y.center, y.scale = y.scale,
                   MIN_EPV = MIN_EPV, #by default 0
                   remove_non_significant = F, alpha = 0.05,
                   FORCE = T, returnData = F, verbose = F)
  }else{
    res_cox <- NA
  }

  #COX
  if(pkg.env$coxSW %in% methods){

    res_coxSW <- coxSW(X = X_train, Y = Y_train,
                       max.variables = ncol(X_train), BACKWARDS = T,
                       x.center = x.center, x.scale = x.scale,
                       #y.center = y.center, y.scale = y.scale,
                       MIN_EPV = MIN_EPV, #by default 0
                       returnData = F, verbose = F)

  }else{
    res_coxSW <- NA
  }

  #COX
  if(pkg.env$coxEN %in% methods){

    if(auto){
      cv.coxEN_res <- cv.coxEN(X = X_train, Y = Y_train,
                               max.variables = ncol(X_train), EN.alpha.list = seq(0,1,0.1), #EN penalization
                               n_run = n_run, k_folds = k_folds,
                               alpha = alpha, remove_non_significant = remove_non_significant, times = times,
                               w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC,
                               MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                               x.scale = x.scale, x.center = x.center,
                               #y.scale = y.scale, y.center = y.center,
                               fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                               pred.attr = pred.attr, pred.method = pred.method, seed = seed)

      ### Optimal number of components
      res_coxEN <- coxEN(X = X_train, Y = data.matrix(Y_train), EN.alpha = cv.coxEN_res$opt.EN.alpha,
                         x.center = x.center, x.scale = x.scale,
                         #y.center = y.center, y.scale = y.scale,
                         MIN_EPV = MIN_EPV, #by default 0
                         remove_non_significant = F, alpha = 0.05, returnData = F)

    }else{

      res_coxEN <- coxEN(X = X_train, Y = Y_train, EN.alpha = EN.alpha,
                         x.center = x.center, x.scale = x.scale,
                         #y.center = y.center, y.scale = y.scale,
                         MIN_EPV = MIN_EPV, #by default 0
                         remove_non_significant = F, alpha = 0.05, returnData = F)

    }

  }else{
    res_coxEN <- NA
  }

  #splsicox-manual
  if(pkg.env$splsicox %in% methods){

    if(auto){
      cv.splsicox_res <- cv.splsicox(X = X_train, Y = data.matrix(Y_train),
                                     max.ncomp = max.ncomp,
                                     n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                     w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     x.scale = x.scale, x.center = x.center,
                                     #y.scale = y.scale, y.center = y.center,
                                     fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                     pred.attr = pred.attr, pred.method = pred.method, seed = seed)

      ### Optimal number of components
      res_splsicox <- splsicox(X = X_train, Y = data.matrix(Y_train),
                             n.comp = cv.splsicox_res$opt.comp,
                             x.center = x.center, x.scale = x.scale,
                             #y.center = y.center, y.scale = y.scale,
                             returnData = F)

    }else{
      res_splsicox <- splsicox(X = X_train, Y = data.matrix(Y_train),
                             n.comp = ncomp,
                             x.center = x.center, x.scale = x.scale,
                             #y.center = y.center, y.scale = y.scale,
                             returnData = F)
    }
  }else{
    res_splsicox <- NA
  }

  #splsdrcox-manual
  if(pkg.env$splsdrcox %in% methods){

    if(auto){
      cv.splsdrcox_res <- cv.splsdrcox(X = X_train, Y = data.matrix(Y_train),
                                     max.ncomp = max.ncomp, eta.list = c(0, 0.25, 0.5, 0.75),
                                     n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                     w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     x.scale = x.scale, x.center = x.center,
                                     #y.scale = y.scale, y.center = y.center,
                                     fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                     pred.attr = pred.attr, pred.method = pred.method, seed = seed)

      res_splsdrcox <- splsdrcox(X = X_train,
                               Y = data.matrix(Y_train),
                               n.comp = cv.splsdrcox_res$opt.comp,
                               eta = cv.splsdrcox_res$opt.eta,
                               x.center = x.center, x.scale = x.scale,
                               #y.scale = y.scale, y.center = y.center,
                               returnData = F)

    }else{
      #solo un spls
      res_splsdrcox <- splsdrcox(X = X_train,
                               Y = data.matrix(Y_train),
                               n.comp = ncomp,
                               eta = eta,
                               x.center = x.center, x.scale = x.scale,
                               #y.scale = y.scale, y.center = y.center,
                               returnData = F)

    }
  }else{
    res_splsdrcox <- NA
  }

  if(pkg.env$splsdrcox_dynamic %in% methods){

    if(auto){
      cv.splsdrcox_dynamic_res <- cv.splsdrcox_dynamic(X = X_train, Y = data.matrix(Y_train),
                                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                                       alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                       w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                       x.scale = x.scale, x.center = x.center,
                                                       #y.scale = y.scale, y.center = y.center,
                                                       vector = vector,
                                                       MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                       EVAL_METHOD = EVAL_METHOD,
                                                       fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                       pred.attr = pred.attr, pred.method = pred.method, seed = seed)

      res_splsdrcox_dynamic <- splsdrcox_dynamic(X = X_train,
                                                 Y = data.matrix(Y_train),
                                                 n.comp = cv.splsdrcox_dynamic_res$opt.comp,
                                                 vector = cv.splsdrcox_dynamic_res$opt.nvar,
                                                 MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                 EVAL_METHOD = EVAL_METHOD,
                                                 x.center = x.center, x.scale = x.scale,
                                                 #y.scale = y.scale, y.center = y.center,
                                                 returnData = F)

    }else{
      #solo un spls
      res_splsdrcox_dynamic <- splsdrcox_dynamic(X = X_train,
                                                 Y = data.matrix(Y_train),
                                                 n.comp = ncomp,
                                                 vector = vector,
                                                 MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                 EVAL_METHOD = EVAL_METHOD,
                                                 x.center = x.center, x.scale = x.scale,
                                                 #y.scale = y.scale, y.center = y.center,
                                                 returnData = F)

    }
  }else{
    res_splsdrcox_dynamic <- NA
  }

  #splsda+cox
  if(pkg.env$splsdacox_dynamic %in% methods){

    event <- Y_train[,"event"]
    if(!is.factor(event)){
      event <- factor(event)
    }

    if(length(levels(event))==1){ #case 100% event [we cannot discriminate between classes]
      res_splsdacox_dynamic = NA
    }else{

      if(auto){
        cv.splsdacox_dynamic_res <- cv.splsdacox_dynamic(X = X_train, Y = data.matrix(Y_train),
                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                       alpha = alpha, remove_non_significant_models = remove_non_significant_models, max.iter = 500,
                                       w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       x.scale = x.scale, x.center = x.center,
                                       #y.scale = y.scale, y.center = y.center,
                                       vector = vector,
                                       MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                       EVAL_METHOD = EVAL_METHOD,
                                       fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                       pred.attr = pred.attr, pred.method = pred.method, seed = seed)

        res_splsdacox_dynamic <- splsdacox_dynamic(X_train, Y_train,
                                 n.comp = cv.splsdacox_dynamic_res$opt.comp,
                                 vector = cv.splsdacox_dynamic_res$opt.nvar,
                                 MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                 EVAL_METHOD = EVAL_METHOD,
                                 x.center = x.center, x.scale = x.scale,
                                 #y.center = y.center, y.scale = y.scale,
                                 max.iter = 500, returnData = F)

      }else{
        res_splsdacox_dynamic <- splsdacox_dynamic(X_train, Y_train,
                                 n.comp = ncomp,
                                 vector = vector,
                                 MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                 EVAL_METHOD = EVAL_METHOD,
                                 x.center = x.center, x.scale = x.scale,
                                 #y.center = y.center, y.scale = y.scale,
                                 max.iter = 500, returnData = F)

      }

    }#more than one Y$event class
  }else{
    res_splsdacox_dynamic <- NA
  }

  #### Save cox_models
  if(!length(res_cox) == 1){
    lst_res[[pkg.env$cox]] <- res_cox
  }else{
    lst_res[[pkg.env$cox]] <- NA
  }

  if(!length(res_coxSW) == 1){
    lst_res[[pkg.env$coxSW]] <- res_coxSW
  }else{
    lst_res[[pkg.env$coxSW]] <- NA
  }

  if(!length(res_coxEN) == 1){
    lst_res[[pkg.env$coxEN]] <- res_coxEN
  }else{
    lst_res[[pkg.env$coxEN]] <- NA
  }

  if(!length(res_splsicox) == 1){
    lst_res[[pkg.env$splsicox]] <- res_splsicox
  }else{
    lst_res[[pkg.env$splsicox]] <- NA
  }

  if(!length(res_splsdrcox) == 1){
    lst_res[[pkg.env$splsdrcox]] <- res_splsdrcox
  }else{
    lst_res[[pkg.env$splsdrcox]] <- NA
  }

  if(!length(res_splsdrcox_dynamic) == 1){
    lst_res[[pkg.env$splsdrcox_dynamic]] <- res_splsdrcox_dynamic
  }else{
    lst_res[[pkg.env$splsdrcox_dynamic]] <- NA
  }

  if(!length(res_splsdacox_dynamic) == 1){
    lst_res[[pkg.env$splsdacox_dynamic]] <- res_splsdacox_dynamic
  }else{
    lst_res[[pkg.env$splsdacox_dynamic]] <- NA
  }

  return(lst_res)
}

boxplot.performance <- function(df, x.var, y.var, x.fill = NULL, x.alpha = NULL, x.lab = NULL, y.lab = NULL, fill.lab = NULL, alpha.lab = NULL, title = NULL, y.limit = NULL, y.limit.exception = NULL, jitter = T, test = "anova", eval_method = "auto", show.median = T, round.median = 3){

  if(!eval_method %in% c("median", "mean", "auto")){
    stop("Eval_method must be one of: 'mean' or 'median'.")
  }

  if(eval_method == "auto"){
    if(!is.null(test) && test %in% c("t.test", "anova")){
      eval_method = "mean"
    }else{
      eval_method = "median"
    }
  }

  df <- df[,unique(c(x.var, y.var, x.fill, x.alpha))]
  df <- df[!is.na(df[,x.var]),]

  #remove NA before get the comparisons
  df <- df[!is.na(df[,y.var]),]
  #drop levels with 0 values
  df <- droplevels.data.frame(df)

  if(!is.null(x.fill)){
    df <- df[!is.na(df[,x.fill]),]
  }

  if(!is.null(x.alpha)){
    df <- df[!is.na(df[,x.alpha]),]
  }

  max <- max(df[!is.na(df[,y.var,drop=T]),y.var,drop=T])

  tests <- c("t.test","wilcox.test","anova","kruskal.test")
  comparisons <-  list()
  cont = 1
  if(!is.null(test)){
    for(i in 1:(length(levels(df[,x.var, drop=T]))-1)){
      for(j in (i+1):length(levels(df[,x.var, drop=T]))){
        comparisons[[cont]] <- c(levels(df[,x.var, drop=T])[i], levels(df[,x.var, drop=T])[j])
        cont = cont + 1
      }
    }
    if(!test %in% tests){
      stop_quietly(paste0("Variables test must be one of the following: ", paste0(tests, collapse = ", ")))
    }
  }

  median.val <- NULL
  for(m in levels(df[,x.var, drop=T])){
    sub_value <- df[df[,x.var, drop=T]==m,y.var,drop=T]
    if(eval_method=="median"){
      median.val <- c(median.val, median(sub_value, na.rm = T))
    }else{
      median.val <- c(median.val, mean(sub_value, na.rm = T))
    }
  }

  if(!is.null(median.val)){
    names(median.val) <- levels(df[,x.var,drop=T])
    median.val <- round(median.val, round.median)

    if(eval_method=="median"){
      x_names <- paste0(levels(df[,x.var,drop=T]), "\nMedian: ", median.val)
    }else{
      x_names <- paste0(levels(df[,x.var,drop=T]), "\nMean: ", median.val)
    }
  }

  if(is.null(x.fill)){
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, fill = x.var, alpha = x.alpha)) +
      geom_boxplot() +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab)) +
      theme(legend.position = "none") +

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "complete")
      }

      if(jitter){
        ggp <- ggp + geom_jitter(color="black", size=1, alpha=0.25, width = 0.2)
      }

  }else{
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, fill = x.fill, alpha = x.alpha)) +
      geom_boxplot(position = position_dodge2(preserve = "single")) +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab)) +
      theme(legend.position = "bottom")

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "complete")
    }

    if(jitter){
      ggp <- ggp + geom_point(position=position_jitterdodge(), color="black", size=1, alpha=0.25)
    }

  }

  if(!is.null(x.alpha)){
    dim_alpha <- length(levels(df[,x.alpha,drop=T]))

    if(!dim_alpha==1){
      s <- 1/dim_alpha
      alpha_values <- seq(1,0+s,-s)

      ggp <- ggp +
        scale_alpha_manual(values=alpha_values) +
        guides(alpha=guide_legend(override.aes=list(fill=grDevices::hcl(c(177,177),74.7,32.5,alpha=alpha_values), colour=NA))) #should be a base R package
    }
  }

  if(!is.null(y.limit) & !y.var %in% y.limit.exception){
    ggp <- ggp + scale_y_continuous(limits = y.limit, n.breaks = 15)
  }else{
    ggp <- ggp + scale_y_continuous(n.breaks = 15)
  }

  if(!is.null(test)){ #with less than
    ggp <- tryCatch(
      # Specifying expression
      expr = {
        if(test=="anova" | test=="kruskal.test"){
          ggp <- ggp + ggpubr::stat_compare_means(method = test, label.x.npc = "center", label.y = 1.025*max)
          ggp
        }else if(length(unique(unlist(comparisons)))==2){
          ggp <- ggp + ggpubr::stat_compare_means(method = test, label.x.npc = "center", label.y = 1.025*max)
          ggp
        }else{
          #some input is generated but I do not want it to be printed.
          output_txt <- capture.output(ggp <- ggp + ggpubr::stat_compare_means(method = test, comparisons = comparisons))
          ggp
        }
      },
      # Specifying error message
      error = function(e){
        ggp
      },
      # Specifying warning message
      warning = function(e){
        ggp
      }
    )
  }

  if(show.median & !is.null(median.val)){
    ggp <- ggp + scale_x_discrete(labels = x_names)
  }

  if(!is.null(fill.lab)){
    ggp <- ggp + guides(fill=guide_legend(title=fill.lab))
  }

  if(!is.null(alpha.lab)){
    ggp <- ggp + guides(alpha=guide_legend(title=alpha.lab))
  }

  if(!is.null(title)){
    ggp <- ggp + ggtitle(title)
  }

  return(ggp)
}

lineplot.performace <- function(df, x.var = "time", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, point = T){
  if(point){
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, color = x.color)) +
      geom_line(aes_string(group = x.color), size = 1) +
      geom_point() +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_color_conesa(palette = "complete")
    }

  }else{
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, color = x.color)) +
      geom_line(aes_string(group = x.color), size = 1) +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_color_conesa(palette = "complete")
    }

  }

  if(length(df[,x.var]>30)){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  if(!is.null(y.limit)){
    ggp <- ggp + ylim(y.limit)
  }

  return(ggp)
}
