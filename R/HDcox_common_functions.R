#' @importFrom caret nearZeroVar createFolds
#' @importFrom cowplot plot_grid
#' @importFrom future availableCores plan
#' @import ggrepel
#' @import ggplot2
#' @importFrom ggpubr ggarrange annotate_figure
#' @import glmnet
#' @import progress
#' @import purrr
#' @importFrom scattermore geom_scattermore
#' @import stats
#' @import survival
#@importFrom survAUC predErr
#' @importFrom survcomp sbrier.score2proba
#' @import survminer
#' @importFrom tidyr pivot_longer starts_with
#' @import utils
#'
#' @importFrom mixOmics spls plsda block.spls block.splsda tune.spls
#' @import furrr

pkg.env <- new.env(parent = emptyenv())
assign(x = 'model_class', value = "HDcox", pkg.env)
assign(x = 'eval_class', value = "evaluation", pkg.env)

assign(x = 'cox', value = c("cox"), pkg.env)
assign(x = 'coxSW', value = c("coxSW"), pkg.env)
assign(x = 'coxEN', value = c("coxEN"), pkg.env)
assign(x = 'classical_methods', value = c(pkg.env$cox, pkg.env$coxSW, pkg.env$coxEN), pkg.env)

assign(x = 'splsicox', value = c("sPLS-ICOX"), pkg.env)
assign(x = 'splsdrcox', value = c("sPLS-DRCOX"), pkg.env)
assign(x = 'splsdrcox_dynamic', value = c("sPLS-DRCOX-Dynamic"), pkg.env)
assign(x = 'splsdacox_dynamic', value = c("sPLS-DACOX-Dynamic"), pkg.env)
assign(x = 'pls_methods', value = c(pkg.env$splsicox, pkg.env$splsdrcox, pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic), pkg.env)

assign(x = 'sb.splsicox', value = c("SB.sPLS-ICOX"), pkg.env)
assign(x = 'sb.splsdrcox', value = c("SB.sPLS-DRCOX"), pkg.env)
assign(x = 'mb.splsdrcox', value = c("MB.sPLS-DRCOX"), pkg.env)
assign(x = 'mb.splsdacox', value = c("MB.sPLS-DACOX"), pkg.env)
assign(x = 'multiblock_methods', value = c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox), pkg.env)
assign(x = 'all_methods', value = c(pkg.env$classical_methods, pkg.env$pls_methods, pkg.env$multiblock_methods), pkg.env)

assign(x = 'cv.coxEN', value = c("cv.coxEN"), pkg.env)
assign(x = 'classical_cv', value = pkg.env$cv.coxEN, pkg.env)

assign(x = 'cv.splsicox', value = c("cv.sPLS-ICOX"), pkg.env)
assign(x = 'cv.splsdrcox', value = c("cv.sPLS-DRCOX"), pkg.env)
assign(x = 'cv.splsdrcox_dynamic', value = c("cv.sPLS-DRCOX-Dynamic"), pkg.env)
assign(x = 'cv.splsdacox_dynamic', value = c("cv.sPLS-DACOX-Dynamic"), pkg.env)
assign(x = 'pls_cv', value = c(pkg.env$cv.splsicox, pkg.env$cv.splsdrcox, pkg.env$cv.splsdrcox_dynamic, pkg.env$cv.splsdacox_dynamic), pkg.env)

assign(x = 'cv.sb.splsicox', value = c("cv.SB.sPLS-ICOX"), pkg.env)
assign(x = 'cv.sb.splsdrcox', value = c("cv.SB.sPLS-DRCOX"), pkg.env)
assign(x = 'cv.mb.splsdrcox', value = c("cv.MB.sPLS-DRCOX"), pkg.env)
assign(x = 'cv.mb.splsdacox', value = c("cv.MB.sPLS-DACOX"), pkg.env)
assign(x = 'multiblock_cv', value = c(pkg.env$cv.sb.splsicox, pkg.env$cv.sb.splsdrcox, pkg.env$cv.mb.splsdrcox, pkg.env$cv.mb.splsdacox), pkg.env)
assign(x = 'all_cv', value = c(pkg.env$classical_cv, pkg.env$pls_cv, pkg.env$multiblock_cv), pkg.env)

assign(x = 'AUC_survivalROC', value = c("survivalROC"), pkg.env)
assign(x = 'AUC_cenROC', value = c("cenROC"), pkg.env)
assign(x = 'AUC_nsROC', value = c("nsROC"), pkg.env)
assign(x = 'AUC_smoothROCtime_C', value = c("smoothROCtime_C"), pkg.env)
assign(x = 'AUC_smoothROCtime_I', value = c("smoothROCtime_I"), pkg.env)
assign(x = 'AUC_risksetROC', value = c("risksetROC"), pkg.env)
assign(x = 'AUC_evaluators', value = c(pkg.env$AUC_survivalROC, pkg.env$AUC_cenROC, pkg.env$AUC_nsROC, pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I, pkg.env$AUC_risksetROC), pkg.env)

assign(x = 'Brier', value = c("Brier_score"), pkg.env)

#assign(x = 'IllegalChars', value = c("`", "*"), pkg.env)
assign(x = 'IllegalChars', value = c("`"), pkg.env)


#' norm01
#' @description Normalize all values into 0-1 range.
#' @param x matrix or data.frame
norm01 <- function(x){
  if(max(x)-min(x) != 0){
    return((x-min(x))/(max(x)-min(x)))
  }else{
    return(x/length(x))
  }
}

#' print.HDcox
#'
#' @param x HDcox object
#' @param ... further arguments passed to or from other methods.
#'
#' @export

print.HDcox <- function(x, ...){

  if(attr(x, "model") %in% pkg.env$all_methods){

    cat("The method used is ", attr(x, "model"), ".\n\n", sep = "")

    if("removed_variables" %in% names(x) && !is.null(x$removed_variables)){
      cat("A total of ", length(x$nzv), " variables have been removed due to Zero or Near-Zero Variance filter", ".\n\n", sep = "")
    }

    if("removed_variables_cox" %in% names(x) && !is.null(x$removed_variables_cox)){
      cat("A total of ", length(x$removed_variables_cox), " variables have been removed due to non-significance filter inside cox model", ".\n\n", sep = "")
    }

    if(!all(is.null(x$survival_model))){
      cat("Survival model:\n", sep = "")
      tab <- summary(x$survival_model$fit)[[7]]

      print(tab)
    }else{
      cat("Survival model could not be computed.\n\n", sep = "")
    }

  }else if(attr(x, "model") %in% pkg.env$all_cv){

    cat("The cross validation method used is ", attr(x, "model"), ".\n\n", sep = "")

    if(!is.null(x$best_model_info)){
      cat("Best model is:\n\n")
      print(x$best_model_info)
    }else{
      cat("Best model could NOT be obtained. All models computed present problems.\n\n")
    }

  }else if(attr(x, "model") %in% pkg.env$eval_class){

    cat(paste0("Evaluation performed for methods: ", paste0(levels(x$df$method), collapse = ", "), ".\n\n", sep = ""))

    for(m in levels(x$df$method)){
      cat(paste0(m,": \n"))
      for(c in colnames(x$df)){

        if(c=="method"){
          next
        }else if(c=="time"){
          time_vector <- levels(x$df[x$df$method==m,c,drop=T])
          time_vector <- unlist(lapply(time_vector, function(x){gsub("time_", "", x)[[1]]}))
          cat(paste0("\t",c,": ", paste0(time_vector, collapse = ", "), "\n"))
        }else{
          ave <- mean(x$df[x$df$method==m,c,drop=T], na.rm = T)
          cat(paste0("\t",c,": ", round(ave, 5), "\n"))
        }

      }
      cat("\n")
    }

  }

}

#' getEPV
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#'
#' @export

getEPV <- function(X,Y){
  if("event" %in% colnames(Y)){
    EPV <- sum(Y[,"event"]) / ncol(X)
  }else{
    stop("Column event has not been detected in Y matrix.")
  }

  return(EPV)
}

deleteZeroVarianceVariables <- function(data, mustKeep = NULL, names = NULL, info=T, freqCut = 95/5, onlyZero = F){

  if(!is.null(names)){
    df <- data[[names]]
  }else{
    df <- data
  }

  #Zero Var
  nZ <- caret::nearZeroVar(x = df[,!colnames(df) %in% c("time", "event", "status")], saveMetrics = T, freqCut = freqCut) #to check if we have to do some changes in the data

  if(onlyZero){
    td <- rownames(nZ[nZ$zeroVar==T,])
  }else{
    td <- rownames(nZ[nZ$nzv==T,])
  }

  #Do not delete
  if(any(mustKeep %in% td)){
    td <- td[-which(td %in% mustKeep)]
  }

  lstDeleted <- td
  df <- df[,!colnames(df) %in% lstDeleted, drop=F]

  #df deleted
  df_cn_deleted <- NULL
  if(info){
    if(!length(lstDeleted)==0){
      for(cn in lstDeleted){
        df_cn_deleted <- rbind(df_cn_deleted, c(cn))#, getInfo(cn)$Description))
      }
      if(!is.null(df_cn_deleted)){
        df_cn_deleted <- as.data.frame(df_cn_deleted)
        colnames(df_cn_deleted) <- c("Variables")#, "Description")
        rownames(df_cn_deleted) <- NULL
      }
    }
  }
  return(list(filteredData = df, variablesDeleted = df_cn_deleted))
}

#' deleteZeroOrNearZeroVariance
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform into binary variables.
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering (default: NULL).
#' @param freqCut Numeric. Cutoff for the ratio of the most common value to the second most common value (default: 95/5).
#'
#' @export

deleteZeroOrNearZeroVariance <- function(X, remove_near_zero_variance = F, remove_zero_variance = T, toKeep.zv = NULL, freqCut = 95/5){

  auxX <- X

  variablesDeleted <- NULL
  if(remove_near_zero_variance){
    lst.zv <- deleteZeroVarianceVariables(data = auxX, info = T, mustKeep = toKeep.zv, freqCut = freqCut)
    variablesDeleted <- lst.zv$variablesDeleted[,1]
    if(!is.null(variablesDeleted)){
      auxX <- auxX[,!colnames(auxX) %in% variablesDeleted,drop=F]
    }
  }else if(remove_zero_variance){
    lst.zv <- deleteZeroVarianceVariables(data = auxX, info = T, mustKeep = toKeep.zv, onlyZero = T)
    variablesDeleted <- lst.zv$variablesDeleted[,1]
    if(!is.null(variablesDeleted)){
      auxX <- auxX[,!colnames(auxX) %in% variablesDeleted,drop=F]
    }
  }

  return(list(X = auxX, variablesDeleted = variablesDeleted))

}

#### ### ### ### ##
# Other functions #
#### ### ### ### ##

deleteIllegalChars <- function(chr.vector){
  v <- chr.vector
  for(i in pkg.env$IllegalChars){
    v <- unlist(sapply(v, function(x, i){gsub(i, "", x, fixed = T)}, i = i))
  }
  return(v)
}

#only for FORMULAS
transformIllegalChars <- function(cn){
  #### Formula cannot manage -,+,* symbols in cn
  if(!length(cn)>1){
    if(length(grep("-", cn, fixed = T))>0){
      cn <- gsub("-", "_", x = cn, fixed = T)
    }
    if(length(grep("+", cn, fixed = T))>0){
      cn <- gsub("+", ".", x = cn, fixed = T)
    }
    if(length(grep("*", cn, fixed = T))>0){
      cn <- gsub("*", ".star.", x = cn, fixed = T)
    }
  }else{
    new_cn <- NULL
    for(c in cn){
      if(length(grep("-", c, fixed = T))>0){
        c <- gsub("-", "_", x = c, fixed = T)
      }
      if(length(grep("+", c, fixed = T))>0){
        c <- gsub("+", ".", x = c, fixed = T)
      }
      if(length(grep("*", c, fixed = T))>0){
        c <- gsub("*", ".star.", x = c, fixed = T)
      }
      new_cn <- c(new_cn, c)
    }
    cn <- new_cn
  }
  return(cn)
}

checkColnamesIllegalChars <- function(X){
  new_cn_X <- deleteIllegalChars(colnames(X))

  if(length(unique(new_cn_X)) == length(unique(colnames(X)))){
    colnames(X) <- new_cn_X
  }else{
    stop(paste0("When deleting illegal chars, some colnames in X get the same name. Update manually the colnames to avoid the next chars: ", paste0(pkg.env$IllegalChars, collapse = " ")))
  }

  return(X)

}

stop_quietly <- function(s = NULL) {
  if(!is.null(s)){
    message(s)
  }

  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

checkXY.class <- function(X, Y, verbose = F){
  # Check if X and Y are matrices
  if (!is.matrix(X)) {
    if(is.data.frame(X)){
      if(verbose){
        message("X data is not a matrix, applying data.matrix\n")
      }
      X <- data.matrix(X)
    }else{
      stop("X data is not a matrix or a data.frame")
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

  if(any(is.na(Y))){
    stop("Y data has one or more NA values. Remove those patients in Y and X matrices.")
  }

  return(list(X = X, Y = Y))
}

checkY.colnames <- function(Y){
  if(!all(colnames(Y) %in% c("event", "status", "time"))){
    stop_quietly("Y must contain 'event' or 'status' and 'time' columns.")
  }else{
    if(any(colnames(Y) %in% "status")){
      colnames(Y)[colnames(Y)=="status"] <- "event"
    }
  }
}

XY.scale <- function(X, Y, x.center, x.scale, y.center, y.scale){
  xmeans <- NULL
  xsds <- NULL
  ymeans <- NULL
  ysds <- NULL
  # Centering and/or scaling
  if(x.center | x.scale){
    Xh <- scale(X, center = x.center, scale = x.scale)
    if(x.center) xmeans <- attr(Xh, "scaled:center")
    if(x.scale) xsds <- attr(Xh, "scaled:scale")
  }else{
    Xh <- X
  }

  if(y.center | y.scale){
    Yh_time <- scale(Y[,"time", drop=F], center = y.center, scale = y.scale)
    if(y.center) ymeans <- attr(Y, "scaled:center")
    if(y.scale) ysds <- attr(Y, "scaled:scale")
    Yh <- Y
    Yh[,"time"] <- Yh_time
  }else{
    Yh <- Y
  }

  return(list(Xh = Xh, xmeans = xmeans, xsds = xsds, Yh = Yh, ymeans = ymeans, ysds = ysds))
}

check.cv.weights <- function(vector){
  if(sum(vector)!=1){
    stop_quietly("Total sum of weight must sum 1.")
  }
}

check.ncomp <- function(X, max.ncomp, verbose = F){

  if(length(max.ncomp)>1){
    stop("max.ncomp must be a number, not a vector")
  }

  if(max.ncomp > ncol(X)){
    if(verbose){
      message(paste0("Number of components cannot be grater than the number of variables. Updated to ", ncol(X), "."))
    }
    max.ncomp <- ncol(X)
  }
  return(max.ncomp)
}

check.maxPredictors <- function(X, Y, MIN_EPV = 5, max.variables, verbose = F){
  max_n_predictors <- getMaxNPredictors(n.var = ncol(X), Y, MIN_EPV)
  max.variables.new <- max.variables

  if(max_n_predictors==0){
    stop_quietly(paste0("Less than ", MIN_EPV, " events. Stop program."))
  }

  if(max.variables > max_n_predictors){
    if(verbose){
      message(paste0("As we have a low number of events, we should work with a maximum of ", max_n_predictors, " variables/components."))
    }
    max.variables.new = max_n_predictors
  }

  return(max.variables.new)
}

check.maxPredictors.cox <- function(X, Y, MIN_EPV = 5, FORCE){
  max_n_predictors <- getMaxNPredictors(n.var = ncol(X), Y, MIN_EPV)

  if(max_n_predictors==0){
    msg <- paste0("Less than ", MIN_EPV, " events per variable. Stop program.")
    message(msg)
    return(F)
  }

  if(ncol(X) > max_n_predictors){
    if(!FORCE){
      EPV <- getEPV(X, Y)
      msg <- paste0("The ratio of Y events to the number of variables used does not meet the minimum EPV. The X matrix should have a maximum of ", max_n_predictors, " variables, or the MIN_EPV parameter should be lowered to ", EPV,".")
      message(msg)
      return(F)
    }
    return(T)
  }
  return(T)
}

#### ## ### ## ### ##
# INITIAL FUNCTIONS #
#### ## ### ## ### ##

getInfoCoxModel <- function(cox){
  survival_model <- NULL
  survival_model$fit <- cox

  survival_model$AIC <- extractAIC(survival_model$fit)[2]
  survival_model$BIC <- extractAIC(survival_model$fit, k = log(survival_model$fit$n))[2]

  survival_model$lp <- survival_model$fit$linear.predictors
  survival_model$coef <- coef(survival_model$fit)

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #       Prediction of the components         #
  #     as if missing values (model free)      #
  #       For cross-validating the GLM         #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  survival_model$YChapeau <- predict(survival_model$fit, type='expected')
  survival_model$Yresidus <- residuals(survival_model$fit, type="martingale")
  return(survival_model)
}

getMaxNPredictors <- function(n.var, Y, MIN_EPV){
  n_events <- NULL
  if(is.numeric(Y[,"event"])){
    n_events <- sum(Y[,"event"]==1)
  }else if(is.logical(Y[,"event"])){
    n_events <- sum(Y[,"event"]==T)
  }

  EPV <-  floor(n_events / 1:n.var) #EPV

  max_n_predictors <- n.var
  if(any(EPV > MIN_EPV)){
    max_n_predictors <- max(which(EPV>MIN_EPV))
  }else{
    cat("Minimum EPV not reached. You should be less strict or increse the number of events.\n")
    return(0) #Treatment in the function that calls this one
  }

  if(MIN_EPV==0){
    max_n_predictors <-n.var
  }

  return(max_n_predictors)
}

### EVALUATE MODELS

getPvalFromCox <- function(cox){
  p_val <- summary(cox)[[7]][,"Pr(>|z|)"]
  return(p_val)
}

removeNonSignificativeCox <- function(cox, alpha, cox_input, time.value = NULL, event.value = NULL){

  d <- cox_input

  if(!is.null(time.value) & !is.null(event.value)){
    time <- time.value
    event <- event.value
    d <- cbind(d, time)
    d <- cbind(d, event)
  }

  p_val <- getPvalFromCox(cox)
  removed_variables <- NULL

  while(any(p_val>alpha) && length(p_val)>1){
    to_remove <- names(which.max(p_val))
    to_remove <- deleteIllegalChars(to_remove)
    d <- d[,!colnames(d) %in% c(to_remove),drop=F]
    cox <- tryCatch(
      # Specifying expression
      expr = {
        survival::coxph(formula = survival::Surv(time,event) ~ .,
                        data = d,
                        ties = "efron",
                        singular.ok = T,
                        robust = T,
                        nocenter = rep(1, ncol(d)-2), #2 == ncol(Yh)
                        model=T, x = T)
      },
      # Specifying error message
      error = function(e){
        message(paste0("Updating cox model: ", e))
        invisible(gc())
        return(NA)
      }
    )

    removed_variables <- c(removed_variables, to_remove)
    p_val <- getPvalFromCox(cox)
  }

  return(list(cox = cox, removed_variables = removed_variables))

}

getFAST_LP_AUC <- function(fast_mode, comp_index, eta_index = NULL, run, fold, lst_X_test, lst_Y_test, comp_model_lst, times = NULL, lst_linear.predictors, df_results_evals_AUC, pred.method, pred.attr, PARALLEL = F, verbose = F){
  lst_linear.predictors <- NULL
  Y_test_full <- NULL
  lst_resCOMPLETE_LP <- getCOMPLETE_LP(comp_index = comp_index, eta_index = eta_index, run = run, fold = fold,
                                       lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, Y_test_full = Y_test_full,
                                       comp_model_lst = comp_model_lst, lst_linear.predictors = lst_linear.predictors)
  Y_test <- lst_resCOMPLETE_LP$Y_test
  Y_test_full <- lst_resCOMPLETE_LP$Y_test_full
  lst_linear.predictors <- lst_resCOMPLETE_LP$lst_linear.predictors

  Y_test_to_use <- NULL
  if(fast_mode){
    Y_test_to_use <- Y_test
  }else{
    Y_test_to_use <- Y_test_full
  }

  lst_resCOMPLETE_LP_AUC <- getCOMPLETE_LP_AUC(Y_test_full = Y_test_to_use,
                                               lst_linear.predictors = lst_linear.predictors, times = times,
                                               df_results_evals_AUC = df_results_evals_AUC,
                                               pred.method, pred.attr, PARALLEL = PARALLEL, verbose = verbose)

  lst_AUC_values <- lst_resCOMPLETE_LP_AUC$lst_AUC_values
  df_results_evals_AUC <- lst_resCOMPLETE_LP_AUC$df_results_evals_AUC

  return(list(lst_AUC_values = lst_AUC_values, df_results_evals_AUC = df_results_evals_AUC))
}

getCOMPLETE_LP <- function(comp_index, eta_index = NULL, run, fold, lst_X_test, lst_Y_test, Y_test_full, comp_model_lst, lst_linear.predictors){

  X_test <- lst_X_test[[run]][[fold]]
  Y_test <- lst_Y_test[[run]][[fold]]
  Y_test_full <- rbind(Y_test_full, Y_test)

  if(!is.null(eta_index)){
    model <- comp_model_lst[[comp_index]][[eta_index]][[run]][[fold]]
  }else{
    model <- comp_model_lst[[comp_index]][[run]][[fold]]
  }

  cox <- model$survival_model$fit

  #should be just null, na cannot happen anymore
  if(all(is.null(cox)) | all(is.na(cox))){
    lp <- rep(NA, nrow(X_test))
    names(lp) <- rownames(X_test)
    lst_linear.predictors[["fit"]] <- c(lst_linear.predictors[["fit"]], lp)
    lst_linear.predictors[["se.fit"]] <- c(lst_linear.predictors[["se.fit"]], lp)
    return(list(lst_linear.predictors = lst_linear.predictors, X_test = X_test, Y_test = Y_test, Y_test_full = Y_test_full))
  }else{
    ## Get LP for each fold
    X_test_mod <- predict.HDcox(object = model, newdata = X_test)

    lp <- getLinealPredictors(cox = cox, data = X_test_mod)
    lst_linear.predictors[["fit"]] <- c(lst_linear.predictors[["fit"]], lp$fit) #includes colnames
    lst_linear.predictors[["se.fit"]] <- c(lst_linear.predictors[["se.fit"]], lp$se.fit) #includes colnames

    return(list(lst_linear.predictors = lst_linear.predictors, X_test = X_test, Y_test = Y_test, Y_test_full = Y_test_full))
  }
}

getCOMPLETE_BRIER <- function(comp_index, eta_index = NULL, run, fold, lst_X_test, lst_Y_test, comp_model_lst, times, verbose = verbose){

  X_test <- lst_X_test[[run]][[fold]]
  Y_test <- lst_Y_test[[run]][[fold]]

  if(!is.null(eta_index)){
    model <- comp_model_lst[[comp_index]][[eta_index]][[run]][[fold]]
  }else{
    model <- comp_model_lst[[comp_index]][[run]][[fold]]
  }

  cox <- model$survival_model$fit
  cox$naive.var <- NULL

  #should be just null, na cannot happen anymore
  if(all(is.null(cox)) | all(is.na(cox))){
    brier <- NA
    return(list(brier_score = brier, X_test = X_test, Y_test = Y_test))
  }else{
    ## Get LP for each fold
    X_test_mod <- predict.HDcox(object = model, newdata = X_test)

    ### ##
    # SurvMetrics::Brier
    ### ##

    # PROBLEMS
    # The timepoint must be positive for time 0 and Brier Score metric.
    # t != 0
    # Requires that coxph returns ALWAYS the x matrix (more data to process)

    #brier_survMetrics <- SurvMetrics_BRIER(model = model, X_test_mod = X_test_mod, Y_test = Y_test, times = times)

    ### ##
    # survAUC::predErr
    ### ##

    # need two Surv objects and both linear predictors
    # pretty similar results compared to survcomp but this one can select the times
    # but integrated error changes significantly
    # sometimes generate ibrier greater than 1!!!

    #brier <- survAUC_BRIER(model, X_test, Y_test, times = brier_survcomp$times, raw_test = T)

    ### ##
    # survcomp::sbrier.score2proba
    ### ##

    # Usado por: Deep Learning–Based Multi-Omics Integration Robustly Predicts Survival in Liver Cancer

    # Needed - time/event/score for train and test dataframes
    # score is risk score (I do not know if survival probabilities or lp) - with lp is similar to pec
    # results changes depending which score are you selecting

    brier <- SURVCOMP_BRIER(model = model, X_test_mod = X_test_mod, Y_test = Y_test)

    ### ##
    # pec
    ### ##

    # PROBLEMS with some dependencies or functions:
    # Error in Hist(time, event) : could not find function "Hist" - come from "prodlim"
    # requires cox model, surv functions, and the data"
    # Requires that coxph returns ALWAYS the x matrix (more data to process)

    #brier_pec <- PEC_BRIER(model, X_test_mod, Y_test)

    return(list(brier_score = brier, X_test = X_test, Y_test = Y_test))
  }
}

# If more than 2 events, get the maximum and the minimum time for event patients and compute X time points between them (equally distributed)
# else, do it by the censored patients
getTimesVector <- function(Y, max_time_points = 15, ACCURACY = 0.001){

  if(length(Y[Y[,"event"]==1,"time"])>1){
    if(is.integer(Y[,"time"])){
      inter <- max(Y[Y[,"event"]==1,"time"]) - min(Y[Y[,"event"]==1,"time"])
      times <- seq(min(Y[Y[,"event"]==1,"time"]), max(Y[Y[,"event"]==1,"time"]), inter / (max_time_points-1))
      times <- round2any(times, accuracy = 1, f = ceiling)
    }else{
      inter <- max(Y[Y[,"event"]==1,"time"]) - min(Y[Y[,"event"]==1,"time"])
      times <- seq(min(Y[Y[,"event"]==1,"time"]), max(Y[Y[,"event"]==1,"time"]), inter / (max_time_points-1))
      times <- round2any(times, accuracy = ACCURACY, f = ceiling)
    }
  }else{
    if(is.integer(Y[,"time"])){
      inter <- max(Y[Y[,"event"]==0,"time"]) - min(Y[Y[,"event"]==0,"time"])
      times <- seq(min(Y[Y[,"event"]==0,"time"]), max(Y[Y[,"event"]==0,"time"]), inter / (max_time_points-1))
      times <- round2any(times, accuracy = 1, f = ceiling)
    }else{
      inter <- max(Y[Y[,"event"]==0,"time"]) - min(Y[Y[,"event"]==0,"time"])
      times <- seq(min(Y[Y[,"event"]==1,"time"]), max(Y[Y[,"event"]==1,"time"]), inter / (max_time_points-1))
      times <- round2any(times, accuracy = ACCURACY, f = ceiling)
    }
  }

  return(times)
}

getCOMPLETE_LP_AUC <- function(Y_test_full, lst_linear.predictors, df_results_evals_AUC, times = NULL, pred.method, pred.attr, PARALLEL = F, max_time_points = 15, verbose = F){
  #times

  ACCURACY <- 0.001 #!!! think to a best mode to select the best accuracy for each possible time
  if(is.null(times)){
    times <- getTimesVector(Y_test_full, max_time_points = max_time_points, ACCURACY = ACCURACY)
  }else{
    if(length(times)>max_time_points){
      if(verbose){
        message(paste0("More than ", max_time_points, " time points have been selected. CV processes could take a long time... Below ", max_time_points, " time points is recommended.\n\n"))
      }
    }
  }

  t1 <- Sys.time()
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lst_linear.predictors$fit,
                                       Y = Y_test_full, times = times, bestModel = NULL,
                                       method = pred.method, eval = pred.attr, PARALLEL = PARALLEL, verbose = verbose)
  t2 <- Sys.time()
  df_results_evals_AUC <- c(df_results_evals_AUC, lst_AUC_values$AUC)

  lst_AUC_values$computed_time <- difftime(t2,t1,"secs")

  return(list(lst_AUC_values = lst_AUC_values, df_results_evals_AUC = df_results_evals_AUC))
}

#' predict.HDcox
#'
#' @param object HDcox model
#' @param ... additional arguments affecting the predictions produced.
#' @param newdata Numeric matrix or data.frame. New data for explanatory variables (raw data). Qualitative variables must be transform into binary variables.
#'
#' @return Score values for new data using the HDcox model selected.
#' @export

predict.HDcox <- function(object, ..., newdata = NULL){

  model <- object

  if(!isa(model,pkg.env$model_class)){
    stop_quietly("Model must be an object of class HDcox.")
  }

  if(all(is.null(model$survival_model))){
    message("Survival Model not found.")
    return(NULL)
  }

  ## ## ## ## ## ##
  # SB sPLS-DRCOX #
  ## ## ## ## ## ##
  if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
    #in this case, we should update newdata before run the method cause the normalization is performed before
    x.mean <- NULL
    x.sd <- NULL
    if(!is.null(model$X$x.mean)){
      x.mean <- model$X$x.mean
    }

    if(!is.null(model$X$x.sd)){
      x.sd <- model$X$x.sd
    }

    if(!is.null(newdata)){
      if(!is.null(x.mean) | !is.null(x.sd)){
        X_test <- list()
        for(b in names(model$list_spls_models)){
          if(b %in% names(newdata)){

            if(!isa(x.mean, "list")){
              if(is.null(x.mean)){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }else{
              if(is.null(x.mean[[b]])){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }

            if(!isa(x.sd, "list")){
              if(is.null(x.sd)){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }else{
              if(is.null(x.sd[[b]])){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }

            X_test[[b]] <- scale(newdata[[b]][,names(model$X$x.mean[[b]]),drop=F], center = center_value, scale = scale_value)
          }
        }
      }else{
        X_test = data.matrix(newdata)
      }
    }else{
      X_test = NULL
    }

    predicted_scores <- NULL
    cn.merge = NULL
    for(b in names(model$list_spls_models)){
      predicted_scores <- cbind(predicted_scores, predict.HDcox(object = model$list_spls_models[[b]], newdata = X_test[[b]]))
      cn.merge <- c(cn.merge, paste0(colnames(model$list_spls_models[[b]]$X$scores), "_", b))
      #cn.merge <- c(cn.merge, paste0(colnames(model$list_spls_models[[b]]$X$scores)[1:model$list_spls_models[[b]]$n.comp], "_", b))
    }

    #colnames(predicted_scores) <- apply(expand.grid(colnames(model$list_spls_models[[1]]$X$scores[,,drop=F]), names(model$list_spls_models)), 1, paste, collapse="_")
    colnames(predicted_scores) <- cn.merge

    rn <- NULL
    if(!is.null(newdata)){
      rn <- rownames(newdata[[1]])
    }else{
      rn <- rownames(model$X$data[[1]])
    }
    rownames(predicted_scores) <- rn
    return(predicted_scores)
  }

  ## ## ## ## ##
  # SB sPLS-ICOX #
  ## ## ## ## ##
  if(attr(model, "model") %in% pkg.env$sb.splsicox){
    #in this case, we should update newdata before run the method cause the normalization is performed before
    x.mean <- NULL
    x.sd <- NULL
    if(!is.null(model$X$x.mean)){
      x.mean <- model$X$x.mean
    }

    if(!is.null(model$X$x.sd)){
      x.sd <- model$X$x.sd
    }

    if(!is.null(newdata)){
      if(!is.null(x.mean) | !is.null(x.sd)){
        X_test <- list()
        for(b in names(model$list_pls_models)){
          if(b %in% names(newdata)){

            if(!isa(x.mean, "list")){
              if(is.null(x.mean)){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }else{
              if(is.null(x.mean[[b]])){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }

            if(!isa(x.sd, "list")){
              if(is.null(x.sd)){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }else{
              if(is.null(x.sd[[b]])){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }

            X_test[[b]] <- scale(newdata[[b]][,names(model$X$x.mean[[b]]),drop=F], center = center_value, scale = scale_value)
          }
        }
      }else{
        X_test = data.matrix(newdata)
      }
    }else{
      X_test = NULL
    }

    predicted_scores <- NULL
    cn.merge = NULL
    for(b in names(model$list_pls_models)){
      predicted_scores <- cbind(predicted_scores, predict.HDcox(object = model$list_pls_models[[b]], newdata = X_test[[b]]))
      cn.merge <- c(cn.merge, paste0(colnames(model$list_pls_models[[b]]$X$scores)[1:model$list_pls_models[[b]]$n.comp], "_", b))
    }

    #colnames(predicted_scores) <- apply(expand.grid(colnames(model$list_pls_models[[1]]$X$scores[,,drop=F]), names(model$list_pls_models)), 1, paste, collapse="_")
    colnames(predicted_scores) <- cn.merge

    rn <- NULL
    if(!is.null(newdata)){
      rn <- rownames(newdata[[1]])
    }else{
      rn <- rownames(model$X$data[[1]])
    }
    rownames(predicted_scores) <- rn
    return(predicted_scores)
  }

  ## ## ## ## ## ## ##
  # GET MEAN and SD #
  ## ## ## ## ## ## ##
  x.mean <- NULL
  x.sd <- NULL
  if(!is.null(model$X$x.mean)){
    x.mean <- model$X$x.mean
  }

  if(!is.null(model$X$x.sd)){
    x.sd <- model$X$x.sd
  }

  ## ## ## ## #
  # NORM DATA #
  ## ## ## ## #
  if(is.null(newdata)){
    X_test <- model$X$data
  }else{
    #Update test data
    if(!is.null(x.mean) | !is.null(x.sd)){
      if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
        X_test <- list()
        for(b in names(model$X$loadings)){
          if(b %in% names(newdata)){

            if(!isa(x.mean, "list")){
              if(is.null(x.mean)){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }else{
              if(is.null(x.mean[[b]])){
                center_value = FALSE
              }else{
                center_value = x.mean[[b]]
              }
            }

            if(!isa(x.sd, "list")){
              if(is.null(x.sd)){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }else{
              if(is.null(x.sd[[b]])){
                scale_value = FALSE
              }else{
                scale_value = x.sd[[b]]
              }
            }

            X_test[[b]] <- scale(newdata[[b]][,names(model$X$x.mean[[b]]),drop=F], center = center_value, scale = scale_value)
          }
        }
      }else{
        if(is.null(x.mean)){
          center_value = FALSE
        }else{
          center_value = x.mean
        }

        if(is.null(x.sd)){
          scale_value = FALSE
        }else{
          scale_value = x.sd
        }
        X_test <- scale(newdata[,names(model$X$x.mean),drop=F], center = center_value, scale = scale_value)
      }
    }else{
      if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
        if(all(names(newdata) %in% names(model$X$loadings))){
          X_test <- newdata
        }
      }else{
        X_test <- data.matrix(newdata)
      }
    }
  }

  ### TEST DATA - selected variables

  ### PLS METHODS
  if(attr(model, "model") %in% c(pkg.env$splsicox, pkg.env$splsdrcox, pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic)){
    X_test <- X_test[,rownames(model$X$W.star),drop=F] #splsdacox_dynamic can filter nzv, reason that different variables for W.star
  }

  ### CLASSICAL METHODS
  if(attr(model, "model") %in% pkg.env$classical_methods){
    coeff <- names(model$survival_model$fit$coefficients)
    coeff <- deleteIllegalChars(coeff)

    if(!all(coeff %in% colnames(X_test))){
      stop("Not all coefficients are present in X matrix.")
    }

    return(X_test[,coeff,drop=F]) #just return newdata scaled
  }

  ### MB METHODS
  if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

    if(!all(names(X_test) %in% names(model$X$loadings))){
      stop(paste0("New data must be a list with at least one of the following names: ", paste0(names(model$n.varX), collapse = ", ")))
    }

    lst_coeff = list()
    for(b in names(model$X$loadings)){
      lst_coeff[[b]] <- rownames(model$X$W.star[[b]])

      if(b %in% names(X_test)){

        if(!all(lst_coeff[[b]] %in% colnames(X_test[[b]]))){
          stop("Not all coefficients are present in X matrix.")
        }

        X_test[[b]] = X_test[[b]][,lst_coeff[[b]],drop=F]
      }
    }
  }

  ### Estimate tts
  if(attr(model, "model") %in% pkg.env$splsdrcox_dynamic){

    predicted_scores <- X_test %*% model$X$W.star

    #MixOmics normalization... (still needed?)
    # predicted_scores <- matrix(data = sapply(1:ncol(predicted_scores),
    #                                          function(x) {predicted_scores[, x] * apply(model$X$scores, 2,
    #                                                                                     function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(X_test), ncol = ncol(predicted_scores))
    colnames(predicted_scores) <- colnames(model$X$scores)
    rownames(predicted_scores) <- rownames(X_test)

  }else if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

    lst_predicted_scores <- list()
    for(b in names(model$X$loadings)){
      if(b %in% names(X_test)){

        #Now, W.star has 0 for those variables have not be selected in each component
        #lst_predicted_scores[[b]] <- X_test[[b]] %*% model$X$W.star[[b]][lst_coeff[[b]],,drop=F]
        lst_predicted_scores[[b]] <- X_test[[b]] %*% model$X$W.star[[b]]

        #NORMALIZATION HAS TO BE PERFORM FOR MIXOMICS ALGORITHMS
        lst_predicted_scores[[b]] <- matrix(data = sapply(1:ncol(lst_predicted_scores[[b]]),
                                                 function(x) {lst_predicted_scores[[b]][, x] * apply(model$X$scores[[b]], 2,
                                                                                            function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(X_test[[b]]), ncol = ncol(lst_predicted_scores[[b]]))
        colnames(lst_predicted_scores[[b]]) <- colnames(model$X$scores[[b]])
        rownames(lst_predicted_scores[[b]]) <- rownames(X_test[[b]])
      }
    }

    predicted_scores = NULL
    for(i in 1:length(lst_predicted_scores)){
      predicted_scores <- cbind(predicted_scores, lst_predicted_scores[[i]])
    }

    colnames(predicted_scores) <- apply(expand.grid(colnames(model$X$scores[[1]]), names(model$X$loadings)), 1, paste, collapse="_")

  }else{
    # "sPLS-DACOX-Dynamic" does not perform the normalization for mixOmics, just this
    predicted_scores <- X_test %*% model$X$W.star
    colnames(predicted_scores) <- colnames(model$X$scores)
    rownames(predicted_scores) <- rownames(X_test)
  }

  return(predicted_scores)
}

check_AUC_improvement_spls1 <- function(fast_mode, pred.attr, df_results_evals_AUC, comp_index, eta_index, eta.list, n_run, k_folds, lst_comp_AUC, MIN_COMP_TO_CHECK, MIN_AUC, max.ncomp){
  #CHECK AUC EVOLUTION PER COMPONENT
  comp_AUC <- NULL

  if(fast_mode){
    ini <- ((eta_index-1)*n_run*k_folds+1)+((comp_index-1)*length(eta.list)*n_run*k_folds)
    end <- ((eta_index)*n_run*k_folds)+((comp_index-1)*length(eta.list)*n_run*k_folds)
    comp_AUC <- ifelse(pred.attr=="mean",
                       mean(df_results_evals_AUC[ini:end], na.rm = T),
                       median(df_results_evals_AUC[ini:end], na.rm = T))
  }else{
    ini <- ((eta_index-1)*n_run+1)+((comp_index-1)*length(eta.list)*n_run)
    end <- ((eta_index)*n_run)+((comp_index-1)*length(eta.list)*n_run)
    comp_AUC <- ifelse(pred.attr=="mean",
                       mean(df_results_evals_AUC[ini:end], na.rm = T),
                       median(df_results_evals_AUC[ini:end], na.rm = T))
  }

  ind <- as.character(eta.list[eta_index])
  lst_comp_AUC[[ind]] <- c(lst_comp_AUC[[ind]], comp_AUC)

  return(lst_comp_AUC)
}

check_AUC_improvement_spls2 <- function(fast_mode, pred.attr, df_results_evals_AUC, comp_index, eta_index, eta.list, n_run, k_folds, lst_comp_AUC, MIN_COMP_TO_CHECK, MIN_AUC, MIN_AUC_INCREASE, max.ncomp){
  optimal_comp_index <- NULL
  optimal_eta_index <- NULL
  optimal_eta <- NULL
  optimal_auc <- NULL
  optimal_comp_flag <- FALSE

  #For at least MIN_COMP_TO_CHECK
  if(comp_index > MIN_COMP_TO_CHECK){
    #If first comp_index is greater than the minimum AUC
    if(any(unlist(purrr::map(lst_comp_AUC, comp_index-MIN_COMP_TO_CHECK))>MIN_AUC)){
      diff_vector <- NULL
      for(e in 1:length(eta.list)){
        ind <- as.character(eta.list[e])
        #Check if the difference between newest components is lesser than the AUC
        for(c in (comp_index-MIN_COMP_TO_CHECK+1):comp_index){
          diff_vector[[ind]] <- c(diff_vector[[ind]], abs(lst_comp_AUC[[ind]][[c]] - lst_comp_AUC[[ind]][[comp_index-MIN_COMP_TO_CHECK]]))
        }
      }

      #all vectors has a lesser difference than the MIN_AUC
      if(any(unlist(purrr::map(diff_vector, function(x){all(x < MIN_AUC_INCREASE)})))){

        #if any, is the maximum value we can obtain?
        optimal_comp_index <- comp_index-MIN_COMP_TO_CHECK
        optimal_eta_index <- which.max(unlist(purrr::map(lst_comp_AUC, comp_index-MIN_COMP_TO_CHECK)))
        optimal_eta <- eta.list[which.max(unlist(purrr::map(lst_comp_AUC, comp_index-MIN_COMP_TO_CHECK)))]
        optimal_auc <- lst_comp_AUC[[optimal_eta_index]][optimal_comp_index]

        #if real maximum stop
        if(all(unlist(purrr::map(lst_comp_AUC, ~max(.))) - optimal_auc < MIN_AUC_INCREASE)){
          optimal_comp_flag <- TRUE
          txt.plscox <- paste0("\n Not improvement found. Stop in component: ", max.ncomp[[comp_index]])
          message(txt.plscox)
          return(list(optimal_comp_index = optimal_comp_index, optimal_eta_index = optimal_eta_index, optimal_eta = optimal_eta, optimal_auc = optimal_auc, optimal_comp_flag = optimal_comp_flag, lst_comp_AUC = lst_comp_AUC))
        }
      }
    }
  }#comp

  return(list(optimal_comp_index = optimal_comp_index, optimal_eta_index = optimal_eta_index, optimal_eta = optimal_eta, optimal_auc = optimal_auc, optimal_comp_flag = optimal_comp_flag, lst_comp_AUC = lst_comp_AUC))

}

check_AUC_improvement <- function(fast_mode, pred.attr, df_results_evals_AUC, comp_index, n_run, k_folds, lst_comp_AUC, MIN_COMP_TO_CHECK, MIN_AUC, MIN_AUC_INCREASE, max.ncomp, method.train = "pls"){
  #CHECK AUC EVOLUTION PER COMPONENT
  comp_AUC <- NULL
  if(fast_mode){
    ini <- (comp_index-1)*n_run*k_folds+1
    end <- comp_index*n_run*k_folds
    comp_AUC <- ifelse(pred.attr=="mean",
                       mean(df_results_evals_AUC[ini:end], na.rm = T),
                       median(df_results_evals_AUC[ini:end], na.rm = T))
  }else{
    ini <- (comp_index-1)*n_run+1
    end <- comp_index*n_run
    comp_AUC <- ifelse(pred.attr=="mean",
                       mean(df_results_evals_AUC[ini:end], na.rm = T),
                       median(df_results_evals_AUC[ini:end], na.rm = T))
  }

  lst_comp_AUC <- c(lst_comp_AUC, comp_AUC)
  optimal_comp_index <- NULL
  optimal_comp_flag <- FALSE

  #For at least MIN_COMP_TO_CHECK
  if(comp_index > MIN_COMP_TO_CHECK){
    if(!is.na(lst_comp_AUC[[comp_index-MIN_COMP_TO_CHECK]])){ #component must not to be NA
      #If first l is greater than the minimum AUC
      if(lst_comp_AUC[[comp_index-MIN_COMP_TO_CHECK]] > MIN_AUC){
        diff_vector <- NULL
        #Check if the difference between newest components is lesser than the AUC
        for(c in (comp_index-MIN_COMP_TO_CHECK+1):comp_index){
          diff_vector <- c(diff_vector, lst_comp_AUC[c] - lst_comp_AUC[[comp_index-MIN_COMP_TO_CHECK]])
        }
        if(all(diff_vector[!is.na(diff_vector)] < MIN_AUC_INCREASE)){ #checking not NA values (maybe is not the better option)
          #all components of study increses less than the minimum - stop processing new components
          optimal_comp_index <- comp_index-MIN_COMP_TO_CHECK
          optimal_comp_flag <- TRUE

          txt.coxEn <- paste0("\n Not improvement found. Stop in lambda: ", max.ncomp[[comp_index]])
          txt.plscox <- paste0("\n Not improvement found. Stop in component: ", max.ncomp[[comp_index]])

          message(ifelse(method.train==pkg.env$coxEN,txt.coxEn, txt.plscox))
          return(list(optimal_comp_index = optimal_comp_index, optimal_comp_flag = optimal_comp_flag, lst_comp_AUC = lst_comp_AUC))
        }
      }
    }
  }

  return(list(optimal_comp_index = optimal_comp_index, optimal_comp_flag = optimal_comp_flag, lst_comp_AUC = lst_comp_AUC))
}

getAUC_RUN_AND_COMP <- function(mode = "AUC", fast_mode, max.ncomp, n_run,
                                df_results_evals, optimal_comp_flag, optimal_comp_index,
                                MIN_COMP_TO_CHECK, lst_AUC_component, df_results_evals_run,
                                df_results_evals_comp, method.train = "other"){

  if(is.null(optimal_comp_flag)){
    optimal_comp_flag <- F
  }

  if(is.null(optimal_comp_index)){
    optimal_comp_index <- 0 #need a value
  }

  if(length(max.ncomp)==1){
    max.ncomp <- 1:max.ncomp
  }

  if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdrcox_dynamic)){ #num.var is factor, convert to numeric
    df_results_evals$n.var <- as.numeric(as.character(df_results_evals$n.var))
    #df_results_evals = df_results_evals[,!colnames(df_results_evals) %in% "n.var"]
  }else{
    df_results_evals = df_results_evals
  }

  df_results_evals_run <- NULL
  df_results_evals_comp <- NULL

  ### ###
  # AUC #
  ### ###
  if(!fast_mode){
    for(l in unique(df_results_evals$n.comps)){
      l.index <- which(l == unique(df_results_evals$n.comps))
      # EVAL PER RUN
      eval_aux.run <- NULL
      for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){
        aux.run <- df_results_evals[which(df_results_evals$n.comps==l & df_results_evals$runs==r),!colnames(df_results_evals) %in% c("fold")]

        ### ### ### ### ###
        # DYNAMIC METHODS #
        ### ### ### ### ###
        if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
          eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
          eval_aux.r <- as.data.frame(t(eval_aux.r))

          eval_aux.nvar.r <- aux.run[,colnames(aux.run) %in% c("n.var")]
          max_val <- max(table(eval_aux.nvar.r))
          names_max <- names(table(eval_aux.nvar.r))[table(eval_aux.nvar.r)==max_val]

          #same max value, takes best c-index (no AUC bc if method complete no-exits)
          #mean c-index before take it
          if(length(names_max)>1){
            sub_aux <- aux.run[,colnames(aux.run) %in% c("n.var", "c_index")]
            sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
            sub_aux$n.var <- factor(sub_aux$n.var)

            best_n_var = NULL
            best_n_var_c_index = 0
            for(lv in levels(sub_aux$n.var)){
              aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$c_index, rm.na=T)
              if(aux_c > best_n_var_c_index){
                best_n_var_c_index = aux_c
                best_n_var = lv
              }
            }

            names_max <- best_n_var
          }

          eval_aux.nvar.r <- names_max
          eval_aux.r$n.var <- as.numeric(eval_aux.nvar.r)
        }else{
          eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
        }

        if(optimal_comp_flag & l.index > (optimal_comp_index+MIN_COMP_TO_CHECK)){
          eval_aux.r[[mode]] <- NA
        }else{
          eval_aux.r[[mode]] <- lst_AUC_component[[l.index]][[r]][[mode]]
        }
        eval_aux.run <- rbind(eval_aux.run, eval_aux.r)
      }

      df_results_evals_run <- rbind(df_results_evals_run, eval_aux.run)

      # EVAL PER COMPONENT
      aux.l <- df_results_evals[which(df_results_evals$n.comps==l),!colnames(df_results_evals) %in% c("fold", "runs")]

      if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
        eval_aux <- apply(aux.l, 2, function(x){mean(x, na.rm = T)})
        eval_aux <- as.data.frame(t(eval_aux))
        eval_aux.nvar <- aux.l[,colnames(aux.l) %in% c("n.var")]
        max_val <- max(table(eval_aux.nvar))
        names_max <- names(table(eval_aux.nvar))[table(eval_aux.nvar)==max_val]

        #same max value, takes best c-index (no AUC bc if method complete no-exits)
        if(length(names_max)>1){
          sub_aux <- aux.run[,colnames(aux.run) %in% c("n.var", "c_index")]
          sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
          sub_aux$n.var <- factor(sub_aux$n.var)

          best_n_var = NULL
          best_n_var_c_index = 0
          for(lv in levels(sub_aux$n.var)){
            aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$c_index, rm.na=T)
            if(aux_c > best_n_var_c_index){
              best_n_var_c_index = aux_c
              best_n_var = lv
            }
          }

          names_max <- best_n_var
        }

        eval_aux.nvar <- names_max
        eval_aux$n.var <- as.numeric(eval_aux.nvar)
      }else{
        eval_aux <- apply(aux.l, 2, function(x){mean(x, na.rm = T)})
        m.freq_var <- as.numeric(names(table(aux.l$n.var)[table(aux.l$n.var) == max(table(aux.l$n.var))]))
        eval_aux[["n.var"]] <- min(m.freq_var) #in case of same quantity, get lower variables
      }

      #show NA values for components where AUC has not been computed
      if(mode == "AUC"){
        AUC_mean <- NULL
        if(optimal_comp_flag & l.index > (optimal_comp_index+MIN_COMP_TO_CHECK)){
          AUC_mean <- NA #IF GREATER THAN OPTIMAL, NA
        }else{
          AUC_v <- NULL
          for(r in 1:n_run){
            AUC_v <- c(AUC_v, c(lst_AUC_component[[l.index]][[r]][[mode]])) #MEAN FOR ALL COMPONENTS
          }
          AUC_mean <- mean(AUC_v, na.rm = T)
        }
        eval_aux[[mode]] <- AUC_mean
      }

      df_results_evals_comp <- rbind(df_results_evals_comp, eval_aux)
    }

  ### ### ### ###
  # AUC + BRIER #
  ### ### ### ###
  }else{
    for(l in unique(df_results_evals$n.comps)){
      l.index <- which(l == unique(df_results_evals$n.comps))
      eval_aux.run <- NULL

      ### ### ### ### ###
      # DYNAMIC METHODS #
      ### ### ### ### ###
      if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

        eval_aux.r <- NULL
        for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){
          aux.run <- df_results_evals[which(df_results_evals$n.comps==l & df_results_evals$runs==r),!colnames(df_results_evals) %in% c("fold")]
          eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})

          ## SELECT A SPECIFIC NUMBER OF VARIABLES (NOT MEAN)
          eval_aux.nvar.r <- aux.run[,colnames(aux.run) %in% c("n.var"),drop=F]
          max_val <- max(table(eval_aux.nvar.r))
          names_max <- names(table(eval_aux.nvar.r))[table(eval_aux.nvar.r)==max_val]

          # same max value, takes best BRIER
          # as AUC could be not compute
          if(length(names_max)>1){
            sub_aux <- aux.run[,colnames(aux.run) %in% c("n.var", "BRIER")]
            sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
            sub_aux$n.var <- factor(sub_aux$n.var)

            best_n_var = NULL
            # best_n_var_c_index = 0 #c-index is higher better
            best_n_var_c_index = 1
            for(lv in levels(sub_aux$n.var)){
              aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$BRIER, rm.na=T)
              # if(aux_c > best_n_var_c_index){ #c-index is higher better
              if(aux_c < best_n_var_c_index){
                best_n_var_c_index = aux_c
                best_n_var = lv
              }
            }
            names_max <- best_n_var
          }

          eval_aux.r[["n.var"]] <- as.numeric(names_max)
          df_results_evals_run <- rbind(df_results_evals_run, eval_aux.r)
        }

      ### ### ### ### #
      # OTHER METHODS #
      ### ### ### ### #
      }else{
        eval_aux.r <- NULL
        for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){
          aux.run <- df_results_evals[which(df_results_evals$n.comps==l & df_results_evals$runs==r),!colnames(df_results_evals) %in% c("fold")]
          eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
          df_results_evals_run <- rbind(df_results_evals_run, eval_aux.r)
        }
      }

      # EVAL PER COMPONENT
      aux.l <- df_results_evals[which(df_results_evals$n.comps==l),!colnames(df_results_evals) %in% c("fold", "runs")]

      ### ### ### ### ###
      # DYNAMIC METHODS #
      ### ### ### ### ###
      if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
        aux <- df_results_evals[which(df_results_evals$n.comps==l),!colnames(df_results_evals) %in% c("fold", "runs")]
        eval_aux <- apply(aux, 2, function(x){mean(x, na.rm = T)})
        ## SELECT A SPECIFIC NUMBER OF VARIABLES (NOT MEAN)
        eval_aux.nvar <- aux[,colnames(aux) %in% c("n.var"),drop=F]
        max_val <- max(table(eval_aux.nvar))
        names_max <- names(table(eval_aux.nvar))[table(eval_aux.nvar)==max_val]

        # same max value, takes best BRIER
        # as AUC could be not compute
        if(length(names_max)>1){
          sub_aux <- aux[,colnames(aux) %in% c("n.var", "BRIER")]
          sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
          sub_aux$n.var <- factor(sub_aux$n.var)

          best_n_var = NULL
          # best_n_var_c_index = 0
          best_n_var_c_index = 1
          for(lv in levels(sub_aux$n.var)){
            aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$BRIER, rm.na=T)
            # if(aux_c > best_n_var_c_index){
            if(aux_c < best_n_var_c_index){
              best_n_var_c_index = aux_c
              best_n_var = lv
            }
          }

          names_max <- best_n_var
        }

        eval_aux[["n.var"]] <- as.numeric(names_max)

      ### ### ### ### #
      # OTHER METHODS #
      ### ### ### ### #
      }else{
        aux <- df_results_evals[which(df_results_evals$n.comps==l),!colnames(df_results_evals) %in% c("fold", "runs")]
        eval_aux <- apply(aux, 2, function(x){mean(x, na.rm = T)})
      }

      # show NA values for components where AUC has not been computed
      if(mode == "AUC"){
        AUC_mean <- NULL
        if(optimal_comp_flag & l.index > (optimal_comp_index+MIN_COMP_TO_CHECK)){
          AUC_mean <- NA #IF GREATER THAN OPTIMAL, NA
        }else{
          AUC_v <- NULL
          auc_per_run <- df_results_evals_run[which(df_results_evals_run[,"n.comps"] == l),][,mode]
          AUC_mean <- mean(auc_per_run, na.rm = T)
        }
        eval_aux[[mode]] <- AUC_mean
      }

      df_results_evals_comp <- rbind(df_results_evals_comp, eval_aux)

    }
  }

  ## clean NaN values for those folds/runs/comps where AUC has not been compute
  if(any(is.nan(df_results_evals_run[,mode]))){
    df_results_evals_run[which(is.nan(df_results_evals_run[,mode])),mode] <- NA
  }
  if(any(is.nan(df_results_evals_comp[,mode]))){
    df_results_evals_comp[which(is.nan(df_results_evals_comp[,mode])),mode] <- NA
  }

  if(method.train==pkg.env$splsdrcox_dynamic){ #num.var is text
    rownames(df_results_evals_run) <- NULL
    colnames(df_results_evals_run) <- c("n.comps", "runs", "n.var", "AIC", "c_index", mode)
    df_results_evals_run <- as.data.frame(df_results_evals_run)

    rownames(df_results_evals_comp) <- NULL
    colnames(df_results_evals_comp) <- c("n.comps", "n.var", "AIC", "c_index", mode)
    df_results_evals_comp <- as.data.frame(df_results_evals_comp)
  }else{
    rownames(df_results_evals_run) <- NULL
    colnames(df_results_evals_run) <- c("n.comps", "runs", "n.var", "AIC", "c_index", mode)
    df_results_evals_run <- as.data.frame(df_results_evals_run)

    rownames(df_results_evals_comp) <- NULL
    colnames(df_results_evals_comp) <- c("n.comps", "n.var", "AIC", "c_index", mode)
    df_results_evals_comp <- as.data.frame(df_results_evals_comp)
  }

  if(method.train %in% c(pkg.env$splsdacox_dynamic,pkg.env$sb.splsicox,pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
    df_results_evals_run$n.var <- factor(df_results_evals_run$n.var)
    df_results_evals_comp$n.var <- factor(df_results_evals_comp$n.var)
  }

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run))
}

getAUC_RUN_AND_COMP_sPLS <- function(mode = "AUC", fast_mode, max.ncomp, eta.list, n_run, df_results_evals, optimal_comp_flag, optimal_comp_index,
                                     MIN_COMP_TO_CHECK, lst_AUC_component, df_results_evals_run, df_results_evals_comp, method.train){

  if(is.null(optimal_comp_flag)){
    optimal_comp_flag <- F
  }

  if(is.null(optimal_comp_index)){
    optimal_comp_index <- 0 #need a value
  }

  if(length(max.ncomp)==1){
    max.ncomp <- 1:max.ncomp
  }

  if(method.train==pkg.env$splsdrcox_dynamic){ #num.var is text
    df_results_evals = df_results_evals[,!colnames(df_results_evals) %in% "n.var"]
  }else{
    df_results_evals = df_results_evals
  }

  if(!fast_mode){ #AUC
    for(l in 1:length(max.ncomp)){
      for(e in 1:length(eta.list)){
        # EVAL PER RUN
        eval_aux.run <- NULL
        for(r in 1:n_run){
          aux.run <- df_results_evals[which(df_results_evals$n.comps==max.ncomp[[l]] & df_results_evals$eta==eta.list[[e]] & df_results_evals$runs==r),!colnames(df_results_evals) %in% c("fold")]

          #could happen cause some times the models compute lesser number of components
          if(nrow(aux.run)==0){
            eval_aux.r <- rep(NA,ncol(df_results_evals))
            eval_aux.run <- rbind(eval_aux.run, eval_aux.r)
            next
          }

          if(method.train %in% pkg.env$sb.splsdrcox){
            eval_aux.r <- apply(aux.run[,!colnames(aux.run) %in% c("n.var")], 2, function(x){mean(x, na.rm = T)})
            eval_aux.r <- as.data.frame(t(eval_aux.r))

            eval_aux.nvar.r <- aux.run[,colnames(aux.run) %in% c("n.var")]
            max_val <- max(table(eval_aux.nvar.r))
            names_max <- names(table(eval_aux.nvar.r))[table(eval_aux.nvar.r)==max_val]

            #same max value, takes best c-index (no AUC bc if method complete no-exits)
            #mean c-index before take it
            if(length(names_max)>1){
              sub_aux <- aux.run[,colnames(aux.run) %in% c("n.var", "c_index")]
              sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
              sub_aux$n.var <- factor(sub_aux$n.var)

              best_n_var = NULL
              best_n_var_c_index = 0
              for(lv in levels(sub_aux$n.var)){
                aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$c_index, rm.na=T)
                if(aux_c > best_n_var_c_index){
                  best_n_var_c_index = aux_c
                  best_n_var = lv
                }
              }

              names_max <- best_n_var
            }

            eval_aux.nvar.r <- names_max
            eval_aux.r$n.var <- eval_aux.nvar.r
            eval_aux.r <- eval_aux.r[,c(1,2,3,6,4,5)]
          }else{
            eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
          }

          if(optimal_comp_flag & l > (optimal_comp_index+MIN_COMP_TO_CHECK)){
            eval_aux.r[[mode]] <- NA
          }else{
            eval_aux.r[[mode]] <- lst_AUC_component[[l]][[e]][[r]][[mode]]
          }
          eval_aux.run <- rbind(eval_aux.run, eval_aux.r)
        }

        if(all(is.na(eval_aux.run))){
          next
        }

        df_results_evals_run <- rbind(df_results_evals_run, eval_aux.run)

        # EVAL PER COMPONENT
        aux.run <- df_results_evals[which(df_results_evals$n.comps==max.ncomp[[l]] & df_results_evals$eta==eta.list[[e]]),!colnames(df_results_evals) %in% c("fold", "runs")]

        if(method.train %in% pkg.env$sb.splsdrcox){
          eval_aux.r <- apply(aux.run[,!colnames(aux.run) %in% c("n.var")], 2, function(x){mean(x, na.rm = T)})
          eval_aux.r <- as.data.frame(t(eval_aux.r))

          eval_aux.nvar.r <- aux.run[,colnames(aux.run) %in% c("n.var")]
          max_val <- max(table(eval_aux.nvar.r))
          names_max <- names(table(eval_aux.nvar.r))[table(eval_aux.nvar.r)==max_val]

          #same max value, takes best c-index (no AUC bc if method complete no-exits)
          #mean c-index before take it
          if(length(names_max)>1){
            sub_aux <- aux.run[,colnames(aux.run) %in% c("n.var", "c_index")]
            sub_aux <- sub_aux[sub_aux$n.var %in% names_max,]
            sub_aux$n.var <- factor(sub_aux$n.var)

            best_n_var = NULL
            best_n_var_c_index = 0
            for(lv in levels(sub_aux$n.var)){
              aux_c <- mean(sub_aux[sub_aux$n.var==lv,]$c_index, rm.na=T)
              if(aux_c > best_n_var_c_index){
                best_n_var_c_index = aux_c
                best_n_var = lv
              }
            }

            names_max <- best_n_var
          }

          eval_aux.nvar.r <- names_max
          eval_aux.r$n.var <- eval_aux.nvar.r
          eval_aux.r <- eval_aux.r[,c(1,2,5,3,4)]
        }else{
          eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
        }

        AUC_mean <- NULL

        if(optimal_comp_flag & l > (optimal_comp_index+MIN_COMP_TO_CHECK)){
          AUC_mean <- NA #IF GREATER THAN OPTIMAL, NA
        }else{
          AUC_v <- NULL
          for(r in 1:n_run){
            AUC_v <- c(AUC_v, c(lst_AUC_component[[l]][[e]][[r]][[mode]])) #MEAN FOR ALL COMPONENTS
          }
          AUC_mean <- mean(AUC_v, na.rm = T)
        }

        eval_aux.r[[mode]] <- AUC_mean
        df_results_evals_comp <- rbind(df_results_evals_comp, eval_aux.r)
      }
    }
  }else{ #AUC / BRIER
    for(l in 1:length(max.ncomp)){
      for(e in 1:length(eta.list)){
        if(method.train %in% c(pkg.env$splsicox, pkg.env$splsdrcox)){
          # EVAL PER COMPONENT
          aux <- df_results_evals[which(df_results_evals$n.comps==max.ncomp[[l]] & df_results_evals$eta==eta.list[[e]]),!colnames(df_results_evals) %in% c("fold", "runs")]
          eval_aux <- apply(aux, 2, function(x){mean(x, na.rm = T)})
          df_results_evals_comp <- rbind(df_results_evals_comp, eval_aux)

          # EVAL PER RUN
          eval_aux.r <- NULL
          for(r in 1:n_run){
            aux.run <- df_results_evals[which(df_results_evals$n.comps==max.ncomp[[l]] & df_results_evals$eta==eta.list[[e]] & df_results_evals$runs==r),!colnames(df_results_evals) %in% c("fold")]
            eval_aux.r <- apply(aux.run, 2, function(x){mean(x, na.rm = T)})
            df_results_evals_run <- rbind(df_results_evals_run, eval_aux.r)
          }
        }else if(method.train %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$sb.splsicox, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
          # we need to manage n.var as factor
        }else{
          stop(paste0(method.train, " is not a HDcox algorithm."))
        }
      }
    }
  }

  rownames(df_results_evals_run) <- NULL
  colnames(df_results_evals_run) <- c("n.comps", "eta", "runs", "n.var", "AIC", "c_index", mode)
  df_results_evals_run <- as.data.frame(df_results_evals_run)

  rownames(df_results_evals_comp) <- NULL
  colnames(df_results_evals_comp) <- c("n.comps", "eta", "n.var", "AIC", "c_index", mode)
  df_results_evals_comp <- as.data.frame(df_results_evals_comp)

  if(method.train %in% pkg.env$sb.splsdrcox){
    df_results_evals_run$n.var <- factor(df_results_evals_run$n.var)
    df_results_evals_comp$n.var <- factor(df_results_evals_comp$n.var)
  }

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run))
}

get_EVAL_PLOTS <- function(fast_mode, best_model_info, w_AUC, w_BRIER, max.ncomp, eta.list = NULL, df_results_evals_fold, df_results_evals_run, df_results_evals_comp, x.text = "Component", colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", colname_BRIER = "BRIER"){

  df_results_evals_comp_aux <- df_results_evals_comp

  if(length(max.ncomp)==1){
    max.ncomp <- 1:max.ncomp
  }

  sd_vector <- NULL

  #First AIC - C_INDEX - BRIER at FOLD LEVEL ALWAYS
  if(!is.null(eta.list)){
    for(l in 1:length(max.ncomp)){
      for(e in 1:length(eta.list)){
        vector <- c("AIC.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]] & df_results_evals_fold$eta==eta.list[[e]],colname_AIC]),
                    "c_index.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]] & df_results_evals_fold$eta==eta.list[[e]],colname_c_index]),
                    "BRIER.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]] & df_results_evals_fold$eta==eta.list[[e]],colname_BRIER]))

        sd_vector <- rbind(sd_vector, vector)

      }
    }
  }else{
    for(l in 1:length(max.ncomp)){
      vector <- c("AIC.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]],colname_AIC]),
                  "c_index.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]],colname_c_index]),
                  "BRIER.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]],colname_BRIER]))

      sd_vector <- rbind(sd_vector, vector)
    }
  }

  sd_vector_AUC <- NULL
  # AUC - fast_mode
  if(w_AUC!=0){
    if(!is.null(eta.list)){
      if(fast_mode){
        for(l in 1:length(max.ncomp)){
          for(e in 1:length(eta.list)){
            sd_vector_AUC <- rbind(sd_vector_AUC, "AUC.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]] & df_results_evals_fold$eta==eta.list[[e]],colname_AUC]))
          }
        }
      }else{
        for(l in 1:length(max.ncomp)){
          for(e in 1:length(eta.list)){
            sd_vector_AUC <- rbind(sd_vector_AUC, "AUC.sd" = sd(df_results_evals_run[df_results_evals_run$n.comps==max.ncomp[[l]] & df_results_evals_run$eta==eta.list[[e]],colname_AUC]))
          }
        }
      }
    }else{
      if(fast_mode){
        for(l in 1:length(max.ncomp)){
          sd_vector_AUC <- rbind(sd_vector_AUC, "AUC.sd" = sd(df_results_evals_fold[df_results_evals_fold$n.comps==max.ncomp[[l]],colname_AUC]))
        }
      }else{
        for(l in 1:length(max.ncomp)){
          sd_vector_AUC <- rbind(sd_vector_AUC, "AUC.sd" = sd(df_results_evals_run[df_results_evals_run$n.comps==max.ncomp[[l]],colname_AUC]))
        }
      }
    }
  }

  if(!is.null(sd_vector_AUC)){
    colnames(sd_vector_AUC) <- "AUC.sd"
    sd_vector <- cbind(sd_vector, sd_vector_AUC)
  }

  df_results_evals_comp_aux <- cbind(df_results_evals_comp, sd_vector)
  df_results_evals_comp_aux$n.comps <- factor(df_results_evals_comp_aux$n.comps, levels = unique(df_results_evals_comp_aux$n.comps))

  ggp_AUC <- NULL
  ggp_BRIER <- NULL
  ggp_c_index <- NULL
  ggp_AIC <- NULL
  if(!is.null(eta.list)){

    penalty_index <- which(colnames(df_results_evals_comp_aux) %in% c("eta", "spv_penalty"))
    penalty_name <- colnames(df_results_evals_comp_aux)[[penalty_index]]

    df_results_evals_comp_aux[[penalty_index]] <- factor(df_results_evals_comp_aux[[penalty_index]], levels = unique(df_results_evals_comp_aux[[penalty_index]]))
    if(w_AUC!=0){
      ggp_AUC <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_AUC, y.var.sd = "AUC.sd", x.color = penalty_name, best_component = best_model_info$n.comps, x.text = x.text, best_eta = best_model_info[[penalty_index]])
    }

    ggp_BRIER <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_BRIER, y.var.sd = "BRIER.sd", x.color = penalty_name, best_component = best_model_info$n.comps, x.text = x.text, best_eta = best_model_info[[penalty_index]])
    ggp_c_index <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_c_index, y.var.sd = "c_index.sd", x.color = penalty_name, best_component = best_model_info$n.comps, x.text = x.text, best_eta = best_model_info[[penalty_index]])
    ggp_AIC <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_AIC, y.var.sd = "AIC.sd", x.color = penalty_name, best_component = best_model_info$n.comps, x.text = x.text, best_eta = best_model_info[[penalty_index]])
  }else{
    if(w_AUC!=0){
      ggp_AUC <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_AUC, y.var.sd = "AUC.sd", best_component = best_model_info$n.comps, x.text = x.text)
    }

    ggp_BRIER <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_BRIER, y.var.sd = "BRIER.sd", best_component = best_model_info$n.comps, x.text = x.text)
    ggp_c_index <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_c_index, y.var.sd = "c_index.sd", best_component = best_model_info$n.comps, x.text = x.text)
    ggp_AIC <- evalplot_errorbar(df = df_results_evals_comp_aux, x.var = "n.comps", y.var = colname_AIC, y.var.sd = "AIC.sd", best_component = best_model_info$n.comps, x.text = x.text)
  }

  rownames(df_results_evals_comp_aux) <- NULL
  return(list(df_results_evals_comp = df_results_evals_comp_aux, ggp_AUC = ggp_AUC, ggp_BRIER = ggp_BRIER, ggp_c_index = ggp_c_index, ggp_AIC = ggp_AIC))
}

#### ### ### ### ### ##
# AUC Other functions #
#### ### ### ### ### ##
get_COX_evaluation_AIC_CINDEX <- function(comp_model_lst, max.ncomp, eta.list = NULL,
                                          n_run, k_folds, total_models,
                                          remove_non_significant_models, alpha = 0.05, verbose = F){

  if(length(max.ncomp)==1){
    max.ncomp <- 1:max.ncomp
  }

  df_results_evals <- NULL
  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso

  message(paste0("Evaluating COX models (AIC and C-Index)..."))
  pb$tick(0)
  if(is.null(eta.list)){
    for(comp in 1:length(max.ncomp)){
      for(r in 1:n_run){
        for(f in 1:k_folds){
          model <- comp_model_lst[[comp]][[r]][[f]]

          if(all(is.na(model)) || all(is.null(model$survival_model))){
            pb$tick()
            next
          }

          cox <- model$survival_model$fit

          # exception - EN no model
          if(all(is.null(cox)) || all(is.na(cox))){
            pb$tick()
            df_results_evals <- rbind(df_results_evals, cbind(max.ncomp[[comp]], r, f, 0, NA, NA))
            next
          }

          if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
            n_var <- unlist(purrr::map(model$n.varX, ~length(.[[1]])))
            n_var <- paste0(n_var, collapse = "_")
          }else if(attr(model, "model") == pkg.env$sb.splsicox){
            n_var <- purrr::map(model$list_pls_models, ~nrow(.$X$loadings))
            n_var <- paste0(n_var, collapse = "_") #VAR FOR SB.spls IS THE MAX NUMBER OF VARIABLES (PER BLOCK)
          }else if(attr(model, "model") == pkg.env$sb.splsdrcox){
            n_var <- purrr::map(model$list_spls_models, ~sum(rowSums(.$X$weightings!=0)>0)) #this have to be checked !!!
            n_var <- paste0(n_var, collapse = "_")
          }else if(attr(model, "model") %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic)){
            n_var <- unique(apply(model$X$weightings, 2, function(x){sum(x!=0)}))
            #n_var <- paste0(n_var, collapse = "_")
          }else if(attr(model, "model") %in% pkg.env$splsicox){
            n_var <- nrow(model$X$loadings)
          }else if(attr(model, "model") %in% pkg.env$classical_methods){
            n_var <- length(model$survival_model$coef) #COX, COXSW and COXEN
          }

          df <- as.data.frame(summary(cox)[[7]])

          #delete models with non-significant components
          if(remove_non_significant_models == T){
            if(any(df$`Pr(>|z|)`>alpha)){
              pb$tick()
              if(verbose){
                message(paste0("\nModel - Comp: ", comp ,", Run: ",r,", Fold: ",f," - Is a non-significant model"))
              }
              next
            }
          }

          aic <- stats::extractAIC(cox, k=2)[2] #k=2 <- AIC, [2] AIC Value
          c_index <- survival::concordance(cox)$concordance

          df_results_evals <- rbind(df_results_evals, cbind(max.ncomp[[comp]], r, f, n_var, aic, c_index))
          pb$tick()

        }#fold
      }#run
    }#component

    if(!is.null(df_results_evals)){
      colnames(df_results_evals) <- c("n.comps", "runs", "fold", "n.var", "AIC", "c_index")
      df_results_evals <- as.data.frame(df_results_evals)
    }

  }else{

    for(comp in 1:length(max.ncomp)){
      for(e in 1:length(eta.list)){
        for(r in 1:n_run){
          for(f in 1:k_folds){

            model <- comp_model_lst[[comp]][[e]][[r]][[f]]

            if(all(is.na(model)) || all(is.null(model$survival_model))){
              pb$tick()
              next
            }

            cox <- model$survival_model$fit

            #exception
            if(all(is.null(cox)) || all(is.na(cox))){
              pb$tick()
              next
            }

            eta <- model$eta
            if(attr(model, "model") == pkg.env$sb.splsdrcox){
              n_var <- purrr::map(model$list_spls_models, ~length(unique(unlist(.$var_by_component))))
              n_var <- paste0(n_var, collapse = "_") #VAR FOR SB.spls IS THE MAX NUMBER OF VARIABLES (PER BLOCK)
            }else{
              n_var <- nrow(model$X$loadings)
            }

            df <- as.data.frame(summary(cox)[[7]])

            #delete models with non-significant components
            if(remove_non_significant_models == T){
              if(any(df$`Pr(>|z|)`>alpha)){
                pb$tick()
                next
              }
            }

            aic <- stats::extractAIC(cox, k=2)[2] #k=2 <- AIC, [2] AIC Value
            c_index <- survival::concordance(cox)$concordance

            df_results_evals <- rbind(df_results_evals, cbind(max.ncomp[[comp]], eta.list[[e]], r, f, n_var, aic, c_index))
            pb$tick()

          }#fold
        }#run
      }#eta
    }#component

    if(!all(is.null(df_results_evals))){
      colnames(df_results_evals) <- c("n.comps", "eta","runs", "fold", "n.var", "AIC", "c_index")
      df_results_evals <- as.data.frame(df_results_evals)
    }

  }

  if(attr(model, "model") %in% c(pkg.env$splsdrcox_dynamic, pkg.env$multiblock_methods) && !all(is.null(df_results_evals))){
    df_results_evals <- as.data.frame(df_results_evals)
    for(cn in colnames(df_results_evals)){
      if(cn=="n.var"){
        df_results_evals[,cn] <- factor(df_results_evals[,cn])
      }else{
        df_results_evals[,cn] <- as.numeric(df_results_evals[,cn])
      }
    }
  }

  return(df_results_evals)

}

#BRIER is a FOLD LEVEL ALWAYS
get_COX_evaluation_BRIER <- function(comp_model_lst,
                                     lst_X_test, lst_Y_test,
                                     df_results_evals, times = NULL,
                                     pred.method, pred.attr,
                                     max.ncomp, n_run, k_folds,
                                     w_BRIER,
                                     MIN_AUC_INCREASE, MIN_AUC, MIN_COMP_TO_CHECK,
                                     method.train, PARALLEL = F, verbose = F){

  fast_mode = T #fold level in BRIER

  if(length(max.ncomp)==1 & !method.train==pkg.env$coxEN){
    max.ncomp <- 1:max.ncomp
  }

  lst_BRIER_component = NULL
  lst_comp_BRIER <- NULL

  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_BRIER <- NULL #fast mode

  optimal_comp_index <- NULL
  optimal_eta <- NULL
  optimal_eta_index <- NULL
  optimal_comp_flag <- FALSE

  total_models <- nrow(df_results_evals)

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso
  pb$tick(0)

  message(paste0("Evaluating prediction acuracy with Brier Score..."))

  # EVAL AUC FOR EACH FOLD
  for(l in unique(df_results_evals$n.comps)){

    l.index <- which(l == unique(df_results_evals$n.comps))
    lst_BRIER_component_run <- NULL

    for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){

      lst_BRIER_component_folds <- NULL

      for(f in unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r,]$fold)){

        # non-significant models could be filtered, check if the model exist in df_results_evals
        if(nrow(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r & df_results_evals$fold==f,])==0){
          pb$tick()
          next
        }

        lst_BRIER <- getCOMPLETE_BRIER(comp_index = l.index, eta_index = NULL, run = r, fold = f,
                                       lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                       comp_model_lst = comp_model_lst, times = times,
                                       verbose = verbose)

        lst_BRIER_values <- lst_BRIER$brier_score
        lst_BRIER_component_folds[[f]] <- lst_BRIER_values$ierror
        df_results_evals_BRIER <- c(df_results_evals_BRIER, lst_BRIER_values$ierror)

        pb$tick()

      } #fold
      if(!is.null(lst_BRIER_component_folds)){
        names(lst_BRIER_component_folds) <- paste0("fold_",unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r,]$fold))
      }
      lst_BRIER_component_run[[r]] <- lst_BRIER_component_folds
    } #run

    if(!is.null(lst_BRIER_component_run)){
      names(lst_BRIER_component_run) <- paste0("run_",unique(df_results_evals[df_results_evals$n.comps==l,]$runs))
    }

    lst_BRIER_component[[l.index]] <- lst_BRIER_component_run

    #CHECK AUC EVOLUTION PER COMPONENT
    # lst_checkImprovement <- check_AUC_improvement(fast_mode = T, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_BRIER,
    #                                               comp_index = l.index, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_BRIER,
    #                                               MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp, method.train = method.train)
    # optimal_comp_index <- lst_checkImprovement$optimal_comp_index
    # optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
    # lst_comp_BRIER <- lst_checkImprovement$lst_comp_AUC

    # if(optimal_comp_flag){
    #   break
    # }

  } #lambda

  #### ### ### ###
  # GET RESULTS #
  #### ### ### ###
  txt <- NULL
  if(method.train==pkg.env$coxEN){
    txt <- "lambda_"
  }else{
    txt <- "comp_"
  }

  # if(optimal_comp_flag){
  #   #names(lst_AUC_component) <- paste0(txt,max.ncomp[1:min(max(max.ncomp),(optimal_comp_index+MIN_COMP_TO_CHECK))])
  #   names(lst_AUC_component) <- paste0(txt,unique(df_results_evals$n.comps)[optimal_comp_index:(optimal_comp_index+MIN_COMP_TO_CHECK)])
  # }else{
  #   names(lst_AUC_component) <- paste0(txt,max.ncomp)
  # }

  names(lst_BRIER_component) <- paste0(txt,max.ncomp)

  # if(fast_mode){ #AUC per FOLD - add info to df_results_evals
  #   if(optimal_comp_flag){
  #     df_results_evals$AUC <- c(df_results_evals_AUC, rep(NA, nrow(df_results_evals)-length(df_results_evals_AUC)))
  #   }else{
  #     df_results_evals$AUC <- df_results_evals_AUC
  #   }
  # }

  df_results_evals$BRIER <- df_results_evals_BRIER

  #### ### ### ### ### ### ### ###
  # TABLES FOR RUN AND FOLD LEVEL #
  #### ### ### ### ### ### ### ###

  #AUC per RUN AND COMP
  optimal_comp_index <- NULL
  lst_AUC_RUN_COMP <- getAUC_RUN_AND_COMP(mode = "BRIER", fast_mode = fast_mode, max.ncomp = max.ncomp,
                                          n_run = n_run, df_results_evals = df_results_evals,
                                          optimal_comp_flag = optimal_comp_flag,
                                          optimal_comp_index = optimal_comp_index, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          lst_AUC_component = lst_BRIER_component,
                                          df_results_evals_run = df_results_evals_run,
                                          df_results_evals_comp = df_results_evals_comp,
                                          method.train = method.train)

  df_results_evals_run <- lst_AUC_RUN_COMP$df_results_evals_run
  df_results_evals_comp <- lst_AUC_RUN_COMP$df_results_evals_comp

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run, df_results_evals_fold = df_results_evals,
              optimal_comp_index = optimal_comp_index, optimal_comp_flag = optimal_comp_flag))
}

get_COX_evaluation_BRIER_sPLS <- function(comp_model_lst,
                                          lst_X_test, lst_Y_test,
                                          df_results_evals, times = NULL,
                                          pred.method, pred.attr,
                                          max.ncomp, eta.list, n_run, k_folds,
                                          w_BRIER,
                                          MIN_AUC_INCREASE, MIN_AUC, MIN_COMP_TO_CHECK, method.train, PARALLEL = F, verbose = F){

  fast_mode = T #fold level in BRIER

  if(length(max.ncomp)==1 & !method.train==pkg.env$coxEN){
    max.ncomp <- 1:max.ncomp
  }

  lst_BRIER_component = NULL
  lst_comp_BRIER <- NULL

  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_BRIER <- NULL #fast mode

  optimal_comp_index <- NULL
  optimal_eta <- NULL
  optimal_eta_index <- NULL
  optimal_comp_flag <- FALSE

  total_models <- nrow(df_results_evals)

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso
  pb$tick(0)

  message(paste0("Evaluating prediction acuracy with Brier Score..."))

  # EVAL AUC FOR EACH FOLD
  for(l in unique(df_results_evals$n.comps)){

    l.index <- which(l == unique(df_results_evals$n.comps))
    lst_BRIER_component_eta <- NULL

    for(e in 1:length(eta.list)){

      lst_BRIER_component_run <- NULL

      for(r in unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$eta==eta.list[[e]],]$runs)){

        lst_BRIER_component_folds <- NULL

        for(f in unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$eta==eta.list[[e]] & df_results_evals$runs==r,]$fold)){

          # non-significant models could be filtered, check if the model exist in df_results_evals
          if(nrow(df_results_evals[df_results_evals$n.comps==l & df_results_evals$eta==eta.list[[e]] & df_results_evals$runs==r & df_results_evals$fold==f,])==0){
            pb$tick()
            next
          }

          lst_BRIER <- getCOMPLETE_BRIER(comp_index = l.index, eta_index = e, run = r, fold = f,
                                         lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                         comp_model_lst = comp_model_lst, times = times,
                                         verbose = verbose)

          lst_BRIER_values <- lst_BRIER$brier_score
          lst_BRIER_component_folds[[f]] <- lst_BRIER_values$ierror
          df_results_evals_BRIER <- c(df_results_evals_BRIER, lst_BRIER_values$ierror)

          pb$tick()

        } #fold
        if(!is.null(lst_BRIER_component_folds)){
          names(lst_BRIER_component_folds) <- paste0("fold_",unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$eta==eta.list[[e]] & df_results_evals$runs==r,]$fold))
        }
        lst_BRIER_component_run[[r]] <- lst_BRIER_component_folds
      } #run

      if(!is.null(lst_BRIER_component_run)){
        names(lst_BRIER_component_run) <- paste0("run_",unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$eta==eta.list[[e]],]$runs))
      }

      lst_BRIER_component_eta[[e]] <- lst_BRIER_component_run

      #CHECK AUC EVOLUTION PER COMPONENT
      # lst_checkImprovement <- check_AUC_improvement(fast_mode = T, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_BRIER,
      #                                               comp_index = l.index, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_BRIER,
      #                                               MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp, method.train = method.train)
      # optimal_comp_index <- lst_checkImprovement$optimal_comp_index
      # optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
      # lst_comp_BRIER <- lst_checkImprovement$lst_comp_AUC

      # if(optimal_comp_flag){
      #   break
      # }

    } #eta

    if(!is.null(lst_BRIER_component_eta)){
      names(lst_BRIER_component_eta) <- paste0("eta_",unique(df_results_evals[df_results_evals$n.comps==l,]$eta))
    }

    lst_BRIER_component[[l.index]] <- lst_BRIER_component_eta

  } #component

  #### ### ### ###
  # GET RESULTS #
  #### ### ### ###
  txt <- "comp_"

  # if(optimal_comp_flag){
  #   #names(lst_AUC_component) <- paste0(txt,max.ncomp[1:min(max(max.ncomp),(optimal_comp_index+MIN_COMP_TO_CHECK))])
  #   names(lst_AUC_component) <- paste0(txt,unique(df_results_evals$n.comps)[optimal_comp_index:(optimal_comp_index+MIN_COMP_TO_CHECK)])
  # }else{
  #   names(lst_AUC_component) <- paste0(txt,max.ncomp)
  # }

  names(lst_BRIER_component) <- paste0(txt,max.ncomp)

  # if(fast_mode){ #AUC per FOLD - add info to df_results_evals
  #   if(optimal_comp_flag){
  #     df_results_evals$AUC <- c(df_results_evals_AUC, rep(NA, nrow(df_results_evals)-length(df_results_evals_AUC)))
  #   }else{
  #     df_results_evals$AUC <- df_results_evals_AUC
  #   }
  # }

  df_results_evals$BRIER <- df_results_evals_BRIER

  #### ### ### ### ### ### ### ###
  # TABLES FOR RUN AND FOLD LEVEL #
  #### ### ### ### ### ### ### ###

  #AUC per RUN AND COMP
  optimal_comp_index <- NULL
  lst_AUC_RUN_COMP <- getAUC_RUN_AND_COMP_sPLS(mode = "BRIER", fast_mode = fast_mode,
                                               max.ncomp = max.ncomp, eta.list = eta.list,
                                               n_run = n_run, df_results_evals = df_results_evals,
                                               optimal_comp_flag = optimal_comp_flag,
                                               optimal_comp_index = optimal_comp_index, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                               lst_AUC_component = lst_BRIER_component,
                                               df_results_evals_run = df_results_evals_run,
                                               df_results_evals_comp = df_results_evals_comp,
                                               method.train = method.train)

  df_results_evals_run <- lst_AUC_RUN_COMP$df_results_evals_run
  df_results_evals_comp <- lst_AUC_RUN_COMP$df_results_evals_comp

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run, df_results_evals_fold = df_results_evals,
              optimal_comp_index = optimal_comp_index, optimal_comp_flag = optimal_comp_flag))
}

get_COX_evaluation_AUC <- function(comp_model_lst,
                                   lst_X_test, lst_Y_test,
                                   df_results_evals, times = NULL,
                                   fast_mode, pred.method, pred.attr,
                                   max.ncomp, n_run, k_folds,
                                   w_AUC,
                                   MIN_AUC_INCREASE, MIN_AUC, MIN_COMP_TO_CHECK, method.train, PARALLEL = F, verbose = F){

  if(length(max.ncomp)==1 & !method.train==pkg.env$coxEN){
    max.ncomp <- 1:max.ncomp
  }

  lst_AUC_component = NULL
  lst_comp_AUC <- NULL

  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_AUC <- NULL #fast mode

  optimal_comp_index <- NULL
  optimal_eta <- NULL
  optimal_eta_index <- NULL
  optimal_comp_flag <- FALSE

  total_models <- ifelse(!fast_mode, nrow(unique(df_results_evals[,c("n.comps", "runs")])), nrow(df_results_evals))

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso
  pb$tick(0)

  message(paste0("Evaluating prediction acuracy with ", pred.method ," algorithm...", ifelse(fast_mode, " [FAST_MODE]", " [BEST_MODE]")))

  if(fast_mode){ # EVAL AUC FOR EACH FOLD

    for(l in unique(df_results_evals$n.comps)){

      l.index <- which(l == unique(df_results_evals$n.comps))
      lst_AUC_component_run <- NULL

      for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){

        lst_AUC_component_folds <- NULL

        for(f in unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r,]$fold)){

          # non-significant models could be filtered, check if the model exist in df_results_evals
          if(nrow(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r & df_results_evals$fold==f,])==0){
            pb$tick()
            next
          }

          lst_FAST_LP_AUC <- getFAST_LP_AUC(fast_mode = fast_mode, comp_index = l.index, run = r, fold = f,
                                            lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, times = times,
                                            comp_model_lst = comp_model_lst, lst_linear.predictors = lst_linear.predictors,
                                            df_results_evals_AUC = df_results_evals_AUC,
                                            pred.method = pred.method, pred.attr = pred.attr, PARALLEL = F, verbose = verbose)

          lst_AUC_values <- lst_FAST_LP_AUC$lst_AUC_values
          df_results_evals_AUC <- lst_FAST_LP_AUC$df_results_evals_AUC
          lst_AUC_component_folds[[f]] <- lst_AUC_values

          pb$tick()

        } #fold
        if(!is.null(lst_AUC_component_folds)){
          names(lst_AUC_component_folds) <- paste0("fold_",unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r,]$fold))
        }
        lst_AUC_component_run[[r]] <- lst_AUC_component_folds
      } #run
      if(!is.null(lst_AUC_component_run)){
        names(lst_AUC_component_run) <- paste0("run_",unique(df_results_evals[df_results_evals$n.comps==l,]$runs))
      }
      lst_AUC_component[[l.index]] <- lst_AUC_component_run

      #CHECK AUC EVOLUTION PER COMPONENT
      lst_checkImprovement <- check_AUC_improvement(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                    comp_index = l.index, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                    MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp, method.train = method.train)
      optimal_comp_index <- lst_checkImprovement$optimal_comp_index
      optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
      lst_comp_AUC <- lst_checkImprovement$lst_comp_AUC

      if(optimal_comp_flag){
        break
      }

    } #lambda

  }else{ # COMPLETE MODE

    for(l in unique(df_results_evals$n.comps)){
      l.index <- which(l == unique(df_results_evals$n.comps))
      lst_AUC_component_run <- NULL

      for(r in unique(df_results_evals[df_results_evals$n.comps==l,]$runs)){

        Y_test_full <- NULL
        lst_linear.predictors <- NULL

        for(f in unique(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r,]$fold)){

          # non-significant models could be filtered, check if the model exist in df_results_evals
          if(nrow(df_results_evals[df_results_evals$n.comps==l & df_results_evals$runs==r & df_results_evals$fold==f,])==0){
            next
          }

          #comp_index is the index of l for coxEN
          # method updates automaticatly the lst of linear predictors addind each fold
          lst_COMPLETE_LP <- getCOMPLETE_LP(comp_index = l.index, run = r, fold = f,
                                            lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, Y_test_full = Y_test_full,
                                            comp_model_lst = comp_model_lst, lst_linear.predictors = lst_linear.predictors)

          Y_test_full <- lst_COMPLETE_LP$Y_test_full
          lst_linear.predictors <- lst_COMPLETE_LP$lst_linear.predictors

        } #fold

        if(is.null(lst_linear.predictors)){
          pb$tick()
          next #no models computed
        }

        lst_resCOMPLETE_LP_AUC <- getCOMPLETE_LP_AUC(Y_test_full = Y_test_full, lst_linear.predictors = lst_linear.predictors,
                                                     times = times, df_results_evals_AUC = df_results_evals_AUC,
                                                     pred.method = pred.method, pred.attr = pred.attr, PARALLEL = PARALLEL, verbose = verbose)
        lst_AUC_values <- lst_resCOMPLETE_LP_AUC$lst_AUC_values
        df_results_evals_AUC <- lst_resCOMPLETE_LP_AUC$df_results_evals_AUC

        if(all(is.na(lst_AUC_values$AUC))){
          lst_AUC_values$AUC <- NA
        }

        lst_AUC_component_run[[r]] <- lst_AUC_values
        pb$tick()

      } #run

      names(lst_AUC_component_run) <- paste0("run_",unique(df_results_evals[df_results_evals$n.comps==l,]$runs))
      lst_AUC_component[[l.index]] <- lst_AUC_component_run

      #CHECK AUC EVOLUTION
      #CHECK AUC EVOLUTION PER COMPONENT
      lst_checkImprovement <- check_AUC_improvement(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                    comp_index = l.index, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                    MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp, method.train = method.train)
      optimal_comp_index <- lst_checkImprovement$optimal_comp_index
      optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
      lst_comp_AUC <- lst_checkImprovement$lst_comp_AUC

      if(optimal_comp_flag){
        break
      }

    } #component
  } #complete mode


  #### ### ### ###
  # GET RESULTS #
  #### ### ### ###
  txt <- NULL
  if(method.train==pkg.env$coxEN){
    txt <- "lambda_"
  }else{
    txt <- "comp_"
  }

  if(optimal_comp_flag){
    #names(lst_AUC_component) <- paste0(txt,max.ncomp[1:min(max(max.ncomp),(optimal_comp_index+MIN_COMP_TO_CHECK))])
    names(lst_AUC_component) <- paste0(txt,unique(df_results_evals$n.comps)[1:(optimal_comp_index+MIN_COMP_TO_CHECK)])
  }else{
    names(lst_AUC_component) <- paste0(txt,max.ncomp)
  }

  # AUC per FOLD - add info to df_results_evals
  # otherwise, the information is store in another site
  if(fast_mode){
    if(optimal_comp_flag){
      df_results_evals$AUC <- c(df_results_evals_AUC, rep(NA, nrow(df_results_evals)-length(df_results_evals_AUC)))
    }else{
      df_results_evals$AUC <- df_results_evals_AUC
    }
  }

  #### ### ### ### ### ### ### ###
  # TABLES FOR RUN AND FOLD LEVEL #
  #### ### ### ### ### ### ### ###

  #AUC per RUN AND COMP
  lst_AUC_RUN_COMP <- getAUC_RUN_AND_COMP(mode = "AUC", fast_mode = fast_mode, max.ncomp = max.ncomp, n_run = n_run, df_results_evals = df_results_evals,
                                          optimal_comp_flag = optimal_comp_flag, optimal_comp_index = optimal_comp_index, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          lst_AUC_component = lst_AUC_component, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp, method.train = method.train)

  df_results_evals_run <- lst_AUC_RUN_COMP$df_results_evals_run
  df_results_evals_comp <- lst_AUC_RUN_COMP$df_results_evals_comp

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run, df_results_evals_fold = df_results_evals,
              optimal_comp_index = optimal_comp_index, optimal_comp_flag = optimal_comp_flag))
}

get_COX_evaluation_AUC_sPLS <- function(comp_model_lst,
                                        lst_X_test, lst_Y_test,
                                        df_results_evals, times = NULL,
                                        fast_mode, pred.method, pred.attr,
                                        max.ncomp, eta.list, n_run, k_folds,
                                        w_AUC,
                                        MIN_AUC_INCREASE, MIN_AUC, MIN_COMP_TO_CHECK, method.train,
                                        PARALLEL = F, verbose = F){

  if(length(max.ncomp)==1){
    max.ncomp <- 1:max.ncomp
  }

  lst_AUC_component = NULL
  lst_comp_AUC <- NULL

  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_AUC <- NULL #fast mode

  optimal_comp_index <- NULL
  optimal_eta <- NULL
  optimal_eta_index <- NULL
  optimal_comp_flag <- FALSE

  total_models <- ifelse(!fast_mode, nrow(unique(df_results_evals[,c("n.comps", "runs", "eta")])), nrow(df_results_evals))

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso
  pb$tick(0)

  message(paste0("Evaluating prediction acuracy with ", pred.method ," algorithm...", ifelse(fast_mode, " [FAST_MODE]", " [BEST_MODE]")))

  if(fast_mode){ # EVAL AUC FOR EACH FOLD

    for(l in 1:length(max.ncomp)){

      lst_AUC_eta_results <- NULL

      for(e in 1:length(eta.list)){

        lst_AUC_component_run <- NULL

        for(r in 1:n_run){

          lst_AUC_component_folds <- NULL

          for(f in 1:k_folds){

            lst_FAST_LP_AUC <- getFAST_LP_AUC(fast_mode = fast_mode, comp_index = l, eta_index = e, run = r, fold = f,
                                              lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, times = times,
                                              comp_model_lst = comp_model_lst, lst_linear.predictors = lst_linear.predictors,
                                              df_results_evals_AUC = df_results_evals_AUC,
                                              pred.method = pred.method, pred.attr = pred.attr, PARALLEL = PARALLEL, verbose = verbose)

            lst_AUC_values <- lst_FAST_LP_AUC$lst_AUC_values
            df_results_evals_AUC <- lst_FAST_LP_AUC$df_results_evals_AUC
            lst_AUC_component_folds[[f]] <- lst_AUC_values

            pb$tick()

          } #fold
          names(lst_AUC_component_folds) <- paste0("fold_",1:k_folds)
          lst_AUC_component_run[[r]] <- lst_AUC_component_folds
        } #run
        names(lst_AUC_component_run) <- paste0("run_",1:n_run)
        lst_AUC_eta_results[[e]] <- lst_AUC_component_run

        #CHECK AUC EVOLUTION PER COMPONENT
        lst_comp_AUC <- check_AUC_improvement_spls1(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                            comp_index = l, eta_index = e, eta.list = eta.list, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                            MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, max.ncomp = max.ncomp)

      } #eta
      names(lst_AUC_eta_results) <- paste0("eta_", eta.list)
      lst_AUC_component[[l]] <- lst_AUC_eta_results

      #CHECK AUC EVOLUTION PER COMPONENT
      lst_checkImprovement <- check_AUC_improvement_spls2(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                          comp_index = l, eta_index = e, eta.list = eta.list, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                          MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp)
      optimal_comp_index <- lst_checkImprovement$optimal_comp_index
      optimal_eta_index <- lst_checkImprovement$optimal_eta_index
      optimal_eta <- lst_checkImprovement$optimal_eta
      optimal_auc <- lst_checkImprovement$optimal_auc
      optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
      lst_comp_AUC <- lst_checkImprovement$lst_comp_AUC

      #For at least MIN_COMP_TO_CHECK

      if(optimal_comp_flag){
        break
      }

    } #component

  }else{ # COMPLETE MODE

    for(l in 1:length(max.ncomp)){

      lst_AUC_eta_results <- NULL

      for(e in 1:length(eta.list)){

        lst_AUC_component_run <- NULL

        for(r in 1:n_run){

          Y_test_full <- NULL
          lst_linear.predictors <- NULL

          for(f in 1:k_folds){

            lst_COMPLETE_LP <- getCOMPLETE_LP(comp_index = l, eta_index = e, run = r, fold = f,
                                              lst_X_test = lst_X_test, lst_Y_test = lst_Y_test, Y_test_full = Y_test_full,
                                              comp_model_lst = comp_model_lst, lst_linear.predictors = lst_linear.predictors)
            Y_test_full <- lst_COMPLETE_LP$Y_test_full
            lst_linear.predictors <- lst_COMPLETE_LP$lst_linear.predictors

          } #fold

          lst_resCOMPLETE_LP_AUC <- getCOMPLETE_LP_AUC(Y_test_full = Y_test_full, lst_linear.predictors = lst_linear.predictors,
                                                       times = times, df_results_evals_AUC = df_results_evals_AUC, pred.method = pred.method, pred.attr = pred.attr, PARALLEL = PARALLEL, verbose = verbose)
          lst_AUC_values <- lst_resCOMPLETE_LP_AUC$lst_AUC_values
          df_results_evals_AUC <- lst_resCOMPLETE_LP_AUC$df_results_evals_AUC

          lst_AUC_component_run[[r]] <- lst_AUC_values
          pb$tick()

        } #run

        names(lst_AUC_component_run) <- paste0("run_",1:n_run)
        lst_AUC_eta_results[[e]] <- lst_AUC_component_run

        #CHECK AUC EVOLUTION PER COMPONENT
        lst_comp_AUC <- check_AUC_improvement_spls1(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                    comp_index = l, eta_index = e, eta.list = eta.list, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                    MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, max.ncomp = max.ncomp)

      } #eta

      names(lst_AUC_eta_results) <- paste0("eta_", eta.list)
      lst_AUC_component[[l]] <- lst_AUC_eta_results

      #CHECK AUC EVOLUTION PER COMPONENT
      lst_checkImprovement <- check_AUC_improvement_spls2(fast_mode = fast_mode, pred.attr = pred.attr, df_results_evals_AUC = df_results_evals_AUC,
                                                          comp_index = l, eta_index = e, eta.list = eta.list, n_run = n_run, k_folds = k_folds, lst_comp_AUC = lst_comp_AUC,
                                                          MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK, MIN_AUC = MIN_AUC, MIN_AUC_INCREASE = MIN_AUC_INCREASE, max.ncomp = max.ncomp)
      optimal_comp_index <- lst_checkImprovement$optimal_comp_index
      optimal_eta_index <- lst_checkImprovement$optimal_eta_index
      optimal_eta <- lst_checkImprovement$optimal_eta
      optimal_auc <- lst_checkImprovement$optimal_auc
      optimal_comp_flag <- lst_checkImprovement$optimal_comp_flag
      lst_comp_AUC <- lst_checkImprovement$lst_comp_AUC

      #For at least MIN_COMP_TO_CHECK

      if(optimal_comp_flag){
        break
      }

    } #comp
  } #complete mode

  #### ### ### ###
  # GET RESULTS #
  #### ### ### ###
  txt <- NULL
  if(method.train==pkg.env$coxEN){
    txt <- "lambda_"
  }else{
    txt <- "comp_"
  }

  if(optimal_comp_flag){
    names(lst_AUC_component) <- paste0(txt,max.ncomp[1:min(max(max.ncomp),(optimal_comp_index+MIN_COMP_TO_CHECK))])
  }else{
    names(lst_AUC_component) <- paste0(txt,max.ncomp)
  }

  if(fast_mode){ #AUC per FOLD - add info to df_results_evals
    if(optimal_comp_flag){
      df_results_evals$AUC <- c(df_results_evals_AUC, rep(NA, nrow(df_results_evals)-length(df_results_evals_AUC)))
    }else{
      df_results_evals$AUC <- df_results_evals_AUC
    }
  }

  #### ### ### ### ### ### ### ###
  # TABLES FOR RUN AND FOLD LEVEL #
  #### ### ### ### ### ### ### ###

  #AUC per RUN AND COMP
  st_AUC_RUN_COMP <- getAUC_RUN_AND_COMP_sPLS(mode = "AUC", fast_mode = fast_mode, max.ncomp = max.ncomp, n_run = n_run, eta.list = eta.list, df_results_evals = df_results_evals,
                                              optimal_comp_flag = optimal_comp_flag, optimal_comp_index = optimal_comp_index, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                              lst_AUC_component = lst_AUC_component, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                              method.train = method.train)

  df_results_evals_run <- st_AUC_RUN_COMP$df_results_evals_run
  df_results_evals_comp <- st_AUC_RUN_COMP$df_results_evals_comp

  if(method.train==pkg.env$splsicox){
    colnames(df_results_evals)[which(colnames(df_results_evals)=="eta")] <- "spv_penalty"
    colnames(df_results_evals_run)[which(colnames(df_results_evals_run)=="eta")] <- "spv_penalty"
    colnames(df_results_evals_comp)[which(colnames(df_results_evals_comp)=="eta")] <- "spv_penalty"
  }

  return(list(df_results_evals_comp = df_results_evals_comp, df_results_evals_run = df_results_evals_run, df_results_evals_fold = df_results_evals,
              optimal_comp_index = optimal_comp_index, optimal_eta_index = optimal_eta_index, optimal_eta = optimal_eta, optimal_comp_flag = optimal_comp_flag))
}

getSubModel <- function(model, comp, remove_non_significant){
  res <- model

  if(all(is.na(res))){
    return(NA)
  }

  #X
  if("loadings" %in% names(res$X) && !all(is.na(res$X$loadings))){
    res$X$loadings <- res$X$loadings[,1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("weightings" %in% names(res$X) && !all(is.na(res$X$weightings))){
    res$X$weightings <- res$X$weightings[,1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("weightings_norm" %in% names(res$X) && !all(is.na(res$X$weightings_norm))){
    res$X$weightings_norm <- res$X$weightings_norm[,1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("W.star" %in% names(res$X) && !all(is.na(res$X$W.star))){
    res$X$W.star <- res$X$W.star[,1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("scores" %in% names(res$X) && !all(is.na(res$X$scores))){
    res$X$scores <- res$X$scores[,1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("E" %in% names(res$X) && !all(is.na(res$X$E))){
    res$X$E <- res$X$E[1:min(ncol(res$X$loadings),comp), drop=F]
  }

  if("var_by_component" %in% names(res) && !all(is.na(res$var_by_component))){
    res$var_by_component <- res$var_by_component[1:min(ncol(res$X$loadings),comp), drop=F]
  }

  res$n.comp <- comp

  #survival_model
  if(!all(is.na(res$X$scores)) & !all(is.na(res$Y$data))){
    cox_model <- cox(X = res$X$scores, Y = res$Y$data, x.center = F, x.scale = F, y.center = F, y.scale = F, remove_non_significant = remove_non_significant, FORCE = T)
    survival_model <- cox_model$survival_model
    res$survival_model <- survival_model
  }

  return(res)
}

getSubModel.mb <- function(model, comp, remove_non_significant){
  res <- model

  if(all(is.na(res))){
    return(NA)
  }

  #X
  if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
    t1 <- Sys.time()
    data <- NULL
    col_names <- NULL
    for(b in names(model$list_spls_models)){
      res$list_spls_models[[b]]$X$loadings <- res$list_spls_models[[b]]$X$loadings[,1:min(ncol(res$list_spls_models[[b]]$X$loadings),comp), drop=F]
      res$list_spls_models[[b]]$X$weightings <- res$list_spls_models[[b]]$X$weightings[,1:min(ncol(res$list_spls_models[[b]]$X$weightings),comp), drop=F]
      res$list_spls_models[[b]]$X$W.star <- res$list_spls_models[[b]]$X$W.star[,1:min(ncol(res$list_spls_models[[b]]$X$W.star),comp),drop=F]
      res$list_spls_models[[b]]$X$scores <- res$list_spls_models[[b]]$X$scores[,1:min(ncol(res$list_spls_models[[b]]$X$scores),comp), drop=F]

      res$list_spls_models[[b]]$var_by_component <- res$list_spls_models[[b]]$var_by_component[1:min(ncol(res$list_spls_models[[b]]$X$scores),comp)]
      res$list_spls_models[[b]]$n.comp <- comp

      #survival_model
      data <- as.data.frame(cbind(data, res$list_spls_models[[b]]$X$scores[,,drop=F]))
      col_names <- c(col_names, paste0(colnames(res$list_spls_models[[b]]$X$scores), "_", b))
    }

    colnames(data) <- col_names
    #survival_model
    cox_model <- cox(X = data, Y = res$Y$data, x.center = F, x.scale = F, y.center = F, y.scale = F, remove_non_significant = remove_non_significant, FORCE = T)
    survival_model <- cox_model$survival_model

    res$survival_model <- survival_model
    res$n.comp <- comp
    t2 <- Sys.time()
    res$time <- difftime(t2,t1,units = "mins")

  }else if(attr(model, "model") %in% pkg.env$sb.splsicox){
    t1 <- Sys.time()
    data <- NULL
    col_names <- NULL
    for(b in names(model$list_pls_models)){
      res$list_pls_models[[b]]$X$loadings <- res$list_pls_models[[b]]$X$loadings[,1:min(ncol(res$list_pls_models[[b]]$X$loadings),comp), drop=F]
      res$list_pls_models[[b]]$X$weightings <- res$list_pls_models[[b]]$X$weightings[,1:min(ncol(res$list_pls_models[[b]]$X$weightings),comp), drop=F]
      res$list_pls_models[[b]]$X$W.star <- res$list_pls_models[[b]]$X$W.star[,1:min(ncol(res$list_pls_models[[b]]$X$W.star),comp),drop=F]
      res$list_pls_models[[b]]$X$scores <- res$list_pls_models[[b]]$X$scores[,1:min(ncol(res$list_pls_models[[b]]$X$scores),comp), drop=F]

      res$list_pls_models[[b]]$var_by_component <- res$list_pls_models[[b]]$var_by_component[1:min(ncol(res$list_pls_models[[b]]$X$scores),comp)]
      res$list_pls_models[[b]]$n.comp <- comp

      #survival_model
      data <- as.data.frame(cbind(data, res$list_pls_models[[b]]$X$scores[,,drop=F]))
      col_names <- c(col_names, paste0(colnames(res$list_pls_models[[b]]$X$scores), "_", b))
    }

    colnames(data) <- col_names
    #survival_model
    cox_model <- cox(X = data, Y = res$Y$data, x.center = F, x.scale = F, y.center = F, y.scale = F, remove_non_significant = remove_non_significant, FORCE = T)
    survival_model <- cox_model$survival_model

    res$survival_model <- survival_model
    res$n.comp <- comp
    t2 <- Sys.time()
    res$time <- difftime(t2,t1,units = "mins")
  }else if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){ #revisar!!!
    t1 <- Sys.time()
    for(b in names(model$mb.model$X)){
      res$X$loadings[[b]] <- res$X$loadings[[b]][,1:min(ncol(res$X$loadings[[b]]),comp), drop=F]
      res$X$weightings[[b]] <- res$X$weightings[[b]][,1:min(ncol(res$X$weightings[[b]]),comp), drop=F]
      res$X$W.star[[b]] <- res$X$W.star[[b]][,1:min(ncol(res$X$W.star[[b]]),comp),drop=F]
      res$X$scores[[b]] <- res$X$scores[[b]][,1:min(ncol(res$X$scores[[b]]),comp), drop=F]
    }

    #survival_model
    data <- as.data.frame(res$X$scores[[1]][,,drop=F])
    for(b in names(res$X$scores)[2:length(res$X$scores)]){
      data <- cbind(data, as.data.frame(res$X$scores[[b]][,,drop=F]))
    }

    colnames(data) <- apply(expand.grid(colnames(res$X$scores[[1]]), names(res$X$scores)), 1, paste, collapse="_")
    #survival_model
    cox_model <- cox(X = data, Y = res$Y$data, x.center = F, x.scale = F, y.center = F, y.scale = F, remove_non_significant = remove_non_significant, FORCE = T)
    survival_model <- cox_model$survival_model

    res$survival_model <- survival_model
    res$n.comp <- comp #revisar!!!
    t2 <- Sys.time()
    res$time <- difftime(t2,t1,units = "mins")
  }

  return(res)
}

get_HDCOX_models2.0 <- function(method = "sPLS-ICOX",
                                lst_X_train, lst_Y_train, vector = NULL,
                                max.ncomp, eta.list = NULL, EN.alpha.list = NULL, max.variables = 15,
                                n_run, k_folds,
                                MIN_NVAR = 10, MAX_NVAR = 10000, MIN_AUC_INCREASE = 0.01, n.cut_points = 5, EVAL_METHOD = "AUC",
                                x.center, x.scale, y.center, y.scale,
                                remove_near_zero_variance = F, remove_zero_variance = F,  toKeep.zv = NULL,
                                remove_non_significant = F,
                                alpha = 0.05, max.iter = 500, returnData = F,
                                total_models, MIN_EPV = 0, tol = 1e-15, PARALLEL = F, verbose = F){

  comp_model_lst <- list()
  fold_list <- list()
  run_list <- list()
  eta_model_lst <- NULL
  info <- NULL # for sPLS

  ## CHECK METHOD
  if(is.null(eta.list) & is.null(EN.alpha.list) & !method %in% c(pkg.env$splsdacox_dynamic, pkg.env$splsdrcox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
    stop_quietly("Method must be one of 'sPLS-DACOX-Dynamic', 'MB.sPLS-DACOX' or 'sPLS-DRCOX-Dynamic' if 'eta.list' and 'EN.alpha.list' is NULL.")
  }else if(!is.null(eta.list) & is.null(EN.alpha.list)  & !method %in% c(pkg.env$splsicox, pkg.env$sb.splsicox, pkg.env$splsdrcox, pkg.env$sb.splsdrcox)){
    stop_quietly("Method must be 'sPLS-ICOX', 'SB.sPLS-ICOX', 'sPLS-DRCOX' or 'SB.sPLS-DRCOX' or if 'eta.list' is not NULL.")
  }else if(!is.null(EN.alpha.list) & !method %in% c(pkg.env$coxEN)){
    stop_quietly("Method must be 'coxEN' if 'EN.alpha.list' is not NULL.")
  }

  pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
  pb <- progress::progress_bar$new(format = pb_text,
                                   total = total_models,
                                   complete = "=",   # Caracteres de las iteraciones finalizadas
                                   incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                   current = ">",    # Caracter actual
                                   clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                   width = 100)      # Ancho de la barra de progreso

  message(paste0("Training all possible models for ", method, "..."))
  pb$tick(0)

  #### ### ### ### ### ### #
  # UPDATING GLOBALS SIZE #
  #### ### ### ### ### ### #
  MB = 4000
  bytes = MB*1024^2
  options(future.globals.maxSize = bytes)

  #### ### ### ### #
  # COMP-REP-FOLDS #
  #### ### ### ### #
  if(method %in% c(pkg.env$splsdrcox_dynamic, pkg.env$splsdacox_dynamic, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

    #function to compute all models at the same time - just last component
    lst_inputs <- list()
    cont = 1
    lst_names = NULL
    for(i in max.ncomp){
      for(r in 1:n_run){
        for(f in 1:k_folds){
          lst_inputs[[cont]] = list()
          lst_inputs[[cont]]$comp = i
          lst_inputs[[cont]]$run <- r
          lst_inputs[[cont]]$fold <- f
          lst_names <- c(lst_names, paste0(i, "_", r, "_", f))
          cont = cont + 1
        }
      }
    }

    names(lst_inputs) <- lst_names

    ## Sometimes, when you compute the model with the higher number of dimensions it fail and should be working with one lesser component,
    ## when you compute the model by iterations it is easy to do, but know the model itself has to compute its better number of components
    ## if using one fail (use one lesser) !!!! HAVE TO BE IMPLEMENTED !!!!

    if(PARALLEL){
      n_cores <- max(future::availableCores() - 1, 1)

      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(lst_inputs), n_cores))
      }else{
        future::plan("multisession", workers = min(length(lst_inputs), n_cores))
      }

      if(method==pkg.env$splsdacox_dynamic){
        lst_all_models <- furrr::future_map(lst_inputs, ~splsdacox_dynamic(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                           Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                           n.comp = .$comp,
                                                                           x.center = x.center, x.scale = x.scale,
                                                                           y.center = y.center, y.scale = y.scale,
                                                                           remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                           remove_non_significant = remove_non_significant,
                                                                           vector = vector,
                                                                           MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                                           MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                           EVAL_METHOD = EVAL_METHOD, tol = tol, alpha = alpha,
                                                                           MIN_EPV = MIN_EPV, returnData = returnData, verbose = verbose), .options = furrr_options(seed = TRUE))

      }else if(method==pkg.env$splsdrcox_dynamic){
        lst_all_models <- furrr::future_map(lst_inputs, ~splsdrcox_dynamic(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                            Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                            n.comp = .$comp,
                                                                            x.center = x.center, x.scale = x.scale,
                                                                            y.center = y.center, y.scale = y.scale,
                                                                            remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                            remove_non_significant = remove_non_significant,
                                                                            vector = vector,
                                                                            MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                                            MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                            EVAL_METHOD = EVAL_METHOD, tol = tol, alpha = alpha,
                                                                            MIN_EPV = MIN_EPV, returnData = returnData), .options = furrr_options(seed = TRUE))

      }else if(method==pkg.env$mb.splsdacox){
        lst_all_models <- furrr::future_map(lst_inputs, ~mb.splsdacox(X = lst_X_train[[.$run]][[.$fold]],
                                                                     Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                     n.comp = .$comp, vector = vector,
                                                                     MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                     x.center = x.center, x.scale = x.scale,
                                                                     y.center = y.center, y.scale = y.scale,
                                                                     remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                     remove_non_significant = remove_non_significant,  tol = tol,  alpha = alpha, max.iter = max.iter,
                                                                     MIN_EPV = MIN_EPV, returnData = returnData, verbose = verbose), .options = furrr_options(seed = TRUE))

      }else if(method==pkg.env$mb.splsdrcox){
        lst_all_models <- furrr::future_map(lst_inputs, ~mb.splsdrcox(X = lst_X_train[[.$run]][[.$fold]],
                                                                              Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                 n.comp = .$comp, vector = vector,
                                                                 MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                 x.center = x.center, x.scale = x.scale,
                                                                 y.center = y.center, y.scale = y.scale,
                                                                 remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                 remove_non_significant = remove_non_significant, tol = tol, alpha = alpha,
                                                                 MIN_EPV = MIN_EPV, returnData = returnData, verbose = verbose), .options = furrr_options(seed = TRUE))
      }
      future::plan("sequential")

    }else{

      if(method==pkg.env$splsdacox_dynamic){
        lst_all_models <- purrr::map(lst_inputs, ~splsdacox_dynamic(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                     Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                     n.comp = .$comp,
                                                                     x.center = x.center, x.scale = x.scale,
                                                                     y.center = y.center, y.scale = y.scale,
                                                                     MIN_EPV = MIN_EPV, remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                     remove_non_significant = remove_non_significant,
                                                                     vector = vector,
                                                                     MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                     EVAL_METHOD = EVAL_METHOD, tol = tol, alpha = alpha,
                                                                     returnData = returnData))

      }else if(method==pkg.env$splsdrcox_dynamic){
        lst_all_models <- purrr::map(lst_inputs, ~splsdrcox_dynamic(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                    Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                    n.comp = .$comp,
                                                                    x.center = x.center, x.scale = x.scale,
                                                                    y.center = y.center, y.scale = y.scale,
                                                                    remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                    remove_non_significant = remove_non_significant,
                                                                    vector = vector,
                                                                    MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                                    MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                                    EVAL_METHOD = EVAL_METHOD, tol = tol, alpha = alpha,
                                                                    MIN_EPV = MIN_EPV, returnData = returnData))

      }else if(method==pkg.env$mb.splsdacox){
        lst_all_models <- purrr::map(lst_inputs, ~mb.splsdacox(X = lst_X_train[[.$run]][[.$fold]],
                                                              Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                              n.comp = .$comp, vector = vector,
                                                              MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                              x.center = x.center, x.scale = x.scale,
                                                              y.center = y.center, y.scale = y.scale,
                                                              remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                              remove_non_significant = remove_non_significant,  tol = tol,  alpha = alpha, max.iter = max.iter,
                                                              MIN_EPV = MIN_EPV, returnData = returnData, verbose = verbose))

      }else if(method==pkg.env$mb.splsdrcox){
        lst_all_models <- purrr::map(lst_inputs, ~mb.splsdrcox(X = lst_X_train[[.$run]][[.$fold]],
                                                               Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                               n.comp = .$comp, vector = vector,
                                                               MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                               x.center = x.center, x.scale = x.scale,
                                                               y.center = y.center, y.scale = y.scale,
                                                               remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                               remove_non_significant = remove_non_significant, tol = tol, alpha = alpha,
                                                               MIN_EPV = MIN_EPV, returnData = returnData))
      }

    }

    rm(lst_inputs)

    ## We need to return a list of lists: COMP->REP->FOLDS
    cont_problems = 0
    comp_model_lst <- list()
    for(c in max.ncomp){
      run_model_lst <- list()
      for(r in 1:n_run){
        fold_model_lst <- list()
        for(f in 1:k_folds){
          name <- paste0(c, "_", r, "_", f)
          if(all(is.na(lst_all_models[[name]])) || all(is.null(lst_all_models[[name]]$survival_model$fit))){
            cont_problems = cont_problems + 1
          }
          fold_model_lst[[f]] = lst_all_models[[name]]
        }
        names(fold_model_lst) <- paste0("fold_",1:k_folds)
        run_model_lst[[r]] <- fold_model_lst
      }
      names(run_model_lst) <- paste0("run_",1:n_run)
      comp_model_lst[[c]] <- run_model_lst
    }
    names(comp_model_lst) <- paste0("comp_",1:max.ncomp)

    ## Before compute all intermediate models, check if all problems
    if(cont_problems == n_run * k_folds){
      # if(verbose){
      #   message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))
      # }
      return(NULL)
    }

    ## We need to fill models from 1:max.ncomp (it uses max.ncomp to fill the others, if it is NULL, method fail!!!)
    pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
    pb <- progress::progress_bar$new(format = pb_text,
                                     total = (max.ncomp-1) * n_run * k_folds,
                                     complete = "=",   # Caracteres de las iteraciones finalizadas
                                     incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                     current = ">",    # Caracter actual
                                     clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                     width = 100)      # Ancho de la barra de progreso

    message(paste0("Creating sub-models for ", method, "..."))
    pb$tick(0)

    for(comp in 1:(max.ncomp-1)){
      run_model_lst <- list()
      for(r in 1:n_run){
        fold_model_lst <- list()
        for(f in 1:k_folds){
          if(method %in% c(pkg.env$sb.splsicox, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
            fold_model <- getSubModel.mb(model = comp_model_lst[[max.ncomp]][[r]][[f]], comp = comp, remove_non_significant = remove_non_significant)
          }else{
            fold_model <- getSubModel(model = comp_model_lst[[max.ncomp]][[r]][[f]], comp = comp, remove_non_significant = remove_non_significant)
          }
          fold_model_lst[[f]] <- fold_model
          pb$tick()
        }
        names(fold_model_lst) <- paste0("fold_",1:k_folds)
        run_model_lst[[r]] <- fold_model_lst

      }
      names(run_model_lst) <- paste0("run_",1:n_run)
      comp_model_lst[[comp]] <- run_model_lst
    }

    #### ### ### ### ### ### #
    # UPDATING GLOBALS SIZE #
    #### ### ### ### ### ### #
    MB = 500
    bytes = MB*1024^2
    options(future.globals.maxSize = bytes)

    return(comp_model_lst)
  }

  #### ### ### ### ##
  # ALPHA-REP-FOLDS #
  #### ### ### ### ##
  if(method %in% c(pkg.env$coxEN)){

    #function to compute all models at the same time
    lst_inputs <- list()
    lst_names = NULL
    cont = 1
    for(i in 1:length(EN.alpha.list)){
      for(r in 1:n_run){
        for(f in 1:k_folds){
          lst_inputs[[cont]] = list()
          lst_inputs[[cont]]$alpha_index = i
          lst_inputs[[cont]]$run <- r
          lst_inputs[[cont]]$fold <- f
          lst_names <- c(lst_names, paste0(i, "_", r, "_", f))
          cont = cont + 1
        }
      }
    }

    names(lst_inputs) <- lst_names

    if(PARALLEL){
      n_cores <- max(future::availableCores() - 1, 1)

      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(lst_inputs), n_cores))
      }else{
        future::plan("multisession", workers = min(length(lst_inputs), n_cores))
      }

      if(method==pkg.env$coxEN){
        lst_all_models <- furrr::future_map(lst_inputs, ~coxEN(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                               Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                               EN.alpha = EN.alpha.list[[.$alpha_index]],
                                                               max.variables = max.variables,
                                                               x.center = x.center, x.scale = x.scale,
                                                               y.center = y.center, y.scale = y.scale,
                                                               remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                               remove_non_significant = remove_non_significant,
                                                               alpha = alpha, MIN_EPV = MIN_EPV, verbose = verbose,
                                                               returnData = returnData), .options = furrr_options(seed = TRUE))

        #test with for:
        for(i in lst_inputs){
          coxEN(X = data.matrix(lst_X_train[[lst_inputs[[1]]$run]][[lst_inputs[[1]]$fold]]),
                Y = data.matrix(lst_Y_train[[lst_inputs[[1]]$run]][[lst_inputs[[1]]$fold]]),
                EN.alpha = EN.alpha.list[[lst_inputs[[1]]$alpha_index]],
                max.variables = max.variables,
                x.center = x.center, x.scale = x.scale,
                y.center = y.center, y.scale = y.scale,
                remove_non_significant = remove_non_significant,
                remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                alpha = alpha, MIN_EPV = MIN_EPV, verbose = verbose,
                returnData = F)
        }

      }

      future::plan("sequential")

    }else{

      if(method==pkg.env$coxEN){
        lst_all_models <- purrr::map(lst_inputs, ~coxEN(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                        Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                        EN.alpha = EN.alpha.list[[.$alpha_index]],
                                                        max.variables = max.variables,
                                                        x.center = x.center, x.scale = x.scale,
                                                        y.center = y.center, y.scale = y.scale,
                                                        remove_non_significant = remove_non_significant,
                                                        remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                        alpha = alpha, MIN_EPV = MIN_EPV, verbose = verbose,
                                                        returnData = returnData))
      }

    }

    rm(lst_inputs)

    ## We need to return a list of lists: ALPHA->REP->FOLDS
    comp_model_lst <- list()
    cont_problems = 0
    for(c in 1:length(EN.alpha.list)){
      run_model_lst <- list()
      for(r in 1:n_run){
        fold_model_lst <- list()
        for(f in 1:k_folds){
          name <- paste0(c, "_", r, "_", f)
          if(all(is.na(lst_all_models[[name]])) || all(is.null(lst_all_models[[name]]$survival_model$fit))){
            cont_problems = cont_problems + 1
          }
          fold_model_lst[[f]] = lst_all_models[[name]]
        }
        names(fold_model_lst) <- paste0("fold_",1:k_folds)
        run_model_lst[[r]] <- fold_model_lst
      }
      names(run_model_lst) <- paste0("run_",1:n_run)
      comp_model_lst[[c]] <- run_model_lst
    }
    names(comp_model_lst) <- paste0("alpha_",EN.alpha.list)

    ## Before compute all intermediate models, check if all problems
    if(cont_problems == n_run * k_folds){
      if(verbose){
        message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))
      }
      return(NULL)
    }

    #### ### ### ### ### ### #
    # UPDATING GLOBALS SIZE #
    #### ### ### ### ### ### #
    MB = 500
    bytes = MB*1024^2
    options(future.globals.maxSize = bytes)

    return(comp_model_lst)
  }

  #### ### ### ### ### ### #
  # COMP-VECTOR-REP-FOLDS #
  #### ### ### ### ### ### #
  if(method %in% c(pkg.env$splsicox, pkg.env$sb.splsicox, pkg.env$splsdrcox, pkg.env$sb.splsdrcox)){

    #function to compute all models at the same time
    lst_inputs <- list()
    cont = 1
    lst_names = NULL
    for(c in max.ncomp){ #computing all components it is the same as computing by iterations
      for(e in 1:length(eta.list)){
        for(r in 1:n_run){
          for(f in 1:k_folds){
            lst_inputs[[cont]] = list()
            lst_inputs[[cont]]$comp <- c
            lst_inputs[[cont]]$eta_index <- e
            lst_inputs[[cont]]$run <- r
            lst_inputs[[cont]]$fold <- f

            lst_names <- c(lst_names, paste0(c, "_", e, "_", r, "_", f))
            cont = cont + 1
          }
        }
      }
    }

    names(lst_inputs) <- lst_names

    ## Sometimes, when you compute the model with the higher number of dimensions it fail and should be working with one lesser component,
    ## when you compute the model by iterations it is easy to do, but know the model itself has to compute its better number of components
    ## if using one fail (use one lesser) !!!! HAVE TO BE IMPLEMENTED !!!!

    ## splsdrcox.modelPerComponent returns a list of models (per each component computed till max)
    ## so in this case we will have:
    ## [[max.comp]][[rep]][[fold]][[others_components_1:max.comp]]

    t1 <- Sys.time()
    if(PARALLEL){
      n_cores <- max(future::availableCores() - 1, 1)

      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(lst_inputs), n_cores))
      }else{
        future::plan("multisession", workers = min(length(lst_inputs), n_cores))
      }

      if(method==pkg.env$splsicox){
        lst_all_models <- furrr::future_map(lst_inputs, ~splsicox(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                  Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                  n.comp = .$comp, spv_penalty = eta.list[[.$eta_index]],
                                                                  x.center = x.center, x.scale = x.scale,
                                                                  y.center = y.center, y.scale = y.scale,
                                                                  remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                  remove_non_significant = remove_non_significant, tol = tol, alpha = alpha,
                                                                  MIN_EPV = MIN_EPV, returnData = returnData, verbose = verbose), .options = furrr_options(seed = TRUE))
      }else if(method==pkg.env$sb.splsicox){
        lst_all_models <- furrr::future_map(lst_inputs, ~sb.splsicox(X = lst_X_train[[.$run]][[.$fold]],
                                                                     Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                     n.comp = .$comp, spv_penalty = eta.list[[.$eta_index]],
                                                                     x.center = x.center, x.scale = x.scale,
                                                                     y.center = y.center, y.scale = y.scale,
                                                                     remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                     remove_non_significant = remove_non_significant, tol = tol, alpha = alpha,
                                                                     MIN_EPV = MIN_EPV, returnData = returnData), .options = furrr_options(seed = TRUE))
      }else if(method==pkg.env$splsdrcox){
        lst_all_models <- furrr::future_map(lst_inputs, ~splsdrcox(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                                  Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                  n.comp = .$comp, eta = eta.list[[.$eta_index]],
                                                                  x.center = x.center, x.scale = x.scale,
                                                                  y.center = y.center, y.scale = y.scale,
                                                                  remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                  remove_non_significant = remove_non_significant,
                                                                  verbose = verbose, alpha = alpha, tol = tol,
                                                                  MIN_EPV = MIN_EPV, returnData = returnData), .options = furrr_options(seed = TRUE))
      }else if(method==pkg.env$sb.splsdrcox){
        lst_all_models <- furrr::future_map(lst_inputs, ~sb.splsdrcox(X = lst_X_train[[.$run]][[.$fold]],
                                                                     Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                                  n.comp = .$comp, eta = eta.list[[.$eta_index]],
                                                                  x.center = x.center, x.scale = x.scale,
                                                                  y.center = y.center, y.scale = y.scale,
                                                                  remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                                  remove_non_significant = remove_non_significant,
                                                                  verbose = verbose, alpha = alpha, tol = tol,
                                                                  MIN_EPV = MIN_EPV, returnData = FreturnData), .options = furrr_options(seed = TRUE))
      }

      future::plan("sequential")

    }else{

      if(method==pkg.env$splsicox){
        lst_all_models <- purrr::map(lst_inputs, ~splsicox(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                           Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                           n.comp = .$comp, spv_penalty = eta.list[[.$eta_index]],
                                                           x.center = x.center, x.scale = x.scale,
                                                           y.center = y.center, y.scale = y.scale,
                                                           MIN_EPV = MIN_EPV, remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                           remove_non_significant = remove_non_significant, tol = tol, alpha = alpha,
                                                           returnData = returnData))
      }else if(method==pkg.env$sb.splsicox){
        lst_all_models <- purrr::map(lst_inputs, ~sb.splsicox(X = lst_X_train[[.$run]][[.$fold]],
                                                              Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                              n.comp = .$comp, spv_penalty = eta.list[[.$eta_index]],
                                                              x.center = x.center, x.scale = x.scale,
                                                              y.center = y.center, y.scale = y.scale,
                                                              remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                              remove_non_significant = remove_non_significant, tol = tol,
                                                              MIN_EPV = MIN_EPV, returnData = returnData), .options = furrr_options(seed = TRUE))
      }else if(method==pkg.env$splsdrcox){
        lst_all_models <- purrr::map(lst_inputs, ~splsdrcox(X = data.matrix(lst_X_train[[.$run]][[.$fold]]),
                                                           Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                           n.comp = .$comp, eta = eta.list[[.$eta_index]],
                                                           x.center = x.center, x.scale = x.scale,
                                                           y.center = y.center, y.scale = y.scale,
                                                           remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                           remove_non_significant = remove_non_significant,
                                                           verbose = verbose, alpha = alpha, tol = tol,
                                                           MIN_EPV = MIN_EPV, returnData = returnData))
      }else if(method==pkg.env$sb.splsdrcox){
        lst_all_models <- purrr::map(lst_inputs, ~sb.splsdrcox(X = lst_X_train[[.$run]][[.$fold]],
                                                              Y = data.matrix(lst_Y_train[[.$run]][[.$fold]]),
                                                              n.comp = .$comp, eta = eta.list[[.$eta_index]],
                                                              x.center = x.center, x.scale = x.scale,
                                                              y.center = y.center, y.scale = y.scale,
                                                              remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                              remove_non_significant = remove_non_significant,
                                                              verbose = verbose, alpha = alpha, tol = tol,
                                                              MIN_EPV = MIN_EPV, returnData = returnData))
      }

    }
    t2 <- Sys.time()
    t2-t1

    rm(lst_inputs)

    ## We need to return a list of lists: COMP->ETA->REP->FOLDS
    comp_model_lst <- list()
    cont_problems = 0
    for(c in max.ncomp){
      eta_model_lst <- list()
      for(e in 1:length(eta.list)){
        run_model_lst <- list()
        for(r in 1:n_run){
          fold_model_lst <- list()
          for(f in 1:k_folds){
            name <- paste0(c, "_", e, "_", r, "_", f)
            if(all(is.na(lst_all_models[[name]])) || all(is.null(lst_all_models[[name]]$survival_model$fit))){
              cont_problems = cont_problems + 1
            }
            fold_model_lst[[f]] = lst_all_models[[name]]
          }
          names(fold_model_lst) <- paste0("fold_", 1:k_folds)
          run_model_lst[[r]] <- fold_model_lst
        }
        names(run_model_lst) <- paste0("run_", 1:n_run)
        eta_model_lst[[e]] <- run_model_lst
      }
      if(method %in% c(pkg.env$splsicox, pkg.env$sb.splsicox)){
        names(eta_model_lst) <- paste0("eta_", eta.list) #eta de momento
      }else if(method %in% c(pkg.env$splsdrcox, pkg.env$sb.splsdrcox)){
        names(eta_model_lst) <- paste0("eta_", eta.list)
      }

      comp_model_lst[[c]] <- eta_model_lst
    }
    names(comp_model_lst) <- paste0("comp_", 1:max.ncomp)

    info <- "No info"

    ## Before compute all intermediate models, check if all problems
    if(cont_problems == n_run * k_folds * length(eta.list)){
      if(verbose){
        message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))
      }
      return(NULL)
    }

    ## We need to fill models from 1:max.ncomp (it uses max.ncomp to fill the others, if it is NULL, method fail!!!)
    for(comp in 1:(max.ncomp-1)){
      eta_model_lst <- list()
      for(e in 1:length(eta.list)){
        run_model_lst <- list()
        for(r in 1:n_run){
          fold_model_lst <- list()
          for(f in 1:k_folds){
            if(method %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
              fold_model <- getSubModel.mb(model = comp_model_lst[[max.ncomp]][[e]][[r]][[f]], comp = comp, remove_non_significant = remove_non_significant)
            }else{
              fold_model <- getSubModel(model = comp_model_lst[[max.ncomp]][[e]][[r]][[f]], comp = comp, remove_non_significant = remove_non_significant)
            }
            fold_model_lst[[f]] <- fold_model
          }
          names(fold_model_lst) <- paste0("fold_", 1:k_folds)
          run_model_lst[[r]] <- fold_model_lst
        }
        names(run_model_lst) <- paste0("run_", 1:n_run)
        eta_model_lst[[e]] <- run_model_lst
      }
      if(method %in% c(pkg.env$splsicox, pkg.env$sb.splsicox)){
        names(eta_model_lst) <- paste0("eta_", eta.list) #eta de momento
      }else if(method %in% c(pkg.env$splsdrcox, pkg.env$sb.splsdrcox)){
        names(eta_model_lst) <- paste0("eta_", eta.list)
      }
      comp_model_lst[[comp]] <- eta_model_lst
    }

    #### ### ### ### ### ### #
    # UPDATING GLOBALS SIZE #
    #### ### ### ### ### ### #
    MB = 500
    bytes = MB*1024^2
    options(future.globals.maxSize = bytes)

    return(list(comp_model_lst = comp_model_lst, info = info))

  }
}

#### ### ### #
# EVALUATION #
#### ### ### #

getBetasFromCOXNET <- function(fit, bestFitModelIndex){
  if(isa(fit,"cv.glmnet")){
    aux = fit$glmnet.fit
  }
  if(all(isa(fit,c("coxnet", "glmnet")))){
    aux = fit
  }
  betas <- aux$beta[,bestFitModelIndex]
  return(betas)
}

getLinealPredictors <- function(cox, data, bestModel = NULL, center=T, scale=F){

  if(all(is.na(cox)) | all(is.null(cox))){
    return(rep(NA, nrow(data)))
  }

  #Linear predictors by hand [los resultados de la curva roc son independientes de si se centran o no los predictores]
  c <- class(cox)
  if(length(c)>1)
    c <- c[1]

  if(c=="coxph"){
    lst_cn <- attr(cox$terms, which = "term.labels")
    lst_cn <- unlist(lapply(lst_cn, function(x){gsub('`',"",x)})) #remove forbiddent characters {`}
    d <- as.data.frame(data[,lst_cn,drop=F])

    lp <- predict(object = cox, newdata = d,
                  type="lp",
                  se.fit=T, na.action=na.pass, reference=c("strata"))

    # lp <- tryCatch({ #si warning, stop
    #   apply(data[,!colnames(data) %in% c("time","event","status"),drop=F], 1, function(x) sum(t((x-cox$means)*cox$coefficients))) #X*betas but they are centered
    #
    #   #save the error
    #   warning=function(cond) {
    #     # Choose a return value in case of error
    #     stop("The length for data and cox variables are not the same. Update the input matrix.")
    #   }
    # })

    #lp <- apply(data[,!colnames(data) %in% c("time","event","status"),drop=F], 1, function(x) sum(t((x-cox$means)*cox$coefficients))) #X*betas but they are centered
  }else if(c=="coxnet"){
    betas <- getBetasFromCOXNET(cox, bestModel)
    lp <- apply(data[,!colnames(data) %in% c("time","event","status"),drop=F], 1, function(x) sum(t(x*betas)))
  }

  if(center & c!="coxph"){
    lp.center <- scale(lp, center = center, scale = scale)
    lp.center <- as.numeric(lp.center)
    names(lp.center) <- names(lp)
    return(lp.center)
  }else{
    return(lp)
  }
}

getVectorCuts <- function(vector, cut_points, verbose = F){

  times = NULL
  if(length(vector)==1){
    return(vector)#if just one value, return that value
  }else{
    times <- vector
  }

  if(cut_points <= 1){
    if(verbose){
      message("Cutpoints must be greater than 1 and less than the vector length.")
    }
    cut_points = 2
  }else if(cut_points>length(times)){
    if(verbose){
      message("Cutpoints must be greater than 1 and less than the vector length.")
    }
    cut_points = length(times)
  }

  diff <- round((length(times)-1) / (cut_points-1))
  res <- min(times) + diff*0:(cut_points-1)
  res[1] = min(times)
  res[length(res)] = max(times)

  return(res)

}

# getVectorOfTime <- function(Y, max_time_points = 15, verbose = F){
#   times = NULL
#   if(is.null(times)){
#     #a maximum of max_time_points time points will be selected
#     if(verbose){
#       message(paste0("\t A maximum of ", max_time_points, " time points will be used."))
#     }
#
#     max_points = max_time_points
#     times <- 1:as.numeric(max(Y[Y[,"event"]==1,"time"]))
#     times <- times[times %% (length(times) / max_points) < 1]
#     if(verbose){
#       message(paste0("\t Time point selected are: ", paste0(times, collapse = ", "), ".\n"))
#       message("\t For specific time points, use the argument 'times'.\n\n")
#     }
#
#   }else{
#     if(length(times)>20){
#       if(verbose){
#         message("\t More than 20 time points have been selected. CV processes could take a long time... Between 5-15 are recommended.\n")
#       }
#     }
#
#     if(max(times)>max(Y[,"time"])){
#       times = times[times <= max(Y[,"time"])]
#       if(verbose){
#         message(paste0("\t It has been selected a vector of times greater than the maximum Y time (censored) event. Time vector updated to: ", paste0(times, collapse = ", "),"\n"))
#       }
#     }
#   }
#
#   return(times)
# }

checkLibraryEvaluator <- function(pred.method){
  #Check evaluator installed:
  if(!pred.method == pkg.env$AUC_cenROC){
    if(!pred.method %in% c(pkg.env$AUC_survivalROC, pkg.env$AUC_nsROC, pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I, pkg.env$AUC_risksetROC)){
      stop(paste0("A non-valid method has been selected. Select one of: ", paste0(c(pkg.env$AUC_cenROC, pkg.env$AUC_survivalROC, pkg.env$AUC_nsROC, pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I, pkg.env$AUC_risksetROC), collapse = ", ")))
    }else if(pred.method == pkg.env$AUC_survivalROC & !requireNamespace("survivalROC", quietly = TRUE)){
      stop("Library 'survivalROC' is required to evaluate with the selected method.")
    }else if(pred.method == pkg.env$AUC_nsROC & !requireNamespace("nsROC", quietly = TRUE)){
      stop("Library 'nsROC' is required to evaluate with the selected method.")
    }else if(pred.method %in% c(pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I) & !requireNamespace("smoothROCtime", quietly = TRUE)){
      stop("Library 'smoothROCtime' is required to evaluate with the selected method.")
    }else if(pred.method == pkg.env$AUC_risksetROC & !requireNamespace("risksetROC", quietly = TRUE)){
      stop("Library 'risksetROC' is required to evaluate with the selected method.")
    }
  }
}

#### ### ### ### ### ### ### ### ### ### ### ### ### ## #
## Eval all models by the pred.methods the user defined #
#### ### ### ### ### ### ### ### ### ### ### ### ### ## #

checkAtLeastTwoEvents <- function(X_test, Y_test){

  if(!isa(X_test, "list")){
    rn_X <- rownames(X_test)

    if(!all(rn_X %in% rownames(Y_test))){
      stop("Rownames of X_test must be in Y_test")
    }

    sub_Y <- Y_test[rn_X,,drop=F]

    if(sum(sub_Y[,"event"])<2){
      stop("To evaluate a model, at least two events are mandatory in TEST data set.")
    }

  }else{ # MULTIOMIC APPROACH
    for(block in names(X_test)){
      rn_X <- rownames(X_test[[block]])

      if(!all(rn_X %in% rownames(Y_test))){
        stop(paste0("Rownames of X_test and block ", block," must be in Y_test"))
      }

      sub_Y <- Y_test[rn_X,,drop=F]

      if(sum(sub_Y[,"event"])<2){
        stop("To evaluate a model, at least two events are mandatory in TEST data set.")
      }
    }
  }
}

#' eval_HDcox_models
#' @description Evaluate multiple HDcox models simultaneously.
#'
#' @param lst_models List of HDcox models. Each object of the list must be named.
#' @param X_test Numeric matrix or data.frame. Explanatory variables for test data (raw format). Qualitative variables must be transform into binary variables.
#' @param Y_test Numeric matrix or data.frame. Response variables for test data. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance. Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I" (default: "cenROC").
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model (default: 15).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param progress_bar Logical. If progress_bar = TRUE, progress bar is shown (default = TRUE).
#'
#' @export

## Eval all models by the pred.methods the user defined
eval_HDcox_models <- function(lst_models, X_test, Y_test, pred.method, pred.attr = "mean", times = NULL, PARALLEL = F, max_time_points = 15, verbose = F, progress_bar = T){

  #Check at least two events in total
  checkAtLeastTwoEvents(X_test, Y_test)

  #Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  if(all(is.null(lst_models))){
    stop("List of model is NULL")
  }

  if(!pred.attr %in% c("mean", "median")){
    stop("pred.attr parameter must be one of: 'mean' or 'median'")
  }

  t1 <- Sys.time()

  if(verbose){
    message(paste0("\nEvaluating with ", pred.method, "...\n"))
  }

  #TEST DATA
  if(is.null(times)){
    times <- getTimesVector(Y_test, max_time_points)
  }

  #MULTIBLOCK
  if(isa(X_test, "list")){
    X_test_ori  <- purrr::map(X_test, ~data.matrix(.))
    Y_test  <- Y_test
  }else{
    X_test_ori  <- data.matrix(X_test)
    Y_test  <- Y_test
  }

  #### ### ### #
  # Evaluation #
  #### ### ### #

  #models not NAs
  names_lst_models <- unlist(lapply(lst_models, function(x){!all(is.na(x))}))
  names_lst_models <- names(names_lst_models)[names_lst_models==T]

  lst_models <- lst_models[names_lst_models]

  if(progress_bar){
    total_models <- length(lst_models)
    pb_text <- "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated remaining time: :eta]"
    pb <- progress::progress_bar$new(format = pb_text,
                                     total = total_models,
                                     complete = "=",   # Caracteres de las iteraciones finalizadas
                                     incomplete = "-", # Caracteres de las iteraciones no finalizadas
                                     current = ">",    # Caracter actual
                                     clear = FALSE,    # Si TRUE, borra la barra cuando termine
                                     width = 100)      # Ancho de la barra de progreso
    pb$tick(0)
  }

  lst_eval <- list()
  #For smoothROC parallel per method
  if(PARALLEL & pred.method %in% c(pkg.env$AUC_smoothROCtime_C, pkg.env$AUC_smoothROCtime_I, pkg.env$AUC_survivalROC, pkg.env$AUC_nsROC)){

    n_cores <- max(future::availableCores() - 1, 1)

    if(.Platform$OS.type == "unix") {
      future::plan("multicore", workers = min(length(lst_models), n_cores))
    }else{
      future::plan("multisession", workers = min(length(lst_models), n_cores))
    }

    lst_eval <- furrr::future_map(lst_models, ~evaluation_list_HDcox(model = ., X_test = X_test, Y_test = Y_test, pred.method = pred.method, pred.attr = pred.attr, times = times, PARALLEL = F, verbose = verbose, progress_bar = progress_bar))
    future::plan("sequential")
  }else{
    lst_eval <- purrr::map(lst_models, ~evaluation_list_HDcox(model = ., X_test = X_test, Y_test = Y_test, pred.method = pred.method, pred.attr = pred.attr, times = times, PARALLEL = F, verbose = verbose, progress_bar = progress_bar))
  }

  names(lst_eval) <- names(lst_models)

  lst_AUC <- list()
  df <- NULL
  for(m in names(lst_eval)){
    lst_AUC[[m]] <- lst_eval[[m]]$lst_AUC_values

    #if AUC_values is NA, we cannot access to lst_AUC_values$AUC.vector
    if(!all(is.na(lst_eval[[m]]$lst_AUC_values))){
      df <- rbind(df, c(m, lst_eval[[m]]$model_time, lst_eval[[m]]$comp.time, lst_eval[[m]]$aic.cox, lst_eval[[m]]$c_index.cox, lst_eval[[m]]$lst_AUC_values$AUC.vector))
    }else{
      df <- rbind(df, c(m, lst_eval[[m]]$model_time, lst_eval[[m]]$comp.time, lst_eval[[m]]$aic.cox, lst_eval[[m]]$c_index.cox, rep(NA, length(times))))
    }

    df <- as.data.frame(df)

  }

  final_times <- times #all the same

  if(is.null(df) || (ncol(df) != 5+length(final_times) & all(is.na(df[,2])))){
    df <- as.data.frame(matrix(data = NA, nrow = length(names_lst_models), ncol = 5+length(final_times)))
    df[,1] <- names_lst_models
  }
  # }else if(ncol(df) != 1+length(final_times) & all(is.na(df[,2]))){
  #   df_na <- as.data.frame(matrix(data = NA, nrow = nrow(df), ncol = length(final_times)))
  #   df <- cbind(df[,1,drop=F], df_na)
  # }

  if(all(is.na(df[,2])) & ncol(df) < (5+length(final_times))){
    colnames(df) <- c("method", "training.time","evaluating.time", "AIC", "c.index", "AUC")
    df <- as.data.frame(df)
    new_df <- tidyr::pivot_longer(df, cols = starts_with("time_"), names_to = "time", values_to = "AUC",)
    new_df$time <- factor(new_df$time, levels = unique(new_df$time))
  }else{
    colnames(df) <- c("method", "training.time","evaluating.time", "AIC", "c.index", paste0("time_",final_times))
    df <- as.data.frame(df)
    df$method <- factor(df$method, levels = unique(df$method))
    df[,!colnames(df) %in% "method"] <- apply(df[,!colnames(df) %in% "method"], 2, as.numeric)

    new_df <- tidyr::pivot_longer(df, cols = starts_with("time_"), names_to = "time", values_to = "AUC",)
    new_df$time <- factor(new_df$time, levels = unique(new_df$time))
  }

  #Look for problems !!!!
  #Could be cause by "sample is too sparse to find TD for time 3.25 and method 'cenROC'"
  #Means SAMPLE have only two values for linear predictors (only one variable and its binary), so AUC cannot be compute for that model.

  for(m in unique(new_df$method)){
    if(all(is.na(new_df[new_df$method==m,]$AUC))){
      message(paste0("Problems computing AUC metric for '", pred.method, "' evaluator and model '", m, "'"))
    }
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  if(verbose){
    message(paste0("\nTime for ", pred.method, ": ", as.character(round(time, 5))))
  }

  return(evaluation_HDcox_class(list(df = new_df, lst_AUC = lst_AUC, time = time)))
}

evaluation_list_HDcox <- function(model, X_test, Y_test, pred.method, pred.attr = "mean", times = NULL, PARALLEL = F, verbose = F, progress_bar = F){

  t3 <- Sys.time()

  #atomic vector
  if(all(is.na(model))){
    return(list(model_time = NA, comp.time = NA, aic.cox = NA, c_index.cox = NA, lst_AUC_values = NA))
  }

  #NULL in HDcox object (NA no anymore)
  if(all(is.null(model$survival_model))){
    return(list(model_time = NA, comp.time = NA, aic.cox = NA, c_index.cox = NA, lst_AUC_values = NA))
  }

  cox <- model$survival_model$fit
  aic.cox <- stats::extractAIC(cox, k=2)[2] #k=2 <- AIC, [2] AIC Value
  c_index.cox <- survival::concordance(cox)$concordance

  #linear predictors
  if(isa(X_test, "list") & !attr(model, "model") %in% pkg.env$multiblock_methods){ #mix between multiblock and all PLS - Special case
    which_block = purrr::map(names(X_test), ~length(grep(., m, fixed = F))>0)
    names(which_block) <- names(X_test)
    X_test_mod <- predict.HDcox(object = model, newdata = X_test[[names(which_block)[which_block==T]]])
  }else{
    X_test_mod <- predict.HDcox(object = model, newdata = X_test) #all multiblock or all PLS - Ok
  }

  lp <- getLinealPredictors(cox = cox, data = X_test_mod)
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Y_test, times = times, bestModel = NULL, eval = pred.attr, method = pred.method, PARALLEL = PARALLEL, verbose = verbose)
  #lst_AUC[[m]] <- lst_AUC_values

  t4 <- Sys.time()
  comp.time <- difftime(t4,t3,units = "mins")

  #df <- rbind(df, c(m, model$time, comp.time, aic.cox, c_index.cox, lst_AUC_values$AUC.vector))

  return(list(model_time = model$time, comp.time = comp.time, aic.cox = aic.cox, c_index.cox = c_index.cox, lst_AUC_values = lst_AUC_values))
}

evaluation_HDcox_class = function(object, ...) {
  model = structure(object, class = pkg.env$model_class,
                    model = pkg.env$eval_class)
  return(model)
}

#' cox.prediction
#'
#' @param model HDcox model.
#' @param new_data Numeric matrix or data.frame. New explanatory variables (raw data). Qualitative variables must be transform into binary variables.
#' @param time Numeric. Time point where the AUC will be evaluated (default: NULL).
#' @param type Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").
#' @param method Character. Prediction method. It can be compute by using the cox model "cox" or by using W.star "W.star" (default: "cox"). (not implemented for MB approaches)!!!
#'
#' @return Return the lp or other metric for the patient
#' @export

cox.prediction <- function(model, new_data, time = NULL, type = "lp", method = "cox"){
  #could be obtain by predicting scores or by computing W*

  if(method == "cox"){
    scores <- predict.HDcox(object = model, newdata = new_data) #X must be original X data

    if(type %in% c("expected", "survival") & is.null(time)){
      stop("For survivial or expected prediction, you must provided a specific time of study.")
    }

    if(type %in% c("expected", "survival")){
      df <- as.data.frame(cbind(scores, data.frame(time = time, event = 0)))
    }else{
      df <- as.data.frame(scores)
    }

    if(all(is.null(model$survival_model))){
      stop("Survival model not found.")
    }

    pred.value <- predict(object = model$survival_model$fit, newdata = df, type = type)

    return(pred.value)

  }else if(method == "W.star"){

    beta <- as.matrix(model$survival_model$fit$coefficients)
    W.star <- model$X$W.star
    coefficients <- W.star %*% beta

    #norm patient
    if(!is.null(model$X$x.mean) & !is.null(model$X$x.sd)){
      norm_patient <- scale(new_data, center = model$X$x.mean, scale = model$X$x.sd)
    }else if(!is.null(model$X$x.mean)){
      norm_patient <- scale(new_data, center = model$X$x.mean, scale = F)
    }else if(!is.null(model$X$x.sd)){
      norm_patient <- scale(new_data, center = F, scale = model$X$x.sd)
    }else{
      norm_patient <- new_data
    }

    pred.value <- norm_patient[,rownames(coefficients)] %*% coefficients
    pred.value <- pred.value[[1]]
    names(pred.value) <- rownames(new_data)

  }

}

getBestVector <- function(Xh, DR_coxph = NULL, Yh, n.comp, max.iter, vector, MIN_AUC_INCREASE, MIN_NVAR = 10, MAX_NVAR = 10000, cut_points = 5, EVAL_METHOD = "AUC", EVAL_EVALUATOR = "cenROC", PARALLEL = F, mode = "spls", times = NULL, max_time_points = 15, verbose = F){

  if(!mode %in% c("spls", "splsda")){
    stop("Mode must be one of: 'spls' or 'splsda'")
  }

  if(!EVAL_METHOD %in% c("AUC", "BRIER", "c_index")){
    stop("Evaluation method must be one of: 'AUC', 'BRIER' or 'c_index'")
  }

  max_ncol <- ncol(Xh)

  if(is.null(vector)){
    vector <- getVectorCuts(vector = c(min(MIN_NVAR, max_ncol):min(max_ncol, MAX_NVAR)), cut_points = cut_points, verbose = verbose)
  }else{
    #check if each value is less than the ncol
    if(is.numeric(vector)){
      vector <- vector[vector<=ncol(Xh)]
      message(paste0("The initial vector is: ", paste0(vector, collapse = ", "), "\n"))
    }else{
      message("Your vector must be a numeric vector. A starting vector is created:")
      vector <- getVectorCuts(vector = c(min(MIN_NVAR, max_ncol):min(max_ncol, MAX_NVAR)), cut_points = cut_points, verbose = verbose)
      message(paste0("The initial vector is: ", paste0(vector, collapse = ", ")))
    }
  }

  if(verbose){
    message(paste0("Original vector: "))
    message(paste0("Values: ", paste0(vector, collapse = ", "), "\n"))
  }

  count = 1
  var_exp = NULL

  #if n_col is minimum than MIN_NVAR, values could be the same, so delete duplicates
  vector <- unique(vector)

  if(PARALLEL){
    n_cores <- future::availableCores() - 1

    if(.Platform$OS.type == "unix") {
      future::plan("multicore", workers = min(length(vector), n_cores))
    }else{
      future::plan("multisession", workers = min(length(vector), n_cores))
    }

    t1 <- Sys.time()
    if(mode %in% "spls"){
      lst_cox_value <- furrr::future_map(vector, ~getCIndex_AUC_CoxModel_spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
    }else{
      lst_cox_value <- furrr::future_map(vector, ~getCIndex_AUC_CoxModel_splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
    }
    t2 <- Sys.time()
    future::plan("sequential")
  }else{
    t1 <- Sys.time()
    if(mode %in% "spls"){
      lst_cox_value <- purrr::map(vector, ~getCIndex_AUC_CoxModel_spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
    }else{
      lst_cox_value <- purrr::map(vector, ~getCIndex_AUC_CoxModel_splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, EVAL_EVALUATOR = EVAL_EVALUATOR, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
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

  rownames(df_cox_value) <- vector

  index <- which.max(df_cox_value) #MAX CONCORDANCE/AUC

  #Select best vector
  keepX <- vector[[index]]
  FLAG = T
  cont = 0

  best_c_index <- df_cox_value[index]
  best_keepX <- keepX

  if(verbose){
    message(paste0("First selection: \n"), paste0(paste0("Value ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), "Pred. Value: ", round(best_c_index, 4), "\n")
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

    new_vector <- NULL

    #before_vector - always go two sides
    for(b in best_keepX){
      aux <- b
      index <- which(aux_vector < aux)
      if(length(index)==0){
        value = aux_vector[which(aux_vector == aux)]
      }else{
        value = round(mean(c(aux, aux_vector[index][[length(aux_vector[index])]]))) #last smaller
      }
      new_vector <- unique(c(new_vector, aux, value))
    }

    #next_vector - always go two sides
    for(b in best_keepX){
      aux <- best_keepX
      index <- which(aux_vector > aux)
      if(length(index)==0){
        value = aux_vector[which(aux_vector == aux)]
      }else{
        value = round(mean(c(aux, aux_vector[index][[1]]))) #first greater
      }
      new_vector <- unique(c(new_vector, aux, value))
    }

    new_vector <- new_vector[!new_vector %in% names(p_val)[!is.na(p_val)]]

    if(verbose){
      message(paste0("Testing: \n"), paste0("Value ", names(best_keepX), ": ", new_vector, "\n"))
    }

    new_vector <- unlist(new_vector)

    all_comb <- expand.grid(new_vector)

    ### update all vector
    aux_vector <- unique(c(aux_vector, new_vector))
    aux_vector <- aux_vector[order(aux_vector)]
    names(aux_vector) <- aux_vector

    p_val <- c(p_val, rep(NA, length(new_vector)))
    names(p_val)[(length(p_val)-length(new_vector)+1):length(p_val)] <- new_vector
    p_val <- p_val[order(as.numeric(names(p_val)))]

    ### OTHER KEEP_VECTOR
    vector_aux <- new_vector
    vector_aux <- vector_aux[order(vector_aux)]

    ## Compute next c_index
    if(PARALLEL){
      n_cores <- future::availableCores() - 1

      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = min(length(vector_aux), n_cores))
      }else{
        future::plan("multisession", workers = min(length(vector_aux), n_cores))
      }

      t1 <- Sys.time()
      if(mode %in% "spls"){
        lst_cox_value <- furrr::future_map(vector_aux, ~getCIndex_AUC_CoxModel_spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
      }else{
        lst_cox_value <- furrr::future_map(vector_aux, ~getCIndex_AUC_CoxModel_splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
      }
      t2 <- Sys.time()
      future::plan("sequential")
    }else{
      t1 <- Sys.time()
      if(mode %in% "spls"){
        lst_cox_value <- purrr::map(vector_aux, ~getCIndex_AUC_CoxModel_spls(Xh = Xh, DR_coxph_ori = DR_coxph, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
      }else{
        lst_cox_value <- purrr::map(vector_aux, ~getCIndex_AUC_CoxModel_splsda(Xh = Xh, Yh = Yh, n.comp = n.comp, keepX = ., scale = F, near.zero.var = F, max.iter = max.iter, times = times, max_time_points = max_time_points), .progress = F)
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

    rownames(df_cox_value_aux) <- vector_aux
    #index <- which.max(rowSums(df_cox_value_aux)) #MAX VAR_MEDIA
    #index <- which.max(df_cox_value_aux[,"Y"]) #MAX Y?
    index <- which.max(df_cox_value_aux) #MAX CONCORDANCE
    best_c_index_aux <- df_cox_value_aux[index]

    p_val[rownames(df_cox_value_aux)] <- df_cox_value_aux

    if(best_c_index >= best_c_index_aux | best_c_index_aux-best_c_index <= MIN_AUC_INCREASE){
      FLAG = F
      if(verbose){
        message(paste0("End: \n"), paste0(paste0("Value ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value: ", round(best_c_index, 4), "\n"))
      }
    }else{
      best_c_index <- best_c_index_aux
      best_keepX <- vector_aux[[index]]
      if(verbose){
        message(paste0("New Vector: \n"), paste0(paste0("Value ", names(best_keepX), ": ", unlist(purrr::map(best_keepX, ~unique(.)))), "\n"), paste0("Pred. Value: ", round(best_c_index_aux, 4), "n"))
      }
    }
  }

  keepX <- best_keepX
  return(list(best.keepX = keepX, p_val = p_val))
}

getCIndex_AUC_CoxModel_spls <- function(Xh, DR_coxph_ori, Yh, n.comp, keepX, scale = F, near.zero.var = F, EVAL_EVALUATOR = "cenROC", max.iter = 200, verbose = F, times = NULL, max_time_points = 15){
  model <- mixOmics::spls(X = Xh, Y = DR_coxph_ori, ncomp = n.comp, keepX = rep(keepX, n.comp), scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDR = model$variates

  d <- as.data.frame(tt_mbsplsDR[[1]][,,drop=F]) #X

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
                      model=T, x = T)
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

getCIndex_AUC_CoxModel_splsda <- function(Xh, Yh, n.comp, keepX, scale = F, near.zero.var = F, EVAL_EVALUATOR = "cenROC", max.iter = 100, verbose = F, times = NULL, max_time_points = 15){
  model <- mixOmics::splsda(X = Xh, Y = Yh[,"event"], ncomp = n.comp, keepX = rep(keepX, n.comp), scale = scale, near.zero.var = near.zero.var, max.iter = max.iter)
  tt_mbsplsDA = model$variates

  d <- as.data.frame(tt_mbsplsDA[[1]][,,drop=F]) #X

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
                      model=T, x = T)
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

  #C-index and AUC
  lp <- getLinealPredictors(cox = cox_model$fit, data = d)
  lst_AUC_values <- getAUC_from_LP_2.0(linear.predictors = lp, Y = Yh, times = times, bestModel = NULL, eval = "mean", method = EVAL_EVALUATOR, PARALLEL = FALSE, verbose = verbose)

  #BRIER
  #lst_BRIER_values <- survAUC_BRIER_LP(lp = lp$fit, Y = Yh, lp_new = lp$fit, Y_test = Yh)
  lst_BRIER_values <- SURVCOMP_BRIER_LP(lp_train = lp$fit, Y_train = Yh, lp_test = lp$fit, Y_test = Yh)

  return(list("c_index" = cox_model$fit$concordance["concordance"], "AUC" = lst_AUC_values$AUC, "BRIER" = lst_BRIER_values$ierror))
}
