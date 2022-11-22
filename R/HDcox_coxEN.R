#### ### ##
# METHODS #
#### ### ##

#' coxEN
#' @description Performs a cox elastic net model (based on glmnet R package).
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param EN.alpha Numeric. Elasticnet mixing parameter. EN.alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param returnData Logical. Return original and normalized X and Y matrices.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "coxEN". The class contains the following elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
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
#' \code{EN.alpha}: EN.alpha selected
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{selected_variables}: Variables selected to the cox model.
#'
#' \code{removed_variables}: Variables removed by EN penalty.
#'
#' \code{opt.lambda}: Optimal lambda used in EN.
#'
#' \code{convergence_issue}: Whether convergence issue is detected
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export

coxEN <- function(X, Y,
                  EN.alpha = 0.5,
                  x.center = TRUE, x.scale = FALSE,
                  y.center = FALSE, y.scale = FALSE,
                  remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                  remove_non_significant = F, alpha = 0.05,
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

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  checkY.colnames(Y)

  #### SCALING
  lst_scale <- XY.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  ####MAX PREDICTORS
  max_variables <- check.maxPredictors(X, Y, MIN_EPV, ncol(X))

  #### REMOVING event in time 0 patients
  if(any(Yh[,"time"]==0)){
    if(verbose){
      message("Some patients get the event at time 0. Those patients will be removed for the analysis")
    }
    pat_names <- rownames(Yh)[Yh[,"time"]==0]
    Yh <- Yh[!rownames(Yh) %in% pat_names,]
    Xh <- Xh[rownames(Xh) %in% rownames(Yh),]
  }

  #### MAX PREDICTORS
  max_n_predictors <- getMaxNPredictors(n.var = ncol(X), Y, MIN_EPV)

  #### INITIALISING VARIABLES
  best_cox <- NULL
  best_lambda <- NULL
  selected_variables <- NULL
  problem <- FALSE #convergence issues

  #A small detail in the Cox model: if death times are tied with censored times, we assume the censored times occurred
  #just before the death times in computing the Breslow approximation; if users prefer the usual convention of after,
  #they can add a small number to all censoring times to achieve this effect.
  for(t in unique(Yh[,"time"])){
    set_pat <- Yh[Yh[,"time"]==t,,drop=F]
    if(nrow(set_pat)>1 & length(unique(set_pat[,"event"]))>1){
      names_censored <- names(which(set_pat[,"event"]==0))
      Yh[names_censored,"time"] <- Yh[names_censored,"time"] + 0.0001
    }
  }

  EN_cox <- tryCatch(
    # Specifying expression
    expr = {
      glmnet::glmnet(x = Xh, y = survival::Surv(time = Yh[,"time"], event = Yh[,"event"]), family = "cox", alpha = EN.alpha, standardize = F, pmax = max_n_predictors, nlambda=200)
    },
    # Specifying error message
    error = function(e){
      message(paste0("coxEN: ", e))
      invisible(gc())
      return(NA)
    },
    warning = function(e){
      #message(paste0(e$message))
      if(verbose){
        message("Model probably has a convergence issue...")
      }
      defaultW <- getOption("warn")
      options(warn = -1)
      res <- glmnet::glmnet(x = Xh, y = survival::Surv(time = Yh[,"time"], event = Yh[,"event"]), family = "cox", EN.alpha = EN.alpha, standardize = F, pmax = max_n_predictors, nlambda=200)
      options(warn = defaultW)
      list(res = res, problem = T)
    }
  )

  #only if problems
  if(all(class(EN_cox) %in% c("list"))){
    problem = EN_cox$problem
    EN_cox = EN_cox$res
  }

  if(all(class(EN_cox) == c("coxnet", "glmnet"))){
    best_lambda <- EN_cox$lambda[which.max(EN_cox$dev.ratio)]

    coef.matrix <- as.matrix(coef(EN_cox, s = best_lambda))
    selected_variables <- rownames(coef.matrix)[which(coef.matrix != 0)]

    d <- as.data.frame(cbind(Xh[,selected_variables,drop=F], Yh)) #data

    best_cox <- tryCatch(
      # Specifying expression
      expr = {
        survival::coxph(formula = survival::Surv(time,event) ~ .,
                        data = d,
                        ties = "efron",
                        singular.ok = T,
                        robust = T,
                        nocenter = rep(1, ncol(Xh)),
                        model=T)
      },
      # Specifying error message
      error = function(e){
        message(paste0("COX: ", e))
        invisible(gc())
        return(NA)
      }
    )

    # best_cox <- tryCatch(
    #   expr = {
    #   # Specifying expression
    #   d <- as.data.frame(Xh[,selected_variables,drop=F]) #data
    #   cox(X = d, Y = Yh, FORCE = T,
    #       x.center = x.center, x.scale = x.scale,
    #       y.center = y.center, y.scale = y.scale)
    #   },
    #   # Specifying error message
    #   error = function(e){
    #     message(paste0("COX: ", e))
    #     invisible(gc())
    #     return(NA)
    #   }
    # )
  }

  if(is.na(best_cox) || (problem & all(best_cox$linear.predictors==0))){
    best_cox = NA
    best_lambda = NA
    selected_variables = NA

    func_call <- match.call()

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")

    survival_model <- NULL
    removed_variables <- NULL

    invisible(gc())
    return(coxEN_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                            Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                            survival_model = survival_model,
                            EN.alpha = EN.alpha,
                            #alpha = alpha,
                            call = func_call,
                            X_input = if(returnData) X_original else NA,
                            Y_input = if(returnData) Y_original else NA,
                            nzv = variablesDeleted,
                            selected_variables = selected_variables,
                            removed_variables = removed_variables,
                            opt.lambda = best_lambda,
                            convergence_issue = problem,
                            time = time)))
  }

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one
  removed_variables <- NULL

  if(remove_non_significant){
    p_val <- summary(best_cox)[[7]][,"Pr(>|z|)"]
    while(any(p_val>alpha)){
      to_remove <- names(which.max(p_val))
      to_remove <- deleteIllegalChars(to_remove)
      d <- d[,!colnames(d) %in% c(to_remove)]
      best_cox <- tryCatch(
        # Specifying expression
        expr = {
          survival::coxph(formula = survival::Surv(time,event) ~ .,
                          data = d,
                          ties = "efron",
                          singular.ok = T,
                          robust = T,
                          nocenter = rep(1, ncol(d)-ncol(Yh)),
                          model=T)
        },
        # Specifying error message
        error = function(e){
          message(paste0("coxEN: ", e))
          invisible(gc())
          return(NA)
        }
      )

      removed_variables <- c(removed_variables, to_remove)
      p_val <- summary(best_cox)[[7]][,"Pr(>|z|)"]
    }
  }

  if(class(best_cox)[[1]]=="coxph.null"){
    survival_model <- NULL
  }else if(class(best_cox)=="coxph"){
    survival_model <- getInfoCoxModel(best_cox)
  }else{
    survival_model <- NULL
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(coxEN_class(list(X = list("data" = if(returnData) Xh else NA, "x.mean" = xmeans, "x.sd" = xsds),
                        Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                        survival_model = survival_model,
                        opt.lambda = best_lambda,
                        EN.alpha = EN.alpha,
                        selected_variables = selected_variables,
                        call = func_call,
                        X_input = if(returnData) X_original else NA,
                        Y_input = if(returnData) Y_original else NA,
                        convergence_issue = problem,
                        alpha = alpha,
                        removed_variables = removed_variables,
                        nzv = variablesDeleted,
                        time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Cross validation COXEN
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param EN.alpha.list Numeric vector. Elasticnet mixing parameter values to test in cross validation. EN.alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
#' @param n_run Number. Number of runs for cross validation.
#' @param k_folds Number. Number of folds for cross validation.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant variables in final cox model will be removed until all variables are significant (forward selection).
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If MIN_EPV is not meet, the value will be change automatically.
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
#' @return Instance of class "HDcox" and model "cv.coxEN".
#' @export

cv.coxEN <- function(X, Y,
                     EN.alpha.list = seq(0,1,0.1),
                     n_run = 10, k_folds = 10,
                     x.center = TRUE, x.scale = FALSE,
                     y.center = FALSE, y.scale = FALSE,
                     remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                     remove_non_significant = F, alpha = 0.05,
                     max.variables = 15,
                     w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
                     MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                     pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
                     MIN_EPV = 5, return_models = F,
                     PARALLEL = F, verbose = F, seed = 123){

  t1 <- Sys.time()

  ############
  # WARNINGS #
  ############

  #### REQUIREMENTS
  checkY.colnames(Y)
  check.cv.weights(c(w_AIC, w_c.index, w_AUC))
  max.variables <- check.ncomp(X, max.variables)
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### MAX PREDICTORS
  max.variables <- check.maxPredictors(X, Y, MIN_EPV, max.variables, verbose = verbose)

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
  total_models <- k_folds * n_run * length(EN.alpha.list)

  comp_model_lst <- get_HDCOX_models2.0(method = pkg.env$coxEN,
                                     lst_X_train = lst_X_train, lst_Y_train = lst_Y_train, eta.list = NULL, max.ncomp = NULL,
                                     EN.alpha.list = EN.alpha.list, n_run = n_run, k_folds = k_folds,
                                     x.center = x.center, x.scale = x.scale,
                                     y.center = y.center, y.scale = y.scale,
                                     remove_non_significant = remove_non_significant,
                                     remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL,
                                     alpha = alpha,
                                     total_models = total_models, MIN_EPV = MIN_EPV,
                                     PARALLEL = PARALLEL, verbose = verbose)

  # comp_model_lst <- get_HDCOX_models(method = pkg.env$coxEN,
  #                                    lst_X_train = lst_X_train, lst_Y_train = lst_Y_train, eta.list = NULL,
  #                                    EN.alpha.list = EN.alpha.list, n_run = n_run, k_folds = k_folds,
  #                                    x.center = x.center, x.scale = x.scale,
  #                                    y.center = y.center, y.scale = y.scale,
  #                                    remove_non_significant = remove_non_significant,
  #                                    remove_near_zero_variance = remove_near_zero_variance,
  #                                    alpha = alpha,
  #                                    total_models = total_models, MIN_EPV = MIN_EPV)

  count_problems = 0
  count_null = 0
  count = 0
  for(l in names(comp_model_lst)){
    for(r in names(comp_model_lst[[l]])){
      for(f in names(comp_model_lst[[l]][[r]])){
        if(is.null(comp_model_lst[[l]][[r]][[f]]$convergence_issue)){
          count_null = count_null + 1
        }else if(comp_model_lst[[l]][[r]][[f]]$convergence_issue){
          count_problems = count_problems + 1
        }else{
          count = count + 1
        }
      }
    }
  }

  if(count_problems>0){
    message(paste0("There are a total of ", count_problems, " out of ", total_models, " models that could present convergence issues.\n"))
  }

  ##########################
  # BEST MODEL FOR CV DATA #
  ##########################
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst,
                                                    max.ncomp = EN.alpha.list, eta.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = F) #deletion at variable level for EN

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.coxEN_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.EN.alpha = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, time = time)))
    }else{
      return(cv.coxEN_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.EN.alpha = NULL, plot_AUC = NULL, plot_c_index = NULL, plot_AIC = NULL, time = time)))
    }
  }

  ##################
  # EVALUATING AUC #
  ##################
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_comp_flag <- NULL

  if(w_AUC!=0){
    total_models <- ifelse(!fast_mode, n_run * length(EN.alpha.list), k_folds * n_run * length(EN.alpha.list))

    #As we are measuring just one evaluator and one method - PARALLEL=F
    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     lst_X_test = lst_X_test, lst_Y_test = lst_Y_test,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = EN.alpha.list, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, total_models = total_models, method.train = pkg.env$coxEN, PARALLEL = F)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
    optimal_comp_index <- lst_df$optimal_comp_index
    optimal_comp_flag <- lst_df$optimal_comp_flag
  }else{
    df_results_evals_fold <- df_results_evals
  }

  ##############
  # BEST MODEL #
  ##############

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_c.index, w_AUC, colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC")

  flag_no_models = F

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==EN.alpha.list[[optimal_comp_index]],, drop=F][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{

    if(all(is.nan(df_results_evals_comp[,"score"]))){
      #message("No models computed. All of them have problems in convergency. Probably due to a high number of variables.")
      best_model_info <- df_results_evals_comp[1,,drop=F]
      flag_no_models = T
    }else{
      best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = T)),, drop=F][1,]
      best_model_info <- as.data.frame(best_model_info)
    }
  }

  ########
  # PLOT #
  ########
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, max.ncomp = EN.alpha.list,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "c_index", colname_AUC = "AUC", x.text = "EN.alpha Penalization")
  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp
  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_c_index <- lst_EVAL_PLOTS$ggp_c_index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  #########################
  # CHANGE 1s COLUMN_NAME #
  #########################
  colnames(best_model_info)[1] <- "eta"
  colnames(df_results_evals_fold)[1] <- "eta"
  colnames(df_results_evals_run)[1] <- "eta"
  colnames(df_results_evals_comp)[1] <- "eta"
  colnames(best_model_info)[1] <- "eta"
  colnames(best_model_info)[1] <- "eta"

  ##########
  # RETURN #
  ##########

  if(!flag_no_models){
    message(paste0("Best model obtained."))
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  if(return_models){
    return(cv.coxEN_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.EN.alpha = best_model_info$eta, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, time = time)))
  }else{
    return(cv.coxEN_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.EN.alpha = best_model_info$eta, plot_AUC = ggp_AUC, plot_c_index = ggp_c_index, plot_AIC = ggp_AIC, time = time)))
  }
}

### ## ##
# CLASS #
### ## ##

coxEN_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$coxEN)
  return(model)
}

cv.coxEN_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$cv.coxEN)
  return(model)
}
