## ---- include = FALSE---------------------------------------------------------
dpi = 125

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=dpi, 
  fig.retina=1, 
  fig.width=1440/dpi, #4:3 FHD
  fig.height=1080/dpi, 
  out.width="100%",
  crop = NULL,
  warning = T, 
  error = T
)

rm(dpi)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ConesaLab/Coxmos")

## ----setup, eval=FALSE, results = "hide"--------------------------------------
#  # load Coxmos
#  library(Coxmos)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("RColorConesa")
#  library(RColorConesa)

## ----load data----------------------------------------------------------------
# load dataset
data("X_proteomic")
data("Y_proteomic")

X <- X_proteomic
Y <- Y_proteomic

rm(X_proteomic, Y_proteomic)

## ----data dimensions, echo = FALSE--------------------------------------------
knitr::kable(X[1:5,1:5])
knitr::kable(Y[1:5,])

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(dim(X), col.names = "X")
knitr::kable(dim(Y), col.names = "Y")

## -----------------------------------------------------------------------------
ggp_density.event <- plot_events(Y = Y, 
                                 categories = c("Censored","Death"), #name for FALSE/0 (Censored) and TRUE/1 (Event)
                                 y.text = "Number of observations", 
                                 roundTo = 0.5, 
                                 max.breaks = 15)

## ----fig.small = T------------------------------------------------------------
ggp_density.event$plot

## -----------------------------------------------------------------------------
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = .7, # 70% train
                                          list = FALSE,
                                          times = 1)

## -----------------------------------------------------------------------------
X_train <- X[index_train,] #106x369
Y_train <- Y[index_train,]
X_test <- X[-index_train,] #44x369
Y_test <- Y[-index_train,]

## ---- eval=FALSE, message=T---------------------------------------------------
#  # classical approach
#  cox_model <- cox(X = X_train, Y = Y_train,
#                   x.center = T, x.scale = F,
#                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
#                   remove_non_significant = F, alpha = 0.05,
#                   MIN_EPV = 5, FORCE = F, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
EPV <- getEPV(X_train, Y_train)

## -----------------------------------------------------------------------------
EPV

## ---- warning=F, eval=FALSE---------------------------------------------------
#  # run cv.coxEN
#  cv.coxen_res <- cv.coxEN(X = X_train, Y = Y_train,
#                           EN.alpha.list = c(0.1, 0.5, 0.9),
#                           max.variables = ncol(X_train),
#                           n_run = 2, k_folds = 5,
#                           x.center = T, x.scale = F,
#                           remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                           remove_variance_at_fold_level = F,
#                           remove_non_significant = F, alpha = 0.05,
#                           w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
#                           MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                           pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                           MIN_EPV = 5, return_models = F,
#                           returnData = F,
#                           PARALLEL = F, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.coxen_res #0.18min.

## ---- warning=F---------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0.5, #cv.coxen_res$opt.EN.alpha,
                     max.variables = 8, #cv.coxen_res$opt.nvar,
                     x.center = T, x.scale = F,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = F, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## ---- warning=F---------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0.5, #cv.coxen_res$opt.EN.alpha
                     max.variables = 8, #cv.coxen_res$opt.nvar
                     x.center = T, x.scale = F,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = T, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## -----------------------------------------------------------------------------
coxen_model$nsv

## ---- warning=F, eval=FALSE---------------------------------------------------
#  # run cv.plsicox
#  cv.splsicox_res <- cv.splsicox(X = X_train, Y = Y_train,
#                                 max.ncomp = 2, spv_penalty.list = c(0.5, 0.9),
#                                 n_run = 2, k_folds = 5,
#                                 x.center = T, x.scale = F,
#                                 remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                 remove_variance_at_fold_level = F,
#                                 remove_non_significant_models = F, alpha = 0.05,
#                                 w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL, max_time_points = 15,
#                                 MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                 pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                 MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                 PARALLEL = F, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsicox_res #1.13min.

## ---- fig.small=T, warning=F, eval=FALSE--------------------------------------
#  # plot cv.plsicox
#  cv.splsicox_res$plot_AUC

## -----------------------------------------------------------------------------
splsicox_model <- splsicox(X = X_train, Y = Y_train, 
                           n.comp = 1, #cv.splsicox_res$opt.comp, 
                           spv_penalty = 0.9, #cv.splsicox_res$opt.spv_penalty,
                           x.center = T, x.scale = F,
                           remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                           remove_non_significant = T,
                           MIN_EPV = 5, returnData = T, verbose = F)

splsicox_model

## ---- warning==FALSE, eval=FALSE----------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdrcox_res <- cv.splsdrcox(X = X_train, Y = Y_train,
#                                   max.ncomp = 2, eta.list = c(0.5, 0.9),
#                                   n_run = 2, k_folds = 5,
#                                   x.center = T, x.scale = F,
#                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                   remove_non_significant_models = F, alpha = 0.05,
#                                   w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0, times = NULL,
#                                   MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                   pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                   MIN_EPV = 5, return_models = F,
#                                   PARALLEL = F, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsdrcox_res #0.17min

## ---- fig.small=T, warning=F, eval=FALSE--------------------------------------
#  # plot cv.plsicox
#  cv.splsdrcox_res$plot_AUC

## -----------------------------------------------------------------------------
splsdrcox_model <- splsdrcox(X = X_train, Y = Y_train, 
                             n.comp = 2, #cv.splsdrcox_res$opt.comp, 
                             eta = 0.5, #cv.splsdrcox_res$opt.eta,
                             x.center = T, x.scale = F,
                             remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                             remove_non_significant = T,
                             MIN_EPV = 5, returnData = T, verbose = F)

splsdrcox_model

## ---- warning=FALSE, eval=FALSE-----------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdrcox_dynamic_res <- cv.splsdrcox_dynamic(X = X_train, Y = Y_train,
#                                                   max.ncomp = 2, vector = NULL,
#                                                   MIN_NVAR = 10, MAX_NVAR = 400,
#                                                   n.cut_points = 10, EVAL_METHOD = "AUC",
#                                                   n_run = 2, k_folds = 5,
#                                                   x.center = T, x.scale = F,
#                                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                                   remove_non_significant_models = F, alpha = 0.05,
#                                                   remove_variance_at_fold_level = F, remove_non_significant = F,
#                                                   w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0,
#                                                   times = NULL, max_time_points = 15, returnData = F,
#                                                   MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                                   pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                                   MIN_EPV = 5, return_models = F,
#                                                   PARALLEL = F, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsdrcox_dynamic_res #0.7mins

## -----------------------------------------------------------------------------
splsdrcox_dynamic_model <- splsdrcox_dynamic(X = X_train, Y = Y_train, 
                                             n.comp = 2, #cv.splsdrcox_dynamic_res$opt.comp,
                                             vector = 369, #cv.splsdrcox_dynamic_res$opt.nvar,
                                             x.center = T, x.scale = F,
                                             remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                             MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                             MIN_AUC_INCREASE = 0.01,
                                             EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                                             remove_non_significant = T, 
                                             MIN_EPV = 5, returnData = T, verbose = F)

splsdrcox_dynamic_model

## ---- warning=FALSE, eval=FALSE-----------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdacox_dynamic_res <- cv.splsdacox_dynamic(X = X_train, Y = Y_train,
#                                                   max.ncomp = 2, vector = NULL,
#                                                   MIN_NVAR = 10, MAX_NVAR = 400,
#                                                   n.cut_points = 10, EVAL_METHOD = "AUC",
#                                                   n_run = 2, k_folds = 5,
#                                                   x.center = T, x.scale = F,
#                                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                                   remove_variance_at_fold_level = F, remove_non_significant = F,
#                                                   remove_non_significant_models = F, alpha = 0.05,
#                                                   w_AIC = 0, w_c.index = 0, w_AUC = 1, w_BRIER = 0,
#                                                   times = NULL, max_time_points = 15, returnData = F,
#                                                   MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                                   pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                                   MIN_EPV = 5, return_models = F, max.iter = 200,
#                                                   PARALLEL = F, verbose = F, seed = 123)

## ---- eval=FALSE, eval=FALSE--------------------------------------------------
#  cv.splsdacox_dynamic_res #0.7min

## -----------------------------------------------------------------------------
splsdacox_dynamic_model <- splsdacox_dynamic(X = X_train, Y = Y_train, 
                                             n.comp = 2, #cv.splsdacox_dynamic_res$opt.comp, 
                                             vector = 330, #cv.splsdacox_dynamic_res$opt.nvar,
                                             x.center = T, x.scale = F,
                                             remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                             MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                             MIN_AUC_INCREASE = 0.01,
                                             EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                                             remove_non_significant = T, 
                                             MIN_EPV = 5, returnData = T, verbose = F)

splsdacox_dynamic_model

## -----------------------------------------------------------------------------
lst_models <- list("COX-EN" = coxen_model,
                   "sPLS-ICOX" = splsicox_model,
                   "sPLS-DRCOX" = splsdrcox_model,
                   "sPLS-DRCOX-Dynamic" = splsdrcox_dynamic_model,
                   "sPLS-DACOX-Dynamic" = splsdacox_dynamic_model)

eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                  X_test = X_test, Y_test = Y_test, 
                                  pred.method = "cenROC",
                                  pred.attr = "mean",
                                  times = NULL, max_time_points = 15, 
                                  PARALLEL = F)

## ---- eval=FALSE--------------------------------------------------------------
#  lst_evaluators <- c(cenROC = "cenROC", risksetROC = "risksetROC")
#  
#  eval_results <- purrr::map(lst_evaluators, ~eval_Coxmos_models(lst_models = lst_models,
#                                                                X_test = X_test, Y_test = Y_test,
#                                                                pred.method = .,
#                                                                pred.attr = "mean",
#                                                                times = NULL,
#                                                                max_time_points = 15,
#                                                                PARALLEL = F))

## -----------------------------------------------------------------------------
eval_results$cenROC

## ---- warning=F---------------------------------------------------------------
lst_eval_results <- plot_evaluation(eval_results$cenROC, evaluation = "AUC", pred.attr = "mean")
lst_eval_results_BRIER <- plot_evaluation(eval_results$cenROC, evaluation = "Brier", pred.attr = "mean")

## ---- fig.small=T, warning=F--------------------------------------------------
lst_eval_results$lst_plots$lineplot.mean
lst_eval_results$lst_plot_comparisons$anova

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test

## -----------------------------------------------------------------------------
lst_models_time <- list(#cv.coxen_res,
                        coxen_model,
                        #cv.splsicox_res,
                        splsicox_model,
                        #cv.splsdrcox_res,
                        splsdrcox_model,
                        #cv.splsdrcox_dynamic_res,
                        splsdrcox_dynamic_model,
                        #cv.splsdacox_dynamic_res,
                        splsdacox_dynamic_model, 
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.list(lst_models_time)

## ---- fig.small=T, warning=F--------------------------------------------------
ggp_time

## -----------------------------------------------------------------------------
#lst_forest_plot <- plot_forest.list(lst_models)
lst_forest_plot <- plot_forest(lst_models$`sPLS-DRCOX`)

## ---- fig.small=T, warning=F--------------------------------------------------
#lst_forest_plot$`sPLS-DRCOX`
lst_forest_plot

## -----------------------------------------------------------------------------
#lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)
lst_ph_ggplot <- plot_proportionalHazard(lst_models$`sPLS-DRCOX`)

## ---- fig.small=T, warning=F--------------------------------------------------
#lst_ph_ggplot$`sPLS-DRCOX`
lst_ph_ggplot

## -----------------------------------------------------------------------------
#density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")
density.plots.lp <- plot_cox.event(lst_models$`sPLS-DRCOX`, type = "lp")

## ---- fig.small=T, warning=F--------------------------------------------------
# density.plots.lp$`sPLS-DRCOX`$plot.density
# density.plots.lp$`sPLS-DRCOX`$plot.histogram

density.plots.lp$plot.density
density.plots.lp$plot.histogram

## -----------------------------------------------------------------------------
ggp_scores <- plot_PLS_Coxmos(model = lst_models$`sPLS-DRCOX`, 
                             comp = c(1,2), mode = "scores")

## ---- fig.small=T, warning=F--------------------------------------------------
ggp_scores$plot

## -----------------------------------------------------------------------------
ggp_loadings <- plot_PLS_Coxmos(model = lst_models$`sPLS-DRCOX`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10) #length from 0,0

## ---- fig.small=T, warning=F--------------------------------------------------
ggp_loadings$plot

## -----------------------------------------------------------------------------
ggp_biplot <- plot_PLS_Coxmos(model = lst_models$`sPLS-DRCOX`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T,
                             overlaps = 20)

## ---- fig.small=T, warning=F--------------------------------------------------
ggp_biplot$plot

## ---- warning=F---------------------------------------------------------------
variable_auc_results <- eval_Coxmos_model_per_variable(model = lst_models$`sPLS-DRCOX`, 
                                                       X_test = lst_models$`sPLS-DRCOX`$X_input, 
                                                       Y_test = lst_models$`sPLS-DRCOX`$Y_input,
                                                       pred.method = "cenROC", pred.attr = "mean",
                                                       times = NULL, max_time_points = 15,
                                                       PARALLEL = FALSE)

variable_auc_plot_train <- plot_evaluation(variable_auc_results, evaluation = "AUC")

## ---- fig.small=T, warning=F--------------------------------------------------
variable_auc_plot_train$lst_plots$lineplot.mean

## -----------------------------------------------------------------------------
# ggp.simulated_beta <- plot_pseudobeta.list(lst_models = lst_models, 
#                                            error.bar = T, onlySig = T, alpha = 0.05, 
#                                            zero.rm = T, auto.limits = T, top = 20,
#                                            show_percentage = T, size_percentage = 2, verbose = F)

ggp.simulated_beta <- plot_pseudobeta(model = lst_models$`sPLS-DRCOX`, 
                                      error.bar = T, onlySig = T, alpha = 0.05, 
                                      zero.rm = T, auto.limits = T, top = 20,
                                      show_percentage = T, size_percentage = 2)

## ---- fig.small=T, warning=F--------------------------------------------------
#ggp.simulated_beta$`sPLS-DRCOX`$plot
ggp.simulated_beta$plot

## -----------------------------------------------------------------------------
# LST_KM_RES_LP <- getAutoKM.list(type = "LP",
#                                 lst_models = lst_models,
#                                 comp = 1:10,
#                                 top = 10,
#                                 ori_data = T,
#                                 BREAKTIME = NULL,
#                                 only_sig = T, alpha = 0.05)

LST_KM_RES_LP <- getAutoKM(type = "LP",
                           model = lst_models$`sPLS-DRCOX`,
                           comp = 1:10,
                           top = 10,
                           ori_data = T,
                           BREAKTIME = NULL,
                           only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=F--------------------------------------------------
#LST_KM_RES_LP$`sPLS-DRCOX`$LST_PLOTS$LP
LST_KM_RES_LP$LST_PLOTS$LP

## -----------------------------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_LP)
# LST_KM_TEST_LP <- getTestKM.list(lst_models = lst_models, 
#                                  X_test = X_test, Y_test = Y_test, 
#                                  type = "LP",
#                                  BREAKTIME = NULL, n.breaks = 20,
#                                  lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM(model = lst_models$`sPLS-DRCOX`, 
                            X_test = X_test, Y_test = Y_test, 
                            type = "LP",
                            BREAKTIME = NULL, n.breaks = 20,
                            cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
#LST_KM_TEST_LP$`sPLS-DRCOX`
LST_KM_TEST_LP

## -----------------------------------------------------------------------------
# LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
#                                   lst_models = lst_models,
#                                   comp = 1:10,
#                                   top = 10,
#                                   ori_data = T,
#                                   BREAKTIME = NULL,
#                                   n.breaks = 20,
#                                   only_sig = T, alpha = 0.05)

LST_KM_RES_COMP <- getAutoKM(type = "COMP",
                             model = lst_models$`sPLS-DRCOX`,
                             comp = 1:10,
                             top = 10,
                             ori_data = T,
                             BREAKTIME = NULL,
                             n.breaks = 20,
                             only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=F--------------------------------------------------
# LST_KM_RES_COMP$`sPLS-DRCOX`$LST_PLOTS$comp_1
# LST_KM_RES_COMP$`sPLS-DRCOX`$LST_PLOTS$comp_2

LST_KM_RES_COMP$LST_PLOTS$comp_1
LST_KM_RES_COMP$LST_PLOTS$comp_2

## -----------------------------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_COMP)
# LST_KM_TEST_COMP <- getTestKM.list(lst_models = lst_models, 
#                                    X_test = X_test, Y_test = Y_test, 
#                                    type = "COMP",
#                                    BREAKTIME = NULL, n.breaks = 20,
#                                    lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM(model = lst_models$`sPLS-DRCOX`, 
                              X_test = X_test, Y_test = Y_test, 
                              type = "COMP",
                              BREAKTIME = NULL, n.breaks = 20,
                              cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
# all patients could be categorize into the same group
# LST_KM_TEST_COMP$`sPLS-DRCOX`$comp_1
# LST_KM_TEST_COMP$`sPLS-DRCOX`$comp_2

LST_KM_TEST_COMP$comp_1
LST_KM_TEST_COMP$comp_2

## -----------------------------------------------------------------------------
# LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
#                                  lst_models = lst_models,
#                                  comp = 1:10, #select how many components you want to compute for the pseudo beta
#                                  top = 10,
#                                  ori_data = T, #original data selected
#                                  BREAKTIME = NULL,
#                                  only_sig = T, alpha = 0.05)

LST_KM_RES_VAR <- getAutoKM(type = "VAR",
                            model = lst_models$`sPLS-DRCOX`,
                            comp = 1:10, #select how many components you want to compute for the pseudo beta
                            top = 10,
                            ori_data = T, #original data selected
                            BREAKTIME = NULL,
                            only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=F--------------------------------------------------
# LST_KM_RES_VAR$`sPLS-DRCOX`$LST_PLOTS$`840`
# LST_KM_RES_VAR$`sPLS-DRCOX`$LST_PLOTS$`3897`

LST_KM_RES_VAR$LST_PLOTS$`840`
LST_KM_RES_VAR$LST_PLOTS$`3897`

## -----------------------------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_VAR)
# LST_KM_TEST_VAR <- getTestKM.list(lst_models = lst_models, 
#                                   X_test = X_test, Y_test = Y_test, 
#                                   type = "VAR", ori_data = T,
#                                   BREAKTIME = NULL, n.breaks = 20,
#                                   lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM(model = lst_models$`sPLS-DRCOX`, 
                             X_test = X_test, Y_test = Y_test, 
                             type = "VAR", ori_data = T,
                             BREAKTIME = NULL, n.breaks = 20,
                             cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
# LST_KM_TEST_VAR$`sPLS-DRCOX`$`840`
# LST_KM_TEST_VAR$`sPLS-DRCOX`$`3897`

LST_KM_TEST_VAR$`840`
LST_KM_TEST_VAR$`3897`

## -----------------------------------------------------------------------------
new_pat <- X_test[1,,drop=F]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat),])

## ---- warning=F---------------------------------------------------------------
# ggp.simulated_beta_newObs <- plot_pseudobeta_newObservation.list(lst_models = lst_models, 
#                                                              new_observation = new_pat,
#                                                              error.bar = T, onlySig = T, alpha = 0.05,
#                                                              zero.rm = T, auto.limits = T, show.betas = T, top = 20)

ggp.simulated_beta_newObs <- plot_pseudobeta_newObservation(model = lst_models$`sPLS-DRCOX`,
                                                        new_observation = new_pat,
                                                        error.bar = T, onlySig = T, alpha = 0.05,
                                                        zero.rm = T, auto.limits = T, show.betas = T, top = 20)

## ---- fig.small=T, warning=F--------------------------------------------------
#ggp.simulated_beta_newObs$`sPLS-DRCOX`$plot

ggp.simulated_beta_newObs$plot

## -----------------------------------------------------------------------------
pat_density <- plot_observation.eventDensity(observation = new_pat,
                                             model = lst_models$`sPLS-DRCOX`,
                                             time = NULL, 
                                             type = "lp")

## ---- fig.small=T, warning=F--------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_observation.eventHistogram(observation = new_pat, 
                                                 model = lst_models$`sPLS-DRCOX`, 
                                                 time = NULL, 
                                                 type = "lp")

## ---- fig.small=T, warning=F--------------------------------------------------
pat_histogram

## -----------------------------------------------------------------------------
knitr::kable(Y_test[1:5,])

## ---- warning=F---------------------------------------------------------------
# lst_cox.comparison <- plot_LP.multipleObservations.list(lst_models = lst_models,
#                                                     new_data = X_test[1:5,],
#                                                     error.bar = T, zero.rm = T, onlySig = T, 
#                                                     alpha = 0.05, top = 10)

lst_cox.comparison <- plot_LP.multipleObservations(model = lst_models$`sPLS-DRCOX`,
                                                   new_data = X_test[1:5,],
                                                   error.bar = T, zero.rm = T, onlySig = T, 
                                                   alpha = 0.05, top = 10)

## ---- fig.small=T, warning=F--------------------------------------------------
# lst_cox.comparison$`sPLS-DRCOX`$plot

lst_cox.comparison$plot

