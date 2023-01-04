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

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ConesaLab/HDcox", build_vignettes = TRUE)

## ----setup, results = "hide"--------------------------------------------------
# load HDcox
library(HDcox)

## ---- eval = FALSE------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ConesaLab/RColorConesa")

## -----------------------------------------------------------------------------
library(RColorConesa)
#theme_set(theme_colorConesa()) #under development

## -----------------------------------------------------------------------------
# load Tasic dataset
data("X_small_data_E.MTAB.386")
data("Y_small_data_E.MTAB.386")

X <- X_small_data_E.MTAB.386
Y <- Y_small_data_E.MTAB.386

rm(X_small_data_E.MTAB.386, Y_small_data_E.MTAB.386)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(X[1:5,1:5])

knitr::kable(Y[1:5,])

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(dim(X))

## -----------------------------------------------------------------------------
ggp_density.event <- plot_events(Y = Y, roundTo = 0.25, categories = c("Censored","Death")) #name for F and T

## ----fig.small = T------------------------------------------------------------
ggp_density.event$plot

## -----------------------------------------------------------------------------
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = .7, #70% train
                                          list = FALSE,
                                          times = 1)

X_train <- X[index_train,] #1103
Y_train <- Y[index_train,]
X_test <- X[-index_train,] #472
Y_test <- Y[-index_train,]

## ---- eval = FALSE, message=T, error=F----------------------------------------
#  # classical approach
#  cox_model <- cox(X = X_train, Y = Y_train, x.center = T, x.scale = F)

## -----------------------------------------------------------------------------
EPV <- sum(Y_train$event==1) / ncol(X_train)
EPV

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.plsicox
#  cv.plsicox_res <- cv.plsicox(X = X_train, Y = Y_train,
#                               max.ncomp =  4,
#                               n_run = 2, k_folds = 10,
#                               x.scale = T, x.center = T,
#                               y.center = F, y.scale = F,
#                               remove_near_zero_variance = F, remove_zero_variance = F, toKeep.zv = NULL,
#                               remove_non_significant_models = F, alpha = 0.05,
#                               w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
#                               MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                               pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                               MIN_EPV = 5, return_models = F,
#                               PARALLEL = T, verbose = F, seed = 123)
#  cv.plsicox_res #1min 8s.

## ---- eval = FALSE, fig.small=T-----------------------------------------------
#  # plot cv.plsicox
#  cv.plsicox_res$plot_AUC

## -----------------------------------------------------------------------------
plsicox_model <- plsicox(X = X_train, Y = Y_train, 
                         n.comp = 3, #n.comp = cv.plsicox_res$opt.comp
                         x.center = T, x.scale = F)

plsicox_model

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdrcox_res <- cv.splsdrcox(X = X_train, Y = Y_train,
#                                   max.ncomp = 4, eta.list = seq(0,0.5,0.25), #penalty
#                                   n_run = 2, k_folds = 10,
#                                   x.scale = T,
#                                   remove_near_zero_variance = F, remove_zero_variance = F,
#                                   PARALLEL = T, verbose = F)
#  
#  cv.splsdrcox_res #2min 40s.

## -----------------------------------------------------------------------------
splsdrcox_model <- splsdrcox(X = X_train, Y = Y_train, 
                             n.comp = 2, eta = 0.25, #n.comp = cv.splsdrcox_res$opt.comp, eta = cv.splsdrcox_res$opt.eta
                             x.center = T, x.scale = F)

splsdrcox_model

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdacox_res <- cv.splsdacox_mixOmics(X = X_train, Y = Y_train,
#                                          max.ncomp = 4,  #penalty
#                                          n_run = 2, k_folds = 10,
#                                          x.scale = T,
#                                          remove_near_zero_variance = F, remove_zero_variance = F,
#                                          PARALLEL = T, verbose = F)
#  
#  cv.splsdacox_res #2min

## -----------------------------------------------------------------------------
splsdacox_model <- splsdacox_mixOmics(X = X_train, Y = Y_train, 
                                    n.comp = 3, #cv.splsdacox_res$opt.comp
                                    x.center = T, x.scale = F)

splsdacox_model

## -----------------------------------------------------------------------------
lst_models <- list("PLS-ICOX" = plsicox_model,
                   "SPLS-DRCOX" = splsdrcox_model,
                   "PLS-DACOX" = splsdacox_model)

eval_results <- eval_models4.0(lst_models = lst_models,
                               X_test = X_test, Y_test = Y_test, 
                               pred.method = "cenROC",
                               pred.attr = "mean",
                               times = seq(1,4,0.5), max_time_points = 15, 
                               PARALLEL = T)

# lst_evaluators <- c(cenROC = "cenROC", 
#                     risksetROC = "risksetROC")
# 
# eval_results <- purrr::map(lst_evaluators, ~eval_models4.0(lst_models = lst_models,
#                                                            X_test = X_test, Y_test = Y_test, 
#                                                            pred.method = .,
#                                                            pred.attr = "mean",
#                                                            times = seq(1,4,0.5), max_time_points = 15, 
#                                                            PARALLEL = T))

## -----------------------------------------------------------------------------
eval_results
#eval_results$cenROC

## -----------------------------------------------------------------------------
lst_eval_results <- plot_evaluation(eval_results)
#lst_eval_results <- plot_evaluation.list(eval_results)

## ---- fig.small=T-------------------------------------------------------------
lst_eval_results$lst_plots$lineplot.mean
lst_eval_results$lst_plot_comparisons$t.test

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test

## -----------------------------------------------------------------------------
lst_models_time <- list(plsicox_model,
                        splsdrcox_model,
                        splsdacox_model,
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.models(lst_models_time)

## ---- fig.small=T-------------------------------------------------------------
ggp_time

## -----------------------------------------------------------------------------
lst_forest_plot <- purrr::map(lst_models, ~survminer::ggforest(.$survival_model$fit, 
                                                               data = .$survival_model$fit$model))

## ---- fig.small=T-------------------------------------------------------------
lst_forest_plot$`SPLS-DRCOX`

## -----------------------------------------------------------------------------
density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")

## ---- fig.small=T-------------------------------------------------------------
density.plots.lp$`SPLS-DRCOX`$plot.density
density.plots.lp$`SPLS-DRCOX`$plot.histogram

## -----------------------------------------------------------------------------
lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)

## ---- fig.small=T-------------------------------------------------------------
lst_ph_ggplot$`SPLS-DRCOX`

## -----------------------------------------------------------------------------
ggp_scores <- plot_PLS_HDcox(model = lst_models$`SPLS-DRCOX`, 
                             comp = c(1,2), mode = "scores")

## ---- fig.small=T-------------------------------------------------------------
ggp_scores$plot

## -----------------------------------------------------------------------------
ggp_loadings <- plot_PLS_HDcox(model = lst_models$`SPLS-DRCOX`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10) #length from 0,0

## ---- fig.small=T-------------------------------------------------------------
ggp_loadings$plot

## -----------------------------------------------------------------------------
ggp_biplot <- plot_PLS_HDcox(model = lst_models$`SPLS-DRCOX`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T)

## ---- fig.small=T-------------------------------------------------------------
ggp_biplot$plot

## -----------------------------------------------------------------------------
ggp.simulated_beta <- plot_pseudobeta.list(lst_models = lst_models, 
                                           error.bar = T, onlySig = T, alpha = 0.05, 
                                           zero.rm = T, auto.limits = T, top = 20,
                                           show_percentage = T, size_percentage = 3)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta$`SPLS-DRCOX`$plot

## -----------------------------------------------------------------------------
LST_KM_RES_LP <- getAutoKM.list(type = "LP",
                                lst_models = lst_models,
                                comp = 1:4,
                                top = 10,
                                ori_data = T,
                                BREAKTIME = NULL,
                                only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_LP$`SPLS-DRCOX`$LST_PLOTS$LP

## -----------------------------------------------------------------------------
LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
                                  lst_models = lst_models,
                                  comp = 1:4,
                                  top = 10,
                                  ori_data = T,
                                  BREAKTIME = NULL,
                                  only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_COMP$`SPLS-DRCOX`$LST_PLOTS$comp_1
LST_KM_RES_COMP$`SPLS-DRCOX`$LST_PLOTS$comp_2

## -----------------------------------------------------------------------------
LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
                                 lst_models = lst_models,
                                 comp = 1:4,
                                 top = 10,
                                 ori_data = T,
                                 BREAKTIME = NULL,
                                 only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_VAR$`SPLS-DRCOX`$LST_PLOTS$POSTN
LST_KM_RES_VAR$`SPLS-DRCOX`$LST_PLOTS$SIRT5

## -----------------------------------------------------------------------------
new_pat <- X_test[1,,drop=F]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat),])

## -----------------------------------------------------------------------------
ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient.list(lst_models = lst_models, 
                                                             new_pat = new_pat,
                                                             error.bar = T, onlySig = T, alpha = 0.05,
                                                             zero.rm = T, auto.limits = T, show.betas = T, top = 20)

# ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient(model = lst_models$`SPLS-DRCOX`, 
#                                                         new_pat = new_pat,
#                                                         error.bar = T, onlySig = T, alpha = 0.05,
#                                                         zero.rm = T, auto.limits = T, show.betas = T, top = 20)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta_newPat$`SPLS-DRCOX`$plot

## -----------------------------------------------------------------------------
pat_density <- plot_patient.eventDensity(patient = new_pat, 
                                         time = NULL, 
                                         model = lst_models$`SPLS-DRCOX`, 
                                         type = "lp")

## ---- fig.small=T-------------------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_patient.eventHistogram(patient = new_pat, 
                                             time = NULL, 
                                             model = lst_models$`SPLS-DRCOX`, 
                                             type = "lp")

## ---- fig.small=T-------------------------------------------------------------
pat_histogram

## ---- eval=F------------------------------------------------------------------
#  #plot_divergent.biplot - for num and qual variables

## -----------------------------------------------------------------------------
knitr::kable(Y_test[1:5,])

## -----------------------------------------------------------------------------
lst_cox.comparison <- plot_LP.multiplePatients.list(lst_models = lst_models, 
                                     df.pat = X_test[1:5,], 
                                     error.bar = T, zero.rm = T, onlySig = T, alpha = 0.05, top = 5)

## ---- fig.small=T-------------------------------------------------------------
lst_cox.comparison$`SPLS-DRCOX`$plot

