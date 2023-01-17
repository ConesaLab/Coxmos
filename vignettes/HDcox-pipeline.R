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
#  devtools::install_github("ConesaLab/HDcox", build_vignettes = TRUE)

## ----setup, eval=FALSE, results = "hide"--------------------------------------
#  # load HDcox
#  library(HDcox)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ConesaLab/RColorConesa")

## -----------------------------------------------------------------------------
library(RColorConesa)
#theme_set(theme_colorConesa()) #under development

## -----------------------------------------------------------------------------
# load Tasic dataset
data("X_miRNA_glioblastoma")
data("Y_miRNA_glioblastoma")

X <- X_miRNA_glioblastoma
Y <- Y_miRNA_glioblastoma

rm(X_miRNA_glioblastoma, Y_miRNA_glioblastoma)

## ---- echo = FALSE------------------------------------------------------------
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
set.seed(321)
index_train <- caret::createDataPartition(Y$event,
                                          p = .8, # 80% train
                                          list = FALSE,
                                          times = 1)

## -----------------------------------------------------------------------------
X_train <- X[index_train,] #443x534
Y_train <- Y[index_train,]
X_test <- X[-index_train,] #109x534
Y_test <- Y[-index_train,]

## ---- eval=FALSE, message=T, error=F------------------------------------------
#  # classical approach
#  cox_model <- cox(X = X_train, Y = Y_train,
#                   x.center = T, x.scale = F,
#                   y.center = F, y.scale = F,
#                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
#                   remove_non_significant = F, alpha = 0.05,
#                   MIN_EPV = 5, FORCE = F, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
EPV <- getEPV(X_train, Y_train)

## -----------------------------------------------------------------------------
EPV

## ---- eval=FALSE, warning=F---------------------------------------------------
#  # run cv.coxEN
#  cv.coxen_res <- cv.coxEN(X = X_train, Y = Y_train,
#                           EN.alpha.list = seq(0,1,0.1),
#                           max.variables = ncol(X_train),
#                           n_run = 2, k_folds = 10,
#                           x.center = T, x.scale = F,
#                           y.center = F, y.scale = F,
#                           remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                           remove_non_significant = F, alpha = 0.05,
#                           w_AIC = 0, w_c.index = 0, w_AUC = 1, times = NULL,
#                           MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                           pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                           MIN_EPV = 5, return_models = F,
#                           PARALLEL = T, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.coxen_res #2.7min.

## -----------------------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0, #cv.coxen_res$opt.EN.alpha
                     max.variables = 47, #cv.coxen_res$opt.nvar
                     x.center = T, x.scale = F, 
                     y.center = F, y.scale = F, 
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = F, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## -----------------------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0, #cv.coxen_res$opt.EN.alpha
                     max.variables = 47, #cv.coxen_res$opt.nvar
                     x.center = T, x.scale = F, 
                     y.center = F, y.scale = F, 
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = T, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## ---- eval=FALSE, message=F---------------------------------------------------
#  # run cv.plsicox
#  cv.plsicox_res <- cv.plsicox(X = X_train, Y = Y_train,
#                               max.ncomp = 10,
#                               n_run = 2, k_folds = 10,
#                               x.center = T, x.scale = F,
#                               y.center = F, y.scale = F,
#                               remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                               remove_non_significant_models = F, alpha = 0.05,
#                               w_AIC = 0, w_c.index = 0, w_AUC = 1, times = NULL,
#                               MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                               pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                               MIN_EPV = 5, return_models = F,
#                               PARALLEL = T, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.plsicox_res #5.45min.

## ---- eval=FALSE, fig.small=T-------------------------------------------------
#  # plot cv.plsicox
#  cv.plsicox_res$plot_AUC

## -----------------------------------------------------------------------------
plsicox_model <- plsicox(X = X_train, Y = Y_train, 
                         n.comp = 3, #n.comp = cv.plsicox_res$opt.comp
                         x.center = T, x.scale = F,
                         y.center = F, y.scale = F,
                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                         tol = 500, 
                         MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
plsicox_model

## ---- eval=FALSE, message=F---------------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdrcox_res <- cv.splsdrcox(X = X_train, Y = Y_train,
#                                   max.ncomp = 10, eta.list = seq(0,0.9,0.25), #penalty
#                                   n_run = 2, k_folds = 10,
#                                   x.center = T, x.scale = F,
#                                   y.center = F, y.scale = F,
#                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                   remove_non_significant_models = F, alpha = 0.05,
#                                   w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
#                                   MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                   pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                   MIN_EPV = 5, return_models = F,
#                                   PARALLEL = T, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsdrcox_res #10.4min

## ---- eval=FALSE, fig.small=T-------------------------------------------------
#  # plot cv.plsicox
#  cv.splsdrcox_res$plot_AUC

## -----------------------------------------------------------------------------
splsdrcox_model <- splsdrcox(X = X_train, Y = Y_train, 
                             n.comp = 1, eta = 0.75, #n.comp = cv.splsdrcox_res$opt.comp, eta = cv.splsdrcox_res$opt.eta
                             x.center = T, x.scale = F,
                             y.center = F, y.scale = F,
                             remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                             MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
splsdrcox_model

## ---- eval=FALSE--------------------------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdrcox_mo_res <- cv.splsdrcox_mixOmics(X = X_train, Y = Y_train,
#                                               max.ncomp = 10, vector = NULL,
#                                               MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 10, EVAL_METHOD = "cenROC",
#                                               n_run = 2, k_folds = 10,
#                                               x.center = T, x.scale = F,
#                                               y.center = F, y.scale = F,
#                                               remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                               remove_non_significant_models = F, alpha = 0.05,
#                                               w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
#                                               MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                               pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                               MIN_EPV = 5, return_models = F,
#                                               PARALLEL = T, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsdrcox_mo_res #3.45mins

## -----------------------------------------------------------------------------
splsdrcox_mo_model <- splsdrcox_mixOmics(X = X_train, Y = Y_train, 
                                         n.comp = 2, vector = 534,
                                         x.center = T, x.scale = F,
                                         y.center = F, y.scale = F,
                                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                         MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                         MIN_AUC_INCREASE = 0.01,
                                         EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                                         MIN_EPV = 5, returnData = T, verbose = F)

splsdrcox_mo_model

## ---- eval=FALSE--------------------------------------------------------------
#  # run cv.splsdrcox
#  cv.splsdacox_res <- cv.splsdacox_mixOmics(X = X_train, Y = Y_train,
#                                            max.ncomp = 10,  vector = NULL,
#                                            MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 10, EVAL_METHOD = "cenROC",
#                                            n_run = 2, k_folds = 10,
#                                            x.center = T, x.scale = F,
#                                            y.center = F, y.scale = F,
#                                            remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                            remove_non_significant_models = F, alpha = 0.05,
#                                            w_AIC = 0,  w_c.index = 0, w_AUC = 1, times = NULL,
#                                            MIN_AUC_INCREASE = 0.05, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                            pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                            MIN_EPV = 5, return_models = F,
#                                            PARALLEL = T, verbose = F, seed = 123)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.splsdacox_res #4min

## -----------------------------------------------------------------------------
splsdacox_mo_model <- splsdacox_mixOmics(X = X_train, Y = Y_train, 
                                         n.comp = 6, vector = 184,
                                         x.center = T, x.scale = F,
                                         y.center = F, y.scale = F,
                                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                         MIN_NVAR = 10, MAX_NVAR = 1000, n.cut_points = 5,
                                         MIN_AUC_INCREASE = 0.01,
                                         EVAL_METHOD = "AUC", pred.method = "cenROC", max.iter = 200,
                                         MIN_EPV = 5, returnData = T, verbose = F)

splsdacox_mo_model

## -----------------------------------------------------------------------------
lst_models <- list("COX-EN" = coxen_model,
                   "PLS-ICOX" = plsicox_model,
                   "sPLS-DRCOX" = splsdrcox_model,
                   "sPLS-DRCOX-MixOmics" = splsdrcox_mo_model,
                   "sPLS-DACOX-MixOmics" = splsdacox_mo_model)

eval_results <- eval_models4.0(lst_models = lst_models,
                               X_test = X_test, Y_test = Y_test, 
                               pred.method = "cenROC",
                               pred.attr = "mean",
                               times = NULL, max_time_points = 15, 
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
lst_eval_results$lst_plots$lineplot
lst_eval_results$lst_plot_comparisons$anova

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test

## -----------------------------------------------------------------------------
lst_models_time <- list(coxen_model,
                        plsicox_model,
                        splsdrcox_model,
                        splsdrcox_mo_model,
                        splsdacox_mo_model)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.list(lst_models_time)

## ---- fig.small=T-------------------------------------------------------------
ggp_time

## -----------------------------------------------------------------------------
lst_forest_plot <- plot_forest.list(lst_models)

## ---- fig.small=T-------------------------------------------------------------
lst_forest_plot$`sPLS-DRCOX-MixOmics`

## -----------------------------------------------------------------------------
lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)

## ---- fig.small=T-------------------------------------------------------------
lst_ph_ggplot$`sPLS-DRCOX-MixOmics`

## -----------------------------------------------------------------------------
density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")

## ---- fig.small=T-------------------------------------------------------------
density.plots.lp$`sPLS-DRCOX-MixOmics`$plot.density
density.plots.lp$`sPLS-DRCOX-MixOmics`$plot.histogram

## -----------------------------------------------------------------------------
ggp_scores <- plot_PLS_HDcox(model = lst_models$`sPLS-DRCOX-MixOmics`, 
                             comp = c(1,2), mode = "scores")

## ---- fig.small=T-------------------------------------------------------------
ggp_scores$plot

## -----------------------------------------------------------------------------
ggp_loadings <- plot_PLS_HDcox(model = lst_models$`sPLS-DRCOX-MixOmics`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10) #length from 0,0

## ---- fig.small=T-------------------------------------------------------------
ggp_loadings$plot

## -----------------------------------------------------------------------------
ggp_biplot <- plot_PLS_HDcox(model = lst_models$`sPLS-DRCOX-MixOmics`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T,
                             overlaps = 20)

## ---- fig.small=T-------------------------------------------------------------
ggp_biplot$plot

## -----------------------------------------------------------------------------
ggp.simulated_beta <- plot_pseudobeta.list(lst_models = lst_models, 
                                           error.bar = T, onlySig = T, alpha = 0.05, 
                                           zero.rm = T, auto.limits = T, top = 20,
                                           show_percentage = T, size_percentage = 3, verbose = F)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta$`sPLS-DRCOX-MixOmics`$plot

## -----------------------------------------------------------------------------
LST_KM_RES_LP <- getAutoKM.list(type = "LP",
                                lst_models = lst_models,
                                comp = 1:4,
                                top = 10,
                                ori_data = T,
                                BREAKTIME = NULL,
                                only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_LP$`sPLS-DRCOX-MixOmics`$LST_PLOTS$LP

## -----------------------------------------------------------------------------
LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
                                  lst_models = lst_models,
                                  comp = 1:4,
                                  top = 10,
                                  ori_data = T,
                                  BREAKTIME = NULL,
                                  n.breaks = 20,
                                  only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_COMP$`sPLS-DRCOX-MixOmics`$LST_PLOTS$comp_1
LST_KM_RES_COMP$`sPLS-DRCOX-MixOmics`$LST_PLOTS$comp_2

## -----------------------------------------------------------------------------
LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
                                 lst_models = lst_models,
                                 comp = 1:10, #select how many components you want to compute for the pseudo beta
                                 top = 10,
                                 ori_data = T,
                                 BREAKTIME = NULL,
                                 only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_VAR$`sPLS-DRCOX-MixOmics`$LST_PLOTS$hsa_miR_222
LST_KM_RES_VAR$`sPLS-DRCOX-MixOmics`$LST_PLOTS$hsa_miR_182

## -----------------------------------------------------------------------------
new_pat <- X_test[1,,drop=F]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat),])

## -----------------------------------------------------------------------------
ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient.list(lst_models = lst_models, 
                                                             new_pat = new_pat,
                                                             error.bar = T, onlySig = T, alpha = 0.05,
                                                             zero.rm = T, auto.limits = T, show.betas = T, top = 20)

# ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient(model = lst_models$`sPLS-DRCOX-MixOmics`, 
#                                                         new_pat = new_pat,
#                                                         error.bar = T, onlySig = T, alpha = 0.05,
#                                                         zero.rm = T, auto.limits = T, show.betas = T, top = 20)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta_newPat$`sPLS-DRCOX-MixOmics`$plot

## -----------------------------------------------------------------------------
pat_density <- plot_patient.eventDensity(patient = new_pat, 
                                         time = NULL, 
                                         model = lst_models$`sPLS-DRCOX-MixOmics`, 
                                         type = "lp")

## ---- fig.small=T-------------------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_patient.eventHistogram(patient = new_pat, 
                                             time = NULL, 
                                             model = lst_models$`sPLS-DRCOX-MixOmics`, 
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
lst_cox.comparison$`sPLS-DRCOX-MixOmics`$plot

