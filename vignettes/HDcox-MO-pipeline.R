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
# load dataset
data("X_small_mo.data_Glioblastoma", package = "HDcox")
data("Y_small_mo.data_Glioblastoma", package = "HDcox")

X <- X_small_mo.data_Glioblastoma
Y <- Y_small_mo.data_Glioblastoma

rm(X_small_mo.data_Glioblastoma, Y_small_mo.data_Glioblastoma)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(X$genes[1:5,1:5]);knitr::kable(X$miRNA[1:5,1:5]);knitr::kable(X$clinical[1:5,1:5])

knitr::kable(Y[1:5,])

## -----------------------------------------------------------------------------
ggp_density.event <- plot_events(Y = Y, roundTo = 150, categories = c("Censored","Death")) #name for F and T

## ----fig.small = T------------------------------------------------------------
ggp_density.event$plot

## -----------------------------------------------------------------------------
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = .7, #70% train
                                          list = FALSE,
                                          times = 1)

X_train <- list()
X_test <- list()
for(omic in names(X)){
  X_train[[omic]] <- X[[omic]][index_train,,drop=F]
  X_test[[omic]] <- X[[omic]][-index_train,,drop=F]
}

Y_train <- Y[index_train,]
Y_test <- Y[-index_train,]

## -----------------------------------------------------------------------------
for(b in names(X_train)){
  EPV <- sum(Y_train$event==1) / ncol(X_train[[b]])
  message(paste0("EPV = ", round(EPV, 2), ", for block ", b))
}


## ---- eval = T, message=F-----------------------------------------------------
x.center = c(genes = T, miRNA = T, clinical = T) #if vector, must be named
x.scale = c(genes = F, miRNA = F, clinical = T) #if vector, must be named

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.sb.plsicox
#  cv.sb.plsicox_res <- cv.sb.plsicox(X = X_train, Y = Y_train,
#                                     max.ncomp = 4,
#                                     n_run = 2, k_folds = 10,
#                                     x.center = x.center, x.scale = x.scale,
#                                     remove_near_zero_variance = T,
#                                     PARALLEL = T, verbose = F)
#  
#  cv.sb.plsicox_res #3min 10s.

## ---- eval = FALSE, fig.small=T-----------------------------------------------
#  # plot cv.plsicox
#  cv.sb.plsicox_res$plot_AUC

## -----------------------------------------------------------------------------
sb.plsicox_model <- sb.plsicox(X = X_train,
                               Y = Y_train,
                               n.comp = 4, #cv.sb.plsicox_res$opt.comp
                               x.center = x.center, x.scale = x.scale,
                               remove_near_zero_variance = T,
                               returnData = T, verbose = F)

sb.plsicox_model

## -----------------------------------------------------------------------------
sb.plsicox_model <- sb.plsicox(X = X_train,
                               Y = Y_train,
                               n.comp = 4, #cv.sb.plsicox_res$opt.comp
                               x.center = x.center, x.scale = x.scale,
                               remove_near_zero_variance = T,
                               remove_non_significant = T,
                               returnData = T, verbose = F)

sb.plsicox_model

## ---- eval=FALSE, message=F---------------------------------------------------
#  # run cv.sb.plsicox
#  fast.sb.plsicox_model <- fast.cv.sb.plsicox(X = X_train, Y = Y_train,
#                                               max.ncomp = 4,
#                                               n_run = 2, k_folds = 10,
#                                               x.center = x.center, x.scale = x.scale,
#                                               remove_near_zero_variance = T,
#                                               remove_non_significant = T,
#                                               PARALLEL = T, verbose = F)
#  
#  fast.sb.plsicox_model #6min 7s.

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.sb.plsicox
#  cv.sb.splsdrcox_res <- cv.sb.splsdrcox(X = X_train, Y = Y_train,
#                                         max.ncomp = 4, eta.list = c(0.25,0.5,0.75),
#                                         n_run = 2, k_folds = 10,
#                                         x.center = x.center, x.scale = x.scale,
#                                         remove_near_zero_variance = T,
#                                         remove_non_significant = T,
#                                         PARALLEL = T, verbose = F)
#  
#  cv.sb.splsdrcox_res #5min

## -----------------------------------------------------------------------------
sb.splsdrcox_model <- sb.splsdrcox(X = X_train, 
                                   Y = Y_train, 
                                   n.comp = 1, eta = 0.75, #n.comp = cv.splsdrcox_res$opt.comp, eta = cv.splsdrcox_res$opt.eta
                                   x.center = x.center, x.scale = x.scale,
                                   remove_near_zero_variance = T,
                                   remove_non_significant = T,
                                   returnData = T, verbose = F)

sb.splsdrcox_model

## ---- eval=FALSE, message=F---------------------------------------------------
#  # run cv.sb.plsicox
#  fast.sb.splsdrcox_model <- fast.cv.sb.splsdrcox(X = X_train, Y = Y_train,
#                                                  max.ncomp = 4, eta.list = c(0.25,0.5,0.75),
#                                                  n_run = 2, k_folds = 10,
#                                                  x.center = x.center, x.scale = x.scale,
#                                                  remove_near_zero_variance = T,
#                                                  remove_non_significant = T,
#                                                  PARALLEL = T, verbose = F)
#  
#  fast.sb.splsdrcox_model #7.5min

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.splsdrcox
#  cv.mb.splsdrcox_res <- cv.mb.splsdrcox(X = X_train, Y = Y_train,
#                                         max.ncomp = 4, vector = NULL, #NULL - autodetection
#                                         n_run = 2, k_folds = 10,
#                                         x.center = x.center, x.scale = x.scale,
#                                         remove_near_zero_variance = T,
#                                         remove_zero_variance = T,
#                                         PARALLEL = T, verbose = F)
#  
#  cv.mb.splsdrcox_res #2min

## -----------------------------------------------------------------------------
mb.splsdrcox_model <- mb.splsdrcox(X = X_train, Y = Y_train, 
                                        n.comp = 4, #cv.mb.splsdrcox_res$opt.comp
                                        vector = list(genes = 10, miRNA = 10, clinical = 10), #cv.mb.splsdrcox_res$opt.nvar
                                        x.center = x.center, x.scale = x.scale, 
                                        remove_near_zero_variance = T, 
                                        remove_zero_variance = T,
                                        verbose = F)

mb.splsdrcox_model

## ---- eval = FALSE, message=F-------------------------------------------------
#  # run cv.splsdrcox
#  cv.mb.splsdacox_res <- cv.mb.splsdacox(X = X_train, Y = Y_train,
#                                     max.ncomp = 4, vector = NULL, #NULL - autodetection
#                                     n_run = 2, k_folds = 10,
#                                     x.center = x.center, x.scale = x.scale,
#                                     remove_near_zero_variance = T,
#                                     remove_zero_variance = T,
#                                     PARALLEL = T, verbose = F)
#  
#  cv.mb.splsdacox_res #2min

## -----------------------------------------------------------------------------
mb.splsdacox_model <- mb.splsdacox(X = X_train, Y = Y_train, 
                                        n.comp = 3, #cv.mb.splsdacox_res$opt.comp
                                        vector = list(genes = 10, miRNA = 255, clinical = 10), #cv.mb.splsdacox_res$opt.nvar
                                        x.center = x.center, x.scale = x.scale, 
                                        remove_near_zero_variance = T, 
                                        remove_zero_variance = T,
                                        verbose = F)

mb.splsdacox_model

## -----------------------------------------------------------------------------
lst_models <- list("SB.PLS-ICOX" = sb.plsicox_model,
                   "SB.sPLS-DRCOX" = sb.splsdrcox_model,
                   "MB.sPLS-DRCOX" = mb.splsdrcox_model,
                   "MB.sPLS-DACOX" = mb.splsdacox_model)

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
lst_eval_results$lst_plots$lineplot.mean
lst_eval_results$lst_plot_comparisons$t.test

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test

## -----------------------------------------------------------------------------
lst_models_time <- list(sb.plsicox_model,
                        sb.splsdrcox_model,
                        mb.splsdrcox_model,
                        mb.splsdacox_model,
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.models(lst_models_time)

## ---- fig.small=T-------------------------------------------------------------
ggp_time

## -----------------------------------------------------------------------------
lst_forest_plot <- purrr::map(lst_models, ~survminer::ggforest(.$survival_model$fit, 
                                                               data = .$survival_model$fit$model))

## ---- fig.small=T-------------------------------------------------------------
lst_forest_plot$`SB.sPLS-DRCOX`

## -----------------------------------------------------------------------------
density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")

## ---- fig.small=T-------------------------------------------------------------
density.plots.lp$`SB.sPLS-DRCOX`$plot.density
density.plots.lp$`SB.sPLS-DRCOX`$plot.histogram

## -----------------------------------------------------------------------------
lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)

## ---- fig.small=T-------------------------------------------------------------
lst_ph_ggplot$`SB.sPLS-DRCOX`

## -----------------------------------------------------------------------------
ggp_scores <- plot_PLS_HDcox(model = lst_models$`SB.sPLS-DRCOX`, 
                             comp = c(1,2), mode = "scores")

## ---- fig.small=T-------------------------------------------------------------
ggp_scores$plot_block

## -----------------------------------------------------------------------------
ggp_loadings <- plot_PLS_HDcox(model = lst_models$`SB.sPLS-DRCOX`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10) #length from 0,0

## ---- fig.small=T-------------------------------------------------------------
ggp_loadings$plot_block

## -----------------------------------------------------------------------------
ggp_biplot <- plot_PLS_HDcox(model = lst_models$`SB.sPLS-DRCOX`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T)

## ---- fig.small=T-------------------------------------------------------------
ggp_biplot$plot_block

## -----------------------------------------------------------------------------
ggp.simulated_beta <- plot_pseudobeta.list(lst_models = lst_models, 
                                           error.bar = T, onlySig = T, alpha = 0.05, 
                                           zero.rm = T, auto.limits = T, top = 20,
                                           show_percentage = T, size_percentage = 3)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta$`SB.sPLS-DRCOX`$plot

## -----------------------------------------------------------------------------
LST_KM_RES_LP <- getAutoKM.list(type = "LP",
                                lst_models = lst_models,
                                comp = 1:4,
                                top = 10,
                                ori_data = T,
                                BREAKTIME = NULL,
                                only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_LP$`SB.sPLS-DRCOX`$LST_PLOTS$LP

## -----------------------------------------------------------------------------
LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
                                  lst_models = lst_models,
                                  comp = 1:4,
                                  top = 10,
                                  ori_data = T,
                                  BREAKTIME = NULL,
                                  only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_COMP$`SB.sPLS-DRCOX`$LST_PLOTS$genes$comp_1

## -----------------------------------------------------------------------------
LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
                                 lst_models = lst_models,
                                 comp = 1:4,
                                 top = 10,
                                 ori_data = T,
                                 BREAKTIME = NULL,
                                 only_sig = T, alpha = 0.05)

## ---- fig.small=T-------------------------------------------------------------
LST_KM_RES_VAR$`SB.sPLS-DRCOX`$LST_PLOTS$genes$MOXD1
LST_KM_RES_VAR$`SB.sPLS-DRCOX`$LST_PLOTS$miRNA$hsa_miR_148a

## -----------------------------------------------------------------------------
new_pat <- list()
for(b in names(X_test)){
  new_pat[[b]] <- X_test[[b]][1,,drop=F]
}


## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat$genes),])

## -----------------------------------------------------------------------------
ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient.list(lst_models = lst_models, 
                                                             new_pat = new_pat,
                                                             error.bar = T, onlySig = T, alpha = 0.05,
                                                             zero.rm = T, auto.limits = T, show.betas = T, top = 20)

# ggp.simulated_beta_newPat <- plot_pseudobeta_newPatient(model = lst_models$`MB.sPLS-DACOX`, 
#                                                         new_pat = new_pat,
#                                                         error.bar = T, onlySig = T, alpha = 0.05,
#                                                         zero.rm = T, auto.limits = T, show.betas = T, top = 20)

## ---- fig.small=T-------------------------------------------------------------
ggp.simulated_beta_newPat$`SB.sPLS-DRCOX`$plot

## -----------------------------------------------------------------------------
pat_density <- plot_patient.eventDensity(patient = new_pat, 
                                         time = NULL, 
                                         model = lst_models$`SB.sPLS-DRCOX`, 
                                         type = "lp")

## ---- fig.small=T-------------------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_patient.eventHistogram(patient = new_pat, 
                                             time = NULL, 
                                             model = lst_models$`SB.sPLS-DRCOX`, 
                                             type = "lp")

## ---- fig.small=T-------------------------------------------------------------
pat_histogram

## -----------------------------------------------------------------------------
sub_X_test <- list()
for(b in names(X_test)){
  sub_X_test[[b]] <- X_test[[b]][1:5,]
}


## -----------------------------------------------------------------------------
knitr::kable(Y_test[1:5,])

## -----------------------------------------------------------------------------
lst_cox.comparison <- plot_LP.multiplePatients.list(lst_models = lst_models, 
                                     df.pat = sub_X_test, 
                                     error.bar = T, zero.rm = T, onlySig = T, alpha = 0.05, top = 5)

# lst_cox.comparison <- plot_LP.multiplePatients(model = lst_models$`SB.PLS-ICOX`, 
#                                      df.pat = sub_X_test, 
#                                      error.bar = T, zero.rm = T, onlySig = T, alpha = 0.05, top = 5)

## ---- fig.small=T-------------------------------------------------------------
lst_cox.comparison$`SB.sPLS-DRCOX`$plot

