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
#  devtools::install_github("ConesaLab/HDcox")

## ----setup, eval=FALSE, results = "hide"--------------------------------------
#  # load HDcox
#  library(HDcox)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ConesaLab/RColorConesa")

## -----------------------------------------------------------------------------
library(RColorConesa)
#theme_set(theme_colorConesa()) #under development

## ----load data----------------------------------------------------------------
# load dataset
data("X_miRNA_glioblastoma")
data("Y_miRNA_glioblastoma")

X <- X_miRNA_glioblastoma
Y <- Y_miRNA_glioblastoma

rm(X_miRNA_glioblastoma, Y_miRNA_glioblastoma)

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
X_train <- X[index_train,] #443x534
Y_train <- Y[index_train,]
X_test <- X[-index_train,] #109x534
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

## ---- eval=FALSE, warning=F---------------------------------------------------
#  # run cv.coxEN
#  cv.coxen_res <- cv.coxEN(X = X_train, Y = Y_train,
#                           EN.alpha.list = seq(0.1,1,0.1),
#                           max.variables = ncol(X_train),
#                           n_run = 2, k_folds = 10,
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
#  cv.coxen_res #3min.

## -----------------------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0.9, #cv.coxen_res$opt.EN.alpha
                     max.variables = 22, #cv.coxen_res$opt.nvar
                     x.center = T, x.scale = F,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = F, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## -----------------------------------------------------------------------------
coxen_model <- coxEN(X = X_train, Y = Y_train, 
                     EN.alpha = 0.9, #cv.coxen_res$opt.EN.alpha
                     max.variables = 22, #cv.coxen_res$opt.nvar
                     x.center = T, x.scale = F,
                     remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                     remove_non_significant = T, alpha = 0.05, 
                     MIN_EPV = 5, returnData = T, verbose = F)

## -----------------------------------------------------------------------------
coxen_model

## -----------------------------------------------------------------------------
coxen_model$removed_variables_cox

## ---- eval=FALSE--------------------------------------------------------------
#  # run cv.plsicox
#  cv.splsicox_res <- cv.splsicox(X = X_train, Y = Y_train,
#                                 max.ncomp = 5, spv_penalty.list = seq(0.1,1,0.3),
#                                 n_run = 2, k_folds = 10,
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
#  cv.splsicox_res #8min.

## ---- eval=FALSE, fig.small=T, warning=F--------------------------------------
#  # plot cv.plsicox
#  cv.splsicox_res$plot_AUC

