#Script for Garnatxa cluster

#libraries
library(ggplot2) #done
library(svglite) #done
library(hrbrthemes) #done
library(caret) #done
library(Matrix) #done

library(purrr) #done
library(furrr) #done
library(ggpubr) #done
library(plyr) #done
library(tidyr) #done
library(dplyr) #done
library(survival) #done
library(survminer) #done

library(glmnet) #done
library(mixOmics) #done
library(progress) #done

library(survivalROC) #done
library(smoothROCtime) #done
library(risksetROC) #done
library(nsROC) #done
library(cenROC) #from GitHub #devtools

library(stats) #done
library(RColorConesa) #from GitHub #devtools #usethis and gert install.packages('hrbrthemes', repos='http://cran.us.r-project.org')

#Prepare functions
#Load HDcox functions
loadHDcox <- function(path){
  file <- paste0(path,"HDcox.R")
  source(file, echo = F)
}

#Load ggplot theme
loadGgplotTheme <- function(path){
  file <- paste0(path,"ggplot_theme.R")
  source(file, echo = F)
}

loadFunctions <- function(path){
  loadHDcox(path)
  loadGgplotTheme(path)
}

#path <- "D:/Pedro/Mega/Doctorado/Otros proyectos/"
path = "/home/salguero/Proyectos/"
loadFunctions(path)

#load data
#load("D:/Pedro/Mega/Doctorado/Otros proyectos/data_GBM.mo.RData")
load("data_GBM.mo.RData") #X and Y matrixes

DEBUG = F

if(DEBUG){
  X.expr <- X.expr[,1:(min(50, ncol(X.expr)))] # SELECTING ONLY 2000 GENES
  X.mirna <- X.mirna[,1:(min(55, ncol(X.mirna)))] # SELECTING ONLY 2000 MIRNA
  X.methyl <- X.methyl[,1:(min(60, ncol(X.methyl)))] # SELECTING ONLY 2000 METHYL
}

#Parameters
if(DEBUG){
  n_run = 1
  k_folds = 10
}else{
  n_run = 5
  k_folds = 10
}

fast_mode = F
pred.method = "cenROC"

todaydate <- format(Sys.time(), '%Y-%m-%d')
txt_folder <- paste0("TCGA_MO_TUNE_",ifelse(fast_mode, "FAST_", "COMPLETE_"), pred.method, "_runs_", n_run, "_folds_", k_folds)
folder <- paste0(txt_folder,"_",todaydate,"/")

dir.create(folder)

####################################################
# Set Train and Test - BLOCKS - INTERSECT PATIENTS #
####################################################
Y.full <- Y.expr
Y.full <- rbind(Y.full, Y.mirna[!rownames(Y.mirna) %in% rownames(Y.full),])
Y.full <- rbind(Y.full, Y.methyl[!rownames(Y.methyl) %in% rownames(Y.full),])
Y.full <- rbind(Y.full, Y.clinical[!rownames(Y.clinical) %in% rownames(Y.full),])

X = list(genes = X.expr, miRNA = X.mirna, clinical = X.clinical)#, methyl = X.methyl)
all_X <- purrr::map(X, ~rownames(.))

#rownames 1st omic and in the others, including Y
X1_NONAs = X[[1]]
final_names <- rownames(X1_NONAs)
for(i in 2:length(X)){
  final_names <- intersect(all_X[[i]], final_names)
}
final_names <- intersect(rownames(Y.full), final_names)

#update X and Y
for(i in 1:length(X)){
  X[[i]] = X[[i]][final_names,,drop=F]
}

Y.omics <- Y.full[final_names,,drop=F]

#Set Train and Test
set.seed(123)
index_train.mo <- caret::createDataPartition(Y.omics$event,
                                             p = .7, #70% train
                                             list = FALSE, 
                                             times = 1)[,1]

X_train.mo <- list()
X_test.mo <- list()
for(omic in names(X)){
  X_train.mo[[omic]] <- X[[omic]][index_train.mo,,drop=F]
  X_test.mo[[omic]] <- X[[omic]][-index_train.mo,,drop=F]
}

Y_train.mo <- Y.omics[index_train.mo,]
Y_test.mo <- Y.omics[-index_train.mo,]

#Algorithms Parameters
if(DEBUG){
  max.ncomp = 4
}else{
  max.ncomp = 10
}

x.center = c(genes = T, miRNA = T, clinical = T) #if vector, must be named
x.scale = c(genes = F, miRNA = F, clinical = T) #if vector, must be named
y.center = F
y.scale = F

#Weights Parameters
w_AIC = 0
w_c.index = 0
w_AUC = 1

#Eval stop detection
MIN_AUC_INCREASE = 0.01
MIN_AUC = 0.75
MIN_COMP_TO_CHECK = 3

#others
return_models = F
MIN_EPV = 5
pred.attr = "mean"
seed = 123
alpha = 0.05
remove_non_significant = T #cox and coxEN variables
remove_non_significant_models = F  #plscox methods
remove_near_zero_variance = T
times = NULL
PARALLEL = T
remove_near_zero_variance = T
toKeep.zv = NULL
max.iter = 500

options(future.globals.maxSize = 1600 * 1024^2) #1.6Gb for each worker because we need more memmory due size of X

#######
# COX #
#######
# aux_folder = paste0(folder, "cox_plot/")
# dir.create(aux_folder)
# best_cox <- cox(X = X_train, Y = Y_train, 
#                 x.center = x.center, x.scale = x.scale, 
#                 y.center = y.center, y.scale = y.scale, 
#                 remove_non_significant = remove_non_significant)
# save(list = c("best_cox"), file = paste0(aux_folder, "cox.RData"))

############
# COX - SW #
############
# aux_folder = paste0(folder, "coxSW_plot/")
# dir.create(aux_folder)
# best_coxSW <- coxSW(X = X_train, Y = Y_train,
#                     x.center = x.center, x.scale = x.scale, 
#                     y.center = y.center, y.scale = y.scale,
#                     initialModel = "NULL", keepVariables = NULL, 
#                     SIG_COX = 0.05, SIG_ENT = 0.1, SIG_OUT = 0.15, SIG_PH = 0.05, check_PH = F, 
#                     BACKWARDS = T, verbose = F)
# save(list = c("best_coxSW"), file = paste0(aux_folder, "coxSW.RData"))

###########
# COX-GLM #
###########
# aux_folder = paste0(folder, "coxEN_plot/")
# dir.create(aux_folder)
# 
# cv.coxEN_res <- cv.coxEN(X = X_train, Y = Y_train,
#                          max.variables = ncol(X_train), EN.alpha.list = seq(0,1,0.1), #EN penalization
#                          n_run = n_run, k_folds = k_folds, 
#                          alpha = alpha, remove_non_significant = remove_non_significant, times = times,
#                          w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, 
#                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
#                          x.scale = x.scale, x.center = x.center, 
#                          y.scale = y.scale, y.center = y.center, 
#                          fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
#                          pred.attr = pred.attr, pred.method = pred.method, seed = seed)
# 
# save_ggplot.svg(cv.coxEN_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.coxEN_res_AUC.svg")
# save_ggplot.svg(cv.coxEN_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.coxEN_res_c_index.svg")
# save_ggplot.svg(cv.coxEN_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.coxEN_res_AIC.svg")
# 
# best_coxEN <- coxEN(X = X_train, Y = data.matrix(Y_train), EN.alpha = cv.coxEN_res$opt.EN.alpha, 
#                     x.center = x.center, x.scale = x.scale, 
#                     y.center = y.center, y.scale = y.scale, 
#                     MIN_EPV = MIN_EPV, remove_non_significant = F, alpha = 0.05)
# 
# save(list = c("cv.coxEN_res", "best_coxEN"), file = paste0(aux_folder, "coxEN.RData"))

###########
# PLSRCOX #
###########
# aux_folder = paste0(folder, "plsRcox_plot/")
# dir.create(aux_folder)
# 
# cv.plsRcox_res <- cv.plsRcox(X = X_train, Y = data.matrix(Y_train),
#                           max.ncomp = max.ncomp, 
#                           n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
#                           w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
#                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
#                           x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center, 
#                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
#                           pred.attr = pred.attr, pred.method = pred.method, seed = seed)
# 
# save_ggplot.svg(cv.plsRcox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.plsRcox_res_AUC.svg")
# save_ggplot.svg(cv.plsRcox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.plsRcox_res_c_index.svg")
# save_ggplot.svg(cv.plsRcox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.plsRcox_res_AIC.svg")
# 
# best_plsRcox <- plsRcox(X = X_train, Y = data.matrix(Y_train), 
#                         n.comp = cv.plsRcox_res$opt.comp, 
#                         x.center = x.center, x.scale = x.scale, 
#                         y.center = y.center, y.scale = y.scale)
# 
# save(list = c("cv.plsRcox_res", "best_plsRcox"), file = paste0(aux_folder, "plsRcox.RData"))

##############
# MO-PLSRCOX #
##############
aux_folder = paste0(folder, "mo.sb_plsRcox_plot/")
dir.create(aux_folder)

cv.sb.plsRcox_res <- cv.sb.plsRcox(X = X_train.mo, Y = data.matrix(Y_train.mo),
                                     max.ncomp = max.ncomp,
                                     n_run = n_run, k_folds = k_folds, alpha = alpha,
                                     remove_near_zero_variance = remove_near_zero_variance,
                                     remove_non_significant_models = remove_non_significant_models,
                                     remove_non_significant = remove_non_significant,
                                     w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                                     fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                     pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

save_ggplot.svg(cv.sb.plsRcox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.sb.plsRcox_res_AUC.svg")
save_ggplot.svg(cv.sb.plsRcox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.sb.plsRcox_res_c_index.svg")
save_ggplot.svg(cv.sb.plsRcox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.sb.plsRcox_res_AIC.svg")

best_sb.plsRcox <- sb.plsRcox(X = X_train.mo,
                              Y = data.matrix(Y_train.mo),
                              n.comp = cv.sb.plsRcox_res$opt.comp, #1
                              x.center = x.center, x.scale = x.scale,
                              y.scale = y.scale, y.center = y.center,
                              remove_non_significant = remove_non_significant, #cox variables
                              remove_near_zero_variance = remove_near_zero_variance, toKeep.zv = NULL,
                              returnData = T, verbose = F)

best_sb.fast.plsRcox <- fast.cv.sb.plsRcox(X = X_train.mo, Y = data.matrix(Y_train.mo),
                                           max.ncomp = max.ncomp,
                                           n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                           w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                           x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                           pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                           returnData = T, verbose = T, PARALLEL = PARALLEL)

############
# sPLSRCOX #
############
# aux_folder = paste0(folder, "splsRcox_plot/")
# dir.create(aux_folder)
# 
# cv.splsRcox_genes_res <- cv.splsRcox(X = X_train.mo$genes, Y = data.matrix(Y_train.mo),
#                            max.ncomp = max.ncomp, eta.list = c(0, 0.25, 0.5, 0.75),
#                            n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
#                            w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
#                            MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
#                            x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
#                            fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
#                            pred.attr = pred.attr, pred.method = pred.method, seed = seed)
# 
# save_ggplot.svg(cv.splsRcox_genes_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsRcox_genes_res_AUC.svg")
# save_ggplot.svg(cv.splsRcox_genes_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsRcox_genes_res_c_index.svg")
# save_ggplot.svg(cv.splsRcox_genes_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsRcox_res_genes_AIC.svg")
# 
# best_splsRcox_genes <- splsRcox(X = X_train.mo$genes,
#                           Y = data.matrix(Y_train.mo),
#                           n.comp = cv.splsRcox_genes_res$opt.comp,
#                           eta = cv.splsRcox_genes_res$opt.eta,
#                           x.center = x.center, x.scale = x.scale, 
#                           y.scale = y.scale, y.center = y.center)

###############
# MO-sPLSRCOX #
###############
aux_folder = paste0(folder, "mo.sb_splsRcox_plot/")
dir.create(aux_folder)

cv.sb.splsRcox_res <- cv.sb.splsRcox(X = X_train.mo, Y = data.matrix(Y_train.mo),
                               max.ncomp = max.ncomp, eta.list = c(0, 0.25, 0.5, 0.75),
                               n_run = n_run, k_folds = k_folds, alpha = alpha, 
                               remove_near_zero_variance = remove_near_zero_variance, 
                               remove_non_significant_models = remove_non_significant_models, 
                               remove_non_significant = remove_non_significant,
                               w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                               MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                               x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                               fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                               pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

save_ggplot.svg(cv.sb.splsRcox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.sb.splsRcox_res_AUC.svg")
save_ggplot.svg(cv.sb.splsRcox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.sb.splsRcox_res_c_index.svg")
save_ggplot.svg(cv.sb.splsRcox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.sb.splsRcox_res_AIC.svg")

best_sb.splsRcox <- sb.splsRcox(X = X_train.mo,
                                Y = data.matrix(Y_train.mo),
                                n.comp = cv.sb.splsRcox_res$opt.comp, #2
                                eta = cv.sb.splsRcox_res$opt.eta, #0.5
                                x.center = x.center, x.scale = x.scale,
                                y.scale = y.scale, y.center = y.center,
                                remove_non_significant = remove_non_significant, #cox variables
                                remove_near_zero_variance = remove_near_zero_variance, toKeep.zv = NULL,
                                returnData = T, verbose = F)

best_sb.fast.splsRcox <- fast.cv.sb.splsRcox(X = X_train.mo, Y = data.matrix(Y_train.mo),
                                          max.ncomp = max.ncomp, eta.list = c(0, 0.25, 0.5, 0.75),
                                          n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                          w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                          MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                          x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                                          fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                          pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                          returnData = T, verbose = T, PARALLEL = PARALLEL)

#####################
# sPLSRCOX_mixOmics #
#####################
# aux_folder = paste0(folder, "splsRcox_mixOmics_plot/")
# dir.create(aux_folder)
# 
# test.keepX <- NULL
# max <- ncol(X_train)
# min <- 5
# if(max > 20){
#   test.keepX <- c(test.keepX, seq(1, max, ceiling(max/20)))
# }else{
#   test.keepX <- c(test.keepX, seq(1, max, 1))
# }
# 
# cv.splsRcox_mixOmics_res <- cv.splsRcox_mixOmics(X = X_train, Y = data.matrix(Y_train),
#                                              max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds, 
#                                              n_run_mixomics = 3, k_folds_mixomics = 4, test.keepX = test.keepX,
#                                              alpha = alpha, remove_non_significant_models = remove_non_significant_models,
#                                              w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
#                                              MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
#                                              x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
#                                              fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
#                                              pred.attr = pred.attr, pred.method = pred.method, seed = seed)
# 
# save_ggplot.svg(cv.splsRcox_mixOmics_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsRcox_mixOmics_res_AUC.svg")
# save_ggplot.svg(cv.splsRcox_mixOmics_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsRcox_mixOmics_res_c_index.svg")
# save_ggplot.svg(cv.splsRcox_mixOmics_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsRcox_mixOmics_res_AIC.svg")
# 
# best_splsRcox_mixOmics <- splsRcox_mixOmics(X = X_train,
#                                             Y = data.matrix(Y_train),
#                                             n.comp = cv.splsRcox_mixOmics_res$opt.comp,
#                                             x.center = x.center, x.scale = x.scale, 
#                                             y.scale = y.scale, y.center = y.center,
#                                             n_run_mixomics = 5, k_folds_mixomics = 10, test.keepX = test.keepX)

########################
# MO-sPLSRCOX-MixOmics #
########################

aux_folder = paste0(folder, "mo.mb_splsRcox_mixOmics_plot/")
dir.create(aux_folder)

cv.mb.splsRcox_mixOmics_res <- cv.mb.splsRcox_mixOmics(X = X_train.mo, Y = data.matrix(Y_train.mo),
                                               max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds, 
                                               vector = NULL, #c(5, 30, 100, 200),
                                               alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                               w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
                                               MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                               x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                                               fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                               pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

save_ggplot.svg(cv.mb.splsRcox_mixOmics_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.mb.splsRcox_mixOmics_res_AUC.svg")
save_ggplot.svg(cv.mb.splsRcox_mixOmics_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.mb.splsRcox_mixOmics_res_c_index.svg")
save_ggplot.svg(cv.mb.splsRcox_mixOmics_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.mb.splsRcox_mixOmics_res_AIC.svg")

best_mb.splsRcox_mixOmics <- mb.splsRcox_mixOmics(X = X_train.mo,
                                         Y = data.matrix(Y_train.mo),
                                         n.comp = cv.mb.splsRcox_mixOmics_res$opt.comp,
                                         vector = cv.mb.splsRcox_mixOmics_res$opt.nvar,
                                         x.center = x.center, x.scale = x.scale, 
                                         y.scale = y.scale, y.center = y.center,
                                         remove_non_significant = remove_non_significant, #cox variables
                                         remove_near_zero_variance = remove_near_zero_variance, toKeep.zv = NULL,
                                         returnData = T, verbose = T)

#CHECK PREDICTIONS
# X_test_mod <- predict.Xscores.HDcox(model = best_sb.splsRcox, newdata = X_test.mo)
# X_test_mod1 <- predict.Xscores.HDcox(model = best_sb.splsRcox$list_spls_models$genes, newdata = scale(X_test.mo$genes, center = best_sb.splsRcox$X$x.mean$genes, scale = F))
# 
# nd_scaled <- list()
# for(b in names(best_mb.splsRcox_mixOmics$X$data)){
#   nd_scaled[[b]] <- scale(X_test.mo[[b]], center = best_mb.splsRcox_mixOmics$X$x.mean[[b]], scale = F)
# }
# X_test_mod <- predict.Xscores.HDcox(model = best_mb.splsRcox_mixOmics, newdata = X_test.mo)
# X_test_mod1 <- predict(best_mb.splsRcox_mixOmics$mb.model, newdata = nd_scaled)

#############
# plsDA-COX #
#############

# #RAW DATA AS INPUT
# aux_folder = paste0(folder, "plsdacox_mixOmics_plot/")
# dir.create(aux_folder)
# 
# cv.plsdacox_mixOmics_res <- cv.plsdacox(X = X_train, Y = data.matrix(Y_train),
#                            max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds, 
#                            alpha = alpha, remove_non_significant_models = remove_non_significant_models, max.iter = max.iter,
#                            w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times,
#                            MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
#                            x.scale = x.scale, x.center = x.center, 
#                            y.scale = y.scale, y.center = y.center,
#                            fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
#                            pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)
# 
# save_ggplot.svg(cv.plsdacox_mixOmics_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.plsdacox_res_AUC.svg")
# save_ggplot.svg(cv.plsdacox_mixOmics_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.plsdacox_res_c_index.svg")
# save_ggplot.svg(cv.plsdacox_mixOmics_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.plsdacox_res_AIC.svg")
# 
# best_plsdacox_mixOmics <- plsdacox(X_train, Y_train, 
#                           n.comp = cv.plsdacox_mixOmics_res$opt.comp, 
#                           x.center = x.center, x.scale = x.scale, 
#                           y.center = y.center, y.scale = y.scale, 
#                           max.iter = max.iter)

################
# MB-sPLSDACOX #
################

aux_folder = paste0(folder, "mo.mb_splsda_mixOmics_plot/")
dir.create(aux_folder)

cv.mb.splsdacox_res <- cv.mb.splsda_mixOmics(X = X_train.mo, Y = data.matrix(Y_train.mo),
                                             max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds, 
                                             vector = NULL, #c(5, 30, 100, 200),
                                             alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                             w_AIC = w_AIC, w_c.index = w_c.index, w_AUC = w_AUC, times = times, max.iter = max.iter,
                                             MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                             x.scale = x.scale, x.center = x.center, y.scale = y.scale, y.center = y.center,
                                             fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                             pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

save_ggplot.svg(cv.mb.splsdacox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.mb.plsdacox_res_AUC.svg")
save_ggplot.svg(cv.mb.splsdacox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.mb.plsdacox_res_c_index.svg")
save_ggplot.svg(cv.mb.splsdacox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.mb.plsdacox_res_AIC.svg")

best_mb.splsda_mixOmics <- mb.splsdacox_mixOmics(X = X_train.mo,
                                                 Y = data.matrix(Y_train.mo),
                                                 n.comp = cv.mb.splsdacox_res$opt.comp, #8
                                                 vector = cv.mb.splsdacox_res$opt.nvar,
                                                 x.center = x.center, x.scale = x.scale, 
                                                 y.scale = y.scale, y.center = y.center,
                                                 remove_non_significant = remove_non_significant, #cox variables
                                                 remove_near_zero_variance = remove_near_zero_variance, toKeep.zv = NULL,
                                                 returnData = T, verbose = T)

###############
# SAVE MODELS #
###############
lst_models_full <- list("cv.sb.plsRcox_res" = cv.sb.plsRcox_res,
                   "best_sb.plsRcox" = best_sb.plsRcox,
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox,
                   "cv.sb.splsRcox_res" = cv.sb.splsRcox_res,
                   "best_sb.splsRcox" = best_sb.splsRcox,
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox,
                   "cv.mb.splsRcox_mixOmics_res" = cv.mb.splsRcox_mixOmics_res,
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics,
                   "cv.mb.splsdacox_res" = cv.mb.splsdacox_res,
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics)

lst_models <- list("best_sb.plsRcox" = best_sb.plsRcox,
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox,
                   "best_sb.splsRcox" = best_sb.splsRcox,
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox,
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics,
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics)

save(list = c("lst_models_full"), file = paste0(folder, "lst_models_full.RData"))
save(list = c("lst_models"), file = paste0(folder, "lst_models.RData"))

##############
# EVALUATION #
##############
lst_evaluations <- c("survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I", "risksetROC")
names(lst_evaluations) <- lst_evaluations

# lst_models <- list("sb.splsRcox" = best_sb.splsRcox)

lst_models <- list("sb.plsRcox" = best_sb.plsRcox,
                   "sb.fast.plsRcox" = best_sb.fast.plsRcox,
                   "sb.splsRcox" = best_sb.splsRcox,
                   "sb.fast.splsRcox" = best_sb.fast.splsRcox,
                   "mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics,
                   "mb.splsda_mixOmics" = best_mb.splsda_mixOmics)

eval_results <- purrr::map(lst_evaluations, ~eval_models4.0(lst_models = lst_models,
                                                            X_test = X_test.mo, Y_test = Y_test.mo, pred.method = .,
                                                            pred.attr = pred.attr, 
                                                            times = times, max_time_points = 15, PARALLEL = PARALLEL))

#############
# EVALPLOTS #
#############
evaluation_folder = paste0(folder, "evaluation_plot/")
dir.create(evaluation_folder)

lst_ggp <- plot.multiple.evaluations(eval_results)

save_ggplot_lst.svg(lst_plots = lst_ggp, object_name = "lineplot.mean", folder = evaluation_folder, wide = T, prefix = "eval_")

for(eval_name in names(lst_ggp)){
  for(comp_name in names(lst_ggp[[eval_name]]$comparisons)){
    #save_ggplot.svg(plot = lst_ggp[[eval_name]]$comparisons[[comp_name]], name = paste0("comparison_", eval_name, "_", comp_name), folder = evaluation_folder, wide = T)
    save_ggplot.svg(plot = lst_ggp[[eval_name]]$comparisons[[comp_name]], name = paste0("comparison_", eval_name, "_", comp_name), folder = evaluation_folder, wide = T)
  }
}

save(list = c("eval_results", "lst_ggp"), file = paste0(evaluation_folder, "eval_results.RData"))

########
# TIME #
########
time_folder = paste0(folder, "time_plot/")
dir.create(time_folder)

lst_models <- list("cv.sb.plsRcox_res" = cv.sb.plsRcox_res,
                   "best_sb.plsRcox" = best_sb.plsRcox,
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox,
                   "cv.sb.splsRcox_res" = cv.sb.splsRcox_res,
                   "best_sb.splsRcox" = best_sb.splsRcox,
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox,
                   "cv.mb.splsRcox_mixOmics_res" = cv.mb.splsRcox_mixOmics_res,
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics,
                   "cv.mb.splsdacox_res" = cv.mb.splsdacox_res,
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics,
                   "eval_results" = eval_results)

ggp_time <- models.time.plot(lst_models)

save_ggplot.svg(ggp_time, folder = time_folder, name = "time.svg", wide = T)
save_ggplot.svg(ggp_time, folder = time_folder, name = "time.svg", wide = T)

save(list = c("lst_models"), file = paste0(time_folder, "time_results.RData"))
gc()

##################
## DENSITY PLOTS #
##################

density_folder = paste0(folder, "density_plot/")
dir.create(density_folder)

lst_models <- list("best_sb.plsRcox" = best_sb.plsRcox, 
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox, 
                   "best_sb.splsRcox" = best_sb.splsRcox, 
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox, 
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics, 
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics)

density.plots.lp <- purrr::map(lst_models, ~eventDensityByLP.plot(., type = "lp"))
density.plots.risk <- purrr::map(lst_models, ~eventDensityByLP.plot(., type = "risk"))
density.plots.expected <- purrr::map(lst_models, ~eventDensityByLP.plot(., type = "expected"))
density.plots.survival <- purrr::map(lst_models, ~eventDensityByLP.plot(., type = "survival"))

for(cn in names(density.plots.lp)){
  
  if(cn %in% names(density.plots.lp) && !is.na(density.plots.lp[[cn]])){
    save_ggplot.svg(plot = density.plots.lp[[cn]]$plot.density, folder = density_folder, name = paste0("pred.lp_density_",cn))  
    save_ggplot.svg(plot = density.plots.lp[[cn]]$plot.histogram, folder = density_folder, name = paste0("pred.lp_histogram_",cn))
  }
  
  if(cn %in% names(density.plots.risk) && !is.na(density.plots.risk[[cn]])){
    save_ggplot.svg(plot = density.plots.risk[[cn]]$plot.density, folder = density_folder, name = paste0("pred.risk_density_",cn))  
    save_ggplot.svg(plot = density.plots.risk[[cn]]$plot.histogram, folder = density_folder, name = paste0("pred.risk_histogram_",cn))
  }
  
  if(cn %in% names(density.plots.expected) && !is.na(density.plots.expected[[cn]])){
    save_ggplot.svg(plot = density.plots.expected[[cn]]$plot.density, folder = density_folder, name = paste0("pred.expected_density_",cn))  
    save_ggplot.svg(plot = density.plots.expected[[cn]]$plot.histogram, folder = density_folder, name = paste0("pred.expected_histogram_",cn))
  }
  
  if(cn %in% names(density.plots.survival) && !is.na(density.plots.survival[[cn]])){
    save_ggplot.svg(plot = density.plots.survival[[cn]]$plot.density, folder = density_folder, name = paste0("pred.expected_density_",cn))  
    save_ggplot.svg(plot = density.plots.survival[[cn]]$plot.histogram, folder = density_folder, name = paste0("pred.expected_histogram_",cn))
  }
}

save(list = c("lst_models"), file = paste0(density_folder, "density_results.RData"))

#################
## FOREST PLOTS #
#################
forest_folder = paste0(folder, "forest_plot/")
dir.create(forest_folder)

lst_forest_plot <- purrr::map(lst_models, ~survminer::ggforest(.$survival_model$fit, data = .$survival_model$fit$model))
save_ggplot_lst.svg(lst_plots = lst_forest_plot, folder = forest_folder, wide = T, prefix = "forest_")

####################
# DIAGNOSTIC PLOTS #
####################
ph_folder = paste0(folder, "ph_plot/")
dir.create(ph_folder)

lst_ph_preplot <- purrr::map(lst_models, ~survival::cox.zph(.$survival_model$fit))
lst_ph_plot <- purrr::map(lst_ph_preplot, ~survminer::ggcoxzph(.))
lst_ph_ggplot <- purrr::map2(.x = lst_ph_preplot, .y = lst_ph_plot, ~ggcoxzph2ggplot(.x, .y))

save_ggplot_lst.svg(lst_plots = lst_ph_ggplot, folder = ph_folder, wide = T, prefix = "ph_")

################
## EVENT PLOTS #
################
event_folder = paste0(folder, "event_plot/")
dir.create(event_folder)

ggp_density.event <- density.event.plot(Y = Y_train.mo, round = T, roundTo = 50)
save_ggplot.svg(plot = ggp_density.event, folder = event_folder, name = "density_events")

###################
## PseudoBeta COX #
###################

psbeta_folder = paste0(folder, "pbetacox_plot/")
dir.create(psbeta_folder)

lst_models <- list("best_sb.plsRcox" = best_sb.plsRcox, 
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox, 
                   "best_sb.splsRcox" = best_sb.splsRcox, 
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox, 
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics, 
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics)

ggp.simulated_beta <- purrr::map(lst_models, ~cox.pls.variable.plot(model = ., error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T, auto.limits = T, top = 20))

for(m in names(ggp.simulated_beta)){
  for(b in names(ggp.simulated_beta[[m]]$plot)){
    save_ggplot.svg(plot = ggp.simulated_beta[[m]]$plot[[b]], folder = psbeta_folder, wide = T, name = paste0("pbetacox_", m, "_", b, "_"))
  }
}

################
# Kaplan-Meier #
################

km_folder = paste0(folder, "km_plot/")
dir.create(km_folder)

lst_models <- list("best_sb.plsRcox" = best_sb.plsRcox, 
                   "best_sb.fast.plsRcox" = best_sb.fast.plsRcox, 
                   "best_sb.splsRcox" = best_sb.splsRcox, 
                   "best_sb.fast.splsRcox" = best_sb.fast.splsRcox, 
                   "best_mb.splsRcox_mixOmics" = best_mb.splsRcox_mixOmics, 
                   "best_mb.splsda_mixOmics" = best_mb.splsda_mixOmics)
  

LST_KM_RES <- purrr::map(lst_models, ~getAutoKM(model = .,
                                                comp = 1:max.ncomp,
                                                top = 10,
                                                ori_data = T,
                                                BREAKTIME = NULL,
                                                only_sig = T, alpha = 0.05))

for(m in names(LST_KM_RES)){
  for(b in names(LST_KM_RES[[m]]$LST_PLOTS)){
    save_ggplot_lst.svg(lst_plots = LST_KM_RES[[m]]$LST_PLOTS[[b]], object_name = NULL, folder = km_folder, wide = T, prefix = paste0("km_", m, "_", b, "_"))
  }
}

#################
## SAVE RESULTS #
#################

##Save all results
save.image(paste0(folder,"output.performance_",todaydate,".RData"))

message("The script has successfully finished!")