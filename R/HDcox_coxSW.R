#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' coxSW
#' @description Performs a standard cox stepwise model (based on My.stepwise R package).
#'
#' @param X Numeric matrix. Predictor variables
#' @param Y Numeric matrix. Response variables. It assumes it has two columns named as "time" and "event". For event column, values can be 0/1 or FALSE/TRUE for censored and event samples.
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param y.center Logical. If y.center = TRUE, Y matrix is centered to zero means (default: FALSE).
#' @param y.scale Logical. If y.scale = TRUE, Y matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, remove_near_zero_variance variables will be removed.
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, remove_zero_variance variables will be removed.
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance filtering.
#' @param initialModel Character vector. Name of variables in X to include in the initial model (default: "NULL").
#' @param toKeep.sw Character vector. Name of variables in X to not be deleted by Step-wise selection.
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If MIN_EPV is not meet, the value will be change automatically.
#' @param alpha Numeric. Cutoff for establish significant variables. Below the number are considered as significant (default: 0.05).
#' @param alpha_ENT Numeric. Maximum P-Value for a variable to enter the model (default: 0.10).
#' @param alpha_OUT Numeric. Minimum P-Value for a variable to leave the model (default: 0.15).
#' @param alpha_PH Numeric. Maximum P-Value for a variable to meet Proportional Hazard Assumption (default: 0.05).
#' @param check_PH Logical. If check_PH = TRUE, PH Assumption will be check (default: FALSE).
#' @param boostDeletion Logical. If boostDeletion = TRUE, all variables with a P-Value greater than alpha_OUT will be deleted in each step. Instead of deleting one by one (default: FALSE).
#' @param BACKWARDS Logical. If BACKWARDS = TRUE, backward strategy is perform.
#' @param MIN_EPV Minimum number of Events Per Variable you want reach for the final cox model. Used to restrict the number of variables can appear in cox model. If the minimum is not meet, the model is not computed.
#' @param returnData Logical. Return original and normalized X and Y matrices.
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "HDcox" and model "coxSW". The class contains the following elements:
#'
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(data)}: normalized Y matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix
#'  }
#' \code{survival_model}: List of survival model information
#' \itemize{
#'  \item \code{fit}: coxph object.
#'  \item \code{AIC}: AIC of cox model.
#'  \item \code{BIC}: BIC of cox model.
#'  \item \code{lp}: linear predictors for train data.
#'  \item \code{coef}: Coefficients for cox model.
#'  \item \code{YChapeau}: Y Chapeau residuals.
#'  \item \code{Yresidus}: Y residuals.
#' }
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @export

coxSW <- function(X, Y,
                  x.center = TRUE, x.scale = FALSE,
                  y.center = FALSE, y.scale = FALSE,
                  remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                  remove_non_significant = F,
                  initialModel = "NULL", toKeep.sw = NULL, max.variables = 20,
                  alpha = 0.05, alpha_ENT = 0.1, alpha_OUT = 0.15, alpha_PH = 0.05, check_PH = F, boostDeletion = F, BACKWARDS = F,
                  MIN_EPV = 5, returnData = T, verbose = F){

  #DFCALLS
  REAL_INTERACTIONS <- symbol <- path.saveImage <- NULL


  #function updated from My.stepwise package
  t1 <- Sys.time()

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### REQUIREMENTS - X to MATRIX
  #### in this case we want to keep it as data.frame
  # lst_check <- checkXY.class(X, Y, verbose = verbose)
  # X <- lst_check$X
  # Y <- lst_check$Y

  if(class(X)[1]!="data.frame" & class(X)[1]=="matrix"){
    X <- as.data.frame(X)
  }
  if(class(Y)[1]!="data.frame" & class(Y)[1]=="matrix"){
    Y <- as.data.frame(Y)
  }

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                          remove_near_zero_variance = remove_near_zero_variance,
                                          remove_zero_variance = remove_zero_variance,
                                          toKeep.zv = toKeep.zv,
                                          freqCut = 95/5)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  max.variables <- check.ncomp(X, max.variables, verbose = verbose)

  #### MAX PREDICTORS
  max.variables <- check.maxPredictors(X, Y, MIN_EPV, max.variables, verbose = verbose)

  #### SCALING
  lst_scale <- XY.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh
  data <- cbind(Xh, Yh)

  #### INITIALISING VARIABLES
  data.all <- NULL
  oneToDelete <- c("")
  allVariablesDeletedSW <- c() #must be returned
  lstMeetAssumption = NULL
  cont = 1
  problem_flag = F

  while(!is.null(oneToDelete)){
    ########
    # CONT #
    ########
    if(verbose){
      message(paste0("Iteration: ", cont))
    }
    cont <- cont+1

    ###################
    # INTERESTING VAR #
    ###################
    in.variable = "NULL"

    #######
    #MODEL#
    #######
    model <- tryCatch(
      # Specifying expression
      expr = {
        coxphSW(time = "time", status = "event", variable.list = NULL, in.variable = in.variable, data = data, max.variables = max.variables,
                   sle = alpha_ENT, sls = alpha_OUT, startBackwards = BACKWARDS, REAL_INTERACTIONS = REAL_INTERACTIONS, symbol=symbol,
                   verbose = verbose)
      },
      # Specifying error message
      error = function(e){
        message(paste0("coxSW: ", e))
        invisible(gc())
        return(NA)
      }
    )

    if(!all(is.na(model))){
      variables.df <- as.data.frame(summary(model)[7])[,5,drop=F]
      variables.df$variables <- rownames(variables.df)
      variables.df <- variables.df[,c(2,1)]
      colnames(variables.df) <- c("variables","Pr(>|z|)")

      coxph.sw <- model
    }else{
      problem_flag = T
      break
    }

    ####################
    # CHECK PH AND SIG #
    ####################
    anyProblem = FALSE

    #############
    # ANY NO PH #
    #############
    if(check_PH){
      fit_zph <- survival::cox.zph(coxph.sw, transform = "km", terms = T, global = F) #do not check global value
      fit_zph.aux <- fit_zph
      #any significant?
      if(any(fit_zph$table[,"p"]<=alpha_PH)){
        variablesNotHorizontal <- rownames(fit_zph$table)[which(fit_zph$table[,"p",drop=F]<=alpha_PH)]
        if(check_PH){
          #for each variable with NPH
          for(oneToDelete in variablesNotHorizontal){
            oneToDelete <- gsub('[\t\`]', '', oneToDelete)
            oneToDelete.raw <- names(model$coefficients)[which(startsWith(names(model$coefficients), oneToDelete))]
            if(length(oneToDelete.raw)>1){
              oneToDelete.raw <- names(which.max(sapply(oneToDelete.raw, nchar))) #if more than one word, just take it the longer match
            }

            oneToDelete <- colnames(data)[which(startsWith(oneToDelete, colnames(data)))] #get real name if factor: update factor variables
            if(length(oneToDelete)>1){
              oneToDelete <- names(which.max(sapply(oneToDelete, nchar))) #if more than one word, just take it the longer match
            }
            if(oneToDelete %in% lstMeetAssumption) #if accepted manually next
              next

            #Get plot
            # grDevices::dev.new()
            # plot(fit_zph[which(rownames(fit_zph$table)==oneToDelete)]) #básica
            #ggcoxzph(fit_zph[which(rownames(fit_zph$table)==oneToDelete)]) #de un paquete

            #creando mi propia linea de fit_zph
            oneToDelete.aux <- oneToDelete
            zph_plot <- plotZPH(fit_zph, oneToDelete)
            #zph_plot <- plotZPH(fit_zph.aux, oneToDelete.aux)
            #Get plot

            grDevices::dev.new(width = 12, height = 12)
            print(zph_plot)

            #ask to the user
            question1 <- readline(paste0("Do you consider the PH Assumption is broken for ", oneToDelete, " with a p-value of ", round(fit_zph$table[oneToDelete,"p"],4),
                                         " (beta of ", round(model$coefficients[oneToDelete.raw],4) ,")? (Y/N): "))
            question1 <- tolower(question1)
            while(!question1 %in% c("n", "no", "yes", "y")){
              question1 <- readline(paste0('You have to enter a n/no or y/yes character: '))
            }

            if(regexpr(question1, 'n', ignore.case = TRUE) == 1){
              grDevices::dev.off() #close window
              lstMeetAssumption <- c(lstMeetAssumption, oneToDelete) #add variable to list


              if(!is.null(path.saveImage)){
                grDevices::png(filename=paste0(path.saveImage,"_", oneToDelete, "_OK.png"), units = "in", res = 100, width = 12, height = 12)
                print(zph_plot)
                grDevices::dev.off()
              }

              next #if everything it is ok, check next variable
            } else if (regexpr(question1, 'y', ignore.case = TRUE) == 1){
              grDevices::dev.off() #close window

              if(!is.null(path.saveImage)){
                grDevices::png(filename=paste0(path.saveImage,"_", oneToDelete, "_BAD.png"), units = "in", res = 100, width = 12, height = 12)
                print(zph_plot)
                grDevices::dev.off()
              }

              #delete variable
              if(verbose){
                print(paste0("Deleting variable ", oneToDelete, ": No PH in Cox Model"))
                message(c("Deleting variable ", oneToDelete, ": No PH in Cox Model"))
              }
              allVariablesDeletedSW <- rbind(allVariablesDeletedSW, c(oneToDelete, "No-PH", fit_zph$table[variablesNotHorizontal,"p"]))

              #sometimes the interaction is founded in other direction
              lstDeleted <- deleteVariablesCox(data, td = oneToDelete, interactions = F)
              data <- lstDeleted$cox.train
              anyProblem = TRUE
              break #just one variable
            }
          } #for PH
        }
      } #PH SIG
    }

    ##############
    # PVAL > SIG #
    ##############

    #if pval == NA, near zero variance problem, assing p.val = 1
    if(length(variables.df[is.na(variables.df$`Pr(>|z|)`),"Pr(>|z|)"])>0){
      variables.df[is.na(variables.df$`Pr(>|z|)`),"Pr(>|z|)"] <- 1
    }

    if(!anyProblem & any(variables.df$`Pr(>|z|)`>alpha)){
      lstHighPVAL <- rownames(variables.df)[variables.df$`Pr(>|z|)`>alpha]

      #do not delete toKeep.sw - only if just one with that prefix length(grep)==1
      #but they are categorical and we will find age501 or sex1
      minusLastChart <- sapply(lstHighPVAL, function(x){substr(x, 0, nchar(x)-1)})
      if(any(toKeep.sw %in% minusLastChart)){ #delete just last char

        #check if grep > 1: T then delete variable for candidates
        for(oneKeep in toKeep.sw[toKeep.sw %in% minusLastChart]){
          if(length(grep(oneKeep, minusLastChart))>1)
            lstHighPVAL <- lstHighPVAL[!minusLastChart %in% oneKeep]
        }
      }

      if(!length(lstHighPVAL)==0){ #if any

        if(boostDeletion){
          index <- which(variables.df[lstHighPVAL,]$`Pr(>|z|)` == max(variables.df[lstHighPVAL,]$`Pr(>|z|)`)) #get max p-val
        }else{
          index <- which.max(variables.df[lstHighPVAL,]$`Pr(>|z|)`) #get max p-val
        }

        oneToDelete <- lstHighPVAL[index]
        oneToDelete <- gsub('[\t\`]', '', oneToDelete)

        getPval <- variables.df[variables.df$variables %in% oneToDelete,]$`Pr(>|z|)`
        if(length(getPval)==0){
          getPval <- variables.df[variables.df$variables %in% paste0("`",oneToDelete,"`"),]$`Pr(>|z|)`
        }

        #print(paste0("Deleting variable ", oneToDelete, ": No Significative in Cox Model (", round(getPval,4), ")\n"))
        if(verbose){
          message(paste0("Deleting variable '", oneToDelete, "': No Significative in Cox Model (", round(getPval,4), ")\n"))
        }
        allVariablesDeletedSW <- rbind(allVariablesDeletedSW, cbind(oneToDelete, rep("No-SIG", length(oneToDelete)), getPval))

        #get real name
        if(verbose){
          print(oneToDelete)
        }

        oneToDelete <- getRealNameVarFromCoxModel(data, oneToDelete)
        oneToDelete <- gsub('[\t\`]', '', oneToDelete)

        lstDeleted <- deleteVariablesCox(data, td = oneToDelete, interactions = F)

        data <- lstDeleted$cox.train
        anyProblem = TRUE
      }
    }

    ##########################
    # NUMBER OF VARIABLES >> #
    ##########################

    if(max.variables < length(coxph.sw$coefficients) & !anyProblem){
      lstHighPVAL <- rownames(variables.df)

      #do not delete toKeep.sw - only if just one with that prefix length(grep)==1
      #but they are categorical and we will find age501 or sex1
      minusLastChart <- sapply(lstHighPVAL, function(x){substr(x, 0, nchar(x)-1)})
      if(any(toKeep.sw %in% minusLastChart)){ #delete just last char

        #check if grep > 1: T then delete variable for candidates
        for(oneKeep in toKeep.sw[toKeep.sw %in% minusLastChart]){
          if(length(grep(oneKeep, minusLastChart))>1)
            lstHighPVAL <- lstHighPVAL[!minusLastChart %in% oneKeep]
        }
      }

      if(!length(lstHighPVAL)==0){ #if any

        if(boostDeletion){
          index <- which(variables.df[lstHighPVAL,]$`Pr(>|z|)` == max(variables.df[lstHighPVAL,]$`Pr(>|z|)`)) #get max p-val
        }else{
          index <- which.max(variables.df[lstHighPVAL,]$`Pr(>|z|)`) #get max p-val
        }

        ##################################################################
        # IF index == all dataset, keep the maximum of variables allowed #
        ##################################################################
        #in other case, just delete higher p-value variables
        if(length(index) == length(lstHighPVAL)){
          index <- index[(max.variables+1):length(index)] #delete latest X variables
        }

        oneToDelete <- lstHighPVAL[index]
        oneToDelete <- gsub('[\t\`]', '', oneToDelete)

        getPval <- variables.df[variables.df$variables %in% oneToDelete,]$`Pr(>|z|)`
        if(length(getPval)==0){
          getPval <- variables.df[variables.df$variables %in% paste0("`",oneToDelete,"`"),]$`Pr(>|z|)`
        }

        #print(paste0("Deleting variable ", oneToDelete, ": No Significative in Cox Model (", round(getPval,4), ")\n"))
        if(verbose){
          message(paste0("Deleting variable '", oneToDelete, "': Too many variables in Cox Model (", round(getPval,4), ")\n"))
        }
        allVariablesDeletedSW <- rbind(allVariablesDeletedSW, cbind(oneToDelete, rep("To many variables", length(oneToDelete)), getPval))
        #sometimes the interaction is founded in other direction: Now interactions are normal variables

        #get real name
        if(verbose){
          print(oneToDelete)
        }
        oneToDelete <- getRealNameVarFromCoxModel(data, oneToDelete)
        oneToDelete <- gsub('[\t\`]', '', oneToDelete)

        lstDeleted <- deleteVariablesCox(data, td = oneToDelete, interactions = F)

        data <- lstDeleted$cox.train
        anyProblem = TRUE
      }
    }

    ##################
    # IF NO PROBLEMS #
    ##################

    if(!anyProblem){
      oneToDelete = NULL
      data.all <- data #must be returned
      if(verbose){
        message("Final model Reached!")
      }
    }

  } #while

  if(!problem_flag){
    allVariablesDeletedSW <- as.data.frame(allVariablesDeletedSW)

    if(nrow(allVariablesDeletedSW)>0){
      colnames(allVariablesDeletedSW) <- c("Variable", "Reason", "P-Value")
      rownames(allVariablesDeletedSW) <- allVariablesDeletedSW$Variable
      allVariablesDeletedSW$Variable = NULL
    }
  }else{
    coxph.sw <- NULL
    data.all <- NA #no model was created, so we do not have which variables enter
  }

  removed_variables <- NULL
  if(class(coxph.sw)=="coxph"){

    #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
    if(remove_non_significant){
      lst_rnsc <- removeNonSignificativeCox(cox = coxph.sw, alpha = alpha, cox_input = d)

      coxph.sw <- lst_rnsc$cox
      removed_variables <- lst_rnsc$removed_variables
    }

    survival_model <- getInfoCoxModel(coxph.sw)
  }else{
    survival_model <- NULL
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  invisible(gc())
  return(coxSW_class(list(X = list("data" = if(returnData) data.all else NA, "x.mean" = xmeans, "x.sd" = xsds),
                        Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                        survival_model = survival_model,
                        call = func_call,
                        X_input = if(returnData) X_original else NA,
                        Y_input = if(returnData) Y_original else NA,
                        alpha = alpha,
                        removed_variables_cox = removed_variables,
                        nzv = variablesDeleted,
                        class = pkg.env$coxSW,
                        time = time)))
}

#### ### ##
# METHODS #
#### ### ##

coxphSW <- function(time = "time", status = "event", variable.list = NULL, in.variable, data, max.variables = NULL, sle = 0.10, sls = 0.15, REAL_INTERACTIONS = T, symbol = ":", startBackwards = F, verbose = F){

  if(is.null(variable.list)){
    variable.list <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
  }

  if(verbose){
    message("VERBOSE MODE")
  }

  start_time <- Sys.time()
  finalModel <- stepwise.coxph(Time = time, Status = status, variable.list = variable.list, in.variable = in.variable, data = data, max.variables = max.variables, sle = sle, sls = sls, startBackwards = startBackwards, REAL_INTERACTIONS = REAL_INTERACTIONS, verbose = verbose)
  end_time <- Sys.time()

  if(verbose){
    print(end_time - start_time)
  }

  return(finalModel)
}

stepwise.coxph <- function (Time = NULL, T1 = NULL, T2 = NULL, Status = NULL, variable.list,
                            in.variable = "NULL", data, max.variables = NULL, sle = 0.15, sls = 0.15, REAL_INTERACTIONS = T, symbol = ":", startBackwards = F, verbose = F){
  univar.pvalue <- NULL
  temp.model <- NULL

  if (is.null(T2)) {

    if(startBackwards){
      in.variable <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
      if(!is.null(max.variables) & length(in.variable) > max.variables){
        set.seed(123)
        subcox <- cox(X = data[,colnames(data) %in% c(in.variable), drop=F],
                      Y = data[,colnames(data) %in% c("time", "event", "status")],
                      remove_non_significant = F, MIN_EPV = 3, FORCE = T)

        vars_aux <- summary(subcox$survival_model$fit)[[7]]
        vars_aux <- vars_aux[order(vars_aux[,"Pr(>|z|)"]),]
        vars_aux <- rownames(vars_aux)[1:min(max.variables, length(rownames(vars_aux)))]
        in.variable <- vars_aux #random selection of n_predictors because backwards selection
      }
      initial.model <- survival::coxph(as.formula(paste("Surv(", Time, ", ", Status, ") ~ .")), data = data[,colnames(data) %in% c(in.variable, "time", "event", "status")],
                                       method = "efron", model = T, singular.ok = T)
    }else{
      initial.model <- survival::coxph(as.formula(paste("Surv(", Time, ", ", Status, ") ~ .")), data = data[,colnames(data) %in% c(in.variable, "time", "event", "status")],
                                       method = "efron", model = T, singular.ok = T)
    }

  }
  else if (is.null(Time)) {

    if(startBackwards){
      in.variable <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
      if(!is.null(max.variables)){
        set.seed(123)
        in.variable <- sample(1:length(in.variable), max.variables) #random selection of n_predictors because backwards selection
      }
      initial.model <- survival::coxph(as.formula(paste("Surv(", T1, ", ", T2, ", ", Status, ") ~ .")),
                                       data = data[,colnames(data) %in% c(in.variable, "time", "event", "status")],
                                       method = "efron", model = T, singular.ok = T)
    }else{
      initial.model <- survival::coxph(as.formula(paste("Surv(", T1, ", ", T2, ", ", Status, ") ~ .")),
                                       data = data[,colnames(data) %in% c(in.variable, "time", "event", "status")],
                                       method = "efron", model = T, singular.ok = T)
    }

  }
  if (is.null(initial.model$coefficients)) { #si no tenemos variables por defecto...
    for (i in 1:length(variable.list)) { #por cada variable de estudio...
      utils::flush.console()
      if(is.null(T2)) {

        uni.model <- tryCatch({ #si NA saltará aviso
          survival::coxph(as.formula(paste("Surv(", Time, ", ", Status, ") ~ .")),
                          data = data[,colnames(data) %in% c(variable.list[i], "time", "event", "status")],
                          method = "efron", model = T, singular.ok = T)},
          #save the error
          error=function(cond) {
            # Choose a return value in case of error
            return(NA)
          })

        if(all(is.na(uni.model))){
          print(paste0("Variable ", variable.list[i], " may contain NAs. Fix it to include it to the model."))
          next
        }
      }
      if (is.null(Time)) {

        uni.model <- tryCatch({ #si NA saltará aviso
          survival::coxph(as.formula(paste("Surv(", T1, ", ", T2, ", ", Status, ") ~ .")),
                          data = data[,colnames(data) %in% c(variable.list[i], "time", "event", "status")],
                          method = "efron", model = T, singular.ok = T)},
          #save the error
          error=function(cond) {
            # Choose a return value in case of error
            return(NA)
          })

        if(is.na(uni.model)){
          print(paste0("Variable ", variable.list[i], " may contain NAs. Fix it to include it to the model."))
          next
        }
      }
      univar.pvalue[i] <- summary(uni.model)$coefficients[5] #cogemos el pvalor de cada variable univariante
    }
    variable.list1 <- variable.list[univar.pvalue <= 0.9 & !is.na(univar.pvalue)] #variables con un pvalor menor a 0.9
    univar.pvalue1 <- univar.pvalue[univar.pvalue <= 0.9 & !is.na(univar.pvalue)]
    uni.x <- variable.list1[which.min(univar.pvalue1)] #coger la de minimo pvalor
    if (length(uni.x) > 0) {
      if (is.null(T2)) {
        formula <- as.formula(paste("Surv(", Time, ", ", Status, ") ~ ", uni.x, sep = ""))
        temp.model <- survival::coxph(formula, data = data[,colnames(data) %in% c(uni.x, "time", "event", "status")], method = "efron", model = T, singular.ok = T)
      }
      if (is.null(Time)) {
        formula <- as.formula(paste("Surv(", T1, ", ", T2, ", ", Status, ") ~ ."))
        temp.model <- survival::coxph(formula, data = data[,colnames(data) %in% c(uni.x, "time", "event", "status")], method = "efron", model = T, singular.ok = T)
      }
    }
  }
  else if (!is.null(initial.model$coefficients)) {
    temp.model <- initial.model
  }

  #stepwise algorithm model

  i <- 0
  break.rule <- TRUE
  broke = F
  iter = 0
  last_removed = ""
  while (break.rule) {
    if(verbose){
      message(paste0("break.rule loop: ", iter))
      print(temp.model)
      iter = iter + 1
    }
    i <- i + 1
    if (i == 1) {
      variable.list2 <- setdiff(variable.list, all.vars(temp.model$formula))
    }
    else {
      variable.list2 <- setdiff(variable.list, c(all.vars(temp.model$formula), out.x))
      out.x <- NULL
    }

    #look for new variables if the limit is not reach
    if (length(variable.list2) != 0 & length(all.vars(temp.model$formula)) < max.variables) {
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      for (k in 1:length(variable.list2)) {
        utils::flush.console()

        #this update can generate problems for some variables:
        #Error in coxph.wtest(fit$var[nabeta, nabeta], temp, control$toler.chol) : NA/NaN/Inf in foreign function call (arg 3)
        model <- tryCatch({
          #line to check
          lst <- NULL
          if(REAL_INTERACTIONS){
            if(isInteraction(variable.list2[k])){
              aux <- strsplit(variable.list2[k], symbol)[[1]]

              # if(!aux[1] %in% names(temp.model$coefficients)){
              #   if(verbose)
              #     print("added ", aux[1], " as single part of interaction")
              #   invisible(update(temp.model, as.formula(paste(". ~ . + ", aux[1], sep = ""))))
              # }
              # if(!aux[2] %in% names(temp.model$coefficients)){
              #   if(verbose){
              #     print("added ", aux[2], " as single part of interaction")
              #   }
              #   invisible(update(temp.model, as.formula(paste(". ~ . + ", aux[2], sep = ""))))
              # }

            }
          }

          invisible(update(temp.model, as.formula(paste(". ~ . + ", variable.list2[k], sep = ""))))}, #with invisible function print nothing // does not work
          #save the error
          error=function(cond){
            # Choose a return value in case of error
            return(NA)
          })

        if(class(model)!="coxph"){# an error happend, next variable
          next
        }

        propHazard <- tryCatch({ #esto se usará en el futuro para comprobar si cumple el PH, tiene que ser invertible para calcularse, sino otra variable
          #line to check
          cox.zph(model, transform = "km", terms = T, global = T)},
          #save the error
          error=function(cond) {
            # Choose a return value in case of error
            return(NA)
          })

        if(class(propHazard)!="cox.zph"){# an error happend, next variable
          next
        }
        if(any(is.na(propHazard$table))){
          next
        }

        if (length(model$coefficients) > 1) {
          if (sum(is.na(model$coefficients)) != 0) {
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
          } else {
            anova.pvalue[k] <- anova(temp.model, model)[2,"P(>|Chi|)"]
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients),"Pr(>|z|)"]
          }
        }
      }
      #list of pvalues complete
      variable.list2.1 <- variable.list2[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      anova.pvalue2 <- anova.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      mv.pvalue2 <- mv.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      enter.x <- variable.list2.1[anova.pvalue2 == min(anova.pvalue2, na.rm = TRUE) & anova.pvalue2 <= sle]
      if(last_removed %in% enter.x){enter.x <- enter.x[-which(last_removed==enter.x)]} #last removed want to enter again, do not allow it
      if(verbose){
        message(paste0("ENTER: ",enter.x))
        print(paste0("ENTER: ",enter.x))
      }
      wald.p <- mv.pvalue2[anova.pvalue2 == min(anova.pvalue2, na.rm = TRUE) & anova.pvalue2 <= sle]
      if (length(setdiff(enter.x, NA)) != 0) {
        if (length(enter.x) > 1) {
          enter.x <- enter.x[which.min(wald.p)]
        }
        temp.model <- update(temp.model, as.formula(paste(". ~ . + ", enter.x, sep = "")))
      }
    } else {
      enter.x <- NULL
    }

    if (i == 1 & length(enter.x) == 0) {
      broke = T
      break
    } else {
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients),
                                c(enter.x, in.variable))
      if (length(variable.list3) != 0) {
        anova.pvalue <- NULL
        for (k in 1:length(variable.list3)) {

          model <- tryCatch({
            lst <- NULL
            if(REAL_INTERACTIONS){
              if(isInteraction(variable.list3[k])){
                aux <- strsplit(variable.list3[k], symbol)[[1]]

                # if(!aux[1] %in% names(temp.model$coefficients)){
                #   if(verbose)
                #     print("added ", aux[1], " as single part of interaction")
                #   invisible(update(temp.model, as.formula(paste(". ~ . + ", aux[1], sep = ""))))
                # }
                # if(!aux[2] %in% names(temp.model$coefficients)){
                #   if(verbose)
                #     print("added ", aux[1], " as single part of interaction")
                #   invisible(update(temp.model, as.formula(paste(". ~ . + ", aux[2], sep = ""))))
                # }

              }
            }
            #line to check
            update(temp.model, as.formula(paste(". ~ . - ", variable.list3[k], sep = "")))},
            #save the error
            error=function(cond) {
              # Choose a return value in case of error
              return(NA)
            })

          if(class(model)!="coxph"){# an error happend, next variable
            next
          }

          propHazard <- tryCatch({ #esto se usará en el futuro para comprobar si cumple el PH, tiene que ser invertible para calcularse, sino otra variable
            #line to check
            cox.zph(model, transform = "km", terms = T, global = T)},
            #save the error
            error=function(cond) {
              # Choose a return value in case of error
              return(NA)
            })

          if(class(propHazard)!="cox.zph"){# an error happend, next variable
            message("ph is empty")
            next
          }
          if(any(is.na(propHazard$table))){
            message("ph is NA")
            next
          }

          anova.pvalue[k] <- anova(model, temp.model)[2, "P(>|Chi|)"]
        }

        if(is.null(anova.pvalue)){
          message("anova is null")
          broke = T
          break.rule = F
          next
        }
        out.x <- variable.list3[anova.pvalue == max(anova.pvalue, na.rm = TRUE) & anova.pvalue > sls]
        out.x <- setdiff(out.x, NA)

        if(verbose){
          message(paste0("OUT: ",out.x))
          print(paste0("OUT: ",out.x))
        }

        if (length(out.x) != 0) {
          if (length(out.x) > 1) {
            out.x.1 <- out.x
            for (j in 1:length(out.x)) {
              out.x[j] <- out.x.1[(length(out.x) - j + 1)] #a la inversa
            }
            wald.p <- rep(NA, length(out.x))
            for (j in 1:length(out.x)) {
              wald.p[j] <- summary(temp.model)$coefficients[, "Pr(>|z|)"][rownames(summary(temp.model)$coefficients) == out.x[j]]
            }
          }

          lst_wald.p <- wald.p
          lst_out.x <- out.x
          while(length(lst_out.x)>0){
            out.x <- lst_out.x[which.max(lst_wald.p)]
            out.x.column <- colnames(data)[which(startsWith(out.x, colnames(data)))] #if factor, the name is not the column name (we have to update it)
            if(length(out.x.column)>1){
              out.x.column <- names(which.max(sapply(out.x.column, nchar))) #if more than one word, just take it the longer match
            }

            if(verbose){
              message(paste0("Wald: ",lst_wald.p))
              print(paste0("Wald: ",lst_wald.p))
            }
            if(verbose){
              message(paste0("Selected: ",out.x))
              print(paste0("Selected: ",out.x))
            }

            #check enter==out
            equal = F
            if(length(enter.x)>0){equal <- enter.x == out.x.column}
            if(max(lst_wald.p) > sls & !equal){ #if the PVal is grater than sls then remove, else, althouth the anova wants to remove, dont do it
              last_removed = out.x.column

              if(verbose){
                message("Original")
                print(temp.model)
              }

              temp.model <- update(temp.model, as.formula(paste(". ~ . - ", out.x.column, sep = "")))
              if(verbose){
                message("Modified")
                print(temp.model)
              }

              break #if create new temp.model break the while loop

            }else{
              #not grater than sls or equal to enter.x
              lst_out.x = lst_out.x[-which(lst_out.x==out.x)]
              lst_wald.p = lst_wald.p[-which(lst_wald.p==max(lst_wald.p))]
            }

          }

          if(length(lst_out.x)==0){ #if we end the while, out==NULL, if break the while do nothing
            out.x <- NULL
          }

        }
      }
      else {
        out.x <- NULL
      }
    }

    if ((length(enter.x) + length(out.x)) == 0) {
      final.model <- temp.model
      break.rule <- FALSE
    }
    enter.x <- NULL
  }

  if(broke){
    final.model <- temp.model #the model has not be updated because any variable enter to the model
  }

  return(final.model)
}

#### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Get the real name for the cox model variable that it is used in the matrix/data.frame
# because for categorical variables, cox add the level at the end and we can not do any match.
# We are not paying attention to different order, but yes for multiple interactions +2
#### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# data.frame/matrix  -  var_cox  -  interaction symbol #

getRealNameVarFromCoxModel <- function(x_sw_train, var, symbol = ":"){
  real.symbol <- symbol

  if(symbol=="."){
    symbol = "\\."
  }

  if(!var %in% colnames(x_sw_train)){

    if(isInteraction(var)){#check if each individual_variable is a factor
      n <- strsplit(var, symbol)[[1]]
      in_vector <- sapply(n, function(x){x %in% colnames(x_sw_train)})
      names_updated <- sapply(names(in_vector[in_vector==F]), function(x){x <- substr(x, 0, nchar(x)-1)})
      names(in_vector)[in_vector==F] <- names_updated

      res_var <- paste(names(in_vector), collapse = real.symbol)

    }else{#not interaction
      res_var <- substr(var, 0, nchar(var)-1)
      if(res_var %in% colnames(x_sw_train)){
        return(res_var)
      }

      res_var <- paste0("`",var,"`")
      if(res_var %in% colnames(x_sw_train)){
        return(res_var)
      }
    }
  }else{#var in colnames
    res_var = var
  }

  if(res_var %in% colnames(x_sw_train))
    return(res_var)
  else
    return(NULL)
}

isInteraction <- function(cn, symbol = ":"){
  if(symbol==".")
    symbol = "\\."
  if(length(grep(symbol,cn))>0){
    return(T)
  }else{
    return(F)
  }
}

deleteVariablesCox <- function(x_sw_train, x_sw_test=NULL, td, interactions=F, symbol = ":"){
  real.symbol = symbol

  if(symbol==".")
    symbol = "\\."

  if(length(td)>0){
    index <- match(td,colnames(x_sw_train)) #perfect matchs
    index <- index[!is.na(index)]
    #if it is NOT an interaction, delete all the interactions where the variable appears
    if(interactions){
      subtd <- unlist(lapply(td, function(x){if(!isInteraction(x))return(x)})) #for individual variables...
      inter <- NA
      splited <- lapply(colnames(x_sw_train), function(x){
        if(isInteraction(x))
          strsplit(x, symbol)})
      for(cn in splited){
        if(is.null(cn))
          next
        if(cn[[1]][1] %in% subtd | cn[[1]][2] %in% subtd){
          inter <- c(inter, paste0(cn[[1]][1], real.symbol, cn[[1]][2]))
        }
      }
      inter <- inter[!is.na(inter)]
      index <- unique(c(index, match(inter, colnames(x_sw_train))))
    }
    if(length(index)>0){
      if(class(x_sw_train)[1]=="data.frame"){
        x_sw_train[,index] <- NULL
        if(!is.null(x_sw_test))
          x_sw_test[,index] <- NULL
      }else if(class(x_sw_train)[1]=="matrix"){
        x_sw_train <- x_sw_train[,-index]
        if(!is.null(x_sw_test))
          x_sw_test <- x_sw_test[,-index]
      }
    }
  }
  return(list(cox.train = x_sw_train, cox.test = x_sw_test))
}

#creando mi propia linea de fit_zph
plotZPH <- function(fit_zph, oneToDelete, df=3){

  #DFCALLS
  x <- y <- NULL

  #https://rpubs.com/Cristina_Gil/Regr_no_lineal
  aux <- as.data.frame(cbind(fit_zph$time, fit_zph$y[,oneToDelete]))
  colnames(aux) <- c("x","y")
  ggp <- ggplot(aux, aes(x=x, y=y)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ splines::ns(x, df = df, intercept = T), color = "red",
                se = TRUE, level = 0.95) +
    xlab(label = "Time") + ylab(label = paste0("Beta(t) for ", oneToDelete))
  return(ggp)
}

### ## ##
# CLASS #
### ## ##

coxSW_class = function(cox_model, ...) {
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$coxSW)
  return(model)
}
