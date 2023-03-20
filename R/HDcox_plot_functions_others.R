#' loadingplot.HDcox
#'
#' @param model HDcox model.
#' @param zero.rm Logical. Remove variables equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @export

loadingplot.HDcox <- function(model, zero.rm = T, top = NULL, auto.limits = T){

  #DFCALLS
  variables <- pp <- NULL

  loading_values <- model$X$loadings
  ggp_loading <- NULL
  df <- NULL
  limit_color = 300

  if(auto.limits){
    auto.limits <- round2any(max(abs(loading_values)), accuracy = 0.1, f = ceiling)
  }else{
    auto.limits <- 1
  }

  for(i in 1:ncol(loading_values)){
    df <- as.data.frame(loading_values[,i])
    df <- cbind(df, rownames(loading_values))
    colnames(df) <- c("pp", "variables")

    col_name <- colnames(loading_values)[[i]]

    if(zero.rm){
      df <- df[!abs(df$pp)==0,]
    }

    if(!is.null(top)){
      if(top < nrow(df)){
        aux_df <- df
        aux_df$pp <- abs(aux_df$pp)
        aux_df <- aux_df[order(aux_df$pp, decreasing = T),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = T),]

    ggp <- NULL
    if(nrow(df)>limit_color){
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=pp))
    }else{
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=1))
    }

    ggp <- ggp +
      geom_bar(stat = "identity") +
      guides(color = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      #scale_fill_discrete(name = "New Legend Title") +
      xlab(label = paste0("Variables")) +
      ylab(label = paste0("Loading Value")) +
      ggtitle(paste0(attr(model, "model"), " - ", col_name))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                        mid = "white", midpoint = 0,
                                        high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }else{
      ggp <- ggp + scale_fill_gradient2(low = "blue",
                                        mid = "white", midpoint = 0,
                                        high = "magenta",
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }

    if(nrow(df)>limit_color){

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                           mid = "white", midpoint = 0,
                                           high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }else{
        ggp <- ggp + scale_color_gradient2(low = "blue",
                                           mid = "white", midpoint = 0,
                                           high = "magenta",
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }

    }

    if(auto.limits){
      ggp <- ggp + ylim(c(-1*auto.limits,auto.limits))
    }

    ggp_loading[[i]] = ggp
  }
  names(ggp_loading) <- colnames(loading_values)
  return(ggp_loading)
}

#' loadingplot.fromVector.HDcox
#'
#' @param model HDcox model.
#' @param vector Vector of loading
#' @param zero.rm Logical. Remove variables equal to 0 (default: FALSE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @export

loadingplot.fromVector.HDcox <- function(model, vector, zero.rm = F, top = NULL, auto.limits = T){

  #DFCALLS
  variables <- pp <- NULL

  loading_values <- vector
  ggp_loading <- NULL
  df <- NULL
  limit_color = 300

  if(auto.limits){
    auto.limits <- round2any(max(abs(loading_values)), accuracy = 0.1, f = ceiling)
  }else{
    auto.limits <- 1
  }

  for(i in 1:ncol(loading_values)){
    df <- as.data.frame(loading_values[,i])
    df <- cbind(df, rownames(loading_values))
    colnames(df) <- c("pp", "variables")

    col_name <- colnames(loading_values)[[i]]

    if(zero.rm){
      df <- df[!abs(df$pp)==0,]
    }

    if(!is.null(top)){
      if(top < nrow(df)){
        aux_df <- df
        aux_df$pp <- abs(aux_df$pp)
        aux_df <- aux_df[order(aux_df$pp, decreasing = T),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = T),]

    ggp <- NULL
    if(nrow(df)>limit_color){
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=pp))
    }else{
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=1))
    }

    ggp <- ggp +
      geom_bar(stat = "identity") +
      guides(color = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab(label = paste0("Variables")) +
      ylab(label = paste0("Loading Value")) +
      ggtitle(paste0(attr(model, "model"), " - ", col_name))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                        mid = "white", midpoint = 0,
                                        high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }else{
      ggp <- ggp + scale_fill_gradient2(low = "blue",
                                        mid = "white", midpoint = 0,
                                        high = "red",
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }

    if(nrow(df)>limit_color){

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                           mid = "white", midpoint = 0,
                                           high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }else{
        ggp <- ggp + scale_color_gradient2(low = "blue",
                                           mid = "white", midpoint = 0,
                                           high = "red",
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }

    }

    if(auto.limits){
      ggp <- ggp + ylim(c(-1*auto.limits,auto.limits))
    }

    ggp_loading[[i]] = ggp
  }
  names(ggp_loading) <- colnames(loading_values)
  return(ggp_loading)
}

#' w.starplot.HDcox
#'
#' @param model HDcox model.
#' @param zero.rm Logical. Remove variables equal to 0 (default: FALSE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @export

w.starplot.HDcox <- function(model, zero.rm = F, top = NULL, auto.limits = T){

  #DFCALLS
  variables <- pp <- NULL

  loading_values <- model$X$W.star
  ggp_loading <- NULL
  df <- NULL
  limit_color = 300

  if(auto.limits){
    auto.limits <- round2any(max(abs(loading_values)), accuracy = 0.1, f = ceiling)
  }else{
    auto.limits <- 1
  }

  for(i in 1:ncol(loading_values)){
    df <- as.data.frame(loading_values[,i])
    df <- cbind(df, rownames(loading_values))
    colnames(df) <- c("pp", "variables")

    col_name <- colnames(loading_values)[[i]]

    if(zero.rm){
      df <- df[!abs(df$pp)==0,]
    }

    if(!is.null(top)){
      if(top < nrow(df)){
        aux_df <- df
        aux_df$pp <- abs(aux_df$pp)
        aux_df <- aux_df[order(aux_df$pp, decreasing = T),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = T),]

    ggp <- NULL
    if(nrow(df)>limit_color){
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=pp))
    }else{
      ggp <- ggplot(df, aes(x = reorder(variables, -pp), y = pp, fill=pp, color=1))
    }

    ggp <- ggp +
      geom_bar(stat = "identity") +
      guides(color = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      #scale_fill_discrete(name = "New Legend Title") +
      xlab(label = paste0("Variables")) +
      ylab(label = paste0("W.star Value")) +
      ggtitle(paste0(attr(model, "model"), " - ", col_name))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                        mid = "white", midpoint = 0,
                                        high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }else{
      ggp <- ggp + scale_fill_gradient2(low = "blue",
                                        mid = "white", midpoint = 0,
                                        high = "magenta",
                                        limits = c(-1*auto.limits,auto.limits), name = "Values")
    }

    if(nrow(df)>limit_color){

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                           mid = "white", midpoint = 0,
                                           high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }else{
        ggp <- ggp + scale_color_gradient2(low = "blue",
                                           mid = "white", midpoint = 0,
                                           high = "red",
                                           limits = c(-1*auto.limits,auto.limits), name = "Values")
      }

    }

    if(auto.limits){
      ggp <- ggp + ylim(c(-1*auto.limits,auto.limits))
    }

    ggp_loading[[i]] = ggp
  }
  names(ggp_loading) <- colnames(loading_values)
  return(ggp_loading)
}

#' plot_classicalcox.comparePatients
plot_classicalcox.comparePatients <- function(model, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                              auto.limits = T, top = NULL){

  #DFCALLS
  value <- patients <- NULL

  coefficients <- model$survival_model$coef
  coefficients <- as.data.frame(coefficients)
  colnames(coefficients) <- "value"
  coefficients <- coefficients[order(coefficients$value, decreasing = T),,drop=F]

  if(!is.null(top)){
    if(top < nrow(coefficients)){
      aux_df <- coefficients
      aux_df[,"value"] <- abs(aux_df[,"value",drop=F])
      aux_df <- aux_df[order(aux_df[,"value",drop=T], decreasing = T),,drop=F]
      aux_df <- aux_df[1:top,,drop=F]
      coefficients <- coefficients[rownames(coefficients) %in% rownames(aux_df),,drop=F]
    }
  }

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

  #lp.new_pat_manual <- norm_patient[,rownames(coefficients)] %*% coefficients #predict lp
  lp.new_pat_variable <- apply(norm_patient[,deleteIllegalChars(rownames(coefficients)),drop=F], 1, function(x){
    x * coefficients$value #predict terms
  })

  #Compute LP without top variables
  #can be change for cox.prediction(model = model, new_data = patient, time = time, type = type, method = "cox")
  #for each patient on the data frame

  lp.pats <- norm_patient[,deleteIllegalChars(names(model$survival_model$coef))] %*% model$survival_model$coef
  colnames(lp.pats) <- "linear predictor"

  rownames(lp.new_pat_variable) <- rownames(coefficients)
  lp.new_pat_variable <- rbind(lp.new_pat_variable, lp.pats[,1])
  rownames(lp.new_pat_variable)[nrow(lp.new_pat_variable)] <- "linear predictor"
  lp.new_pat_variable <- as.data.frame(lp.new_pat_variable)
  lp.new_pat_variable$var <- rownames(lp.new_pat_variable)

  lp.new_pat_variable <- tidyr::pivot_longer(lp.new_pat_variable, !var, names_to = "patients", values_to = "value")

  lp.new_pat_variable$var <- factor(lp.new_pat_variable$var, levels = unique(lp.new_pat_variable$var))

  lp.new_pat_variable$lp.flag <- ifelse(lp.new_pat_variable$var == "linear predictor", T, F)
  lp.new_pat_variable$lp.flag <- factor(lp.new_pat_variable$lp.flag)

  lp.new_pat_variable$patients <- factor(lp.new_pat_variable$patients, levels = rownames(new_data))

  accuracy <- 0.1
  auto.limits.flag = T

  df_cox_sd <- summary(model$survival_model$fit)[[7]][,"se(coef)"]

  sd.min <- coefficients - as.data.frame(df_cox_sd[rownames(coefficients)])
  sd.max <- coefficients + as.data.frame(df_cox_sd[rownames(coefficients)])
  auto.limits <- NULL
  if(auto.limits.flag){
    if(!is.null(sd.min) & !is.null(sd.max)){
      auto.limits_min <- round2any(x = max(c(abs(coefficients$value-sd.min$value),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(x = max(c(abs(coefficients$value+sd.max$value),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
      auto.limits <- max(auto.limits_min, auto.limits_max)
    }else{
      auto.limits <- round2any(max(abs(lp.new_pat_variable$value)), accuracy = accuracy, f = ceiling)
    }
  }else{
    auto.limits <- round2any(max(c(abs(sd.max), abs(sd.min), abs(lp.new_pat_variable$value))), accuracy = accuracy, f = ceiling)
  }

  ggp <- ggplot(lp.new_pat_variable[lp.new_pat_variable$lp.flag==F,], aes(x = var, y = value, fill = patients)) +
    geom_bar(stat = "identity", position = "dodge") + xlab(label = "Variables")
  ggp2 <- ggplot(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,], aes(x = var, y = value, fill = patients)) +
    geom_bar(stat = "identity", position = "dodge")
  #guides(color = "none")

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "complete", continuous = F)
    ggp2 <- ggp2 + RColorConesa::scale_fill_conesa(palette = "complete", continuous = F)
  }

  if(!auto.limits.flag){
    #ggp <- ggp + scale_y_continuous(breaks=seq(-1*auto.limits, auto.limits, 0.1))
    ggp <- ggp + scale_y_continuous(n.breaks = 10)
    ggp2 <- ggp2 + scale_y_continuous(n.breaks = 10)
  }else{
    #ggp <- ggp + scale_y_continuous(breaks=seq(-1*auto.limits, auto.limits, 0.1), limits = c(-1*auto.limits, auto.limits))
    ggp <- ggp + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))
    ggp2 <- ggp2 + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))
  }

  if(length(unique(lp.new_pat_variable$var))>15){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggp2 <- ggp2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  res_all.plot <- ggp
  res_lp.plot <- ggp2 + xlab(label = "")

  ggp <- ggp + guides(fill = "none")
  ggp2 <- ggp2 + ylab(label = "") + xlab(label = "")

  pp <- ggpubr::ggarrange(ggp, ggp2, ncol = 2, widths = c(0.8, 0.2), align = "h")

  return(list(plot = pp, var.plot = res_all.plot, lp.plot = res_lp.plot, lp = lp.pats, lp.var = lp.new_pat_variable, norm_patients = norm_patient, patients = new_data))
}

#' prop.vector.density
prop.vector.density <- function(df, breaks = 20){
  min <- min(df$lp)
  max <- max(df$lp)

  inc <- (max - min)/breaks
  for(i in seq(min+inc, max, inc)){
    mmin <- i-inc
    mmax <- i
    prop.between2values(df, mmin, mmax)
  }
}
