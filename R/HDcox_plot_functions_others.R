#' loadingplot.HDcox
#'
#' @param model HDcox model.
#' @param zero.rm Logical. Remove variables equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @export

loadingplot.HDcox <- function(model, zero.rm = T, top = NULL, auto.limits = T){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

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

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

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

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

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

# prop.vector.density
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
