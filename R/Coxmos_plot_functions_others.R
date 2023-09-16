#' loadingplot.Coxmos
#'
#' @description
#' The `loadingplot.Coxmos` function visualizes the loading values of a given Coxmos model. The
#' function produces a series of bar plots for each component's loading values, offering a
#' comprehensive view of the model's variable contributions. The plots can be customized to exclude
#' zero loadings, display only the top variables, and automatically adjust the color scale limits.
#' @details
#' The primary objective of the `loadingplot.Coxmos` function is to facilitate the interpretation of
#' Coxmos models by visualizing the loading values of each component. The function first verifies the
#' class of the provided model to ensure it is a valid Coxmos model.
#'
#' The loading values are extracted from the model and processed based on the user's specifications.
#' If the `zero.rm` parameter is set to TRUE, variables with zero loadings are excluded from the
#' visualization. Additionally, if the `top` parameter is specified, only the top variables, ranked
#' by their absolute loading values, are displayed.
#'
#' The function employs the 'ggplot2' framework for visualization. The color scale of the plots can be
#' automatically adjusted based on the maximum absolute loading value when `auto.limits` is set to
#' TRUE. If the `RColorConesa` package is available, it utilizes its color palettes for enhanced
#' visualization; otherwise, default colors are applied.
#'
#' @param model Coxmos model.
#' @param zero.rm Logical. Remove variables equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @return A list of \code{ggplot2} objects, each representing the loading values for a component of
#' the Coxmos model.
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' X <- X_proteomic[,1:50]
#' Y <- Y_proteomic
#' splsicox.model <- splsicox(X, Y, n.comp = 2, spv_penalty = 0.5, x.center = TRUE, x.scale = TRUE)
#' loadingplot.Coxmos(model = splsicox.model)

loadingplot.Coxmos <- function(model, zero.rm = TRUE, top = NULL, auto.limits = TRUE){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class Coxmos.")
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
        aux_df <- aux_df[order(aux_df$pp, decreasing = TRUE),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = TRUE),]

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

#' loadingplot.fromVector.Coxmos
#'
#' @param model Coxmos model.
#' @param vector Vector of loading
#' @param zero.rm Logical. Remove variables equal to 0 (default: FALSE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).

loadingplot.fromVector.Coxmos <- function(model, vector, zero.rm = FALSE, top = NULL, auto.limits = TRUE){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class Coxmos.")
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
        aux_df <- aux_df[order(aux_df$pp, decreasing = TRUE),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = TRUE),]

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

#' w.starplot.Coxmos
#'
#' @description
#' The `w.starplot.Coxmos` function offers a graphical representation of the W* (W star) values from
#' a given Coxmos model. Through this visualization, users can gain insights into the variable
#' contributions and their significance in the model. The function provides options for customization,
#' allowing users to focus on specific variables, exclude zero values, and adjust the visual limits.
#' @details
#' The `w.starplot.Coxmos` function is tailored to visualize the W* values, which are indicative of
#' the variable contributions in a Coxmos model. Initially, the function checks the class of the
#' provided model to ensure its compatibility with the Coxmos framework.
#'
#' The W* values are extracted from the model and subsequently processed based on user-defined
#' parameters. The `zero.rm` option allows users to exclude variables with zero W* values, ensuring
#' a more concise visualization. If the `top` parameter is specified, the function focuses on
#' displaying only the top-ranked variables based on their absolute W* values.
#'
#' The visualization is constructed using the 'ggplot2' framework. The color scale can be automatically
#' adjusted to the maximum absolute W* value when the `auto.limits` parameter is set to TRUE. The
#' function also checks for the availability of the `RColorConesa` package. If present, it leverages
#' its color palettes for a more refined visualization; in its absence, default color schemes are applied.
#'
#' @param model Coxmos model.
#' @param zero.rm Logical. Remove variables equal to 0 (default: FALSE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#'
#' @return A list of \code{ggplot2} objects, each representing the W* values for a component of
#' the Coxmos model.
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' X <- X_proteomic[,1:50]
#' Y <- Y_proteomic
#' splsicox.model <- splsicox(X, Y, n.comp = 2, spv_penalty = 0.5, x.center = TRUE, x.scale = TRUE)
#' w.starplot.Coxmos(model = splsicox.model)

w.starplot.Coxmos <- function(model, zero.rm = FALSE, top = NULL, auto.limits = TRUE){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class Coxmos.")
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
        aux_df <- aux_df[order(aux_df$pp, decreasing = TRUE),]
        aux_df <- aux_df[1:top,]
        df <- df[df$variables %in% aux_df$variables,]
      }
    }

    df <- df[order(df$pp, decreasing = TRUE),]

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
