#### ### ### #
# SAVE PLOTS #
#### ### ### #

#' save_ggplot
#' @description Allows to save ggplot2 objects in .tiff format based on an specific resolution.
#'
#' @param plot ggplot2 object
#' @param folder folder
#' @param name file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#'
#' @export
#'
#' @examples
#' \dontrun{
#' save_ggplot(plot, folder)
#' }

save_ggplot <- function(plot, folder = NULL, name = "plot", wide = T, quality = "4K", dpi = 80, custom = NULL){
  width=NULL
  height=NULL

  if(!quality %in% c("HD", "FHD", "2K", "4K", "8K")){
    stop("uality must be one of the following options: 'HD', 'FHD', '2K', '4K', '8K'")
  }

  ratios <- c(1.5,1.333333,1.5,2)

  if(wide){
    if(quality == "HD"){
      width = 1280/dpi#4.266667
      height = 720/dpi#2.4
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1920/dpi#6.4
      height = 1080/dpi#3.6
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 2560/dpi#8.533333
      height = 1440/dpi#4.8
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 3840/dpi#12.8
      height = 2160/dpi#7.2
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 7680/dpi#25.6
      height = 4320/dpi#14.4
    }
  }else{
    if(quality == "HD"){
      width = 960/dpi#3.19992
      height = 720/dpi
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1440/dpi#4.79988
      height = 1080/dpi
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 1920/dpi#6.39984
      height = 1440/dpi
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 2880/dpi#9.59976
      height = 2160/dpi
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 5760/dpi#19.19952
      height = 4320/dpi
    }
  }

  if(!is.null(custom)){
    if(length(custom)==2){
      width = custom[1]
      height = custom[2]
    }
  }

  if(!endsWith(name,".tiff")){
    name <- paste0(name, ".tiff")
  }

  #remove illegal characters
  name <- transformIllegalChars(name, except = c("-"))

  if(class(plot)[1] %in% "ggsurvplot"){
    plot_surv = plot$plot
    if("table" %in% names(plot)){
      p2 = plot$table
      plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
    }
    ggsave(plot = plot_surv, filename = paste0(folder,name), width = width, height = height, device='tiff', dpi=dpi)
  }else{
    ggsave(plot = plot, filename = paste0(folder,name), width = width, height = height, device='tiff', dpi=dpi)
  }
}

#' save_ggplot.svg
#' @description Allows to save ggplot2 objects in .svg format based on an specific resolution.
#'
#' @param plot ggplot2 object
#' @param folder folder
#' @param name file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#'
#' @export
#'
#' @examples
#' \dontrun{
#' save_ggplot.svg(plot, folder)
#' }

save_ggplot.svg <- function(plot, folder = NULL, name = "plot", wide = T, quality = "4K", dpi = 80, custom = NULL){
  width=NULL
  height=NULL

  if(!quality %in% c("HD", "FHD", "2K", "4K", "8K")){
    stop("quality must be one of the following options: 'HD', 'FHD', '2K', '4K', '8K'")
  }

  ratios <- c(1.5,1.333333,1.5,2)

  if(wide){
    if(quality == "HD"){
      width = 1280/dpi#4.266667
      height = 720/dpi#2.4
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1920/dpi#6.4
      height = 1080/dpi#3.6
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 2560/dpi#8.533333
      height = 1440/dpi#4.8
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 3840/dpi#12.8
      height = 2160/dpi#7.2
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 7680/dpi#25.6
      height = 4320/dpi#14.4
    }
  }else{
    if(quality == "HD"){
      width = 960/dpi#3.19992
      height = 720/dpi
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1440/dpi#4.79988
      height = 1080/dpi
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 1920/dpi#6.39984
      height = 1440/dpi
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 2880/dpi#9.59976
      height = 2160/dpi
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 5760/dpi#19.19952
      height = 4320/dpi
    }
  }

  if(!is.null(custom)){
    if(length(custom)==2){
      width = custom[1]
      height = custom[2]
    }
  }

  if(!endsWith(name,".svg")){
    name <- paste0(name, ".svg")
  }

  #remove illegal characters
  name <- transformIllegalChars(name, except = c("-"))

  if(class(plot)[1] %in% "ggsurvplot"){
    plot_surv = plot$plot
    if("table" %in% names(plot)){
      p2 = plot$table
      plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
    }
    ggsave(plot = plot_surv, filename = paste0(folder,name), width = width, height = height, device='svg', dpi=dpi)
  }else{
    ggsave(plot = plot, filename = paste0(folder,name), width = width, height = height, device='svg', dpi=dpi)
  }
}

#' save_ggplot_lst
#' @description Allows to save a list of ggplot2 objects in .tiff format based on an specific resolution.
#'
#' @param lst_plots list of ggplot2
#' @param folder folder
#' @param prefix prefix for file name
#' @param suffix sufix for file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#' @param object_name file name
#'
#' @export
#'
#' @examples
#' \dontrun{
#' save_ggplot_lst(plot_lst, folder)
#' }

save_ggplot_lst <- function(lst_plots, folder = NULL, prefix = NULL, suffix = NULL, wide = T, quality = "4K", dpi = 80, custom = NULL, object_name = NULL){
  width=NULL
  height=NULL

  if(!quality %in% c("HD", "FHD", "2K", "4K", "8K")){
    stop("quality must be one of the following options: 'HD', 'FHD', '2K', '4K', '8K'")
  }

  ratios <- c(1.5,1.333333,1.5,2)

  if(wide){
    if(quality == "HD"){
      width = 1280/dpi#4.266667
      height = 720/dpi#2.4
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1920/dpi#6.4
      height = 1080/dpi#3.6
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 2560/dpi#8.533333
      height = 1440/dpi#4.8
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 3840/dpi#12.8
      height = 2160/dpi#7.2
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 7680/dpi#25.6
      height = 4320/dpi#14.4
    }
  }else{
    if(quality == "HD"){
      width = 960/dpi#3.19992
      height = 720/dpi
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1440/dpi#4.79988
      height = 1080/dpi
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 1920/dpi#6.39984
      height = 1440/dpi
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 2880/dpi#9.59976
      height = 2160/dpi
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 5760/dpi#19.19952
      height = 4320/dpi
    }
  }

  if(!is.null(custom)){
    if(length(custom)==2){
      width = custom[1]
      height = custom[2]
    }
  }

  if(!is.null(names(lst_plots))){
    for(cn in names(lst_plots)){

      name <- paste0(prefix,cn,suffix)
      #remove illegal characters
      name <- transformIllegalChars(name, except = c("-"))

      name <- paste0(folder,name)
      if(!endsWith(name,".tiff")){
        name <- paste0(name, ".tiff")
      }

      if(is.null(object_name)){
        if(class(lst_plots[[cn]])[1] %in% "ggsurvplot"){
          plot_surv = lst_plots[[cn]]$plot
          if("table" %in% names(lst_plots[[cn]])){
            p2 = lst_plots[[cn]]$table
            plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
          }
          ggsave(plot = plot_surv, filename = name, width = width, height = height, device='tiff', dpi=dpi)
        }else{
          ggsave(plot = lst_plots[[cn]], filename = name, width = width, height = height, device='tiff', dpi=dpi)
        }
      }else{
        if(class(lst_plots[[cn]][[object_name]])[1] %in% "ggsurvplot"){
          plot_surv = lst_plots[[cn]][[object_name]]$plot
          if("table" %in% names(lst_plots[[cn]][[object_name]])){
            p2 = lst_plots[[cn]][[object_name]]$table
            plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
          }
          ggsave(plot = plot_surv, filename = name, width = width, height = height, device='tiff', dpi=dpi)
        }else{
          ggsave(plot = lst_plots[[cn]][[object_name]], filename = name, width = width, height = height, device='tiff', dpi=dpi)
        }
      }
    }
  }else{
    for(cn in 1:length(lst_plots)){

      name <- paste0(prefix,cn,suffix)
      #remove illegal characters
      name <- transformIllegalChars(name, except = c("-"))

      name <- paste0(folder,name)
      if(!endsWith(name,".tiff")){
        name <- paste0(name, ".tiff")
      }

      if(is.null(object_name)){
        ggsave(plot = lst_plots[[cn]], filename = name, width = width, height = height, device='tiff', dpi=dpi)
      }else{
        ggsave(plot = lst_plots[[cn]][[object_name]], filename = name, width = width, height = height, device='tiff', dpi=dpi)
      }
    }
  }
}

#' save_ggplot_lst.svg
#' @description Allows to save a list of ggplot2 objects in .svg format based on an specific resolution.
#'
#' @param lst_plots list of ggplot2
#' @param folder folder
#' @param prefix prefix for file name
#' @param suffix sufix for file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#' @param object_name file name
#'
#' @export
#'
#' @examples
#' \dontrun{
#' save_ggplot_lst.svg(plot_lst, folder)
#' }

save_ggplot_lst.svg <- function(lst_plots, folder = NULL, prefix = NULL, suffix = NULL, wide = T, quality = "4K", dpi = 80, custom = NULL, object_name = NULL){
  width=NULL
  height=NULL

  if(!quality %in% c("HD", "FHD", "2K", "4K", "8K")){
    stop("quality must be one of the following options: 'HD', 'FHD', '2K', '4K', '8K'")
  }

  ratios <- c(1.5,1.333333,1.5,2)

  if(wide){
    if(quality == "HD"){
      width = 1280/dpi#4.266667
      height = 720/dpi#2.4
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1920/dpi#6.4
      height = 1080/dpi#3.6
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 2560/dpi#8.533333
      height = 1440/dpi#4.8
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 3840/dpi#12.8
      height = 2160/dpi#7.2
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 7680/dpi#25.6
      height = 4320/dpi#14.4
    }
  }else{
    if(quality == "HD"){
      width = 960/dpi#3.19992
      height = 720/dpi
    }else if(quality == "FHD"){
      dpi = dpi * ratios[1]
      width = 1440/dpi#4.79988
      height = 1080/dpi
    }else if(quality == "2K"){
      dpi = dpi * ratios[1] * ratios[2]
      width = 1920/dpi#6.39984
      height = 1440/dpi
    }else if(quality == "4K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3]
      width = 2880/dpi#9.59976
      height = 2160/dpi
    }else if(quality == "8K"){
      dpi = dpi * ratios[1] * ratios[2] * ratios[3] * ratios[4]
      width = 5760/dpi#19.19952
      height = 4320/dpi
    }
  }

  if(!is.null(custom)){
    if(length(custom)==2){
      width = custom[1]
      height = custom[2]
    }
  }

  if(!is.null(names(lst_plots))){
    for(cn in names(lst_plots)){

      name <- paste0(prefix,cn,suffix)
      #remove illegal characters
      name <- transformIllegalChars(name, except = c("-"))

      name <- paste0(folder,name)
      if(!endsWith(name,".svg")){
        name <- paste0(name, ".svg")
      }

      if(is.null(object_name)){

        if(class(lst_plots[[cn]])[1] %in% "ggsurvplot"){
          plot_surv = lst_plots[[cn]]$plot
          if("table" %in% names(lst_plots[[cn]])){
            p2 = lst_plots[[cn]]$table
            plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
          }
          ggsave(plot = plot_surv, filename = name, width = width, height = height, device='svg', dpi=dpi)
        }else{
          ggsave(plot = lst_plots[[cn]], filename = name, width = width, height = height, device='svg', dpi=dpi)
        }

      }else{
        if(class(lst_plots[[cn]][[object_name]])[1] %in% "ggsurvplot"){
          plot_surv = lst_plots[[cn]][[object_name]]$plot
          if("table" %in% names(lst_plots[[cn]][[object_name]])){
            p2 = lst_plots[[cn]][[object_name]]$table
            plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
          }
          ggsave(plot = plot_surv, filename = name, width = width, height = height, device='svg', dpi=dpi)
        }else{
          ggsave(plot = lst_plots[[cn]][[object_name]], filename = name, width = width, height = height, device='svg', dpi=dpi)
        }
      }
    }
  }else{
    for(cn in 1:length(lst_plots)){

      name <- paste0(prefix,cn,suffix)
      #remove illegal characters
      name <- transformIllegalChars(name, except = c("-"))

      name <- paste0(folder,name)
      if(!endsWith(name,".svg")){
        name <- paste0(name, ".svg")
      }

      if(is.null(object_name)){
        ggsave(plot = lst_plots[[cn]], filename = name, width = width, height = height, device='svg', dpi=dpi)
      }else{
        ggsave(plot = lst_plots[[cn]][[object_name]], filename = name, width = width, height = height, device='svg', dpi=dpi)
      }
    }
  }
}

#### ### ### ### #
# TIME CONSUMING #
#### ### ### ### #

#' Time consuming plot.
#' @description Generate a ggplot2 with the consumed time for each model included in "lst_models". The models have to be generated by the HDcox methods (include the variable time inside the model objects).
#'
#' @param lst_models List of HDcox models. Each HDCox object has the attribute time measured in minutes (cross-validation models could be also added to this function).
#' @param x.text Character. X axis title.
#' @param y.text Character. Y axis title. If y.text = NULL, then y.text = "Time (mins)" (default: NULL).
#'
#' @return A ggplot2 bar plot.
#' @export
#'
#' @details For time comparison between best models, use the cross validation objects instead the final models. The training has to be taken into account.
#'
#' @examples
#' \dontrun{
#'   lst_models = list("cox" = cox.model,
#'   "cv.sPLS-ICOX" = cv.splsicox_model,
#'   "sPLS-ICOX" = splsicox_model)
#'   plot_time.list(lst_models, x.text = "Method")
#' }

plot_time.list <- function(lst_models, x.text = "Method", y.text = NULL){

  if(is.null(names(lst_models))){
    names(lst_models) <- unlist(lapply(lst_models, function(x){
      if("time" %in% names(x)){
        attr(x, "model")
      }else if("time" %in% names(x[[1]])){
        attr(x[[1]], "model")
      }
    }))
  }

  lst_times <- list()
  for(m in names(lst_models)){
    if(isa(lst_models[[m]],pkg.env$model_class)){
      lst_times[[m]] <- lst_models[[m]]$time
    }else if(isa(lst_models[[m]][[1]],pkg.env$model_class)){
      eval_sum <- lst_models[[m]][[1]]$time
      if(length(lst_models[[m]])>1){
        for(i in 2:length(lst_models[[m]])){
          eval_sum <- eval_sum + lst_models[[m]][[i]]$time
        }
      }
      lst_times[[m]] <- eval_sum
    }

  }

  total_time <- lst_times[[1]]
  if(length(lst_times)>1){
    for(m in 2:length(lst_times)){
      total_time <- total_time + lst_times[[m]]
    }
  }

  lst_times$Total <- total_time

  df.times <- do.call(rbind.data.frame, lst_times)
  colnames(df.times) <- "times"
  df.times$method <- names(lst_times)
  rownames(df.times) <- NULL

  roundTo = 0
  max.breaks = 10
  if(roundTo == 0){
    #select the decimals of Y
    if(length(grep("\\.", df.times$times))>0){
      ch <- gsub("\\.", "", as.character(format(min(df.times$times)/max.breaks, scientific = F, trim = T)))
      cont = 0
      for(c in 1:nchar(ch)){
        if(substr(ch,c,c) == "0"){
          cont = cont + 1
        }else{
          break
        }
      }
      roundTo = 1*10^-cont
    }else{
      roundTo = 0.1
    }
  }

  breaks_size = round2any(max(df.times$times), roundTo, f = ceiling) / max.breaks

  roundTo = 0
  max.breaks = 10
  if(roundTo == 0){
    #select the decimals of Y
    if(length(grep("\\.", df.times$times))>0){
      ch <- gsub("\\.", "", as.character(format(max(df.times$times)/max.breaks, scientific = F, trim = T)))
      cont = 1
      for(c in 1:nchar(ch)){
        if(substr(ch,c,c) == "0"){
          cont = cont + 1
        }else{
          break
        }
      }
      roundTo = 1*10^-cont
    }else{
      roundTo = 0.1
    }
  }

  breaks_size = round2any(breaks_size, roundTo, f = ceiling)
  breaks = seq(0, max(df.times$times)+breaks_size, by=breaks_size)

  accuracy <- roundTo
  max <- max(breaks)

  df.times$times <- round(df.times$times, digits = 4)
  x.var = "method"
  y.var = "times"
  x.color = "method"
  x.text = x.text
  if(is.null(y.text)){
    y.text = paste0("Time (",attr(lst_times[["Total"]], "units"),")")
  }

  df.times$method <- factor(df.times$method, levels = df.times$method)

  ggp_time <- ggplot(df.times, aes_string(x = x.var, y = y.var, fill = x.color)) +
    geom_bar(stat="identity") +
    scale_y_continuous(breaks = breaks) +
    geom_text(aes_string(label = "times"), vjust = 0, nudge_y = accuracy)

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp_time <- ggp_time + RColorConesa::scale_fill_conesa(palette = "complete")
  }

  if(!is.null(y.text)){
    ggp_time = ggp_time + ylab(label = y.text)
  }
  if(!is.null(x.text)){
    ggp_time = ggp_time + xlab(label = x.text)
  }

  return(ggp_time)
}

#### ### ### ### ###
# MODEL EVALUATION #
#### ### ### ### ###

#' plot_evaluation.list
#' @description Run the function "plot_evaluation" for a list of results. More information in "?plot_evaluation".
#'
#' @param lst_eval_results List (named) of HDcox evaluation results.
#' @param evaluation Character. Perform the evaluation using the "AUC" or "Brier" metric (default: "AUC").
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param y.min Numeric. Minimum Y value for establish the Y axis value. If y.min = NULL, automatic detection is performed (default: NULL).
#' @param type Character. Plot type. Must be one of the following: "both", "line" or "mean". In other case, "both" will be selected (default: "both").
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #lst_eval_results #list of multiple eval_HDcox_models()
#' plot_evaluation.list(lst_eval_results)
#' }

plot_evaluation.list <- function(lst_eval_results, evaluation = "AUC", pred.attr = "mean", y.min = NULL, type = "both"){

  lst_res <- purrr::map(lst_eval_results, ~plot_evaluation(eval_results = .,
                                                           evaluation = evaluation,
                                                           pred.attr = pred.attr,
                                                           y.min = y.min, type = type))

  return(lst_res)

}

#' plot_evaluation
#' @description Perform a evaluation test for the result object provided. It includes test comparison as "t.test", "anova", "wilcoxon", "kruscal-wallis". Furthermore, it generates plots for comparing the results.
#'
#' @param eval_results HDcox evaluation object.
#' @param evaluation Character. Perform the evaluation using the "AUC" or "Brier" metric (default: "AUC").
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following: "mean" or "median" (default: "mean").
#' @param y.min Numeric. Minimum Y value for establish the Y axis value. If y.min = NULL, automatic detection is performed (default: NULL).
#' @param type Character. Plot type. Must be one of the following: "both", "line" or "mean". In other case, "both" will be selected (default: "both").
#' @param legend_title Character. Legend title (default: "Method").
#' @param legend_size_text Numeric. Text size for legend title (default: 12).
#' @param x_axis_size_text Numeric. Text size for x axis (default: 10).
#' @param y_axis_size_text Numeric. Text size for y axis (default: 10).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' eval_results <- eval_HDcox_models(lst_models, X_test, Y_test, pred.method)
#' plot_evaluation(eval_results)
#' }

plot_evaluation <- function(eval_results, evaluation = "AUC", pred.attr = "mean", y.min = NULL, type = "both", legend_title = "Method", legend_size_text = 12, x_axis_size_text = 10, y_axis_size_text = 10){

  if(!evaluation %in% c("AUC", "Brier")){
    message("Evaluation parameter is not 'AUC' or 'Brier'. Changed to 'AUC'.")
    type = "AUC"
  }

  if(!pred.attr %in% c("mean", "median")){
    stop("pred.attr parameter must be one of: 'mean' or 'median'")
  }

  if(!type %in% c("both", "line", "mean")){
    type = "both"
  }

  #select minimum for all evals
  if(is.null(y.min)){
    if(evaluation=="AUC"){
      y.min <- floor(min(eval_results$df$AUC, na.rm = T)*10)/10
    }else{
      y.min <- floor(min(eval_results$df$Brier, na.rm = T)*10)/10
    }

  }

  if(is.infinite(y.min)){
    if(evaluation=="AUC"){
      message("All AUC is NA. Returning NA.")
    }else{
      message("All Brier is NA. Returning NA.")
    }
    return(NA)
  }

  lst_ggp <- list()

  lst_plots <- comboplot.performance2.0(df = eval_results$df,
                                        x.var = ifelse(evaluation=="AUC", "time", "brier_time"),
                                        y.var = evaluation,
                                        y.lab = ifelse(evaluation=="AUC", "AUC", "Brier"),
                                        x.color = "method",
                                        legend_title = legend_title,
                                        y.limit = c(y.min, 1), pred.attr = pred.attr,
                                        legend_size_text = legend_size_text, x_axis_size_text = x_axis_size_text, y_axis_size_text = y_axis_size_text)
  if(type == "both"){
    lst_ggp <- lst_plots
  }else if(type == "line"){
    lst_ggp <- lst_plots$lineplot
  }else if(type == "mean"){
    lst_ggp <- lst_plots$lineplot.mean
  }

  lst_tests <-c("t.test", "anova","wilcox.test", "kruskal.test", "NULL")
  lst_plot_comparisons <- list()
  for(t in 1:length(lst_tests)){

    if(lst_tests[[t]]!="NULL"){
      test_comparations = lst_tests[[t]]
    }else{
      test_comparations = NULL
    }

    plot <- boxplot.performance(df = eval_results$df,
                                x.var = "method",
                                y.var = evaluation,
                                x.fill = "method",
                                x.alpha = NULL,
                                alpha.lab = NULL,
                                x.lab = "Method",
                                y.lab = ifelse(evaluation=="AUC", "AUC", "Brier Score"),
                                fill.lab = NULL,
                                title = paste0("Method Performance"),
                                y.limit = NULL,
                                y.limit.exception = NULL,
                                jitter = F,
                                test = test_comparations,
                                show.median = T,
                                round.median = 3,
                                legend_title = legend_title,
                                legend_size_text = legend_size_text,
                                x_axis_size_text = x_axis_size_text,
                                y_axis_size_text = y_axis_size_text)

    if(lst_tests[[t]] == "NULL"){
      lst_plot_comparisons[["no_test"]] <- plot
    }else{
      lst_plot_comparisons[[lst_tests[[t]]]] <- plot
    }

  }

  table <- NULL
  for(m in unique(eval_results$df$method)){
    for(c in colnames(eval_results$df)){
      if(c=="method" | c=="time" | c=="brier_time"){
        next
      }else{
        vector <- c(m, c,
                    mean(eval_results$df[eval_results$df$method==m,c,drop=T]),
                    median(eval_results$df[eval_results$df$method==m,c,drop=T]),
                    sd(eval_results$df[eval_results$df$method==m,c,drop=T]))
        table <- rbind(table, vector)
      }
    }
  }

  table <- as.data.frame(table)
  rownames(table) <- NULL
  colnames(table) <- c("method","metric","mean","median","sd")

  table$method <- factor(table$method)
  table$metric <- factor(table$metric)
  table$mean <- as.numeric(table$mean)
  table$median <- as.numeric(table$median)
  table$sd <- as.numeric(table$sd)

  return(list("lst_plots" = lst_ggp, "lst_plot_comparisons" = lst_plot_comparisons, df = table))
}

coxweightplot.fromVector.HDcox <- function(model, vector, sd.min = NULL, sd.max = NULL, zero.rm = F, top = NULL, auto.limits = T,
                                           block = NULL, show_percentage = T, size_percentage = 3){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NA)
  }

  #DFCALLS
  variables <- pp <- NULL

  loading_values <- vector
  ggp_loading <- NULL
  lst_top_loadings <- NULL
  lst_all_loadings <- NULL
  df <- NULL
  limit_color = 300

  #accuracy <- ifelse(max(vector)-min(vector) < 0.15, 0.01, 0.1)
  accuracy <- 0.1

  if(auto.limits){
    if(!is.null(sd.min) & !is.null(sd.max)){
      auto.limits_min <- round2any(max(abs(sd.min)), accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(max(abs(sd.max)), accuracy = accuracy, f = ceiling)
      auto.limits <- max(auto.limits_min, auto.limits_max)
    }else{
      auto.limits <- round2any(max(abs(loading_values)), accuracy = accuracy, f = ceiling)
    }
  }else{
    auto.limits <- round2any(max(c(abs(sd.max), abs(sd.min), abs(loading_values))), accuracy = accuracy, f = ceiling)
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

    #mid point 0 - cause we are working with coefficients instead of e^b
    ggp <- ggp +
      geom_bar(stat = "identity") +
      guides(color = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      #scale_fill_discrete(name = "New Legend Title") +
      xlab(label = paste0("Variables")) +
      ylab(label = paste0("Estimate Beta value"))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                        mid = "white", midpoint = 0,
                                        high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                        limits = c(-1*auto.limits,auto.limits), name = "Beta value")
    }else{
      ggp <- ggp + scale_fill_gradient2(low = "blue",
                                        mid = "white", midpoint = 0,
                                        high = "red",
                                        limits = c(-1*auto.limits,auto.limits), name = "Beta value")
    }

    #add total positive and negative values
    risk_t.val = sum(loading_values[loading_values>0,])
    preventive_t.val = sum(loading_values[loading_values<=0,])

    risk_val = sum(df[df$pp>0,]$pp)
    preventive_val = sum(df[df$pp<=0,]$pp)

    perc_risk = risk_t.val / (risk_t.val+abs(preventive_t.val))
    perc_preventive = abs(preventive_t.val) / (risk_t.val+abs(preventive_t.val))

    risk_explained = risk_val/risk_t.val*100
    preventive_explained = preventive_val/preventive_t.val*100

    if(is.nan(risk_explained)){
      risk_explained <- 0
    }
    if(is.nan(preventive_explained)){
      preventive_explained <- 0
    }

    total_explained = risk_explained*perc_risk + preventive_explained*perc_preventive

    if(!is.null(top)){
      if(top < nrow(loading_values)){
        txt.subtitle = paste0("Top ", top, " variables explain a ", round(total_explained, 2), " % of the model.")
      }else{
        #all variables selected
        txt.subtitle = paste0("Variables explain a ", round(total_explained, 2), " % of the model.")
      }

    }else{
      txt.subtitle = paste0("Variables explain a ", round(total_explained, 2), " % of the model.")
    }

    explained_perc = NULL
    for(value in df$pp){
      if(value>0){
        explained_perc = c(explained_perc, value / risk_t.val * perc_risk * 100)
      }else{
        explained_perc = c(explained_perc, abs(value) / abs(preventive_t.val) * perc_preventive * 100)
      }
    }
    df$explained = explained_perc

    df.all <- as.data.frame(loading_values)
    colnames(df.all) <- "value"
    explained_perc = NULL
    for(value in df.all$value){
      if(value>0){
        explained_perc = c(explained_perc, value / risk_t.val * perc_risk * 100)
      }else{
        explained_perc = c(explained_perc, abs(value) / abs(preventive_t.val) * perc_preventive * 100)
      }
    }
    df.all$perc.explained = explained_perc

    if(show_percentage & !is.null(top)){
      df$explained_text <- paste0(round(df$explained, 2), " %")
      ggp <- ggp + geom_text(aes(label = df$explained_text, y = sign(df$pp)*max(pp)*0.025), size = size_percentage)
    }

    if(is.null(block)){
      ggp <- ggp + ggtitle(paste0(attr(model, "model"), " - Survival Weight"), subtitle = txt.subtitle)
    }else{
      ggp <- ggp + ggtitle(paste0(attr(model, "model"), " - Survival Weight [", block, "]"), subtitle = txt.subtitle)
    }

    if(nrow(df)>limit_color){

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                           mid = "white", midpoint = 0,
                                           high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                           limits = c(-1*auto.limits,auto.limits), name = "Beta value")
      }else{
        ggp <- ggp + scale_color_gradient2("blue",
                                           mid = "white", midpoint = 0,
                                           high = "red",
                                           limits = c(-1*auto.limits,auto.limits), name = "Beta value")
      }

    }

    if(auto.limits){
      #ggp <- ggp + scale_y_continuous(breaks=seq(-1*auto.limits, auto.limits, 0.1))
      ggp <- ggp + scale_y_continuous(n.breaks = 10)
    }else{
      #ggp <- ggp + scale_y_continuous(breaks=seq(-1*auto.limits, auto.limits, 0.1), limits = c(-1*auto.limits, auto.limits))
      ggp <- ggp + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))
    }

    if(!is.null(sd.min) & !is.null(sd.max)){
      sd.min <- sd.min[rownames(df),,drop=F]
      sd.max <- sd.max[rownames(df),,drop=F]
      ggp <- ggp + geom_errorbar(aes(ymin=sd.min, ymax=sd.max), width=.35, position=position_dodge(.2))
    }

    if(ncol(loading_values)==1){
      return(list(plot = ggp, top_coefficients = df, coefficients = df.all))
    }

    ggp_loading[[i]] = ggp
    lst_top_loadings[[i]] <- df
    lst_all_loadings[[i]] <- df.all
  }
  names(ggp_loading) <- colnames(loading_values)
  names(lst_top_loadings) <- colnames(loading_values)
  names(lst_all_loadings) <- colnames(loading_values)
  return(list(plot = ggp_loading, top_coefficients = lst_top_loadings, coefficients = lst_all_loadings))
}



evalplot_errorbar <- function(df, x.var, y.var, y.var.sd, x.color = NULL, best_component = NULL, best_eta = NULL, x.text = "Component"){

  line_size = 1.25
  dot_size = 2.5
  error_width = 0.5
  error_pos = 0.15 #0.3
  error_size = 0.75
  best_flag = F

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    color_conesa <- RColorConesa::colorConesa(1)
  }else{
    color_conesa <- "blue"
  }


  if(!is.null(x.color) & !is.null(best_component) & !is.null(best_eta)){
    best_flag = T
    best_df <- df[df[,x.var] == best_component,,drop=F]
    best_df[!best_df[,x.color] == as.character(best_eta),c(y.var, y.var.sd)] <- NA #I need NA because is moved (position_dodge)
    #best_df <- best_df[best_df[,x.color] == as.character(best_eta),]
  }else if(!is.null(best_component)){
    best_flag = T
    best_df <- df[df[,x.var] == best_component,,drop=F]
  }

  #ROUND AUC VALUES - 3 decimal digits
  df[,y.var] <- round2any(df[,y.var], accuracy = 0.001)
  df[,y.var] <- round2any(df[,y.var], accuracy = 0.001)
  best_df[,y.var] <- round2any(best_df[,y.var], accuracy = 0.001)

  if(!is.null(x.color)){
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, color = x.color, group = x.color)) +
      geom_line(aes_string(x = x.var, y = y.var, color = x.color), size = line_size, position=position_dodge(error_pos)) +
      geom_point(aes_string(color = x.color), size = dot_size, position=position_dodge(error_pos)) +
      geom_errorbar(aes(ymin=df[,y.var]-df[,y.var.sd],
                        ymax=df[,y.var]+df[,y.var.sd],
                        x = df[,x.var],
                        color=df[,x.color]),
                    width=error_width,
                    size = error_size,
                    position=position_dodge(error_pos)) +
      scale_x_discrete(x.text, labels = df[,x.var], breaks = df[,x.var])

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_color_conesa(palette = "complete")
    }

  }else{
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var)) +
      geom_line(color = color_conesa, group = x.var, size = line_size) +
      geom_point(color = color_conesa, size = dot_size) +
      geom_errorbar(aes(ymin=df[,y.var]-df[,y.var.sd],
                        ymax=df[,y.var]+df[,y.var.sd],
                        x = df[,x.var]),
                    color=color_conesa,
                    width=error_width,
                    size = error_size,
                    position=position_dodge(error_pos)) +
      scale_x_discrete(x.text, labels = df[,x.var], breaks = df[,x.var])
  }

  if(best_flag){
    if(!is.null(x.color)){
      ggp <- ggp + geom_point(data = best_df, aes_string(x = x.var, y = y.var, color = x.color, group = x.color),
                              position=position_dodge(error_pos), size = dot_size, shape = 23, fill = "white",
                              stroke = 2, show.legend = F)
    }else{
      ggp <- ggp + geom_point(data = best_df, aes_string(x = x.var, y = y.var), position=position_dodge(error_pos),
                              size = dot_size, shape = 23, fill = "white", color = color_conesa,
                              stroke = 2, show.legend = F)
    }
  }

  return(ggp)
}

lineplot.performace2.0 <- function(df, x.var = "time", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, point = T, mean = F, legend_rm = T, legend_title = "Method", legend_size_text = 12, x_axis_size_text = 10, y_axis_size_text = 10){

  MAX_X_ELEMENTS = 20

  if(mean){
    mean_vector = NULL
    for(m in unique(df$method)){
      mean_vector <- c(mean_vector, colMeans(df[df$method==m,y.var,drop=F], na.rm = T))
    }
    names(mean_vector) <- unique(df$method)
    mean_vector <- data.frame(mean_vector)
    mean_vector$method <- rownames(mean_vector)
    mean_vector <- mean_vector[,c(2,1)]
    rownames(mean_vector) <- NULL
  }

  if(point){
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, color = x.color)) +
      geom_line(aes_string(group = x.color), size = 1) +
      geom_point() +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_color_conesa(palette = "complete")
    }

  }else{
    ggp <- ggplot2::ggplot(df, aes_string(x = x.var, y = y.var, color = x.color)) +
      geom_line(aes_string(group = x.color), size = 1) +
      xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
      ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab))

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_color_conesa(palette = "complete")
    }

  }

  if(length(unique(df[,x.var,drop=T]))>MAX_X_ELEMENTS){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_size_text))
  }else{
    ggp <- ggp + theme(axis.text.x = element_text(vjust = 0.5, size = x_axis_size_text))
  }

  ggp <- ggp + theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = y_axis_size_text))

  # if(!is.null(y.limit)){
  #   ggp <- ggp + ylim(y.limit)
  # }

  if(mean){
    ggp <- ggp + geom_hline(data = mean_vector, aes_string(yintercept = mean_vector$mean_vector, color = x.color), size = 1)
  }

  if(legend_rm){
    ggp <- ggp + theme(legend.position = "none")
  }else{
    ggp <- ggp + theme(legend.text=element_text(size = legend_size_text), legend.title = element_text(size=legend_size_text, face = "bold"))
    ggp <- ggp + guides(color=guide_legend(title=legend_title))
  }

  if(!is.null(y.limit)){
    ggp <- ggp + scale_y_continuous(limits = y.limit,
                                    minor_breaks = seq(y.limit[1], y.limit[2], 0.05),
                                    labels = as.character(format(seq(y.limit[1], y.limit[2], 0.05), nsmall = 2)),
                                    breaks = seq(y.limit[1], y.limit[2], 0.05))
  }else{
    ggp <- ggp + scale_y_continuous(minor_breaks = seq(floor(min(df$AUC)*10)/10, ceiling(max(df$AUC)*10)/10, 0.05),
                                    labels = as.character(format(seq(floor(min(df$AUC)*10)/10, ceiling(max(df$AUC)*10)/10, 0.05), nsmall = 2)),
                                    breaks = length(seq(round(min(df$AUC)*10)/10, ceiling(max(df$AUC)*10)/10, 0.05)))
  }


  return(ggp)
}

barplot.mean_performace2.0 <- function(df, x.var = "method", y.var="AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, hide_labels = T, legend_rm = NULL, legend_title = "Method", legend_size_text = 12, x_axis_size_text = 10, y_axis_size_text = 10){

  #DFCALLS
  MAX_X_ELEMENTS = 20
  method <- NULL

  mean_vector = NULL
  for(m in unique(df$method)){
    mean_vector <- c(mean_vector, colMeans(df[df$method==m,y.var,drop=F], na.rm = T))
  }
  names(mean_vector) <- unique(df$method)
  mean_vector <- data.frame(mean_vector)
  mean_vector$method <- rownames(mean_vector)
  mean_vector <- mean_vector[,c(2,1)]
  rownames(mean_vector) <- NULL

  mean_vector <- mean_vector[order(mean_vector$mean_vector, decreasing = T),]
  #mean_vector$method <- factor(mean_vector$method, levels = mean_vector$method)

  ggp <- ggplot2::ggplot(mean_vector, aes(x = reorder(method, -mean_vector), y = mean_vector, fill = method, color = method)) +
    #geom_col(position = "identity", size = 0.5) +
    geom_point(position = "identity", size = 2) +
    xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
    ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab))

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "complete") + RColorConesa::scale_color_conesa(palette = "complete")
  }

  if(legend_rm){
    ggp <- ggp + theme(legend.position = "none")
  }else{
    ggp <- ggp + theme(legend.text=element_text(size = legend_size_text), legend.title = element_text(size=legend_size_text, face = "bold"))
    ggp <- ggp + guides(color=guide_legend(title=legend_title))
  }

  if(length(unique(df[,x.var,drop=T]))>MAX_X_ELEMENTS){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_size_text))
  }else{
    ggp <- ggp + theme(axis.text.x = element_text(vjust = 0.5, size = x_axis_size_text))
  }
  ggp <- ggp + theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = y_axis_size_text))

  if(!is.null(y.limit)){
    ggp <- ggp + coord_cartesian(ylim = y.limit)
  }

  if(hide_labels){
    ggp <- ggp +  ylab("") + xlab("") + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
  }

  return(ggp)
}

point.sd.mean_performace2.0 <- function(df, x.var = "method", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, pred.attr = "mean", hide_labels = T, legend_rm = NULL, legend_title = "Method", legend_size_text = 12, x_axis_size_text = 10, y_axis_size_text = 10){

  #DFCALLS
  MAX_X_ELEMENTS = 20
  method <- NULL

  mean_vector = NULL
  sd_vector = NULL
  for(m in unique(df$method)){
    if(pred.attr %in% "mean"){
      mean_vector <- c(mean_vector, colMeans(df[df$method==m,y.var,drop=F], na.rm = T))
    }else if(pred.attr %in% "median"){
      mean_vector <- c(mean_vector, apply(df[df$method==m,y.var,drop=F], 2, function(x){median(x, na.rm = T)}))
    }
    sd_vector <- c(sd_vector, sd(df[df$method==m,y.var,drop=F][[y.var]], na.rm = T))
  }
  sd_vector[is.na(sd_vector)] <- 0 #if NA is because we do not have sd for that vector of AUC
  names(mean_vector) <- unique(df$method)
  mean_vector <- data.frame(mean_vector)
  mean_vector$method <- rownames(mean_vector)
  mean_vector <- mean_vector[,c(2,1)]
  rownames(mean_vector) <- NULL
  mean_vector$sd <- sd_vector

  min <- round2any(min(mean_vector$mean_vector-mean_vector$sd) * 10, 0.5, floor) / 10
  max <- round2any(max(mean_vector$mean_vector+mean_vector$sd) * 10, 0.5, ceiling) / 10

  mean_vector$method <- factor(x = mean_vector$method, levels = levels(df$method))

  #mean_vector <- mean_vector[order(mean_vector$mean_vector, decreasing = T),]
  #mean_vector$method <- factor(mean_vector$method, levels = mean_vector$method)

  ggp <- ggplot2::ggplot(mean_vector, aes(x = reorder(method, -mean_vector), y = mean_vector, fill = method, color = method)) +
    #geom_col(position = "identity", size = 0.5) +
    geom_point(position = "identity", size = 2.5) +
    xlab(ifelse(is.null(x.lab), x.var, x.lab)) +
    ylab(ifelse(is.null(y.lab),toupper(y.var),y.lab)) +
    geom_errorbar(aes(ymin=mean_vector-sd, ymax=mean_vector+sd), width=.4, size = 1.25,
                  position=position_dodge(.9))

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "complete") + RColorConesa::scale_color_conesa(palette = "complete")
  }

  if(legend_rm){
    ggp <- ggp + theme(legend.position = "none")
  }else{
    ggp <- ggp + theme(legend.text=element_text(size = legend_size_text), legend.title = element_text(size=legend_size_text, face = "bold"))
    ggp <- ggp + guides(color=guide_legend(title=legend_title), fill="none")
  }

  if(length(unique(df[,x.var,drop=T]))>MAX_X_ELEMENTS){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_size_text))
  }else{
    ggp <- ggp + theme(axis.text.x = element_text(vjust = 0.5, size = x_axis_size_text))
  }
  ggp <- ggp + theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = y_axis_size_text))


  if(!is.null(y.limit)){
    ggp <- ggp + coord_cartesian(ylim = y.limit)
  }

  if(hide_labels){
    ggp <- ggp +  ylab("") + xlab("") + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
  }

  #ggp <- ggp + scale_y_continuous(minor_breaks = seq(min, max, 0.5), n.breaks = seq(min, max, 0.5))

  if(is.nan(min)){
    min = 0
  }
  if(is.na(max)){
    max = 1
  }

  ggp <- ggp + theme(panel.grid.major.y = element_blank())
  ggp <- ggp + theme(panel.grid.major.x = element_blank())
  ggp <- ggp + scale_y_continuous(minor_breaks = seq(min, max, 0.05),
                                  #labels = as.character(format(seq(min, max, 0.05), nsmall = 2)),
                                  breaks = length(seq(min, max, 0.05)))
  #
  # ggp <- ggp + theme(panel.grid.minor = element_blank(), )
  # ggp <- ggp + xlab("AUC median per method")

  return(ggp)
}

comboplot.performance2.0 <- function(df, x.var = "time", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, pred.attr = "mean", point = T, mean = F, hide_labels = T, legend_title = "Method", legend_size_text = 12, x_axis_size_text = 10, y_axis_size_text = 10){
  a <- lineplot.performace2.0(df = df, x.var = x.var, y.var = y.var, x.color = x.color, x.lab = x.lab, y.lab = y.lab, y.limit = y.limit, point = point, mean = F, legend_rm = T, legend_title = legend_title, legend_size_text = legend_size_text, x_axis_size_text = x_axis_size_text, y_axis_size_text = y_axis_size_text)
  b <- point.sd.mean_performace2.0(df = df, x.var = x.var, y.var = y.var, x.color = x.color, x.lab = NULL, y.lab = NULL, y.limit = y.limit, pred.attr = pred.attr, hide_labels = T, legend_rm = F, legend_title = legend_title, legend_size_text = legend_size_text, x_axis_size_text = x_axis_size_text, y_axis_size_text = y_axis_size_text)

  pp <- ggpubr::ggarrange(a, b, ncol = 2, widths = c(0.8, 0.2), align = "h")

  a <- lineplot.performace2.0(df, x.var, y.var, x.color, x.lab, y.lab, y.limit, point, mean = F, legend_rm = F, legend_title = legend_title, legend_size_text = legend_size_text, x_axis_size_text = x_axis_size_text, y_axis_size_text = y_axis_size_text)
  return(list(lineplot = a, lineplot.mean = pp))
}

plot_VAR_eval <- function(lst_BV, EVAL_METHOD = "AUC", dot_size = 3){
  values = NULL #just in case
  best_keepX <- lst_BV$best.keepX
  best_keepX <- paste0(unlist(lapply(best_keepX, function(x){x[[1]]})), collapse = "_")
  df.pval <- data.frame(names = factor(names(lst_BV$p_val), levels = names(lst_BV$p_val)), values = lst_BV$p_val)
  if(EVAL_METHOD == "BRIER"){
    df.pval$values <- 1- df.pval$values
  }

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    color_conesa <- RColorConesa::colorConesa(1)
  }else{
    color_conesa <- "blue"
  }

  ggp <- ggplot(df.pval, aes(x = names, y = values)) +
    geom_line(group = 1, color = color_conesa, linewidth = 1.5) + ylab("Pred. Value") + xlab("Number of variables")
  ggp <- ggp + geom_point(data = df.pval[df.pval$names==best_keepX,,drop=F],
                          aes(x = names, y = values), color = color_conesa,
                          size = dot_size, shape = 23, fill = "white",
                          stroke = 2, show.legend = F)

  return(ggp)
}

#### ### ### ### ##
# EVENT PLOTS - Y #
#### ### ### ### ##

#' plot_events
#'
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.breaks Numeric. Maximum number of breaks in X axis (default: 20).
#' @param roundTo Numeric. Value to round time. If roundTo = 0.1, the results will be rounded to the tenths (default: 0.1).
#' @param categories Character vector. Vector of length two to name both categories for censored and non-censored observations (default: c("Censored","Death")).
#' @param y.text Character. Y axis title (default: "Number of observations").
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_events(Y, categories = c("Censored","Death"))
#' }

plot_events <- function(Y, max.breaks = 20, roundTo = 0.1, categories = c("Censored","Death"), y.text = "Number of observations", verbose = F){

  #REQUIREMENTS
  if(length(categories)>2 | length(categories)<2 | !is.character(categories)){
    stop("categories parameter must be a character vector of length two.")
  }

  if(!is.character(y.text) | length(y.text)>1){
    stop("y.text parameter must be a character vector of length one.")
  }

  if(!is.numeric(roundTo)){
    stop("roundTo parameter must be a numeric vector of length one.")
  }

  if(!is.numeric(max.breaks)){
    stop("max.breaks parameter must be a numeric vector of length one.")
  }

  if(roundTo == 0){
    #select the decimals of Y
    if(length(grep("\\.", Y$time))>0){
      roundTo = 1*10^-(nchar(gsub("\\.", "", as.character(Y[,"time"][[1]])))-1)
    }else{
      roundTo = 0.1
    }

  }

  #DFCALLS
  Y <- as.data.frame(Y)
  Category <- Time <- Values <- x.names <- breaks<- NULL

  if(!is.logical(Y[,"event"])){
    if(verbose){
      message("Y matrix must has event column as TRUE, FALSE. as.logical() function has been used.")
    }
    Y[,"event"] <- as.logical(Y[,"event"])
  }

  breaks_size = round2any((max(Y[,"time"]) - min(Y[,"time"])) / (max.breaks+1), roundTo, f = ceiling)
  breaks = seq(min(Y[,"time"]), max(Y[,"time"])+breaks_size, by=breaks_size)
  breaks = round2any(breaks, roundTo, f = floor)
  if(max(breaks)<max(Y[,"time"])){breaks=c(breaks, max(breaks)+breaks_size)}
  x.names <- cut(x = Y[,"time"], breaks = breaks, include.lowest = T)

  Y <- cbind(Y, "time_g" = x.names)

  vt=NULL
  vcategory=NULL
  vvalues=NULL
  for(t in levels(x.names)){
    vt <- c(vt, t, t)
    vcategory <- c(vcategory, categories)
    vvalues<- c(vvalues, sum(Y[Y[,"time_g"]==t, "event"]==F), sum(Y[Y[,"time_g"]==t, "event"]==T))
  }

  dd <- data.frame(Time=vt, Category=vcategory, Values=vvalues)
  dd$Time <- factor(dd$Time, levels = levels(x.names))

  #check last group
  if(all(dd$Values[c(length(dd$Values)-1,length(dd$Values))]==0)){
    dd <- dd[-c(length(dd$Values)-1,length(dd$Values)),]
    dd <- droplevels.data.frame(dd)
  }

  ggp_density <- ggplot(dd, aes(fill=Category, x=Time, y=Values)) +
    #geom_bar(position="stack", stat="identity") +
    geom_bar(stat = "identity") +
    ylab(y.text) +
    scale_y_continuous(n.breaks = 10) +
    guides(fill=guide_legend(title="Class"), color = "none")

  if(length(levels(x.names))>15){
    ggp_density <- ggp_density + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp_density <- ggp_density + RColorConesa::scale_fill_conesa()
  }

  return(list(plot = ggp_density, df = dd))
}

#' plot_divergent.biplot
#' @description Two side plot by a qualitative and a quantitative variable and Y event matrix.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables with "NAMEVAR1" and "NAMEVAR2" variables. "NAMEVAR1" must be a factor variable.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param NAMEVAR1 Character. Factor variable name (must be located in colnames(X) and have to have two levels).
#' @param NAMEVAR2 Character. Numerical variable name (must be located in colnames(X)).
#' @param BREAKTIME Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points. If BREAKTIME = NULL, "n.breaks" is used (default: NULL).
#' @param x.text Character. Title for X axis.
#'
#' @return A ggplot2 two side bar plot. X axis represent the number of samples per each NAMEVAR1 factor levels and Y axis, the X NAMEVAR2 numerical variables categorize in groups of breaks.
#' @export
#'
#' @details For time comparison between best models, use the cross validation objects instead the final models. The training has to be taken into account.
#'
#' @examples
#' \dontrun{
#'   NAMEVAR1 = "sex"
#'   NAMEVAR2 = "age"
#'   plot_divergent.biplot(X, Y, NAMEVAR1, NAMEVAR2, BREAKTIME = 5, x.text = "N. of Patients")
#' }

plot_divergent.biplot <- function(X, Y, NAMEVAR1, NAMEVAR2, BREAKTIME, x.text = "N. of Samples"){
  df<-NULL

  VAR1 <- X[rownames(X), NAMEVAR1] #will be a factor
  if(!is.factor(VAR1)){
    VAR1 <- factor(VAR1)
  }
  VAR2 <- X[rownames(X), NAMEVAR2] #must be numerical

  OUTCOME <- Y[rownames(X),"event"]

  df <- data.frame(VAR1, VAR2, OUTCOME)  # merge by row names (by=0 or by="row.names")
  colnames(df) <- c(NAMEVAR1, NAMEVAR2, "event")

  #add age as category
  df.index <- NULL
  cat <- NULL
  index <- NULL

  BREAKTIME = BREAKTIME
  min <- round2any(min(VAR2), accuracy = BREAKTIME, f = floor)
  max <- round2any(max(VAR2), accuracy = BREAKTIME, f = ceiling)
  for(i in seq(min,max,BREAKTIME)){
    if(i!=max){
      new <- which(df[,NAMEVAR2]>=i & df[,NAMEVAR2]<=(i+BREAKTIME-1))
      index <- c(index, new)
      cat <- c(cat, rep(paste0(i,"-",i+BREAKTIME-1), length(new)))
    }else{
      new <- which(df[,NAMEVAR2]>=i)
      index <- c(index, new)
      cat <- c(cat, rep(paste0(i, "<="), length(new)))
    }
  }
  df.index <- as.data.frame(index)
  df.index$cat <- cat
  df.index <- df.index[order(df.index$index),]

  df$cat <- factor(df.index$cat, levels = unique(cat))
  df[,NAMEVAR1] <- factor(df[,NAMEVAR1])

  value_cat <- NULL
  value_var1 <- NULL
  num_event <- NULL
  name_event <- NULL

  dim <- length(levels(df[,NAMEVAR1])) * length(unique(df[,"event"]))

  for(i in levels(df$cat)){
    value_cat <- c(value_cat,rep(i, dim))
    value_var1 <- c(value_var1, rep(levels(df[,NAMEVAR1]), length(unique(df[,"event"]))))

    for(j in levels(df[,NAMEVAR1])){
      num_event<- c(num_event, sum(df[df$cat==i & df[,NAMEVAR1]==j, "event"]==1))
      name_event <- c(name_event, "Event")
    }
    for(j in levels(df[,NAMEVAR1])){
      num_event<- c(num_event, sum(df[df$cat==i & df[,NAMEVAR1]==j, "event"]==0))
      name_event <- c(name_event, "Censored")
    }
  }

  df.final <- data.frame(value_cat,value_var1,num_event,name_event)
  df.final$value_cat <- factor(df.final$value_cat, levels = unique(cat))

  #to divergent graph we need negative values
  #NAMEVAR1 must be a two length factor
  class2 <- which(df.final$value_var1== levels(df[,NAMEVAR1])[1])
  df.final[class2,]$num_event <- df.final[class2,]$num_event*-1

  breaks_values <- pretty(df.final$num_event)

  real_center_deviation <- abs(mean(breaks_values)) / sum(abs(breaks_values))

  ggp_distribution <- df.final %>%
    ggplot(aes(x = value_cat, y = num_event, fill = name_event))+
    geom_bar(position="dodge", stat="identity")+
    coord_flip() +
    geom_hline(yintercept = 0, color="white") +
    ggtitle(paste0(NAMEVAR1,"_",levels(df[,NAMEVAR1])[1], " vs ", NAMEVAR1,"_",levels(df[,NAMEVAR1])[2])) +
    ylab("N. of Patients") + xlab(paste0(NAMEVAR2)) +
    scale_y_continuous(breaks = breaks_values,
                       labels = abs(breaks_values)) +
    guides(fill=guide_legend(title="Event type")) +
    theme(plot.title = element_text(hjust = 0.5 + round2any(real_center_deviation, 0.01))) +
    xlab(label = x.text)

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp_distribution <- ggp_distribution + RColorConesa::scale_fill_conesa()
  }

  return(ggp_distribution)
}

#### ### ### ### ### ### ###
# PLS PLOTS - HDCOX MODELS #
#### ### ### ### ### ### ###

#' plot_PLS_HDcox
#'
#' @param model HDcox model.
#' @param comp Numeric vector. Vector of length two. Select which components to plot (default: c(1,2)).
#' @param mode Character. Choose one of the following plots: "scores", "loadings" o "biplot" (default: "scores").
#' @param factor Factor. Factor variable to color the observations. If factor = NULL, event will be used (default: NULL).
#' @param legend_title Character. Legend title (default: NULL).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param only_top Logical. If "only_top" = TRUE, then only top/radius loading variables will be shown in loading or biplot graph (default: FALSE).
#' @param radius Numeric. Radius size (loading/scale value) to plot variable names that are greater than the radius value (default: NULL).
#' @param names Logical. Show loading names for top variables or for those that are outside the radius size (default: TRUE).
#' @param colorReverse Logical. Reverse palette colors (default: FALSE).
#' @param text.size Numeric. Text size (default: 2).
#' @param overlaps Numeric. Number of overlaps to show when plotting loading names (default: 10).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_PLS_HDcox(model, comp = c(1,2), mode = "scores")
#' }

plot_PLS_HDcox <- function(model, comp = c(1,2), mode = "scores", factor = NULL, legend_title = NULL, top = NULL, only_top = F, radius = NULL, names = T, colorReverse = F, text.size = 2, overlaps = 10){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NA)
  }

  if(attr(model, "model") %in% pkg.env$pls_methods){
    plot_HDcox.PLS.model(model = model,
                         comp = comp,
                         mode = mode,
                         factor = factor,
                         legend_title = legend_title,
                         top = top, only_top = only_top,
                         radius = radius, names = names,
                         colorReverse = colorReverse, text.size = text.size,
                         overlaps = overlaps)

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){
    plot_HDcox.MB.PLS.model(model = model,
                            comp = comp,
                            mode = mode,
                            factor = factor,
                            legend_title = legend_title,
                            top = top, only_top = only_top,
                            radius = radius, names = names,
                            colorReverse = colorReverse, text.size = text.size,
                            overlaps = overlaps)
  }else{
    stop("Model must be a PLS HDcox model.")
  }

}

#' plot_HDcox.PLS.model
#'
#' @param model HDcox model.
#' @param comp Numeric vector. Vector of length two. Select which components to plot (default: c(1,2)).
#' @param mode Character. Choose one of the following plots: "scores", "loadings" o "biplot" (default: "scores").
#' @param factor Factor. Factor variable to color the observations. If factor = NULL, event will be used (default: NULL).
#' @param legend_title Character. Legend title (default: NULL).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param only_top Logical. If "only_top" = TRUE, then only top/radius loading variables will be shown in loading or biplot graph (default: FALSE).
#' @param radius Numeric. Radius size (loading/scale value) to plot variable names that are greater than the radius value (default: NULL).
#' @param names Logical. Show loading names for top variables or for those that are outside the radius size (default: TRUE).
#' @param colorReverse Logical. Reverse palette colors (default: FALSE).
#' @param text.size Numeric. Text size (default: 2).
#' @param overlaps Numeric. Number of overlaps to show when plotting loading names (default: 10).

plot_HDcox.PLS.model <- function(model, comp = c(1,2), mode = "scores", factor = NULL, legend_title = NULL, top = NULL, only_top = F, radius = NULL, names = T, colorReverse = F, text.size = 2, overlaps = 10){

  MAX_POINTS = 1000
  MAX_LOADINGS = 15
  POINT_SIZE = 3
  POINT_SIZE_LOAD = 1.5 #another scale
  POINT_RES = c(1024, 1024)

  ggp = NULL
  aux.model = model

  if(!is.null(top) & !is.null(radius)){
    message("Only top meassure will be used. Radius and top do not work simultaneously.")
    radius <- NULL
  }

  modes <- c("scores", "loadings", "biplot")
  if(!mode %in% modes){
    stop_quietly(paste0("mode must be one of the following: ", paste0(modes, collapse = ", ")))
  }

  if(!is.null(factor)){
    if(!is.factor(factor) & mode %in% c("scores", "biplot")){
      stop_quietly("Factor must be a factor object.")
    }
  }else{
    factor <- factor(model$Y$data[,"event"])
  }

  if(!isa(aux.model,pkg.env$model_class)){
    stop_quietly("'model' must be a HDcox object.")
  }else if(attr(aux.model, "model") %in% c(pkg.env$multiblock_methods)){
    stop_quietly("For single block models, use the function 'plot_HDcox.MB.PLS.model'")
  }else if(!attr(aux.model, "model") %in% c(pkg.env$pls_methods, pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){
    stop_quietly("'model' must be a HDcox object PLS class ('sPLS-ICOX','sPLS-DRCOX','sPLS-DRCOX-Dynamic', or 'sPLS-DACOX-Dynamic'.")
  }

  if(mode=="scores"){

    if(ncol(aux.model$X$scores)==1){
      message("The model has only 1 component")

      df <- cbind(aux.model$X$scores[,1], aux.model$X$scores[,1])
      colnames(df) <- c("p1", "p2")
    }else{
      df <- as.data.frame(aux.model$X$scores)
    }

    subdata_loading = NULL
    ggp <- ggplot(as.data.frame(df))

    if(nrow(df) > MAX_POINTS){
      ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor), pointsize = POINT_SIZE, pixels = POINT_RES)
    }else{
      ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor))
    }

    ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)
    ggp <- ggp + stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F)

    if("R2" %in% names(model)){
      txt.expression <- paste0("Scores (",attr(aux.model, "model"),") - ")
      r2_1 <- round(model$R2[[comp[1]]], 4)
      r2_2 <- round(model$R2[[comp[2]]], 4)
      r2 <- round(sum(unlist(model$R2)), 4)
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
    }else{
      txt.expression <- paste0("Scores (",attr(aux.model, "model"),")")
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
        xlab(label = paste0("comp_",as.character(comp[1]))) +
        ylab(label = paste0("comp_",as.character(comp[2])))
    }

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp +
        RColorConesa::scale_color_conesa(reverse = colorReverse) +
        RColorConesa::scale_fill_conesa(reverse = colorReverse)
    }

  }else if(mode=="loadings"){

    if(ncol(aux.model$X$loadings)==1){
      message("The model has only 1 component")

      df <- as.data.frame(cbind(aux.model$X$loadings,aux.model$X$loadings))
      colnames(df) <- c("p1", "p2")

    }else{
      df <- as.data.frame(aux.model$X$loadings)
    }

    if(nrow(df)<MAX_LOADINGS){
      subdata_loading <- df
    }else if(!is.null(top)){
      aux_loadings <- apply(df,1,function(x){sqrt(crossprod(as.numeric(x[comp])))})
      aux_loadings <- aux_loadings[order(aux_loadings, decreasing = T)]
      subdata_loading <- df[names(aux_loadings)[1:top],]
    }else if(!is.null(radius)){
      subdata_loading <- df_loading[apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))>radius}),]
    }else{
      subdata_loading <- NULL
    }

    ggp <- ggplot(as.data.frame(df))

    if(nrow(df) > MAX_POINTS){
      ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]]), pointsize = POINT_SIZE, pixels = POINT_RES)
    }else{
      ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]]))
    }

    ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)

    if("R2" %in% names(model)){
      txt.expression <- paste0("Loadings (",attr(aux.model, "model"),") - ")
      r2_1 <- round(model$R2[[comp[1]]], 4)
      r2_2 <- round(model$R2[[comp[2]]], 4)
      r2 <- round(sum(r2_1, r2_2), 4)
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
    }else{
      txt.expression <- paste0("Loadings (",attr(aux.model, "model"),")")
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
        xlab(label = paste0("comp_",as.character(comp[1]))) +
        ylab(label = paste0("comp_",as.character(comp[2])))
    }

    if(names & !is.null(subdata_loading)){
      ggp <- ggp + ggrepel::geom_text_repel(data = subdata_loading, aes(x = subdata_loading[,comp[1]],
                                                                        y = subdata_loading[,comp[2]]),
                                            max.overlaps = getOption("ggrepel.max.overlaps", default = overlaps),
                                            label = rownames(subdata_loading), size=text.size)
    }

    if(!is.null(radius) & !is.null(subdata_loading)){
      if(requireNamespace("ggforce", quietly = TRUE)){
        ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
      }
    }

  }else if(mode=="biplot"){
    if(ncol(aux.model$X$loadings)==1){
      message("The model has only 1 component")

      df <- as.data.frame(cbind(aux.model$X$scores, aux.model$X$scores))
      colnames(df) <- c("p1", "p2")

      df_loading <- as.data.frame(cbind(aux.model$X$loadings[,1], aux.model$X$loadings[,1]))
      max.loadings <- apply(abs(df_loading), 2, max)
      max.scores <- apply(abs(df), 2, max)
    }else{
      df <- as.data.frame(aux.model$X$scores)

      df_loading <- as.data.frame(aux.model$X$loadings)
      max.loadings <- apply(abs(df_loading), 2, max)
      max.scores <- apply(abs(df), 2, max)
    }

    #scale scores to -1,1
    df <- norm01(df[,comp])*2-1
    ggp <- ggplot(as.data.frame(df))

    if(nrow(df) > MAX_POINTS){
      ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor), pointsize = POINT_SIZE, pixels = POINT_RES)
    }else{
      ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor))
    }

    ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)
    ggp <- ggp + stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F)

    if("R2" %in% names(model)){
      txt.expression <- paste0("Biplot (",attr(aux.model, "model"),") - ")
      r2_1 <- round(model$R2[[comp[1]]], 4)
      r2_2 <- round(model$R2[[comp[2]]], 4)
      r2 <- round(sum(r2_1, r2_2), 4)
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
    }else{
      txt.expression <- paste0("Biplot (",attr(aux.model, "model"),")")
      ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
        xlab(label = paste0("comp_",as.character(comp[1]))) +
        ylab(label = paste0("comp_",as.character(comp[2])))
    }

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp +
        RColorConesa::scale_color_conesa(reverse = colorReverse) +
        RColorConesa::scale_fill_conesa(reverse = colorReverse)
    }

    if(nrow(df_loading)<MAX_LOADINGS){
      subdata_loading <- df_loading
    }else if(!is.null(top)){
      aux_loadings <- apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))})
      aux_loadings <- aux_loadings[order(aux_loadings, decreasing = T)]
      subdata_loading <- df_loading[names(aux_loadings)[1:top],]
    }else if(!is.null(radius)){
      subdata_loading <- df_loading[apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))>radius}),]
    }else{
      subdata_loading <- NULL
    }

    #depending on DF instead of df_loadings - ARROWS
    if(any(!is.null(top), !is.null(radius))){

      no_selected_loadings <- df_loading[!rownames(df_loading) %in% rownames(subdata_loading),]
      if(nrow(no_selected_loadings)!=0 & !only_top){
        ggp <- ggp + geom_segment(data = no_selected_loadings, lineend = "butt", linejoin = "mitre", size = 0.2,
                                  aes(x = 0, y = 0, xend = no_selected_loadings[,comp[1]],
                                      yend = no_selected_loadings[,comp[2]]),
                                  arrow = arrow(length = unit(0.1, "cm")))
      }

      ggp <- ggp + geom_segment(data = subdata_loading, lineend = "butt", linejoin = "mitre",
                                size = 0.33, aes(x = 0, y = 0, xend = subdata_loading[,comp[1]],
                                                 yend = subdata_loading[,comp[2]]),
                                arrow = arrow(length = unit(0.1, "cm")))

    }else{
      #show all loadings
      no_selected_loadings <- df_loading[!rownames(df_loading) %in% rownames(subdata_loading),]
      ggp <- ggp + geom_segment(data = no_selected_loadings, lineend = "butt", linejoin = "mitre", size = 0.2,
                                aes(x = 0, y = 0, xend = no_selected_loadings[,comp[1]],
                                    yend = no_selected_loadings[,comp[2]]),
                                arrow = arrow(length = unit(0.1, "cm")))
    }

    if(names & !is.null(subdata_loading)){
      ggp <- ggp + ggrepel::geom_text_repel(data = subdata_loading, aes(x = subdata_loading[,comp[1]],
                                                                        y = subdata_loading[,comp[2]]),
                                            max.overlaps = getOption("ggrepel.max.overlaps", default = overlaps),
                                            label = rownames(subdata_loading), size=text.size)
    }

    if(is.null(top) & !is.null(radius) & nrow(df) < MAX_POINTS){
      ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
    }

  }

  #reorder legend
  if(!is.null(factor) & length(levels(factor))>3){
    ggp <- ggp + guides(color=guide_legend(nrow = ceiling(length(levels(factor))/3), byrow = T))
  }

  return(list(plot = ggp, outliers = rownames(subdata_loading)))
}

#' plot_HDcox.MB.PLS.model
#'
#' @param model HDcox model.
#' @param comp Numeric vector. Vector of length two. Select which components to plot (default: c(1,2)).
#' @param mode Character. Choose one of the following plots: "scores", "loadings" o "biplot" (default: "scores").
#' @param factor Factor. Factor variable to color the observations. If factor = NULL, event will be used (default: NULL).
#' @param legend_title Character. Legend title (default: NULL).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#' @param only_top Logical. If "only_top" = TRUE, then only top/radius loading variables will be shown in loading or biplot graph (default: FALSE).
#' @param radius Numeric. Radius size (loading/scale value) to plot variable names that are greater than the radius value (default: NULL).
#' @param names Logical. Show loading names for top variables or for those that are outside the radius size (default: TRUE).
#' @param colorReverse Logical. Reverse palette colors (default: FALSE).
#' @param text.size Numeric. Text size (default: 2).
#' @param overlaps Numeric. Number of overlaps to show when plotting loading names (default: 10).

plot_HDcox.MB.PLS.model <- function(model, comp = c(1,2), mode = "scores", factor = NULL, legend_title = NULL, top = NULL, only_top = F, radius = NULL, names = T, colorReverse = F, text.size = 2, overlaps = 10){

  MAX_POINTS = 1000
  MAX_LOADINGS = 15
  POINT_SIZE = 3
  POINT_SIZE_LOAD = 1.5 #another scale
  POINT_RES = c(1024, 1024)

  ggp = NULL
  aux.model = model

  if(!is.null(top) & !is.null(radius)){
    message("Only top meassure will be used. Radius and top do not work simultaneously.")
    radius <- NULL
  }

  modes <- c("scores", "loadings", "biplot")
  if(!mode %in% modes){
    stop_quietly(paste0("mode must be one of the following: ", paste0(modes, collapse = ", ")))
  }

  if(!is.null(factor)){
    if(!is.factor(factor) & mode %in% c("scores", "biplot")){
      stop_quietly("Factor must be a factor object.")
    }
  }else{
    factor <- factor(model$Y$data[,"event"])
  }

  if(!isa(aux.model,pkg.env$model_class)){
    stop_quietly("'model' must be a HDcox object.")
  }else if(attr(aux.model, "model") %in% pkg.env$pls_methods){
    stop_quietly("For PLS models, use the function 'plot_HDcox.PLS.model'")
  }else if(!attr(aux.model, "model") %in% pkg.env$multiblock_methods){
    stop_quietly("'model' must be a HDcox object PLS class ('SB.sPLS-ICOX','SB.sPLS-DRCOX','MB.sPLS-DRCOX' or 'MB.sPLS-DACOX').")
  }

  lst_ggp <- list()
  lst_outliers <- list()

  #4 is lst_pls, lst_spls, mb_models...
  for(block in names(aux.model$X$data)){

    lst_ggp[[block]] <- local({

      block <- block

      if(mode=="scores"){

        if(attr(aux.model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
          if(ncol(aux.model[[4]][[block]]$X$scores)==1){
            message("The model has only 1 component")

            df <- cbind(aux.model[[4]][[block]]$X$scores[,1], aux.model[[4]][[block]]$X$scores[,1])
            colnames(df) <- c("p1", "p2")
          }else{
            df <- as.data.frame(aux.model[[4]][[block]]$X$scores)
          }
        }else{ #multiblock
          if(ncol(aux.model$X$scores[[block]])==1){
            message("The model has only 1 component")

            df <- cbind(aux.model$X$scores[[block]][,1], aux.model$X$scores[[block]][,1])
            colnames(df) <- c("p1", "p2")
          }else{
            df <- as.data.frame(aux.model$X$scores[[block]])
          }
        }


        subdata_loading = NULL
        ggp <- ggplot(as.data.frame(df))

        if(nrow(df) > MAX_POINTS){
          ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor), pointsize = POINT_SIZE, pixels = POINT_RES)
        }else{
          ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor))
        }

        ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)
        ggp <- ggp + stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F)

        if("R2" %in% names(model)){
          txt.expression <- paste0("Scores (",attr(aux.model, "model"),") - ", block, " - ")
          r2_1 <- round(model$R2[[comp[1]]], 4)
          r2_2 <- round(model$R2[[comp[2]]], 4)
          r2 <- round(sum(r2_1, r2_2), 4)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
            xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
            ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
        }else{
          txt.expression <- paste0("Scores (",attr(aux.model, "model"),") - ", block)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
            xlab(label = paste0("comp_",as.character(comp[1]))) +
            ylab(label = paste0("comp_",as.character(comp[2])))
        }

        if(requireNamespace("RColorConesa", quietly = TRUE)){
          ggp <- ggp +
            RColorConesa::scale_color_conesa(reverse = colorReverse) +
            RColorConesa::scale_fill_conesa(reverse = colorReverse)
        }

      }else if(mode=="loadings"){

        if(attr(aux.model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
          if(ncol(aux.model[[4]][[block]]$X$loadings)==1){
            message("The model has only 1 component")

            df <- cbind(aux.model[[4]][[block]]$X$loadings[,1], aux.model[[4]][[block]]$X$loadings[,1])
            colnames(df) <- c("p1", "p2")
          }else{
            df <- as.data.frame(aux.model[[4]][[block]]$X$loadings)
          }
        }else{ #multiblock
          if(ncol(aux.model$X$loadings[[block]])==1){
            message("The model has only 1 component")

            df <- cbind(aux.model$X$loadings[[block]][,1], aux.model$X$loadings[[block]][,1])
            colnames(df) <- c("p1", "p2")
          }else{
            df <- as.data.frame(aux.model$X$loadings[[block]])
          }
        }

        if(class(df)[[1]] %in% "matrix"){
          df <- as.data.frame.matrix(df)
        }

        if(nrow(df)<MAX_LOADINGS){
          subdata_loading <- df
        }else if(!is.null(top)){
          aux_loadings <- apply(df,1,function(x){sqrt(crossprod(as.numeric(x[comp])))})
          aux_loadings <- aux_loadings[order(aux_loadings, decreasing = T)]
          subdata_loading <- df[names(aux_loadings)[1:top],]
        }else if(!is.null(radius)){
          subdata_loading <- df_loading[apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))>radius}),]
        }else{
          subdata_loading <- NULL
        }

        ggp <- ggplot(as.data.frame(df))

        if(nrow(df) > MAX_POINTS){
          ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]]), pointsize = POINT_SIZE, pixels = POINT_RES)
        }else{
          ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]]))
        }

        ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)

        if("R2" %in% names(model)){
          txt.expression <- paste0("Loadings (",attr(aux.model, "model"),") - ", block, " - ")
          r2_1 <- round(model$R2[[comp[1]]], 4)
          r2_2 <- round(model$R2[[comp[2]]], 4)
          r2 <- round(sum(r2_1, r2_2), 4)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
            xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
            ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
        }else{
          txt.expression <- paste0("Loadings (",attr(aux.model, "model"),") - ", block)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
            xlab(label = paste0("comp_",as.character(comp[1]))) +
            ylab(label = paste0("comp_",as.character(comp[2])))
        }

        if(names & !is.null(subdata_loading)){
          ggp <- ggp + ggrepel::geom_text_repel(data = subdata_loading, aes(x = subdata_loading[,comp[1]],
                                                                            y = subdata_loading[,comp[2]]),
                                                max.overlaps = getOption("ggrepel.max.overlaps", default = overlaps),
                                                label = rownames(subdata_loading), size=text.size)
        }

        if(!is.null(radius) & !is.null(subdata_loading)){
          if(requireNamespace("ggforce", quietly = TRUE)){
            ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
          }
        }

      }else if(mode=="biplot"){


        if(attr(aux.model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
          if(ncol(aux.model[[4]][[block]]$X$loadings)==1){
            message("The model has only 1 component")

            df <- as.data.frame(cbind(aux.model[[4]][[block]]$X$scores, aux.model[[4]][[block]]$X$scores))
            colnames(df) <- c("p1", "p2")

            df_loading <- as.data.frame(cbind(aux.model[[4]][[block]]$X$loadings[,1], aux.model[[4]][[block]]$X$loadings[,1]))
            max.loadings <- apply(abs(df_loading), 2, max)
            max.scores <- apply(abs(df), 2, max)
          }else{
            df <- as.data.frame(aux.model[[4]][[block]]$X$scores)

            df_loading <- as.data.frame(aux.model[[4]][[block]]$X$loadings)
            max.loadings <- apply(abs(df_loading), 2, max)
            max.scores <- apply(abs(df), 2, max)
          }
        }else{ #multiblock
          if(ncol(aux.model$X$loadings[[block]])==1){
            message("The model has only 1 component")

            df <- as.data.frame(cbind(aux.model$X$scores[[block]], aux.model$X$scores[[block]]))
            colnames(df) <- c("p1", "p2")

            df_loading <- as.data.frame(cbind(aux.model$X$loadings[[block]][,1], aux.model$X$loadings[[block]][,1]))
            max.loadings <- apply(abs(df_loading), 2, max)
            max.scores <- apply(abs(df), 2, max)
          }else{
            df <- as.data.frame(aux.model$X$scores[[block]])

            df_loading <- as.data.frame(aux.model$X$loadings[[block]])
            max.loadings <- apply(abs(df_loading), 2, max)
            max.scores <- apply(abs(df), 2, max)
          }
        }

        #scale scores to -1,1
        df <- norm01(df[,comp])*2-1
        ggp <- ggplot(as.data.frame(df))

        if(nrow(df) > MAX_POINTS){
          ggp <- ggp + scattermore::geom_scattermore(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor), pointsize = POINT_SIZE, pixels = POINT_RES)
        }else{
          ggp <- ggp + geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor))
        }

        ggp <- ggp + labs(color = legend_title) + theme(legend.position="bottom") + coord_fixed(ratio=1)
        ggp <- ggp + stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F)

        if("R2" %in% names(model)){
          txt.expression <- paste0("Biplot (",attr(aux.model, "model"),") - ", block, " - ")
          r2_1 <- round(model$R2[[comp[1]]], 4)
          r2_2 <- round(model$R2[[comp[2]]], 4)
          r2 <- round(sum(r2_1, r2_2), 4)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression) ~R^2 == .(r2))) +
            xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character(r2_1*100), " %)")) +
            ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character(r2_2*100), " %)"))
        }else{
          txt.expression <- paste0("Biplot (",attr(aux.model, "model"),") - ", block)
          ggp <- ggp + ggtitle(label = bquote(.(txt.expression))) +
            xlab(label = paste0("comp_",as.character(comp[1]))) +
            ylab(label = paste0("comp_",as.character(comp[2])))
        }

        if(requireNamespace("RColorConesa", quietly = TRUE)){
          ggp <- ggp +
            RColorConesa::scale_color_conesa(reverse = colorReverse) +
            RColorConesa::scale_fill_conesa(reverse = colorReverse)
        }

        if(nrow(df_loading)<MAX_LOADINGS){
          subdata_loading <- df_loading
        }else if(!is.null(top)){
          aux_loadings <- apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))})
          aux_loadings <- aux_loadings[order(aux_loadings, decreasing = T)]
          subdata_loading <- df_loading[names(aux_loadings)[1:top],]
        }else if(!is.null(radius)){
          subdata_loading <- df_loading[apply(df_loading,1,function(x){sqrt(crossprod(as.numeric(x[comp])))>radius}),]
        }else{
          subdata_loading <- NULL
        }

        #depending on DF instead of df_loadings - ARROWS
        if(any(!is.null(top), !is.null(radius))){

          no_selected_loadings <- df_loading[!rownames(df_loading) %in% rownames(subdata_loading),]
          if(nrow(no_selected_loadings)!=0 & !only_top){
            ggp <- ggp + geom_segment(data = no_selected_loadings, lineend = "butt", linejoin = "mitre", size = 0.2,
                                      aes(x = 0, y = 0, xend = no_selected_loadings[,comp[1]],
                                          yend = no_selected_loadings[,comp[2]]),
                                      arrow = arrow(length = unit(0.1, "cm")))
          }

          ggp <- ggp + geom_segment(data = subdata_loading, lineend = "butt", linejoin = "mitre",
                                    size = 0.33, aes(x = 0, y = 0, xend = subdata_loading[,comp[1]],
                                                     yend = subdata_loading[,comp[2]]),
                                    arrow = arrow(length = unit(0.1, "cm")))

        }else{
          #show all loadings
          no_selected_loadings <- df_loading[!rownames(df_loading) %in% rownames(subdata_loading),]
          ggp <- ggp + geom_segment(data = no_selected_loadings, lineend = "butt", linejoin = "mitre", size = 0.2,
                                    aes(x = 0, y = 0, xend = no_selected_loadings[,comp[1]],
                                        yend = no_selected_loadings[,comp[2]]),
                                    arrow = arrow(length = unit(0.1, "cm")))
        }

        if(names & !is.null(subdata_loading)){
          ggp <- ggp + ggrepel::geom_text_repel(data = subdata_loading, aes(x = subdata_loading[,comp[1]],
                                                                            y = subdata_loading[,comp[2]]),
                                                max.overlaps = getOption("ggrepel.max.overlaps", default = overlaps),
                                                label = rownames(subdata_loading), size=text.size)
        }

        if(is.null(top) & !is.null(radius) & nrow(df) < MAX_POINTS){
          ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
        }

      }

      #reorder legend
      if(!is.null(factor) & length(levels(factor))>3){
        ggp <- ggp + guides(color=guide_legend(nrow = ceiling(length(levels(factor))/3), byrow = T))
      }


      ggp})

  }

  return(list(plot_block = lst_ggp))
}

#### ### ### ### ### ##
# PROPORTIONAL HAZARD #
#### ### ### ### ### ##

#' plot_proportionalHazard.list
#' @description Run the function "plot_proportionalHazard" for a list of models. More information in "?plot_proportionalHazard".
#'
#' @param lst_models List of HDcox models.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_proportionalHazard.list(lst_models)
#' }

plot_proportionalHazard.list <- function(lst_models){

  lst_plots <- purrr::map(lst_models, ~plot_proportionalHazard(model = .))

  return(lst_plots)
}

#' plot_proportionalHazard
#' @description Perform the proportional hazard assumption for a specific HDcox model.
#' It uses the functions "survival::cox.zph" and "survminer::ggcoxzph", but it returns a ggplot2 graph.
#'
#' @param model HDcox model.
#' @return A ggplot2 plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_proportionalHazard(model)
#' }

plot_proportionalHazard <- function(model){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(all(is.null(model$survival_model$fit)) || all(is.na(model$survival_model$fit))){
    message(paste0("Survival model not found for ", attr(model, "model")))
    return(NULL)
  }

  ph_preplot <- survival::cox.zph(model$survival_model$fit)
  ph_plot <- survminer::ggcoxzph(ph_preplot)
  ph_ggplot <- ggcoxzph2ggplot(ph_preplot, ph_plot)
  return(ph_ggplot)
}

ggcoxzph2ggplot <- function(pre.ggcoxzph, ggcoxzph){
  lst_plots <- list()
  for(p in names(ggcoxzph)){
    lst_plots[[p]] <- ggcoxzph[[p]]
  }

  if(length(lst_plots)==1){
    return(lst_plots[[1]])
  }

  global_test.txt <- paste0("Global Schoenfeld Test: ", round(pre.ggcoxzph$table["GLOBAL","p"], digits = 4))

  len <- length(lst_plots)
  p.vector <- my_primeFactors(ifelse(len %% 2 == 1, len+1, len))
  if(length(p.vector)>2){

    while(length(p.vector)>2){
      if(p.vector[1] < p.vector[length(p.vector)]){
        p.vector <- c(p.vector[1] * p.vector[2], p.vector[3:length(p.vector)])
      }else{
        p.vector <- c(p.vector[1:(length(p.vector)-2)], p.vector[length(p.vector)-1] * p.vector[length(p.vector)])
      }
    }

    ncol <- min(p.vector)
    nrow <- max(p.vector)
  }else{
    ncol <- min(p.vector)
    nrow <- max(p.vector)
  }

  ggp <- ggpubr::ggarrange(plotlist = lst_plots, nrow = nrow, ncol = ncol)
  ggp_final <- ggpubr::annotate_figure(ggp, top = global_test.txt)

  return(ggp_final)
}

my_primeFactors <- function(num) {
  current <- num
  ret.vals <- vector()
  x <- 2
  while (x <= num - 1){
    while (current %% x == 0) {
      current <- current / x
      ret.vals <- c(ret.vals, x)
    }
    x <- x + 1
  }
  if (is.logical(ret.vals)) return(num) else return(ret.vals)
}

#### ### ### ##
# FOREST PLOT #
#### ### ### ##

#' plot_forest.list
#' @description Run the function "plot_forest" for a list of models. More information in "?plot_forest".
#'
#' @param lst_models List of HDcox models.
#' @param title Character. Forest plot title (default: "Hazard Ratio").
#' @param cpositions Numeric vector. Relative positions of first three columns in the OX scale (default: c(0.02, 0.22, 0.4)).
#' @param fontsize Numeric. Elative size of annotations in the plot (default: 0.7).
#' @param refLabel Character. Label for reference levels of factor variables (default: "reference").
#' @param noDigits Numeric. Number of digits for estimates and p-values in the plot (default: 2).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_forest.list(lst_models)
#' }

plot_forest.list <- function(lst_models,
                             title = "Hazard Ratio",
                             cpositions = c(0.02, 0.22, 0.4),
                             fontsize = 0.7,
                             refLabel = "reference",
                             noDigits = 2){

  lst_forest_plot <- purrr::map(lst_models, ~plot_forest(model = .,
                                                         title = paste0(title, " - ", .$class), cpositions = cpositions,
                                                         fontsize = fontsize, refLabel = refLabel, noDigits = noDigits))

  return(lst_forest_plot)

}

#' plot_forest
#' @description Forest plot for HDcox models using the R library plot_forest. For more information check ?survminer::ggforest.
#' @param model HDcox model.
#' @param title Character. Forest plot title (default: "Hazard Ratio").
#' @param cpositions Numeric vector. Relative positions of first three columns in the OX scale (default: c(0.02, 0.22, 0.4)).
#' @param fontsize Numeric. Elative size of annotations in the plot (default: 0.7).
#' @param refLabel Character. Label for reference levels of factor variables (default: "reference").
#' @param noDigits Numeric. Number of digits for estimates and p-values in the plot (default: 2).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_forest(model)
#' }

plot_forest <- function(model,
                        title = "Hazard Ratio",
                        cpositions = c(0.02, 0.22, 0.4),
                        fontsize = 0.7,
                        refLabel = "reference",
                        noDigits = 2){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(!attr(model, "model") %in% pkg.env$all_methods){
    stop(paste0("Model must be one of the following HDcox models: ", paste0(pkg.env$all_methods, collapse = ", ")))
  }

  if(all(is.null(model$survival_model$fit)) || all(is.na(model$survival_model$fit)) || all(is.null(model)) || all(is.na(model))){
    message(paste0("Survival model not found for ", attr(model, "model")))
    return(NULL)
  }

  ggp <- survminer::ggforest(model = model$survival_model$fit,
                             data = model$survival_model$fit$model,
                             main = title, cpositions = cpositions, fontsize = fontsize, refLabel = refLabel, noDigits = noDigits)
  return(ggp)
}

#### ### ### ### ### ### ### #
# EVENT DISTRIBUTION - MODEL #
#### ### ### ### ### ### ### #

#' plot_cox.event.list
#' @description Run the function "plot_cox.event" for a list of models. More information in "?plot_cox.event".
#'
#' @param lst_models List of HDcox models.
#' @param type Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").
#' @param n.breaks Numeric. Number of time-break points to compute (default: 20).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cox.event.list(lst_models)
#' }

plot_cox.event.list <- function(lst_models, type = "lp", n.breaks = 20){

  ggp_list <- purrr::map(lst_models, ~plot_cox.event(model = ., type = type, n.breaks = n.breaks))

  return(ggp_list)
}

#' plot_cox.event
#'
#' @param model HDcox model.
#' @param type Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").
#' @param n.breaks Numeric. If BREAKTIME is NULL, "n.breaks" is the number of time-break points to compute (default: 20).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cox.event(model)
#' }

plot_cox.event <- function(model, type = "lp", n.breaks = 20){

  #DFCALLS
  event <- NULL

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  #exits
  if(all(is.null(model$survival_model$fit)) || all(is.na(model$survival_model$fit))){
    message(paste0("Survival model not found for ", attr(model, "model"), "."))
    return(NULL)
  }

  if(type=="survival"){
    lp <- exp(-predict(model$survival_model$fit, type = "expected"))
  }else if(type %in% c("lp", "risk", "expected")){
    lp <- predict(model$survival_model$fit, type = type)
  }else{
    stop_quietly("Type must be one of the follow: 'lp', 'risk', 'expected', 'survival'")
  }
  names(lp) <- rownames(model$survival_model$fit$y)

  df_hr <- cbind(lp, model$Y_input[names(lp),"event"])
  colnames(df_hr) <- c(type, "event")
  df_hr <- as.data.frame(df_hr)
  df_hr$event <- factor(df_hr$event, levels = c(0,1))

  ggp.d <- ggplot(df_hr, aes(x=lp, fill=event)) +
    geom_density(alpha=0.5) +
    ggtitle(attr(model, "model")) +
    xlab(label = type)

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp.d <- ggp.d + RColorConesa::scale_fill_conesa()
  }

  binwidth <- (max(df_hr[,1]) - min(df_hr[,1])) / n.breaks
  breaks <- seq(min(df_hr[,1]), max(df_hr[,1]), binwidth)

  ggp.h <- ggplot(df_hr, aes(x=lp, fill=event, color=event)) +
    geom_histogram(position = "stack", alpha=0.75, breaks = breaks) +
    ggtitle(attr(model, "model")) +
    xlab(label = type)

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp.h <- ggp.h + RColorConesa::scale_fill_conesa() + RColorConesa::scale_color_conesa()
  }

  return(list(df = df_hr, plot.density = ggp.d, plot.histogram = ggp.h))
}

prop.between2values <- function(df, min, max){
  aux.df <- df[min<df$lp & df$lp<=max,]
  count <- table(aux.df$event)
  perc <- round(prop.table(table(aux.df$event))*100,2)

  total_0 <- round(count[[1]] / sum(df$event==levels(df$event)[[1]]) * 100,2)
  total_1 <- round(count[[2]] / sum(df$event==levels(df$event)[[2]]) * 100,2)
  cat(paste0("Between ", min, " and ", max, " there are:\n",
             perc[[1]], " % of censored (",total_0, " % of total censored)\n",
             perc[[2]], " % of events (",total_1, " % of total event)\n\n"))
}

#### ### ### ### ### ### ### ### ##
# EVENT DISTRIBUTION - PREDICTION #
#### ### ### ### ### ### ### ### ##

#' patient.eventDensity
#'
#' @param patient Numeric matrix or data.frame. New explanatory variables (raw data) for one observation. Qualitative variables must be transform into binary variables.
#' @param time Numeric. Time point where the AUC will be evaluated (default: NULL).
#' @param model HDcox model.
#' @param type Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").
#' @param size Numeric. Point size (default: 3).
#' @param color String. R Color.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_patient.eventDensity(patient, time = NULL, model)
#' }

plot_patient.eventDensity <- function(patient, time = NULL, model, type = "lp", size = 3, color = "red"){

  #DFCALLS
  x <- y <- event <- NULL

  pred.value <- cox.prediction(model = model, new_data = patient, time = time, type = type, method = "cox")

  plot <- plot_cox.event(model, type = type)
  plot <- plot$plot.density

  #get density
  density_event <- density(plot$data[plot$data$event==1,1])
  index <- which.min(abs(density_event$x - pred.value))
  y.value_event <- density_event$y[index]

  density_noevent <- density(plot$data[plot$data$event==0,1])
  index <- which.min(abs(density_noevent$x - pred.value))
  y.value_noevent <- density_noevent$y[index]
  y.value <- max(y.value_event, y.value_noevent)

  max <- max(density_event$y) / 10

  df <- data.frame(x = c(pred.value, pred.value), y = c(y.value_noevent, y.value_event), event = factor(c(0,1)))

  plot.new <- plot +
    geom_point(data = df, aes(x = x, y = y, fill = event, color = event), size = size)

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    plot.new <- plot.new + RColorConesa::scale_color_conesa()
  }

  df <- data.frame(x = pred.value, y = y.value + max)

  plot.new <- plot +
    geom_point(data = df, aes(x = x, y = y), inherit.aes = F, color = color, size = size) +
    geom_segment(data = df, aes(x = x, y = 0, xend = x, yend = y), inherit.aes = F, color = color, size = 0.8)

  return(plot.new)
}

#' patient.eventHistogram
#'
#' @param patient Numeric matrix or data.frame. New explanatory variables (raw data) for one observation. Qualitative variables must be transform into binary variables.
#' @param time Numeric. Time point where the AUC will be evaluated (default: NULL).
#' @param model HDcox model.
#' @param type Character. Prediction type: "lp", "risk", "expected" or "survival" (default: "lp").
#' @param size Numeric. Point size (default: 3).
#' @param color String. R Color.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_patient.eventHistogram(patient, time = NULL, model)
#' }

plot_patient.eventHistogram <- function(patient, time = NULL, model, type = "lp", size = 3, color = "red"){

  #DFCALLS
  x <- y <- NULL

  pred.value <- cox.prediction(model = model, new_data = patient, time = time, type = type, method = "cox")

  plot <- plot_cox.event(model, type = type)
  plot <- plot$plot.histogram

  #get histogram
  intervals <- plot$layers[[1]]$stat_params$breaks

  index <- which.min(abs(intervals - pred.value))

  if(pred.value > intervals[index]){
    index <- c(index, index+1)
  }else{
    index <- c(index-1, index)
  }

  #max <- max(density_event$y) / 10

  y.value <- nrow(plot$data[plot$data[,1] >= intervals[index[1]] & plot$data[,1] <= intervals[index[2]],])

  #df <- data.frame(x = pred.value, y = y.value + max)
  df <- data.frame(x = (intervals[index[1]] + intervals[index[2]]) / 2, y = y.value)

  plot.new <- plot +
    geom_point(data = df, aes(x = x, y = y), inherit.aes = F, color = color, size = size) +
    geom_segment(data = df, aes(x = x, y = 0, xend = x, yend = y), inherit.aes = F, color = color, size = 0.8)

  return(plot.new)
}

#### ### ### ### ### ### ### ### ###
# PSEUDOBETA PLOTS - PLSCOX MODELS #
#### ### ### ### ### ### ### ### ###

#' plot_pseudobeta.list
#' @description Run the function "plot_pseudobeta" for a list of models. More information in "?plot_pseudobeta".
#'
#' @param lst_models List of HDcox models.
#' @param error.bar Logical. Show error bar (default: TRUE).
#' @param onlySig Logical. Compute psudobetas using only significant components (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param show_percentage Logical. If show_percentage = TRUE, it shows the contribution percentage for each variable to the full model (default: TRUE).
#' @param size_percentage Numeric. Size of percentage text (default: 3).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pseudobeta.list(lst_models)
#' }

plot_pseudobeta.list <- function(lst_models, error.bar = T, onlySig = F, alpha = 0.05, zero.rm = T, top = NULL, auto.limits = T, show_percentage = T, size_percentage = 3, verbose = F){


  if(all(unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods))){
    sub_lst_models <- lst_models
  }else{
    sub_lst_models <- lst_models[unlist(purrr::map(lst_models, function(x){x$class})) %in% pkg.env$pls_methods]
    if(verbose){
      message(paste0("Model ", paste0(names(lst_models[!unlist(purrr::map(lst_models, function(x){x$class})) %in% pkg.env$pls_methods]), collapse = ", "), " are not based in PLS methodology. Other models computed."))
    }
  }

  lst_plots <- purrr::map(sub_lst_models, ~plot_pseudobeta(model = .,
                                                           error.bar = error.bar,
                                                           onlySig = onlySig, alpha = alpha,
                                                           zero.rm = zero.rm, auto.limits = auto.limits, top = top,
                                                           show_percentage = show_percentage, size_percentage = size_percentage))

  return(lst_plots)
}

#' plot_pseudobeta
#' @description Decompose a PLS-cox model into a original variable pseudo-beta interpretation. As each component has associated a cox coeffitient, and each component has a weight related to the original variable. Cox final formula is expressed in terms of original variables.
#'
#' @param model HDcox model.
#' @param error.bar Logical. Show error bar (default: TRUE).
#' @param onlySig Logical. Compute psudobetas using only significant components (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param show_percentage Logical. If show_percentage = TRUE, it shows the contribution percentage for each variable to the full model (default: TRUE).
#' @param size_percentage Numeric. Size of percentage text (default: 3).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pseudobeta(model)
#' }

plot_pseudobeta <- function(model, error.bar = T, onlySig = F, alpha = 0.05, zero.rm = T, top = NULL, auto.limits = T, show_percentage = T, size_percentage = 3){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(!attr(model, "model") %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)){
    stop("Model must be one of the follow models: 'sPLS-ICOX', 'sPLS-DRCOX', 'sPLS-DRCOX-Dynamic', 'sPLS-DACOX-Dynamic', 'SB.sPLS-ICOX', 'SB.sPLS-DRCOX', 'MB.sPLS-DRCOX', 'MB.sPLS-DACOX'")
  }

  if(all(is.null(model$survival_model))){
    stop("Survival Model not found.")
  }

  df.aux <- as.data.frame(summary(model$survival_model$fit)[[7]])

  if(attr(model, "model") %in% pkg.env$pls_methods){

    if(onlySig){
      rn <- rownames(df.aux)[df.aux$`Pr(>|z|)` <= alpha]
      coefficients <- as.matrix(model$survival_model$fit$coefficients)[rn,,drop=F]
      sd <- df.aux[rn,"se(coef)",drop=F]
      W.star <- model$X$W.star[,rn,drop=F]
    }else{
      coefficients <- as.matrix(model$survival_model$fit$coefficients)
      sd <- df.aux[,"se(coef)",drop=F]
      W.star <- model$X$W.star
    }

    vector <- W.star %*% coefficients

    if(error.bar){
      sd.min <- W.star %*% data.matrix(coefficients-sd)
      sd.max <- W.star %*% data.matrix(coefficients+sd)
    }else{
      sd.min <- NULL
      sd.max <- NULL
    }

    #sort
    vector <- vector[order(vector[,1], decreasing = T),,drop=F]

    if(error.bar){
      sd.min <- sd.min[rownames(vector),,drop=F]
      sd.max <- sd.max[rownames(vector),,drop=F]
    }

    if(all(vector[,1]==0)){
      return(list(beta = vector,
                  plot = NULL,
                  sd.min = sd.min,
                  sd.max = sd.max))
    }

    plot <- coxweightplot.fromVector.HDcox(model = model, vector = vector,
                                           sd.min = sd.min, sd.max = sd.max, auto.limits = auto.limits,
                                           zero.rm = zero.rm, top = top,
                                           show_percentage = show_percentage,
                                           size_percentage = size_percentage)

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){

    if(onlySig){
      rn <- rownames(df.aux)[df.aux$`Pr(>|z|)` <= alpha]
      coefficients <- as.matrix(model$survival_model$fit$coefficients)[rn,,drop=F]
      sd <- df.aux[rn,"se(coef)",drop=F]
      W.star <- list()
      if(attr(model, "model") %in% pkg.env$sb.splsicox){
        for(b in names(model$list_spls_models)){
          W.star[[b]] <- model$list_spls_models[[b]]$X$W.star
        }
      }else if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
        for(b in names(model$list_spls_models)){
          W.star[[b]] <- model$list_spls_models[[b]]$X$W.star
        }
      }else{
        W.star <- model$X$W.star
      }

    }else{
      coefficients <- as.matrix(model$survival_model$fit$coefficients)
      sd <- df.aux[,"se(coef)",drop=F]
      W.star <- list()
      if(attr(model, "model") %in% pkg.env$sb.splsicox){
        for(b in names(model$list_spls_models)){
          W.star[[b]] <- model$list_spls_models[[b]]$X$W.star
        }
      }else if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
        for(b in names(model$list_spls_models)){
          W.star[[b]] <- model$list_spls_models[[b]]$X$W.star
        }
      }else{
        ### IF MODEL COMES FROM MIXOMICS - We should use W* WITHOUT NORMALIZE, normalization is only for predicting new X scores - checked
        W.star <- model$X$W.star
      }
    }

    vector <- list()
    sd.min <- list()
    sd.max <- list()
    plot <- list()

    for(b in names(model$X$data)){
      coeff <- coefficients[grep(b,rownames(coefficients)),,drop=F]
      if(length(coeff)==0){
        next
      }

      components <- unlist(lapply(rownames(coeff), function(x){paste0(strsplit(x, "_")[[1]][1:2], collapse = "_")}))
      vector[[b]] <- W.star[[b]][,components,drop=F] %*% coeff

      if(error.bar){
        sd.min[[b]] <- W.star[[b]][,components,drop=F] %*% data.matrix(coeff-sd[rownames(coeff),,drop=F])
        sd.max[[b]] <- W.star[[b]][,components,drop=F] %*% data.matrix(coeff+sd[rownames(coeff),,drop=F])
      }else{
        sd.min[[b]] <- NULL
        sd.max[[b]] <- NULL
      }

      #sort
      vector[[b]] <- vector[[b]][order(vector[[b]][,1], decreasing = T),,drop=F]

      if(error.bar){
        sd.min[[b]] <- sd.min[[b]][rownames(vector[[b]]),,drop=F]
        sd.max[[b]] <- sd.max[[b]][rownames(vector[[b]]),,drop=F]
      }

      if(all(vector[[b]][,1]==0)){
        plot[[b]] = NULL
      }else{
        plot[[b]] <- coxweightplot.fromVector.HDcox(model = model, vector = vector[[b]],
                                                    sd.min = sd.min[[b]], sd.max = sd.max[[b]], auto.limits = auto.limits,
                                                    zero.rm = zero.rm, top = top, block = b,
                                                    show_percentage = show_percentage,
                                                    size_percentage = size_percentage)
      }

    }

  }

  if(attr(model, "model") %in% pkg.env$pls_methods){
    vector <- plot$coefficients
    return(list(plot = plot$plot,
                beta = vector,
                sd.min = sd.min,
                sd.max = sd.max))
  }else{

    aux_vector <- list()
    aux_plot <- list()
    for(b in names(model$X$data)){
      aux_vector[[b]] <- plot[[b]]$coefficients
      aux_plot[[b]] <- plot[[b]]$plot
    }

    return(list(plot = aux_plot,
                beta = aux_vector,
                sd.min = sd.min,
                sd.max = sd.max))
  }

}

#### ### ### ### ### ### ### ###
# PSEUDOBETA PLOTS - PREDICTION #
#### ### ### ### ### ### ### ###

#' plot_pseudobeta_newPatient.list
#' @description Run the function "plot_pseudobeta_newPatient" for a list of models. More information in "?plot_pseudobeta_newPatient".
#'
#' @param lst_models List of HDcox models.
#' @param new_observation Numeric matrix or data.frame. New explanatory variables (raw data) for one observation. Qualitative variables must be transform into binary variables.
#' @param error.bar Logical. Show error bar (default: TRUE).
#' @param onlySig Logical. Compute psudobetas using only significant components (default: TRUE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param show.betas Logical. Show original betas (default: FALSE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pseudobeta_newPatient.list(lst_models, new_observation)
#' }

plot_pseudobeta_newPatient.list <- function(lst_models, new_observation, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                            top = NULL, auto.limits = T, show.betas = F, verbose = F){

  if(all(unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods))){
    sub_lst_models <- lst_models
  }else{
    sub_lst_models <- lst_models[unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)]
    if(verbose){
      message(paste0("Model ", paste0(names(lst_models[!unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)]), collapse = ", "), " are not based in PLS methodology. Other models computed."))
    }
  }

  lst_plots <- purrr::map(sub_lst_models, ~plot_pseudobeta_newPatient(model = .,
                                                                      new_observation = new_observation,
                                                                      error.bar = error.bar,
                                                                      onlySig = onlySig, alpha = alpha,
                                                                      zero.rm = zero.rm, top = top,
                                                                      auto.limits = auto.limits, show.betas = show.betas))

  return(lst_plots)

}

#' plot_pseudobeta.newPatient
#' @description Generates a ggplot2 to compare the pseudobeta direction from de HDcox model and the values of the new observation.
#'
#' @param model HDcox model.
#' @param new_observation Numeric matrix or data.frame. New explanatory variables (raw data) for one observation. Qualitative variables must be transform into binary variables.
#' @param error.bar Logical. Show error bar (default: TRUE).
#' @param onlySig Logical. Compute psudobetas using only significant components (default: TRUE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).
#' @param top Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If top = NULL, all variables are shown (default: NULL).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param show.betas Logical. Show original betas (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pseudobeta_newPatient(model, new_observation)
#' }

plot_pseudobeta_newPatient <- function(model, new_observation, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                       top = NULL, auto.limits = T, show.betas = F){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(attr(model, "model") %in% pkg.env$pls_methods){
    plot_pseudobeta.newPatient(model = model,
                               new_observation = new_observation,
                               error.bar = error.bar,
                               onlySig = onlySig, alpha = alpha,
                               zero.rm = zero.rm, top = top,
                               auto.limits = auto.limits, show.betas = show.betas)
  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){
    plot_MB.pseudobeta.newPatient(model = model,
                                  new_observation = new_observation,
                                  error.bar = error.bar,
                                  onlySig = onlySig, alpha = alpha,
                                  zero.rm = zero.rm, top = top,
                                  auto.limits = auto.limits, show.betas = show.betas)
  }else{
    stop("Model not belong to any PLS or MB HDcox methods.")
  }
}

plot_pseudobeta.newPatient <- function(model, new_observation, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                       top = NULL, auto.limits = T, show.betas = F){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  #DFCALLS
  lp <- lp.min <- lp.max <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                        alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)
  coefficients <- ggp.simulated_beta$beta

  if(all(coefficients==0)){
    message("No significant variables selected.")
    return(NULL)
  }

  coeff.min <- NULL
  coeff.max <- NULL
  if(error.bar){
    coeff.min <- ggp.simulated_beta$sd.min
    coeff.max <- ggp.simulated_beta$sd.max
  }

  #norm patient
  new_observation <- new_observation[,names(model$X$x.mean),drop=F]

  if(!is.null(model$X$x.mean) & !is.null(model$X$x.sd)){
    norm_patient <- scale(new_observation, center = model$X$x.mean, scale = model$X$x.sd)
  }else if(!is.null(model$X$x.mean)){
    norm_patient <- scale(new_observation, center = model$X$x.mean, scale = F)
  }else if(!is.null(model$X$x.sd)){
    norm_patient <- scale(new_observation, center = F, scale = model$X$x.sd)
  }else{
    norm_patient <- new_observation
  }

  #lp.new_observation_manual <- norm_patient[,rownames(coefficients)] %*% coefficients #predict lp
  lp.new_observation_variable <- as.data.frame(norm_patient[,rownames(coefficients)] * coefficients$value) #predict terms
  colnames(lp.new_observation_variable) <- "value"

  lp.new_observation_variable.min <- NULL
  lp.new_observation_variable.max <- NULL
  if(error.bar){
    lp.new_observation_variable.min <- norm_patient[,rownames(coeff.min)] * coeff.min
    lp.new_observation_variable.max <- norm_patient[,rownames(coeff.max)] * coeff.max
  }

  #filter pat_variables using psudobeta plot (top could be applied)
  lp.new_observation_variable <- lp.new_observation_variable[rownames(ggp.simulated_beta$plot$data),,drop=F]
  lp.new_observation_variable.min <- lp.new_observation_variable.min[rownames(ggp.simulated_beta$plot$data),,drop=F]
  lp.new_observation_variable.max <- lp.new_observation_variable.max[rownames(ggp.simulated_beta$plot$data),,drop=F]

  coefficients <- coefficients[rownames(lp.new_observation_variable),,drop=F]

  #terms
  # df <- as.data.frame(cbind(cbind(ggp.simulated_beta$beta,
  #                                 rep("Beta",nrow(ggp.simulated_beta$beta))),
  #                           rownames(ggp.simulated_beta$beta)))
  # colnames(df) <- c("beta", "type", "var")
  #
  # df$beta <- as.numeric(df$beta)
  # df <- df[order(df$beta, decreasing = T),]
  #
  # df.pat <- cbind(cbind(lp.new_observation_variable,  rep("Patient Linear Predictor", nrow(lp.new_observation_variable))), rownames(lp.new_observation_variable))
  # colnames(df.pat) <- c("beta", "type", "var")
  # df <- rbind(df, df.pat)
  #
  # df$beta <- as.numeric(df$beta)
  # df$var <- factor(df$var, levels = unique(df$var))
  # df$type <- factor(df$type, levels = unique(df$type))

  #terms
  if(error.bar){
    df.pat <- data.frame("lp" = lp.new_observation_variable[,1],
                         "lp.min" = lp.new_observation_variable.min[,1],
                         "lp.max" = lp.new_observation_variable.max[,1],
                         "var" = rownames(lp.new_observation_variable))
  }else{
    df.pat <- data.frame("lp" = lp.new_observation_variable[,1],
                         "lp.min" = 0,
                         "lp.max" = 0,
                         "var" = rownames(lp.new_observation_variable))
  }

  df.pat$lp <- as.numeric(df.pat$lp)
  df.pat$lp.min <- as.numeric(df.pat$lp.min)
  df.pat$lp.max <- as.numeric(df.pat$lp.max)
  df.pat$var <- factor(df.pat$var, levels = unique(df.pat$var))

  accuracy <- 0.1

  #limit based on max value in abs between lower and higher values
  if(show.betas){
    if(error.bar){
      val_min <- as.numeric(min(min(coeff.max), min(df.pat$lp.min)))
      val_max <- as.numeric(max(max(coeff.max), max(df.pat$lp.max)))
      auto.limits_min <- round2any(val_min, accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(val_max, accuracy = accuracy, f = ceiling)
      auto.limits <- max(auto.limits_min, auto.limits_max)
    }else{
      auto.limits <- round2any(max(abs(coefficients), abs(df.pat$lp)), accuracy = accuracy, f = ceiling)
    }
  }else{ #not show.betas
    if(error.bar){
      auto.limits_min <- round2any(max(abs(df.pat$lp.min)), accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(max(abs(df.pat$lp.max)), accuracy = accuracy, f = ceiling)
      auto.limits <- max(auto.limits_min, auto.limits_max)
    }else{
      auto.limits <- round2any(max(abs(df.pat$lp)), accuracy = accuracy, f = ceiling)
    }
  }

  ggp <- ggplot(df.pat, aes(x = var, y = lp, fill = lp, color = 1)) +
    geom_bar(stat = "identity", position = "dodge")

  if(error.bar){
    ggp <- ggp + geom_errorbar(aes(ymin=lp.min, ymax=lp.max), width=.35, position=position_dodge(.2))
  }

  if(!show.betas){
    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                        mid = "white", midpoint = 0,
                                        high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                        limits = c(-1*auto.limits,auto.limits), name = "Beta value")
    }else{
      ggp <- ggp + scale_fill_gradient2(low = "blue",
                                        mid = "white", midpoint = 0,
                                        high = "red",
                                        limits = c(-1*auto.limits,auto.limits), name = "Beta value")
    }
  }

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "warm", continuous = T)
  }

  ggp <- ggp + guides(color= "none")
  ggp <- ggp + ylab(label = "Linear Predictor")
  ggp <- ggp + xlab(label = "Variables")
  ggp <- ggp + ggtitle(label = paste0("Observation - ", rownames(new_observation)))

  if(length(unique(df.pat$var))>15){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  if(show.betas){

    auto.limits.lp <- max(abs(min(df.pat$lp.max)), abs(max(df.pat$lp.max)))
    ggp.aux <- ggp + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits.lp, auto.limits.lp))

    ggp.aux2 <- ggp.simulated_beta$plot
    ggp.aux2 <- ggp.aux2 + guides(fill = "none")
    suppressMessages(
      ggp.aux2 <- ggp.aux2 + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))
    )

    sign.beta <- coefficients$value>0
    names(sign.beta)<-rownames(coefficients)
    sign.pat <- df.pat$lp>0
    same.sign <- sign.beta == sign.pat
    same.sign <- same.sign[rownames(ggp.simulated_beta$plot$data)]

    ggp.aux$mapping$fill[[2]] <- same.sign
    ggp.aux <- ggp.aux + guides(fill = guide_legend(title="Same beta direction:")) + theme(legend.position="left")

    #overwriting fill generates a message
    suppressMessages({
      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp.aux <- ggp.aux + RColorConesa::scale_fill_conesa(reverse = T)
      }else{
        ggp.aux <- ggp.aux + scale_fill_discrete()
      }
    })

    ggp <- ggpubr::ggarrange(ggp.aux, ggp.aux2, ncol = 2, widths = c(0.5, 0.5), align = "h")
  }

  return(list(plot = ggp, lp.var = lp.new_observation_variable, norm_pat = norm_patient, pat = new_observation))

}

plot_MB.pseudobeta.newPatient <- function(model, new_observation, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                          top = NULL, auto.limits = T, show.betas = F){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  #checks
  if(!all(names(new_observation) == names(model$X$data))){
    stop("New patint has to have the same blocks as the model.")
  }

  #DFCALLS
  lp <- lp.min <- lp.max <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                        alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)

  coefficients <- ggp.simulated_beta$beta #list

  coeff.min <- NULL
  coeff.max <- NULL

  if(error.bar){
    coeff.min <- ggp.simulated_beta$sd.min
    coeff.max <- ggp.simulated_beta$sd.max
  }

  #norm patient
  norm_patient <- list()
  lp.new_observation_variable <- list()

  lst_plots <- list()
  lst_lp.var <- list()

  #for each block... that is returned in gg.suimulated_beta...
  for(b in names(model$X$data)[names(model$X$data) %in% names(ggp.simulated_beta$plot)]){

    new_observation[[b]] <- new_observation[[b]][,names(model$X$x.mean[[b]]),drop=F]

    if(!is.null(model$X$x.mean[[b]]) & !is.null(model$X$x.sd[[b]])){
      norm_patient[[b]] <- scale(new_observation[[b]], center = model$X$x.mean[[b]], scale = model$X$x.sd[[b]])
    }else if(!is.null(model$X$x.mean[[b]])){
      norm_patient[[b]] <- scale(new_observation[[b]], center = model$X$x.mean[[b]], scale = F)
    }else if(!is.null(model$X$x.sd[[b]])){
      norm_patient[[b]] <- scale(new_observation[[b]], center = F, scale = model$X$x.sd[[b]])
    }else{
      norm_patient <- new_observation
    }

    lp.new_observation_variable[[b]] <- as.data.frame(norm_patient[[b]][,rownames(coefficients[[b]])] * coefficients[[b]]$value) #predict terms
    colnames(lp.new_observation_variable[[b]]) <- "value"

    lp.new_observation_variable.min <- NULL
    lp.new_observation_variable.max <- NULL

    if(error.bar){
      if(b %in% names(coeff.min)){
        lp.new_observation_variable.min <- norm_patient[[b]][,rownames(coeff.min[[b]])] * coeff.min[[b]]
        lp.new_observation_variable.max <- norm_patient[[b]][,rownames(coeff.max[[b]])] * coeff.max[[b]]
      }
    }

    #filter pat_variables using psudobeta plot (top could be applied)
    lp.new_observation_variable[[b]] <- lp.new_observation_variable[[b]][rownames(ggp.simulated_beta$plot[[b]]$data),,drop=F]
    lp.new_observation_variable.min <- lp.new_observation_variable.min[rownames(ggp.simulated_beta$plot[[b]]$data),,drop=F]
    lp.new_observation_variable.max <- lp.new_observation_variable.max[rownames(ggp.simulated_beta$plot[[b]]$data),,drop=F]

    coefficients[[b]] <- coefficients[[b]][rownames(lp.new_observation_variable[[b]]),,drop=F]

    if(all(coefficients[[b]]==0)){
      message("No significant variables selected.")
      next
    }

    #terms
    if(error.bar){
      df.pat <- data.frame("lp" = lp.new_observation_variable[[b]][,1],
                           "lp.min" = lp.new_observation_variable.min[,1],
                           "lp.max" = lp.new_observation_variable.max[,1],
                           "var" = rownames(lp.new_observation_variable[[b]]))
    }else{
      df.pat <- data.frame("lp" = lp.new_observation_variable[[b]][,1],
                           "lp.min" = 0,
                           "lp.max" = 0,
                           "var" = rownames(lp.new_observation_variable[[b]]))
    }

    df.pat$lp <- as.numeric(df.pat$lp)
    df.pat$lp.min <- as.numeric(df.pat$lp.min)
    df.pat$lp.max <- as.numeric(df.pat$lp.max)
    df.pat$var <- factor(df.pat$var, levels = unique(df.pat$var))

    accuracy <- 0.1

    if(show.betas){
      if(error.bar){
        val_min <- as.numeric(max(abs(coeff.min[[b]]), abs(df.pat$lp.min)))
        val_max <- as.numeric(max(abs(coeff.max[[b]]), abs(df.pat$lp.max)))
        auto.limits_min <- round2any(val_min, accuracy = accuracy, f = ceiling)
        auto.limits_max <- round2any(val_max, accuracy = accuracy, f = ceiling)
        auto.limits <- max(auto.limits_min, auto.limits_max)
      }else{
        auto.limits <- round2any(max(abs(coefficients[[b]]), abs(df.pat$lp)), accuracy = accuracy, f = ceiling)
      }
    }else{ #not show.betas
      if(error.bar){
        auto.limits_min <- round2any(max(abs(df.pat$lp.min)), accuracy = accuracy, f = ceiling)
        auto.limits_max <- round2any(max(abs(df.pat$lp.max)), accuracy = accuracy, f = ceiling)
        auto.limits <- max(auto.limits_min, auto.limits_max)
      }else{
        auto.limits <- round2any(max(abs(df.pat$lp)), accuracy = accuracy, f = ceiling)
      }
    }

    ggp <- ggplot(df.pat, aes(x = var, y = lp, fill = lp, color = 1)) +
      geom_bar(stat = "identity", position = "dodge")

    if(error.bar){
      ggp <- ggp + geom_errorbar(aes(ymin=lp.min, ymax=lp.max), width=.35, position=position_dodge(.2))
    }

    if(!show.betas){
      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_fill_gradient2(low = RColorConesa::getConesaPalettes()$warm["blue"],
                                          mid = "white", midpoint = 0,
                                          high = RColorConesa::getConesaPalettes()$warm["magenta"],
                                          limits = c(-1*auto.limits,auto.limits), name = "Beta value")
      }else{
        ggp <- ggp + scale_fill_gradient2(low = "blue",
                                          mid = "white", midpoint = 0,
                                          high = "red",
                                          limits = c(-1*auto.limits,auto.limits), name = "Beta value")
      }
    }

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp <- ggp + RColorConesa::scale_fill_conesa(palette = "warm", continuous = T)
    }

    ggp <- ggp + guides(color= "none")
    ggp <- ggp + ylab(label = "Linear Predictor")
    ggp <- ggp + xlab(label = "Variables")
    ggp <- ggp + ggtitle(label = paste0("Observation - ", rownames(new_observation[[b]])))

    if(length(unique(df.pat$var))>15){
      ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    }

    if(show.betas){

      auto.limits.lp <- max(abs(min(df.pat$lp.max)), abs(max(df.pat$lp.max)))
      ggp.aux <- ggp + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits.lp, auto.limits.lp))

      ggp.aux2 <- ggp.simulated_beta$plot[[b]]
      ggp.aux2 <- ggp.aux2 + guides(fill = "none")
      suppressMessages(
        ggp.aux2 <- ggp.aux2 + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))
      )

      sign.beta <- coefficients[[b]]$value>0
      names(sign.beta)<-rownames(coefficients[[b]])
      sign.pat <- df.pat$lp>0
      same.sign <- sign.beta == sign.pat
      same.sign <- same.sign[rownames(ggp.simulated_beta$plot[[b]]$data)]

      ggp.aux$mapping$fill[[2]] <- same.sign
      ggp.aux <- ggp.aux + guides(fill = guide_legend(title="Same beta direction:")) + theme(legend.position="left")

      #overwriting fill generates a message
      suppressMessages({
        if(requireNamespace("RColorConesa", quietly = TRUE)){
          ggp.aux <- ggp.aux + RColorConesa::scale_fill_conesa(reverse = T)
        }else{
          ggp.aux <- ggp.aux + scale_fill_discrete()
        }
      })

      ggp <- ggpubr::ggarrange(ggp.aux, ggp.aux2, ncol = 2, widths = c(0.5, 0.5), align = "h")
    }

    lst_plots[[b]] <- ggp
    lst_lp.var[[b]] <- lp.new_observation_variable

  }

  return(list(plot = lst_plots, lp.var = lst_lp.var, norm_pat = norm_patient, pat = new_observation))

}


#### ### ### ###
# KAPLAN MEIER #
#### ### ### ###

#' getAutoKM.list
#' @description Run the function "getAutoKM" for a list of models. More information in "?getAutoKM".
#'
#' @param type Character. Kaplan Meier for complete model linear predictor ("LP"), for PLS components ("COMP") or for original variables ("VAR") (default: LP).
#' @param lst_models List of HDcox models.
#' @param comp Numeric vector. Vector of length two. Select which components to plot (default: c(1,2)).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: 10).
#' @param ori_data Logical. Compute the Kaplan-Meier plot with the raw-data or the normalize-data to compute the best cut-point for splitting the data into two groups. Only used when type = "VAR" (default: TRUE).
#' @param BREAKTIME Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points. If BREAKTIME = NULL, "n.breaks" is used (default: NULL).
#' @param n.breaks Numeric. If BREAKTIME is NULL, "n.breaks" is the number of time-break points to compute (default: 20).
#' @param only_sig Logical. If "only_sig" = TRUE, then only significant log-rank test variables are returned (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param title Character. Kaplan-Meier plot title (default: NULL).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getAutoKM.list(type = "LP", lst_models)
#' }

getAutoKM.list <- function(type = "LP", lst_models, comp = 1:2, top = NULL, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){
  if(!type %in% c("LP", "COMP", "VAR", "LPVAR")){
    stop("Type parameters must be one of the following: LP, COMP, VAR or LPVAR")
  }

  if(type %in% c("LP")){
    lst <- purrr::map(lst_models, ~getLPKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "COMP"){

    if(all(unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods))){
      sub_lst_models <- lst_models
    }else{
      sub_lst_models <- lst_models[unlist(purrr::map(lst_models, function(x){x$class})) %in% pkg.env$pls_methods]
      if(verbose){
        message(paste0("Model ", paste0(names(lst_models[!unlist(purrr::map(lst_models, function(x){x$class})) %in% pkg.env$pls_methods]), collapse = ", "), " are not based in PLS methodology. Other models computed."))
      }
    }

    lst <- purrr::map(sub_lst_models, ~getCompKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "VAR"){
    lst <- purrr::map(lst_models, ~getVarKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "LPVAR"){
    lst <- purrr::map(lst_models, ~getLPVarKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }
  return(lst)
}

#' getAutoKM
#' @description Performs a Kaplan-Meier plot for the HDcox model selected. KM could be based on model Linear Predictor value ("LP"), based on PLS-COX component ("COMP") or at original variable level ("VAR").
#'
#' @param type Character. Kaplan Meier for complete model linear predictor ("LP"), for PLS components ("COMP") or for original variables ("VAR") (default: LP).
#' @param model HDcox model.
#' @param comp Numeric vector. Vector of length two. Select which components to plot (default: c(1,2)).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: 10).
#' @param ori_data Logical. Compute the Kaplan-Meier plot with the raw-data or the normalize-data to compute the best cut-point for splitting the data into two groups. Only used when type = "VAR" (default: TRUE).
#' @param BREAKTIME Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points. If BREAKTIME = NULL, "n.breaks" is used (default: NULL).
#' @param n.breaks Numeric. If BREAKTIME is NULL, "n.breaks" is the number of time-break points to compute (default: 20).
#' @param only_sig Logical. If "only_sig" = TRUE, then only significant log-rank test variables are returned (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param title Character. Kaplan-Meier plot title (default: NULL).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getAutoKM(type = "LP", model)
#' }

getAutoKM <- function(type = "LP", model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){
  if(!type %in% c("LP", "COMP", "VAR")){
    stop("Type parameters must be one of the following: LP, COMP or VAR")
  }

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(length(comp)==1){
    comp <- 1:comp
  }

  if(type == "LP"){
    return(getLPKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "COMP"){
    return(getCompKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "VAR"){
    return(getVarKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "LPVAR"){
    return(getLPVarKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, n.breaks = n.breaks, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }
}

getLPKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  if(length(comp)==1){
    comp <- 1:comp
  }

  if(attr(model, "model") %in% c(pkg.env$classical_methods, pkg.env$pls_methods, pkg.env$multiblock_methods)){

    if(all(is.null(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

  }else{
    if(verbose){
      message("Model not have components or is not a HDcox object.")
    }
    return(NA)
  }

  #select data
  vars_data <- as.data.frame(model$survival_model$fit$linear.predictors)
  rownames(vars_data) <- rownames(model$X$data)
  colnames(vars_data) <- "LP"

  vars_num <- vars_data
  if(all(dim(vars_num)>0)){
    info_logrank_num <- getLogRank_NumVariables(data = vars_num, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
  }else{
    info_logrank_num <- NULL
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / n.breaks
  }

  d <- info_logrank_num$df_numASqual
  rownames(d) <- rownames(model$X$data)
  v_names <- info_logrank_num$df_nvar_lrtest[,1:2]

  if(all(is.null(d) & is.null(v_names))){
    message("Instead of LP Kaplan-Meier curve, Survival function, Hazard Curve and Cumulative Hazard will be returned.")
  }

  LST_SPLOT <- plot_survivalplot.qual(data = d,
                                      sdata = data.frame(model$Y$data),
                                      BREAKTIME = BREAKTIME,
                                      cn_variables = v_names$Variable,
                                      name_data = NULL, title = title)

  return(list(info_logrank_num = info_logrank_num, LST_PLOTS = LST_SPLOT))

}

getCompKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  if(length(comp)==1){
    comp <- 1:comp
  }

  # DFCALLS
  vars <- lst_vars <- info_logrank_qual <- NULL

  if(attr(model, "model") %in% pkg.env$pls_methods){

    if(!all(is.null(model$survival_model))){
      vars <- names(model$survival_model$fit$coefficients)
    }else{
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){

    if(!all(is.null(model$survival_model))){
      for(b in names(model$X$data)){
        if(attr(model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
          lst_vars[[b]] <- colnames(model[[4]][[b]]$X$W.star)
          keep <- which(paste0(lst_vars[[b]],"_",b) %in% names(model$survival_model$fit$coefficients))
          lst_vars[[b]] <- lst_vars[[b]][keep]
        }else{
          lst_vars[[b]] <- colnames(model$X$W.star[[b]])
          keep <- which(paste0(lst_vars[[b]],"_",b) %in% names(model$survival_model$fit$coefficients))
          lst_vars[[b]] <- lst_vars[[b]][keep]
        }
      }
      vars <- names(model$survival_model$fit$coefficients)
    }else{
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

  }else{
    if(verbose){
      message("Model not have components or is not a HDcox object.")
    }
    return(NA)
  }

  #select original or scale data - top X of each component, takes all of them
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    unique_vars <- deleteIllegalChars(unique(unlist(vars)))
    # vars %*% coeff to get component LP
    cn_aux <- colnames(as.data.frame(model$X$scores[rownames(model$X$scores),unique_vars,drop=F]))
    sc_aux <- as.data.frame(model$X$scores[rownames(model$X$scores),unique_vars,drop=F])
    coeff_aux <- model$survival_model$fit$coefficients[cn_aux]
    if(length(names(coeff_aux))>1){
      vars_data <- NULL
      for(cn in colnames(sc_aux)){
        vars_data <- cbind(vars_data, as.matrix(sc_aux[,cn,drop=F]) %*% coeff_aux[cn])
      }
      colnames(vars_data) <- names(unique_vars)
    }else{
      vars_data <- as.matrix(sc_aux) %*% coeff_aux
      colnames(vars_data) <- names(unique_vars)
    }
  }else{
    vars_data <- list()
    for(b in names(model$X$data)){
      # vars %*% coeff to get component LP
      unique_vars <- deleteIllegalChars(unique(unlist(lst_vars[[b]])))
      if(length(unique_vars)==0){next}#no components selected
      if(attr(model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
        cn_aux <- colnames(as.data.frame(model[[4]][[b]]$X$scores[rownames(model[[4]][[b]]$X$scores),unique_vars,drop=F]))
        sc_aux <- as.data.frame(model[[4]][[b]]$X$scores[rownames(model[[4]][[b]]$X$scores),unique_vars,drop=F])
        coeff_aux <- model$survival_model$fit$coefficients[paste0(cn_aux, "_", b)]
        if(length(names(coeff_aux))>1){
          vars_data[[b]] <- NULL
          # if coeff_aux has comp_1_genes and comp_10_genes, both start by comp_1
          # new colnames vector to match
          new_coeff_names <- names(coeff_aux)
          new_coeff_names <- unlist(lapply(new_coeff_names, function(x){paste0(strsplit(x, "_")[[1]][1], "_", strsplit(x, "_")[[1]][2])}))
          for(cn in colnames(sc_aux)){
            idx <- which(new_coeff_names %in% cn)
            vars_data[[b]] <- cbind(vars_data[[b]], as.matrix(sc_aux[,cn,drop=F]) %*% coeff_aux[idx])
          }
          colnames(vars_data[[b]]) <- names(unique_vars)
        }else{
          vars_data[[b]] <- as.matrix(sc_aux) %*% coeff_aux
          colnames(vars_data[[b]]) <- names(unique_vars)
        }
      }else{
        cn_aux <- colnames(as.data.frame(model$X$scores[[b]][rownames(model$X$scores[[b]]),unique_vars,drop=F]))
        sc_aux <- as.data.frame(model$X$scores[[b]][rownames(model$X$scores[[b]]),unique_vars,drop=F])
        coeff_aux <- model$survival_model$fit$coefficients[paste0(cn_aux, "_", b)]
        if(length(names(coeff_aux))>1){
          vars_data[[b]] <- NULL
          # if coeff_aux has comp_1_genes and comp_10_genes, both start by comp_1
          # new colnames vector to match
          new_coeff_names <- names(coeff_aux)
          new_coeff_names <- unlist(lapply(new_coeff_names, function(x){paste0(strsplit(x, "_")[[1]][1], "_", strsplit(x, "_")[[1]][2])}))
          for(cn in colnames(sc_aux)){
            idx <- which(new_coeff_names %in% cn)
            vars_data[[b]] <- cbind(vars_data[[b]], as.matrix(sc_aux[,cn,drop=F]) %*% coeff_aux[idx])
          }
          colnames(vars_data[[b]]) <- names(unique_vars)
        }else{
          vars_data[[b]] <- as.matrix(sc_aux) %*% coeff_aux
          colnames(vars_data[[b]]) <- names(unique_vars)
        }
      }
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    vars_num <- vars_data

    if(all(dim(vars_num)>0)){
      info_logrank_num <- getLogRank_NumVariables(data = vars_num, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
    }else{
      info_logrank_num <- NULL
    }
  }else{
    info_logrank_num <- list()
    vars_num <- list()
    for(b in names(model$X$data)){
      if(!b %in% names(vars_data)){next}
      vars_num[[b]] <- vars_data[[b]]

      if(all(dim(vars_num[[b]]))>0){
        info_logrank_num[[b]] <- getLogRank_NumVariables(data = vars_num[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
      }else{
        info_logrank_num[[b]] <- NULL
      }
    }
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / n.breaks
  }

  ##join data
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
      d <- info_logrank_num$df_numASqual
      v_names <- info_logrank_num$df_nvar_lrtest[,1:2]
  }else{
    v_names <- list()
    d <- list()
    for(b in names(model$X$data)){
      d[[b]] <- info_logrank_num[[b]]$df_numASqual
      v_names[[b]] <- info_logrank_num[[b]]$df_nvar_lrtest[,1:2]
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(only_sig){

      if(length(v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable)==0){
        if(verbose){
          cat("Any variable has a significant log-rank test value. Survival function, Hazard Curve and Cumulative Hazard plots will be returned.")
        }
      }

      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                           sdata = data.frame(model$Y$data),
                                           BREAKTIME = BREAKTIME,
                                           cn_variables = v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable,
                                           name_data = NULL, title = title)
    }else{
      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                           sdata = data.frame(model$Y$data),
                                           BREAKTIME = BREAKTIME,
                                           cn_variables = v_names$Variable,
                                           name_data = NULL, title = title)
    }
  }else{
    LST_SPLOT <- list()
    for(b in names(model$X$data)){
      if(only_sig){
        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                  sdata = data.frame(model$Y$data),
                                                  BREAKTIME = BREAKTIME,
                                                  cn_variables = v_names[[b]][v_names[[b]]$`P-Val (Log Rank)` <= alpha,]$Variable,
                                                  name_data = NULL, title = title)
      }else{
        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                  sdata = data.frame(model$Y$data),
                                                  BREAKTIME = BREAKTIME,
                                                  cn_variables = v_names[[b]]$Variable,
                                                  name_data = NULL, title = title)
      }
    }

  }

  return(list(info_logrank_num = info_logrank_num, LST_PLOTS = LST_SPLOT))

}

getLPVarKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  if(length(comp)==1){
    comp <- 1:comp
  }

  message("LPVAR only implemented for PLS methods. Results are pretty similar to work with ORIGINAL variables.")

  if(attr(model, "model") %in% pkg.env$pls_methods){

    if(all(is.null(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    #selecting pseudo betas
    pseudo_betas <- plot_pseudobeta(model = model,
                                    error.bar = T, onlySig = only_sig, alpha = alpha,
                                    zero.rm = F, auto.limits = F, top = top,
                                    show_percentage = F, size_percentage = 3)
    names_top <- pseudo_betas$plot$data$variables
    pseudo_betas$beta <- pseudo_betas$beta[names_top,]

    pseudo_betas$plot <- NULL
    vars <- rownames(pseudo_betas$beta)

  }else if(attr(model, "model") %in% pkg.env$classical_methods){

    if(all(is.na(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    #in classical methods, select selected variables
    df <- as.data.frame(summary(model$survival_model$fit)[7]$coefficients)
    vars <- rownames(df[order(df$`Pr(>|z|)`, decreasing = F),])[1:min(top, nrow(df))]

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){

    if(all(is.na(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    lst_vars <- list()
    for(b in names(model$X$data)){
      vars <- list()

      if(attr(model, "model") %in% pkg.env$sb.splsicox){
        aux <- model$list_spls_models[[b]]
      }else if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
        aux <- model$list_spls_models[[b]]
      }

      if(attr(model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){

        message("ARREGLAR PARA SB")

      }else if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

        message("ARREGLAR PARA MB")

      }

      # names(vars) <- as.character(1:length(vars))
      # lst_vars[[b]] <- vars

    }

  }

  #select original or scale data - top X of each component, takes all of them
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    unique_vars <- deleteIllegalChars(unique(unlist(vars)))
    if(ori_data){
      vars_data <- as.data.frame(model$X_input[rownames(model$X$data),unique_vars,drop=F])
    }else{
      vars_data <- as.data.frame(model$X$data[,unique_vars,drop=F])
    }

    vars_data <- as.data.frame(scale(vars_data, center = model$X$x.mean[unique_vars], scale = model$X$x.sd[unique_vars]))

    #GET LP_VAR per each patient
    if(attr(model, "model") %in% pkg.env$pls_methods){

      # lp <- model$survival_model$fit$linear.predictors)
      # lp_calculated <- vars_data[,rownames(pseudo_betas$beta)] %*% pseudo_betas$beta$value ## COMPROBATION LP ## !!!!

      aux <- NULL
      for(cn in rownames(pseudo_betas$beta)){
        aux <- cbind(aux, vars_data[,cn,drop=T] * pseudo_betas$beta[cn,]$value)
      }
      aux <- as.data.frame(aux)
      colnames(aux) <- rownames(pseudo_betas$beta)
      vars_data <- aux
    }

  }else{
    vars_data <- list()
    for(b in names(model$X$data)){
      unique_vars <- deleteIllegalChars(unique(unlist(lst_vars[[b]])))
      if(ori_data){
        vars_data[[b]] <- as.data.frame(model$X_input[[b]][rownames(model$X$data[[b]]),unique_vars,drop=F])
      }else{
        vars_data[[b]] <- as.data.frame(model$X$data[[b]][,unique_vars,drop=F])
      }
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(attr(model, "model") %in% pkg.env$pls_methods){
      colnames(vars_data) <- paste0("LP_", colnames(vars_data))
    }

    names_qual <- apply(vars_data, 2, function(x){all(x %in% c(0,1))})
    vars_qual <- vars_data[,names_qual,drop=F]
    vars_num <- vars_data[,!names_qual,drop=F]

    if(all(dim(vars_qual)>0)){
      for(cn in colnames(vars_qual)){vars_qual[,cn] <- factor(vars_qual[,cn], levels = c(0, 1))}
      info_logrank_qual <- getLogRank_QualVariables(data = vars_qual, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL)
    }else{
      info_logrank_qual = NULL
    }

    if(all(dim(vars_num)>0)){
      info_logrank_num <- getLogRank_NumVariables(data = vars_num, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
    }else{
      info_logrank_num <- NULL
    }
  }else{
    info_logrank_qual <- list()
    info_logrank_num <- list()
    vars_qual <- list()
    vars_num <- list()
    for(b in names(model$X$data)){
      names_qual <- apply(vars_data[[b]], 2, function(x){all(x %in% c(0,1))})
      vars_qual[[b]] <- vars_data[[b]][,names_qual,drop=F]
      vars_num[[b]] <- vars_data[[b]][,!names_qual,drop=F]

      if(all(dim(vars_qual[[b]]))>0){
        for(cn in colnames(vars_qual[[b]])){vars_qual[[b]][,cn] <- factor(vars_qual[[b]][,cn], levels = c(0, 1))}
        info_logrank_qual[[b]] <- getLogRank_QualVariables(data = vars_qual[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL)
      }else{
        info_logrank_qual[[b]] = NULL
      }

      if(all(dim(vars_num[[b]]))>0){
        info_logrank_num[[b]] <- getLogRank_NumVariables(data = vars_num[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
      }else{
        info_logrank_num[[b]] <- NULL
      }
    }
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / n.breaks
  }

  ##join data
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(all(dim(vars_qual))>0 & all(dim(vars_num)>0)){
      d <- cbind(vars_qual, info_logrank_num$df_numASqual)
      v_names <- info_logrank_num$df_nvar_lrtest[,1:2]
      v_names <- rbind(v_names, info_logrank_qual)

    }else if(all(dim(vars_qual)>0)){
      d <- vars_qual
      v_names <- info_logrank_qual

    }else{
      d <- info_logrank_num$df_numASqual
      v_names <- info_logrank_num$df_nvar_lrtest[,1:2]
    }
  }else{
    v_names <- list()
    d <- list()
    for(b in names(model$X$data)){
      if(all(dim(vars_qual[[b]]))>0 & all(dim(vars_num[[b]])>0)){
        d[[b]] <- cbind(vars_qual[[b]], info_logrank_num[[b]]$df_numASqual)
        v_names[[b]] <- info_logrank_num[[b]]$df_nvar_lrtest[,1:2]
        v_names[[b]] <- rbind(v_names[[b]], info_logrank_qual[[b]])

      }else if(all(dim(vars_qual[[b]])>0)){
        d[[b]] <- vars_qual[[b]]
        v_names[[b]] <- info_logrank_qual[[b]]

      }else{
        d[[b]] <- info_logrank_num[[b]]$df_numASqual
        v_names[[b]] <- info_logrank_num[[b]]$df_nvar_lrtest[,1:2]
      }
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(only_sig){

      if(length(v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable)==0){
        if(verbose){
          cat("All variables has a non-significant log-rank test value. Survival function, Hazard Curve and Cumulative Hazard plots will be returned.")
        }
      }

      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                          sdata = data.frame(model$Y$data),
                                          BREAKTIME = BREAKTIME,
                                          cn_variables = v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable,
                                          name_data = NULL, title = title)
    }else{
      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                          sdata = data.frame(model$Y$data),
                                          BREAKTIME = BREAKTIME,
                                          cn_variables = v_names$Variable,
                                          name_data = NULL, title = title)
    }
  }else{
    LST_SPLOT <- list()
    for(b in names(model$X$data)){
      if(only_sig){

        if(length(v_names[[b]][v_names[[b]]$`P-Val (Log Rank)` <= alpha,]$Variable)==0){
          if(verbose){
            cat("Any variable has a significant log-rank test value. Survival function, Hazard Curve and Cumulative Hazard plots will be returned.")
          }
        }

        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                 sdata = data.frame(model$Y$data),
                                                 BREAKTIME = BREAKTIME,
                                                 cn_variables = v_names[[b]][v_names[[b]]$`P-Val (Log Rank)` <= alpha,]$Variable,
                                                 name_data = NULL, title = title)
      }else{
        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                 sdata = data.frame(model$Y$data),
                                                 BREAKTIME = BREAKTIME,
                                                 cn_variables = v_names[[b]]$Variable,
                                                 name_data = NULL, title = title)
      }
    }

  }

  return(list(info_logrank_qual = info_logrank_qual, info_logrank_num = info_logrank_num, LST_PLOTS = LST_SPLOT))

}

getVarKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, n.breaks = 20, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  if(length(comp)==1){
    comp <- 1:comp
  }

  if(attr(model, "model") %in% pkg.env$pls_methods){

    if(all(is.null(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    #selecting the variables with a W.star different than 0
    vars_data <- list()
    vars <- list()
    for(c in comp){
      if(ncol(model$X$W.star)>=c){
        rn <- rownames(model$X$W.star[model$X$W.star[,c]!=0,c,drop=F])
        vars[[c]] <- rownames(model$X$W.star[rn,c,drop=F])[order(abs(model$X$W.star[rn,c]), decreasing = T)][1:min(top, length(rn))]
      }else{
        break
      }
    }

  }else if(attr(model, "model") %in% pkg.env$classical_methods){

    if(all(is.na(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    #in classical methods, select selected variables
    df <- as.data.frame(summary(model$survival_model$fit)[7]$coefficients)
    vars <- rownames(df[order(df$`Pr(>|z|)`, decreasing = F),])[1:min(top, nrow(df))]

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){

    if(all(is.na(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

    lst_vars <- list()
    for(b in names(model$X$data)){
      vars <- list()
      vars_data <- list()

      if(attr(model, "model") %in% pkg.env$sb.splsicox){
        aux <- model$list_spls_models[[b]]
      }else if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
        aux <- model$list_spls_models[[b]]
      }

      if(attr(model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){

        for(c in comp){
          if(ncol(aux$X$W.star)>=c){
            rn <- rownames(aux$X$W.star[aux$X$W.star[,c]!=0,c,drop=F])
            vars[[c]] <- rownames(aux$X$W.star[rn,,drop=F])[order(abs(aux$X$W.star[rn,c]), decreasing = T)][1:min(top, length(rn))]
          }else{
            break
          }
        }

      }else if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox)){

        for(c in comp){
          if(ncol(model$X$W.star[[b]])>=c){
            rn <- rownames(model$X$W.star[[b]][model$X$W.star[[b]][,c]!=0,c,drop=F])
            vars[[c]] <- rownames(model$X$W.star[[b]][rn,,drop=F])[order(abs(model$X$W.star[[b]][rn,c]), decreasing = T)][1:min(top, length(rn))]
          }else{
            break
          }
        }

      }

      names(vars) <- as.character(1:length(vars))
      lst_vars[[b]] <- vars
    }

  }

  #select original or scale data - top X of each component, takes all of them
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    unique_vars <- deleteIllegalChars(unique(unlist(vars)))
    if(ori_data){
      vars_data <- as.data.frame(model$X_input[rownames(model$X$data),unique_vars,drop=F])
    }else{
      vars_data <- as.data.frame(model$X$data[,unique_vars,drop=F])
    }
  }else{
    vars_data <- list()
    for(b in names(model$X$data)){
      unique_vars <- deleteIllegalChars(unique(unlist(lst_vars[[b]])))
      if(ori_data){
        vars_data[[b]] <- as.data.frame(model$X_input[[b]][rownames(model$X$data[[b]]),unique_vars,drop=F])
      }else{
        vars_data[[b]] <- as.data.frame(model$X$data[[b]][,unique_vars,drop=F])
      }
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    names_qual <- apply(vars_data, 2, function(x){all(x %in% c(0,1))})
    vars_qual <- vars_data[,names_qual,drop=F]
    vars_num <- vars_data[,!names_qual,drop=F]

    if(all(dim(vars_qual)>0)){
      for(cn in colnames(vars_qual)){vars_qual[,cn] <- factor(vars_qual[,cn], levels = c(0, 1))}
      info_logrank_qual <- getLogRank_QualVariables(data = vars_qual, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL)
    }else{
      info_logrank_qual = NULL
    }

    if(all(dim(vars_num)>0)){
      info_logrank_num <- getLogRank_NumVariables(data = vars_num, sdata = data.frame(model$Y$data),
                                                  VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
    }else{
      info_logrank_num <- NULL
    }
  }else{
    info_logrank_qual <- list()
    info_logrank_num <- list()
    vars_qual <- list()
    vars_num <- list()
    for(b in names(model$X$data)){
      names_qual <- apply(vars_data[[b]], 2, function(x){all(x %in% c(0,1))})
      vars_qual[[b]] <- vars_data[[b]][,names_qual,drop=F]
      vars_num[[b]] <- vars_data[[b]][,!names_qual,drop=F]

      if(all(dim(vars_qual[[b]]))>0){
        for(cn in colnames(vars_qual[[b]])){vars_qual[[b]][,cn] <- factor(vars_qual[[b]][,cn], levels = c(0, 1))}
        info_logrank_qual[[b]] <- getLogRank_QualVariables(data = vars_qual[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL)
      }else{
        info_logrank_qual[[b]] = NULL
      }

      if(all(dim(vars_num[[b]]))>0){
        info_logrank_num[[b]] <- getLogRank_NumVariables(data = vars_num[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
      }else{
        info_logrank_num[[b]] <- NULL
      }
    }
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / n.breaks
  }

  ##join data
  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(all(dim(vars_qual))>0 & all(dim(vars_num)>0)){
      d <- cbind(vars_qual, info_logrank_num$df_numASqual)
      v_names <- info_logrank_num$df_nvar_lrtest[,1:2]
      v_names <- rbind(v_names, info_logrank_qual)

    }else if(all(dim(vars_qual)>0)){
      d <- vars_qual
      v_names <- info_logrank_qual

    }else{
      d <- info_logrank_num$df_numASqual
      v_names <- info_logrank_num$df_nvar_lrtest[,1:2]
    }
  }else{
    v_names <- list()
    d <- list()
    for(b in names(model$X$data)){
      if(all(dim(vars_qual[[b]]))>0 & all(dim(vars_num[[b]])>0)){
        d[[b]] <- cbind(vars_qual[[b]], info_logrank_num[[b]]$df_numASqual)
        v_names[[b]] <- info_logrank_num[[b]]$df_nvar_lrtest[,1:2]
        v_names[[b]] <- rbind(v_names[[b]], info_logrank_qual[[b]])

      }else if(all(dim(vars_qual[[b]])>0)){
        d[[b]] <- vars_qual[[b]]
        v_names[[b]] <- info_logrank_qual[[b]]

      }else{
        d[[b]] <- info_logrank_num[[b]]$df_numASqual
        v_names[[b]] <- info_logrank_num[[b]]$df_nvar_lrtest[,1:2]
      }
    }
  }

  if(!attr(model, "model") %in% pkg.env$multiblock_methods){
    if(only_sig){

      if(length(v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable)==0){
        if(verbose){
          cat("Any variable has a significant log-rank test value. Survival function, Hazard Curve and Cumulative Hazard plots will be returned.")
        }
      }

      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                          sdata = data.frame(model$Y$data),
                                          BREAKTIME = BREAKTIME,
                                          cn_variables = v_names[v_names$`P-Val (Log Rank)` <= alpha,]$Variable,
                                          name_data = NULL, title = title)
    }else{
      LST_SPLOT <- plot_survivalplot.qual(data = d,
                                          sdata = data.frame(model$Y$data),
                                          BREAKTIME = BREAKTIME,
                                          cn_variables = v_names$Variable,
                                          name_data = NULL, title = title)
    }
  }else{
    LST_SPLOT <- list()
    for(b in names(model$X$data)){
      if(only_sig){

        if(verbose & length(v_names[[b]][v_names[[b]]$`P-Val (Log Rank)` <= alpha,]$Variable)==0){
          cat(paste0("Any variable has a significant log-rank test value for block '", b, "'. Survival function, Hazard Curve and Cumulative Hazard plots will be returned."))
        }

        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                 sdata = data.frame(model$Y$data),
                                                 BREAKTIME = BREAKTIME,
                                                 cn_variables = v_names[[b]][v_names[[b]]$`P-Val (Log Rank)` <= alpha,]$Variable,
                                                 name_data = NULL, title = title)
      }else{
        LST_SPLOT[[b]] <- plot_survivalplot.qual(data = d[[b]],
                                                 sdata = data.frame(model$Y$data),
                                                 BREAKTIME = BREAKTIME,
                                                 cn_variables = v_names[[b]]$Variable,
                                                 name_data = NULL, title = title)
      }
    }

  }

  return(list(info_logrank_qual = info_logrank_qual, info_logrank_num = info_logrank_num, LST_PLOTS = LST_SPLOT))

}

getLogRank_QualVariables <- function(data, sdata, VAR_EVENT, name_data = NULL){

  LST_QVAR_SIG = NULL #significant qualitative variables

  if(is.null(name_data)){
    data <- data
  }else{
    data <- data[[name_data]]
  }

  for(cn in colnames(data)){
    if(cn==VAR_EVENT){ #skip outcome variable
      next
    }

    variable <- data[,cn] #select the variable

    tbl <- as.data.frame(sort(table(variable)))
    if(all(dim(tbl)==c(1,1))){
      next #just one factor
    }
    tbl$Rel <- round(tbl$Freq/sum(tbl$Freq), digits = 4)*100

    indexNONA <- which(!is.na(variable))

    aux <- cbind(sdata[indexNONA,], variable[indexNONA])
    colnames(aux)[3] <- cn

    #SA
    f = as.formula(paste0("Surv(time = time, event = event) ~ ", "`",cn,"`"))
    kmsurvival <- tryCatch(
      # Specifying expression
      expr = {
        survminer::surv_fit(formula = f, data = aux)
      },
      # Specifying error message
      error = function(e){
        message(paste0("Problems at variable ", cn, ".\n",e$message),". Try to change the name of the variable.")
        NA
      }
    )

    if(all(is.na(kmsurvival))){
      LST_QVAR_SIG <- rbind(LST_QVAR_SIG, c(cn, NA))
      next
    }else{
      pval <- surv_pvalue(kmsurvival)
      LST_QVAR_SIG <- rbind(LST_QVAR_SIG, c(cn, round(pval$pval,4)))
    }

  }

  LST_QVAR_SIG <- as.data.frame(LST_QVAR_SIG)
  LST_QVAR_SIG[,2] <- as.numeric(LST_QVAR_SIG[,2])

  if(exists("VAR_DESCRIPTION")){
    colnames(LST_QVAR_SIG) <- c("Variable", "P-Val (Log Rank)", "Description")
  }else{
    colnames(LST_QVAR_SIG) <- c("Variable", "P-Val (Log Rank)")
  }

  LST_QVAR_SIG <- LST_QVAR_SIG[order(LST_QVAR_SIG$`P-Val (Log Rank)`),]

  return(LST_QVAR_SIG)
}

getLogRank_NumVariables <- function(data, sdata, VAR_EVENT, name_data = NULL, minProp = 0.1, ROUND_CP = 4){

  if(is.null(name_data)){
    data <- data
  }else{
    data <- data[[name_data]]
  }

  LST_NVAR_SIG = NULL
  df_qualnumvars = NULL
  for(cn in colnames(data)){
    variable <- data[,cn,drop=T]

    auxData <- cbind(sdata, variable)

    cn_ori <- cn
    #### Formula cannot manage -,+,* symbols in cn
    cn <- transformIllegalChars(cn)

    colnames(auxData)[3] <- cn

    # Determine the optimal cutpoint for continuous variables, using the maximally selected rank statistics from the 'maxstat' R package.
    minProp = minProp #we have to establish a minimum number of patients per group in 0-1

    res.cut <- tryCatch(
      expr = {
        survminer::surv_cutpoint(auxData, time="time", event="event", variables = cn, minprop = minProp)
      },
      # Specifying error message
      error = function(e){
        message(paste0("Problems with variable '",cn,"'", ": ", e))
        NA
      }
    )

    if(all(is.na(res.cut))){
      next
    }

    if(res.cut$cutpoint[1,1]<=0){
      cutpoint_value <- round2any(res.cut$cutpoint[1,1], accuracy = 1/(10^ROUND_CP), f = ceiling)
    }else{
      cutpoint_value <- round(res.cut$cutpoint[1,1], ROUND_CP)
    }

    variable <- ifelse(variable>cutpoint_value, paste0("greater than ", cutpoint_value), paste0("lesser/equal than ", cutpoint_value))
    variable <- data.frame(factor(variable))
    colnames(variable) = cn_ori

    if(is.null(df_qualnumvars)){
      #colnames(variable) = cn_ori
      df_qualnumvars <- variable
      colnames(variable) = cn
    }else{
      #colnames(variable) = cn_ori
      df_qualnumvars <- cbind(df_qualnumvars, variable)
      colnames(variable) = cn
    }

    tbl <- as.data.frame(sort(table(variable)))
    tbl$Rel <- round(tbl$Freq/sum(tbl$Freq), digits = 4)*100

    #update of auxData with TRUE/FALSE
    indexNONA <- which(!is.na(variable))

    auxData <- cbind(sdata[indexNONA,], variable[indexNONA,])
    colnames(auxData)[3] <- cn

    #SA
    f = as.formula(paste0("Surv(time = time, event = event) ~ ", "`",cn,"`"))
    kmsurvival <- tryCatch(
      # Specifying expression
      expr = {
        survminer::surv_fit(formula = f, data = auxData)
      },
      # Specifying error message
      error = function(e){
        message(paste0("Problems at variable ", cn, ".\n",e$message),". Try to change the name of the variable.")
        NA
      }
    )

    if(all(is.na(kmsurvival))){
      LST_NVAR_SIG <- rbind(LST_NVAR_SIG, c(cn_ori, NA, NA))
      next
    }else{
      pval <- surv_pvalue(kmsurvival)
      LST_NVAR_SIG <- rbind(LST_NVAR_SIG, c(cn_ori, round(pval$pval,4), cutpoint_value))
    }

  }

  if(!is.null(LST_NVAR_SIG)){
    LST_NVAR_SIG <- as.data.frame(LST_NVAR_SIG)
    LST_NVAR_SIG[,2] <- as.numeric(LST_NVAR_SIG[,2])
    LST_NVAR_SIG[,3] <- as.numeric(LST_NVAR_SIG[,3])

    if(exists("VAR_DESCRIPTION")){
      colnames(LST_NVAR_SIG) <- c("Variable", "P-Val (Log Rank)", "Cutoff", "Description")
    }else{
      colnames(LST_NVAR_SIG) <- c("Variable", "P-Val (Log Rank)", "Cutoff")
    }

    LST_NVAR_SIG <- LST_NVAR_SIG[order(LST_NVAR_SIG$`P-Val (Log Rank)`),]
  }else{
    #any variable have been computed
    message("None of the variables have been selected for computing the Kaplan-Meier plot. The problem could be related to the 'minProp' value. Try to decrease it.")
  }

  return(list(df_numASqual = df_qualnumvars, df_nvar_lrtest = LST_NVAR_SIG))

}

plot_survivalplot.qual <- function(data, sdata, cn_variables, name_data = NULL, BREAKTIME = 5, title = NULL){

  lst_splots <- list()

  if(!length(cn_variables)==0){
    for(cn in cn_variables){
      if(is.null(name_data)){
        if(!cn %in% colnames(data)){
          message(paste0("Variable ", cn, " not found in data."))
          next
        }else{
          aux <- cbind(sdata, data[,cn])
        }
      }else{
        if(!cn %in% colnames(data[[name_data]])){
          message(paste0("Variable ", cn, " not found in data."))
          next
        }else{
          aux <- cbind(sdata, data[[name_data]][,cn])
        }
      }

      #delete NAs
      aux <- aux[!is.na(aux[,3]),]

      cn_ori <- cn
      #### Formula cannot manage -,+,* symbols in cn
      cn <- transformIllegalChars(cn)

      colnames(aux)[3] <- cn

      f = as.formula(paste0("Surv(time = time, event = event) ~ `", cn, "`"))

      kmsurvival <- tryCatch(
        # Specifying expression
        expr = {
          survminer::surv_fit(formula = f, data = aux)
        },
        # Specifying error message
        error = function(e){
          message(paste0("Problems at variable ", cn, ".\n",e$message),". Try to change the name of the variable.")
          NA
        }
      )

      if(all(is.na(kmsurvival))){
        next
      }

      ## change name kmsurvival to original
      # if(cn != cn_ori){
      #   #kmsurvival$strata
      #   aux_strata <- names(kmsurvival$strata)
      #   names_strata <- vapply(aux_strata, function(x) strsplit(x, "=")[[1]], FUN.VALUE = character(2))
      #   names_strata[names_strata==cn] <- cn_ori
      #   new_names <- apply(names_strata, 2, function(x){paste0(x, collapse = "=")})
      #   names(kmsurvival$strata) <- new_names
      #   #kmsurvival$call
      # }

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        colors <- RColorConesa::colorConesa(length(levels(data[,cn_ori])))
        names(colors) <- NULL
      }else{
        colors <- NULL
      }

      # GGSURVPLOT DOES NOT PRINT INTERVALS IF ALL DATA IS NOT SELECTED FOR RIBBON STYLE
      # IF PROBLEMS CHANGE TO STEP STYLE
      kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", palette = colors,
                                      conf.int = TRUE, ggtheme = theme_bw(), legend.labs = levels(aux[,cn]),
                                      conf.int.style = "ribbon",
                                      conf.int.alpha = 0.25,
                                      xlim = c(0, round2any(max(aux$time), 5, ceiling)),
                                      pval = T,
                                      surv.median.line = "hv", # Add medians survival
                                      risk.table = TRUE,
                                      legend.title = cn_ori,
                                      break.time.by = BREAKTIME,
                                      font.caption = 8,
                                      font.x = 10,
                                      font.y = 10,
                                      font.tickslab = 8,
                                      font.legend = 8,
                                      title = title)

      kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
        theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))

      lst_splots[[cn_ori]] <- kmplot
    }
  }else{
    f = as.formula("Surv(time = time, event = event) ~ 1")
    kmsurvival <- survminer::surv_fit(formula = f, data = sdata)

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      colors <- RColorConesa::colorConesa(1)
      names(colors) <- NULL
    } else {
      colors <- NULL
    }

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", palette = colors,
                                    conf.int = TRUE, ggtheme = theme_bw(),
                                    conf.int.style = "ribbon",
                                    conf.int.alpha = 0.25,
                                    xlim = c(0, round2any(max(sdata$time), 5, ceiling)),
                                    pval = T,
                                    surv.median.line = "hv", # Add medians survival
                                    risk.table = TRUE,
                                    title = "Survival Function",
                                    legend = "none",
                                    break.time.by = BREAKTIME,
                                    font.caption = 8,
                                    font.x = 10,
                                    font.y = 10,
                                    font.tickslab = 8,
                                    font.legend = 8)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["SurvivalFunction"]] <- kmplot

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", palette = colors, fun = "event",
                                    conf.int = TRUE, ggtheme = theme_bw(),
                                    conf.int.style = "ribbon",
                                    conf.int.alpha = 0.25,
                                    xlim = c(0, round2any(max(sdata$time), 5, ceiling)),
                                    pval = T,
                                    surv.median.line = "hv", # Add medians survival
                                    risk.table = TRUE,
                                    title = "Hazard Curve",
                                    legend = "none",
                                    break.time.by = BREAKTIME,
                                    font.caption = 8,
                                    font.x = 10,
                                    font.y = 10,
                                    font.tickslab = 8,
                                    font.legend = 8)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["HazardCurve"]] <- kmplot

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", palette = colors, fun = "cumhaz",
                                    conf.int = TRUE, ggtheme = theme_bw(),
                                    conf.int.style = "ribbon",
                                    conf.int.alpha = 0.25,
                                    xlim = c(0, round2any(max(sdata$time), 5, ceiling)),
                                    pval = T,
                                    surv.median.line = "hv", # Add medians survival
                                    risk.table = TRUE,
                                    xlab = "Time (Days)",
                                    ylab = "Cumulative Hazard",
                                    title = "Cumulative Hazard",
                                    legend = "none",
                                    break.time.by = BREAKTIME,
                                    font.caption = 8,
                                    font.x = 10,
                                    font.y = 10,
                                    font.tickslab = 8,
                                    font.legend = 8)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["CumulativeHazard"]] <- kmplot
  }

  return(lst_splots)
}

#### ### ### ### ### ##
# TEST - KAPLAN-MEIER #
#### ### ### ### ### ##

#' getCutoffAutoKM.list
#' @description Run the function "getCutoffAutoKM" for a list of models. More information in "?getCutoffAutoKM".
#'
#' @param lst_results List of lists. Result of getAutoKM.list() function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getCutoffAutoKM.list(lst_results)
#' }

getCutoffAutoKM.list <- function(lst_results){
  LST_RES <- purrr::map(lst_results, ~getCutoffAutoKM(.))
}

#' getCutoffAutoKM
#' @description Gets the cutoff value from the results of getAutoKM() functions.
#'
#' @param result List. Result of getAutoKM() function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getCutoffAutoKM(result)
#' }

getCutoffAutoKM <- function(result){

  if(all(is.null(result)) || all(is.na(result))){
    return(NULL)
  }

  # High dimensional
  if("df_nvar_lrtest" %in% names(result$info_logrank_num)){
    value <- result$info_logrank_num$df_nvar_lrtest$Cutoff
    names(value) <- result$info_logrank_num$df_nvar_lrtest$Variable
  }else{
    # MO
    value <- list()
    cont = 1
    for(b in names(result$info_logrank_num)){
      if(is.null(result$info_logrank_num[[b]]$df_nvar_lrtest)){
        return(NULL)
      }

      value[[cont]] <- result$info_logrank_num[[b]]$df_nvar_lrtest$Cutoff
      names(value[[cont]]) <- paste0(result$info_logrank_num[[b]]$df_nvar_lrtest$Variable, "_", b)
      cont = cont + 1
    }

    value <- unlist(value)

  }

  return(value)
}

#' getTestKM.list
#' @description Run the function "getTestKM" for a list of models. More information in "?getTestKM".
#'
#' @param lst_models List of HDcox model
#' @param X_test Numeric matrix or data.frame. Explanatory variables for test data (raw format). Qualitative variables must be transform into binary variables.
#' @param Y_test Numeric matrix or data.frame. Response variables for test data. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param lst_cutoff Numeric vector. Cutoff vector to split the observations into two groups for each model. Recommended to compute optimal cutoff value with getAutoKM() or getAutoKM.list() functions.
#' @param type Character. Kaplan Meier for complete model linear predictor ("LP"), for PLS components ("COMP") or for original variables ("VAR") (default: LP).
#' @param ori_data Logical. Compute the Kaplan-Meier plot with the raw-data or the normalize-data to compute the best cut-point for splitting the data into two groups. Only used when type = "VAR" (default: TRUE).
#' @param BREAKTIME Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points. If BREAKTIME = NULL, "n.breaks" is used (default: NULL).
#' @param n.breaks Numeric. If BREAKTIME is NULL, "n.breaks" is the number of time-break points to compute (default: 20).
#' @param title Character. Kaplan-Meier plot title (default: NULL).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_LP)
#' getTestKM.list(lst_models, X_test, Y_test, lst_cutoff)
#' }

getTestKM.list <- function(lst_models, X_test, Y_test, lst_cutoff, type = "LP", ori_data = T, BREAKTIME = NULL, n.breaks = 20, title = NULL, verbose = F){

  if(!type %in% c("LP", "COMP", "VAR")){
    stop("Type parameters must be one of the following: LP, COMP or VAR")
  }

  if(type == "COMP"){
    if(all(unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods))){
      sub_lst_models <- lst_models
    }else{
      sub_lst_models <- lst_models[unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)]
      if(verbose){
        message(paste0("Model ", paste0(names(lst_models[!unlist(purrr::map(lst_models, function(x){x$class})) %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)]), collapse = ", "), " are not based in PLS methodology. Other models computed."))
      }
    }
  }else{
    sub_lst_models <- lst_models
  }

  if(!length(sub_lst_models) == length(lst_cutoff) & !length(lst_cutoff) == 1){
    stop("List of models and list of cutoff must have the same length or list of cutoff must be just one value.")
  }

  LST_GGP <- NULL
  if(length(lst_cutoff)==1 && !isa(lst_cutoff, "list")){
    LST_GGP <- purrr::map(sub_lst_models, ~getTestKM(model = .,
                                                 X_test = X_test, Y_test = Y_test,
                                                 cutoff = lst_cutoff, type = type, ori_data = ori_data,
                                                 BREAKTIME = BREAKTIME, n.breaks = n.breaks, title = title))
  }else{
    LST_GGP <- purrr::map2(.x = sub_lst_models, .y = lst_cutoff, ~getTestKM(model = .x,
                                                                            X_test = X_test, Y_test = Y_test,
                                                                            cutoff = .y, type = type, ori_data = ori_data,
                                                                            BREAKTIME = BREAKTIME, n.breaks = n.breaks, title = title))
  }

  return(LST_GGP)

}

#' getTestKM
#' @description Perform a new Kaplan-Meier curve for test data, but using the same cutoff computed in the original model.
#'
#' @param model HDcox model.
#' @param X_test Numeric matrix or data.frame. Explanatory variables for test data (raw format). Qualitative variables must be transform into binary variables.
#' @param Y_test Numeric matrix or data.frame. Response variables for test data. Object must have two columns named as "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event observations.
#' @param cutoff Numeric. Cutoff value to split the observations into two groups. Recommended to compute optimal cutoff value with getAutoKM() function.
#' @param type Character. Kaplan Meier for complete model linear predictor ("LP"), for PLS components ("COMP") or for original variables ("VAR") (default: LP).
#' @param ori_data Logical. Compute the Kaplan-Meier plot with the raw-data or the normalize-data to compute the best cut-point for splitting the data into two groups. Only used when type = "VAR" (default: TRUE).
#' @param BREAKTIME Numeric. Size of time to split the data into "total_time / BREAKTIME + 1" points. If BREAKTIME = NULL, "n.breaks" is used (default: NULL).
#' @param n.breaks Numeric. If BREAKTIME is NULL, "n.breaks" is the number of time-break points to compute (default: 20).
#' @param title Character. Kaplan-Meier plot title (default: NULL).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cutoff <- getCutoffAutoKM(KM_RES_LP)
#' getTestKM(model, X_test, Y_test, cutoff)
#' }

getTestKM <- function(model, X_test, Y_test, cutoff, type = "LP", ori_data = T, BREAKTIME = NULL, n.breaks = 20, title = NULL){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(!type %in% c("LP", "COMP", "VAR")){
    stop("Type parameters must be one of the following: LP, COMP or VAR")
  }

  if(!is.numeric(cutoff)){
    message("cutoff parameter must be numeric. Returning NA")
    return(NA)
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(Y_test[,"time"]) - min(Y_test[,"time"])) / n.breaks
  }

  if(is.null(title)){
    title = attr(model, "model")
  }else{
    title = paste0(attr(model, "model"), " - ", title)
  }

  #create new variable
  if(type=="LP"){
    #predict scores X_test
    test_score <- predict(object = model, newdata = X_test)
    #predict LP using scores
    test_lp <- predict(model$survival_model$fit, newdata = as.data.frame(test_score))

    if(is.na(cutoff)){
      message("Cutoff not found for LP")
      return(NA)
    }

    txt_greater <- paste0("greater than ", cutoff)
    txt_lower <- paste0("lesser/equal than ", cutoff)

    LP <- ifelse(test_lp>cutoff, txt_greater, txt_lower)
    LP <- factor(LP)

    d <- as.data.frame(LP)
    colnames(d) <- type

    ggp <- plot_survivalplot.qual(d,
                                  sdata = data.frame(Y_test),
                                  BREAKTIME = BREAKTIME,
                                  cn_variables = type,
                                  name_data = NULL, title = title)[[type]]

    return(ggp)

  }else if(type=="COMP"){
    lst_test_lp <- NULL
    lst_ggp <- NULL

    #predict scores X_test
    test_score <- predict(model, newdata = X_test)
    test_score <- test_score[,names(model$survival_model$fit$coefficients),drop=F]
    for(cn in names(model$survival_model$fit$coefficients)){
      # check only coef in final model
      if(!cn %in% names(model$survival_model$fit$coefficients)){
        next
      }

      #get LP for individual components
      lst_test_lp[[cn]] <- test_score[,cn,drop=F] %*% model$survival_model$fit$coefficients[cn]
      colnames(lst_test_lp[[cn]]) <- cn

      if(is.na(cutoff[[cn]])){
        message(paste0("Cutoff not found for component: ", cn))
        next
      }

      txt_greater <- paste0("greater than ", cutoff[[cn]])
      txt_lower <- paste0("lesser/equal than ", cutoff[[cn]])

      LP <- ifelse(lst_test_lp[[cn]]>cutoff[[cn]], txt_greater, txt_lower)
      LP <- factor(LP)

      d <- as.data.frame(LP)
      colnames(d) <- cn

      lst_ggp[[cn]] <- plot_survivalplot.qual(d,
                                              sdata = data.frame(Y_test),
                                              BREAKTIME = BREAKTIME,
                                              cn_variables = cn,
                                              name_data = NULL, title = paste0(title," - ",cn))[[cn]]
    }

    return(lst_ggp)

  }else if(type=="VAR"){

    #As deleteIllegalChars() is performed in KM_VAR, run it always for VAR in TEST
    if(!attr(model, "model") %in% pkg.env$multiblock_methods){
      new_cn <- deleteIllegalChars(colnames(X_test))
      colnames(X_test) <- new_cn
    }else if(isa(X_test, "list")){
      for(b in names(X_test)){
        new_cn <- deleteIllegalChars(colnames(X_test[[b]]))
        colnames(X_test[[b]]) <- new_cn
      }
    }

    if(attr(model, "model") %in% c(pkg.env$sb.splsicox, pkg.env$sb.splsdrcox)){
      lst_ggp <- NULL
      ## SB.PLSICOX
      if(attr(model, "model") %in% c(pkg.env$sb.splsicox)){
        for(b in names(model$list_spls_models)){
          new_cutoff <- cutoff[endsWith(names(cutoff), paste0("_",b))]
          names(new_cutoff) <- unlist(lapply(names(new_cutoff), function(x){substr(x, start = 1, stop = nchar(x)-nchar(paste0("_",b)))}))
          lst_ggp[[b]] <- getTestKM(model = model$list_spls_models[[b]], X_test = X_test[[b]], Y_test, new_cutoff, type, ori_data, BREAKTIME, n.breaks, title)
        }
        return(lst_ggp)
      }else{
        ## SB.sPLSDRCOX
        for(b in names(model$list_spls_models)){
          new_cutoff <- cutoff[endsWith(names(cutoff), paste0("_",b))]
          names(new_cutoff) <- unlist(lapply(names(new_cutoff), function(x){substr(x, start = 1, stop = nchar(x)-nchar(paste0("_",b)))}))
          lst_ggp[[b]] <- getTestKM(model$list_spls_models[[b]], X_test[[b]], Y_test, new_cutoff, type, ori_data, BREAKTIME, n.breaks, title)
        }
        return(lst_ggp)
      }
    }else if(attr(model, "model") %in% c(pkg.env$mb.splsdrcox, pkg.env$mb.splsdacox) && isa(X_test, "list")){
      ## MBs.
      lst_ggp <- NULL
      for(b in names(model$mb.model$X)){
        new_cutoff <- cutoff[endsWith(names(cutoff), paste0("_",b))]
        names(new_cutoff) <- unlist(lapply(names(new_cutoff), function(x){substr(x, start = 1, stop = nchar(x)-nchar(paste0("_",b)))}))
        lst_ggp[[b]] <- getTestKM(model, X_test[[b]], Y_test, new_cutoff, type, ori_data, BREAKTIME, n.breaks, title)
      }
      return(lst_ggp)
    }

    X_test <- X_test[,names(cutoff),drop=F]
    lst_ggp <- NULL

    if(!ori_data){
      ori_names <- colnames(X_test)
      c <- F
      if(!all(is.null(model$X$x.mean))){
        c <- model$X$x.mean[ori_names]
      }
      s <- F
      if(!all(is.null(model$X$x.sd))){
        s <- model$X$x.sd[ori_names]
      }

      X_test <- scale(x = X_test, center = c, scale = s)
    }

    for(cn in colnames(X_test)){

        if(is.na(cutoff[[cn]])){
          message(paste0("Cutoff not found for variable: ", cn))
          next
        }

        txt_greater <- paste0("greater than ", cutoff[[cn]])
        txt_lower <- paste0("lesser/equal than ", cutoff[[cn]])

        LP <- ifelse(X_test[,cn]>cutoff[[cn]], txt_greater, txt_lower)
        LP <- factor(LP)

        d <- as.data.frame(LP)
        colnames(d) <- cn

        lst_ggp[[cn]] <- plot_survivalplot.qual(data = d,
                                                sdata = data.frame(Y_test),
                                                BREAKTIME = BREAKTIME,
                                                cn_variables = cn,
                                                name_data = NULL, title = title)[[cn]]
    }

    return(lst_ggp)

  }

}

#### ### ### ### ### ### ### ### #
# PREDICTION - MULTIPLE PATIENTS #
#### ### ### ### ### ### ### ### #

#' plot_LP.multiplePatients.list
#'
#' @description Run the function "plot_LP.multiplePatients" for a list of models. More information in "?plot_LP.multiplePatients".
#'
#' @param lst_models List of HDcox models.
#' @param new_data Numeric matrix or data.frame. New explanatory variables (raw data). Qualitative variables must be transform into binary variables.
#' @param error.bar Logical. Show error bar (default: FALSE).
#' @param onlySig Logical. Compute plot using only significant components (default: TRUE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables equal to 0 (default: TRUE).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_LP.multiplePatients.list(lst_models, new_data)
#' }

plot_LP.multiplePatients.list <- function(lst_models, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                          auto.limits = T, top = NULL){

  lst_plots <- purrr::map(lst_models, ~plot_LP.multiplePatients(model = ., new_data = new_data, error.bar = error.bar, onlySig = onlySig,
                                                                alpha = alpha, zero.rm = zero.rm,
                                                                auto.limits = auto.limits, top = top))

  return(lst_plots)
}

#' plot_LP.multiplePatients
#'
#' @param model HDcox model.
#' @param new_data Numeric matrix or data.frame. New explanatory variables (raw data). Qualitative variables must be transform into binary variables.
#' @param error.bar Logical. Show error bar (default: FALSE).
#' @param onlySig Logical. Compute plot using only significant components (default: TRUE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold (default: 0.05).
#' @param zero.rm Logical. Remove variables equal to 0 (default: TRUE).
#' @param auto.limits Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).
#' @param top Numeric. Show "top" first variables. If top = NULL, all variables are shown (default: NULL).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_LP.multiplePatients(model, new_data)
#' }

plot_LP.multiplePatients <- function(model, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                     auto.limits = T, top = NULL){

  if(!isa(model,pkg.env$model_class)){
    message("Model must be an object of class HDcox.")
    print(model)
    return(NULL)
  }

  if(attr(model, "model") %in% pkg.env$pls_methods){
    plot_cox.comparePatients(model = model,
                             new_data = new_data,
                             error.bar = error.bar,
                             onlySig = onlySig, alpha = alpha,
                             zero.rm = zero.rm, top = top,
                             auto.limits = auto.limits)
  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){
    plot_MB.cox.comparePatients(model = model,
                                new_data = new_data,
                                error.bar = error.bar,
                                onlySig = onlySig, alpha = alpha,
                                zero.rm = zero.rm, top = top,
                                auto.limits = auto.limits)
  }else{ #classical methods
    plot_classicalcox.comparePatients(model = model,
                                      new_data = new_data,
                                      error.bar = error.bar,
                                      onlySig = onlySig, alpha = alpha,
                                      zero.rm = zero.rm, top = top,
                                      auto.limits = auto.limits)
  }
}

plot_classicalcox.comparePatients <- function(model, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                              auto.limits = T, top = NULL){

  #DFCALLS
  value <- patients <- NULL

  coefficients <- model$survival_model$fit$coefficients
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

  lp.pats <- norm_patient[,deleteIllegalChars(names(model$survival_model$fit$coefficients))] %*% model$survival_model$fit$coefficients
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

plot_cox.comparePatients <- function(model, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                     auto.limits = T, top = NULL){

  #DFCALLS
  value <- patients <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                        alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)

  coefficients <- ggp.simulated_beta$beta

  if(all(coefficients==0)){
    message("No significant variables selected.")
    return(NULL)
  }

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

  lp.pats <- norm_patient[,rownames(ggp.simulated_beta$beta)] %*% ggp.simulated_beta$beta$value
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
  sd.min <- ggp.simulated_beta$sd.min[rownames(coefficients),]
  sd.max <- ggp.simulated_beta$sd.max[rownames(coefficients),]
  auto.limits <- NULL
  if(auto.limits.flag){
    if(!is.null(sd.min) & !is.null(sd.max)){
      auto.limits_min <- round2any(x = max(c(abs(coefficients$value-sd.min),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(x = max(c(abs(coefficients$value+sd.max),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
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

plot_MB.cox.comparePatients <- function(model, new_data, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                        auto.limits = T, top = NULL){

  #DFCALLS
  value <- patients <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                        alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)

  lst_coefficients <- ggp.simulated_beta$beta

  lst_plot <- list()
  lst_var.plot <- list()
  lst_lp.plot <- list()
  lst_lp <- list()
  lst_lp.var <- list()
  lst_norm_patients <- list()

  # blocks in ggp.simulated_beta$plot
  for(b in names(model$X$data)[names(model$X$data) %in% names(ggp.simulated_beta$plot)]){
    coefficients <- lst_coefficients[[b]][order(lst_coefficients[[b]]$value, decreasing = T),,drop=F]

    if(all(coefficients==0)){
      message("No significant variables selected.")
      next
    }

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
    if(!is.null(model$X$x.mean[[b]]) & !is.null(model$X$x.sd[[b]])){
      norm_patient <- scale(new_data[[b]][,names(model$X$x.mean[[b]])], center = model$X$x.mean[[b]], scale = model$X$x.sd[[b]])
    }else if(!is.null(model$X$x.mean[[b]])){
      norm_patient <- scale(new_data[[b]][,names(model$X$x.mean[[b]])], center = model$X$x.mean[[b]], scale = F)
    }else if(!is.null(model$X$x.sd[[b]])){
      norm_patient <- scale(new_data[[b]][,names(model$X$x.sd[[b]])], center = F, scale = model$X$x.sd[[b]])
    }else{
      norm_patient <- new_data[[b]]
    }

    #lp.new_pat_manual <- norm_patient[,rownames(coefficients)] %*% coefficients #predict lp
    lp.new_pat_variable <- apply(norm_patient[,deleteIllegalChars(rownames(coefficients)),drop=F], 1, function(x){
      x * coefficients$value #predict terms
    })

    #Compute LP without top variables
    #can be change for cox.prediction(model = model, new_data = patient, time = time, type = type, method = "cox")
    #for each patient on the data frame

    lp.pats <- norm_patient[,rownames(ggp.simulated_beta$beta[[b]])] %*% ggp.simulated_beta$beta[[b]]$value
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

    lp.new_pat_variable$patients <- factor(lp.new_pat_variable$patients, levels = rownames(new_data[[b]]))

    accuracy <- 0.1
    auto.limits.flag = T
    sd.min <- ggp.simulated_beta$sd.min[[b]][rownames(coefficients),]
    sd.max <- ggp.simulated_beta$sd.max[[b]][rownames(coefficients),]
    auto.limits <- NULL
    if(auto.limits.flag){
      if(!is.null(sd.min) & !is.null(sd.max)){
        auto.limits_min <- round2any(x = max(c(abs(coefficients$value-sd.min),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
        auto.limits_max <- round2any(x = max(c(abs(coefficients$value+sd.max),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
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

    lst_plot[[b]] <- pp
    lst_var.plot[[b]] <- res_all.plot
    lst_lp.plot[[b]] <- res_lp.plot
    lst_lp[[b]] <- lp.pats
    lst_lp.var[[b]] <- lp.new_pat_variable
    lst_norm_patients[[b]] <- norm_patient

  }

  return(list(plot = lst_plot, var.plot = lst_var.plot, lp.plot = lst_lp.plot, lp = lst_lp, lp.var = lst_lp.var, norm_patients = lst_norm_patients, patients = new_data))
}
