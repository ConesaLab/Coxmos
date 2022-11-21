#' save_ggplot !!!
#'
#' @param plot ggplot2
#' @param folder folder
#' @param name file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#'
#' @export
#'
save_ggplot <- function(plot, folder = NULL, name = NULL, wide = T, quality = "4K", dpi = 80, custom = NULL){
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

  if(class(plot[1]) == "ggsurvplot"){
    plot_surv = plot$plot
    if("table" %in% names(plot)){
      p2 = plot$table
      plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
    }
    ggsave(plot = plot_surv, filename = name, width = width, height = height, device='tiff', dpi=dpi)
  }else{
    ggsave(plot = plot, filename = paste0(folder,name), width = width, height = height, device='tiff', dpi=dpi)
  }
}

#' save_ggplot.svg !!!
#'
#' @param plot ggplot2
#' @param folder folder
#' @param name file name
#' @param wide 16:9 or 4:3
#' @param quality one of: "HD", "FHD", "2K", "4K", "8K"
#' @param dpi dpi
#' @param custom Custom size. Numeric vector of width and height
#'
#' @export
#'
save_ggplot.svg <- function(plot, folder = NULL, name = NULL, wide = T, quality = "4K", dpi = 80, custom = NULL){
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

  if(class(plot[1]) == "ggsurvplot"){
    plot_surv = plot$plot
    if("table" %in% names(plot)){
      p2 = plot$table
      plot_surv = cowplot::plot_grid(plot_surv,p2,align = "v",ncol =1,rel_heights = c(4,1))
    }
    ggsave(plot = plot_surv, filename = name, width = width, height = height, device='svg', dpi=dpi)
  }else{
    ggsave(plot = plot, filename = paste0(folder,name), width = width, height = height, device='svg', dpi=dpi)
  }
}

#' save_ggplot_lst !!!
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

      name <- paste0(folder,prefix,cn,suffix)
      if(!endsWith(name,".tiff")){
        name <- paste0(name, ".tiff")
      }

      if(is.null(object_name)){
        if(class(lst_plots[[cn]])[1] == "ggsurvplot"){
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
        if(class(lst_plots[[cn]][[object_name]])[1] == "ggsurvplot"){
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

      name <- paste0(folder,prefix,cn,suffix)
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

#' save_ggplot_lst.svg !!!
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

      name <- paste0(folder,prefix,cn,suffix)
      if(!endsWith(name,".svg")){
        name <- paste0(name, ".svg")
      }

      if(is.null(object_name)){

        if(class(lst_plots[[cn]])[1] == "ggsurvplot"){
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
        if(class(lst_plots[[cn]][[object_name]])[1] == "ggsurvplot"){
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

      name <- paste0(folder,prefix,cn,suffix)
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

#' Time consuming plot.
#'
#' @param lst_models A list of HDcox objects. Each HDCox object has the attribute time measured in minutes.
#' @param x.text X axis title.
#' @param y.text X axis title. (default: NULL == "Time (mins)").
#'
#' @return A ggplot2 bar plot for each object in lst_models and the total of all of them.
#' @export
#'
#' @details For time comparison between best models, use the cross validation objects instead the final models. The training has to be taken into account.
#'
#' @examples
#' \dontrun{
#'   lst_models = {"cox" = cox_model, "PLS-ICOX" = cv.plsicox_model, "sPLS-DRCOX" = cv.splsdrcox_model}
#'   plot_time.models(lst_models, x.text = "Method")
#' }

plot_time.models <- function(lst_models, x.text = "Method", y.text = NULL){

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
    if(class(lst_models[[m]])==pkg.env$model_class){
      lst_times[[m]] <- lst_models[[m]]$time
    }else if(class(lst_models[[m]][[1]])==pkg.env$model_class){
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

  lst_times$total <- total_time

  df.times <- do.call(rbind.data.frame, lst_times)
  colnames(df.times) <- "times"
  df.times$method <- names(lst_times)
  rownames(df.times) <- NULL

  max <- round2any(max(df.times$times), 10, f = ceiling)
  divisions <- 10

  if(((max - round(min(df.times$times))) / divisions)>2){
    dist <- round2any(((max - round(min(df.times$times))) / divisions), 10, f = ceiling)
  }else{
    dist <- round2any(((max - round(min(df.times$times))) / divisions), 1, f = ceiling)
  }

  accuracy <- dist * 0.1

  df.times$times <- round(df.times$times, digits = 4)
  x.var = "method"
  y.var = "times"
  x.color = "method"
  x.text = x.text
  if(is.null(y.text)){
    y.text = paste0("Time (",attr(lst_times[["total"]], "units"),")")
  }

  df.times$method <- factor(df.times$method, levels = df.times$method)

  ggp_time <- ggplot(df.times, aes_string(x = x.var, y = y.var, fill = x.color)) +
    geom_bar(stat="identity") +
    scale_y_continuous(breaks = seq(0, max, by = dist)) +
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
    best_df[!best_df[,x.color] == as.character(best_eta),c(y.var, y.var.sd)] <- NA
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
      ggp <- ggp + geom_point(data = best_df, aes_string(x = x.var, y = y.var, color = x.color, group = x.color), position=position_dodge(error_pos), size = dot_size, shape = 23, fill = "white", stroke = 2, show.legend = F)
    }else{
      ggp <- ggp + geom_point(data = best_df, aes_string(x = x.var, y = y.var), position=position_dodge(error_pos), size = dot_size, shape = 23, fill = "white", color = color_conesa, stroke = 2, show.legend = F)
    }
  }

  return(ggp)
}

lineplot.performace2.0 <- function(df, x.var = "time", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, point = T, mean = F, legend_rm = T){

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

  if(length(unique(df[,x.var,drop=T]))>30){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  # if(!is.null(y.limit)){
  #   ggp <- ggp + ylim(y.limit)
  # }

  if(mean){
    ggp <- ggp + geom_hline(data = mean_vector, aes_string(yintercept = mean_vector$mean_vector, color = x.color), size = 1)
  }

  if(legend_rm){
    ggp <- ggp + theme(legend.position = "none")
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

barplot.mean_performace2.0 <- function(df, x.var = "method", y.var="AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, hide_labels = T, legend_rm = NULL){

  #DFCALLS
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
  }

  if(length(unique(df[,x.var,drop=T]))>30){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

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

point.sd.mean_performace2.0 <- function(df, x.var = "method", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, pred.attr = "mean", hide_labels = T, legend_rm = NULL){

  #DFCALLS
  method <- NULL

  mean_vector = NULL
  sd_vector = NULL
  for(m in unique(df$method)){
    if(pred.attr %in% "mean"){
      mean_vector <- c(mean_vector, colMeans(df[df$method==m,y.var,drop=F], na.rm = T))
    }else if(pred.attr %in% "median"){
      mean_vector <- c(mean_vector, apply(df[df$method==m,y.var,drop=F], 2, function(x){median(x, na.rm = T)}))
    }
    sd_vector <- c(sd_vector, sd(df[df$method==m,y.var,drop=F]$AUC, na.rm = T))
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
  }

  if(length(unique(df[,x.var,drop=T]))>30){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

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

comboplot.performance2.0 <- function(df, x.var = "time", y.var = "AUC", x.color = "method", x.lab = NULL, y.lab = NULL, y.limit = NULL, pred.attr = "mean", point = T, mean = F, hide_labels = T){
  a <- lineplot.performace2.0(df = df, x.var = x.var, y.var = y.var, x.color = x.color, x.lab = x.lab, y.lab = y.lab, y.limit = y.limit, point = point, mean = F, legend_rm = T)
  b <- point.sd.mean_performace2.0(df = df, x.var = x.color, x.color = x.color, x.lab = NULL, y.lab = NULL, y.limit = y.limit, pred.attr = pred.attr, hide_labels = T, legend_rm = F)

  pp <- ggpubr::ggarrange(a, b, ncol = 2, widths = c(0.8, 0.2), align = "h")

  a <- lineplot.performace2.0(df, x.var, y.var, x.color, x.lab, y.lab, y.limit, point, mean = F, legend_rm = F)
  return(list(lineplot = a, lineplot.mean = pp))
}

#' plot_evaluation.list
#'
#' @param lst_eval_results List of eval_models4.0 results.
#' @param pred.attr "mean" or "median"
#' @param y.min Minimum Y value to plot. If NULL, automatic detection.
#' @param type type plot. Must be one of the following: "both", "line", "mean". In other case, "both" will be selected.
#'
#' @export

plot_evaluation.list <- function(lst_eval_results, pred.attr = "mean", y.min = NULL, type = "both"){

  lst_res <- purrr::map(lst_eval_results, ~plot_evaluation(eval_results = .,
                                                      pred.attr = pred.attr,
                                                      y.min = y.min, type = type))

  return(lst_res)

}

#' plot_evaluation
#'
#' @param eval_results Eval_models4.0 object
#' @param pred.attr "mean" or "median"
#' @param y.min Minimum Y value to plot. If NULL, automatic detection.
#' @param type type plot. Must be one of the following: "both", "line", "mean". In other case, "both" will be selected.
#'
#' @export

plot_evaluation <- function(eval_results, pred.attr = "mean", y.min = NULL, type = "both"){

  if(!pred.attr %in% c("mean", "median")){
    stop("pred.attr parameter must be one of: 'mean' or 'median'")
  }

  if(!type %in% c("both", "line", "mean")){
    type = "both"
  }

  #select minimum for all evals
  if(is.null(y.min)){
    y.min <- floor(min(eval_results$df$AUC, na.rm = T)*10)/10
  }

  lst_ggp <- list()

  lst_plots <- comboplot.performance2.0(df = eval_results$df,
                                        x.var = "time", y.var = "AUC", x.color = "method",
                                        y.limit = c(y.min, 1), pred.attr = pred.attr)
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
                                y.var = "AUC",
                                x.fill = "method",
                                x.alpha = NULL,
                                x.lab = "Method",
                                y.lab = "AUC",
                                fill.lab = NULL,
                                title = paste0("Method Performance"),
                                y.limit = NULL,
                                y.limit.exception = NULL,
                                jitter = F,
                                test = test_comparations,
                                show.median = T,
                                round.median = 3)

    if(lst_tests[[t]] == "NULL"){
      lst_plot_comparisons[["no_test"]] <- plot
    }else{
      lst_plot_comparisons[[lst_tests[[t]]]] <- plot
    }

  }

  return(list("lst_plots" = lst_ggp, "lst_plot_comparisons" = lst_plot_comparisons))
}

#' loadingplot.HDcox
#'
#' @param model HDcox model
#' @param zero.rm Remove variables equal to 0.
#' @param top Number. Show top variables.
#' @param auto.limits Logical. If TRUE, limits are detected for a better plot.
#'
#' @export

loadingplot.HDcox <- function(model, zero.rm = F, top = NULL, auto.limits = T){

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
#' @param model HDcox model
#' @param vector Vector of loading
#' @param zero.rm Remove variables equal to 0.
#' @param top Number. Show top variables.
#' @param auto.limits Logical. If TRUE, limits are detected for a better plot.
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
#' @param model HDcox model
#' @param zero.rm Remove variables equal to 0.
#' @param top Number. Show top variables.
#' @param auto.limits Logical. If TRUE, limits are detected for a better plot.
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

coxweightplot.fromVector.HDcox <- function(model, vector, sd.min = NULL, sd.max = NULL, zero.rm = F, top = NULL, auto.limits = T, block = NULL){

  #DFCALLS
  variables <- pp <- NULL

  loading_values <- vector
  ggp_loading <- NULL
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

    if(is.null(block)){
      ggp <- ggp + ggtitle(paste0(attr(model, "model"), " - Survival Weight"))
    }else{
      ggp <- ggp + ggtitle(paste0(attr(model, "model"), " - Survival Weight [", block, "]"))
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

    ggp_loading[[i]] = ggp
  }
  names(ggp_loading) <- colnames(loading_values)
  return(ggp_loading)
}

#' plot_pseudobeta.list
#'
#' @param lst_models List of HDcox model
#' @param error.bar Logical. Show error bar
#' @param onlySig Logical. Compute psudobetas using only significant components.
#' @param alpha Significant value (P.value).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0.
#' @param top Plot the top X variables with the higher pseudobetas in absolute value.
#' @param auto.limits Compute Y limit automatically
#'
#' @export

plot_pseudobeta.list <- function(lst_models, error.bar = T, onlySig = F, alpha = 0.05, zero.rm = F, top = NULL, auto.limits = T){

  lst_plots <- purrr::map(lst_models, ~plot_pseudobeta(model = .,
                                                       error.bar = error.bar,
                                                       onlySig = onlySig, alpha = alpha,
                                                       zero.rm = zero.rm, auto.limits = auto.limits, top = top))

  return(lst_plots)
}

#' plot_pseudobeta
#'
#' @param model HDcox model
#' @param error.bar Logical. Show error bar
#' @param onlySig Logical. Compute psudobetas using only significant components.
#' @param alpha Significant value (P.value).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0.
#' @param top Plot the top X variables with the higher pseudobetas in absolute value.
#' @param auto.limits Compute Y limit automatically
#'
#' @export

plot_pseudobeta <- function(model, error.bar = T, onlySig = F, alpha = 0.05, zero.rm = F, top = NULL, auto.limits = T){

  if(!attr(model, "model") %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)){
    stop("Model must be one of the follow models: 'PLS-ICOX', 'sPLS-DRCOX', 'sPLS-DRCOX-MixOmics', 'PLS-DACOX-MixOmics', 'SB.PLS-ICOX', 'SB.sPLS-DRCOX', 'MB.sPLS-DRCOX', 'MB.sPLS-DACOX'")
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
                                           sd.min = sd.min, sd.max = sd.max, auto.limits = T,
                                           zero.rm = zero.rm, top = top)[[1]]

  }else if(attr(model, "model") %in% pkg.env$multiblock_methods){

    if(onlySig){
      rn <- rownames(df.aux)[df.aux$`Pr(>|z|)` <= alpha]
      coefficients <- as.matrix(model$survival_model$fit$coefficients)[rn,,drop=F]
      sd <- df.aux[rn,"se(coef)",drop=F]
      W.star <- list()
      if(attr(model, "model") %in% pkg.env$sb.plsicox){
        for(b in names(model$list_pls_models)){
          W.star[[b]] <- model$list_pls_models[[b]]$X$W.star
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
      if(attr(model, "model") %in% pkg.env$sb.plsicox){
        for(b in names(model$list_pls_models)){
          W.star[[b]] <- model$list_pls_models[[b]]$X$W.star
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
                                               sd.min = sd.min[[b]], sd.max = sd.max[[b]], auto.limits = T,
                                               zero.rm = zero.rm, top = top, block = b)[[1]]
      }

    }

  }

  colnames(vector) <- "coef"

  #it will be a list for mb approaches
  return(list(plot = plot,
              beta = vector,
              sd.min = sd.min,
              sd.max = sd.max))
}

#' plot_pseudobeta.newPatient.list
#'
#' @param lst_models HDcox model
#' @param new_pat Row of new patients
#' @param error.bar Logical. Show error bar
#' @param onlySig Logical. Compute psudobetas using only significant components.
#' @param alpha Significant value (P.value).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0.
#' @param top Plot the top X variables with the higher pseudobetas in absolute value.
#' @param auto.limits Compute Y limit automatically
#' @param show.betas Show original betas
#'
#' @export

plot_pseudobeta.newPatient.list <- function(lst_models, new_pat, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                            top = NULL, auto.limits = T, show.betas = F){

  lst_plots <- purrr::map(lst_models, ~plot_pseudobeta.newPatient(model = .,
                                                                  new_pat = new_pat,
                                                                  error.bar = error.bar,
                                                                  onlySig = onlySig, alpha = alpha,
                                                                  zero.rm = zero.rm, top = top,
                                                                  auto.limits = auto.limits, show.betas = show.betas))

  return(lst_plots)

}

#' plot_pseudobeta.newPatient
#'
#' @param model HDcox model
#' @param new_pat Row of new patients
#' @param error.bar Logical. Show error bar
#' @param onlySig Logical. Compute psudobetas using only significant components.
#' @param alpha Significant value (P.value).
#' @param zero.rm Logical. Remove variables with a pseudobeta equal to 0.
#' @param top Plot the top X variables with the higher pseudobetas in absolute value.
#' @param auto.limits Compute Y limit automatically
#' @param show.betas Show original betas
#'
#' @export

plot_pseudobeta.newPatient <- function(model, new_pat, error.bar = T, onlySig = T, alpha = 0.05, zero.rm = T,
                                       top = NULL, auto.limits = T, show.betas = F){

  #DFCALLS
  lp <- lp.min <- lp.max <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                              alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)
  coefficients <- ggp.simulated_beta$beta

  coeff.min <- NULL
  coeff.max <- NULL
  if(error.bar){
    coeff.min <- ggp.simulated_beta$sd.min
    coeff.max <- ggp.simulated_beta$sd.max
  }

  #norm patient
  new_pat <- new_pat[,names(model$X$x.mean),drop=F]

  if(!is.null(model$X$x.mean) & !is.null(model$X$x.sd)){
    norm_patient <- scale(new_pat, center = model$X$x.mean, scale = model$X$x.sd)
  }else if(!is.null(model$X$x.mean)){
    norm_patient <- scale(new_pat, center = model$X$x.mean, scale = F)
  }else if(!is.null(model$X$x.sd)){
    norm_patient <- scale(new_pat, center = F, scale = model$X$x.sd)
  }else{
    norm_patient <- new_pat
  }

  #lp.new_pat_manual <- norm_patient[,rownames(coefficients)] %*% coefficients #predict lp
  lp.new_pat_variable <- norm_patient[,rownames(coefficients)] * coefficients #predict terms

  lp.new_pat_variable.min <- NULL
  lp.new_pat_variable.max <- NULL
  if(error.bar){
    lp.new_pat_variable.min <- norm_patient[,rownames(coeff.min)] * coeff.min
    lp.new_pat_variable.max <- norm_patient[,rownames(coeff.max)] * coeff.max
  }

  #filter pat_variables using psudobeta plot (top could be applied)
  lp.new_pat_variable <- lp.new_pat_variable[rownames(ggp.simulated_beta$plot$data),,drop=F]
  lp.new_pat_variable.min <- lp.new_pat_variable.min[rownames(ggp.simulated_beta$plot$data),,drop=F]
  lp.new_pat_variable.max <- lp.new_pat_variable.max[rownames(ggp.simulated_beta$plot$data),,drop=F]

  coefficients <- coefficients[rownames(lp.new_pat_variable),,drop=F]

  #terms
  # df <- as.data.frame(cbind(cbind(ggp.simulated_beta$beta,
  #                                 rep("Beta",nrow(ggp.simulated_beta$beta))),
  #                           rownames(ggp.simulated_beta$beta)))
  # colnames(df) <- c("beta", "type", "var")
  #
  # df$beta <- as.numeric(df$beta)
  # df <- df[order(df$beta, decreasing = T),]
  #
  # df.pat <- cbind(cbind(lp.new_pat_variable,  rep("Patient Linear Predictor", nrow(lp.new_pat_variable))), rownames(lp.new_pat_variable))
  # colnames(df.pat) <- c("beta", "type", "var")
  # df <- rbind(df, df.pat)
  #
  # df$beta <- as.numeric(df$beta)
  # df$var <- factor(df$var, levels = unique(df$var))
  # df$type <- factor(df$type, levels = unique(df$type))

  #terms
  if(error.bar){
    df.pat <- data.frame("lp" = lp.new_pat_variable[,1],
                         "lp.min" = lp.new_pat_variable.min[,1],
                         "lp.max" = lp.new_pat_variable.max[,1],
                         "var" = rownames(lp.new_pat_variable))
  }else{
    df.pat <- data.frame("lp" = lp.new_pat_variable[,1],
                         "lp.min" = 0,
                         "lp.max" = 0,
                         "var" = rownames(lp.new_pat_variable))
  }

  df.pat$lp <- as.numeric(df.pat$lp)
  df.pat$lp.min <- as.numeric(df.pat$lp.min)
  df.pat$lp.max <- as.numeric(df.pat$lp.max)
  df.pat$var <- factor(df.pat$var, levels = unique(df.pat$var))

  accuracy <- 0.1

  if(show.betas){
    if(error.bar){
      val_min <- as.numeric(max(abs(coeff.min), abs(df.pat$lp.min)))
      val_max <- as.numeric(max(abs(coeff.max), abs(df.pat$lp.max)))
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
    ggp <- ggp + RColorConesa::scale_color_conesa(palette = "warm", continuous = T)
  }

  ggp <- ggp + guides(color= "none")
  ggp <- ggp + ylab(label = "Linear Predictor")
  ggp <- ggp + xlab(label = "Variables")
  ggp <- ggp + ggtitle(label = paste0("Observation - ", rownames(new_pat)))

  if(length(unique(df.pat$var))>15){
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  if(show.betas){

    ggp.aux <- ggp + scale_y_continuous(n.breaks = 10, limits = c(-1*auto.limits, auto.limits))

    ggp.aux2 <- ggp.simulated_beta$plot
    ggp.aux2 <- ggp.aux2 + guides(fill = "none")

    sign.beta <- coefficients>0
    sign.pat <- df.pat$lp>0
    same.sign <- sign.beta == sign.pat
    same.sign <- same.sign[rownames(ggp.simulated_beta$plot$data),,drop=F]

    ggp.aux$mapping$fill[[2]] <- same.sign
    ggp.aux <- ggp.aux + guides(fill = guide_legend(title="Same beta direction:")) + theme(legend.position="left")

    if(requireNamespace("RColorConesa", quietly = TRUE)){
      ggp.aux <- ggp.aux + RColorConesa::scale_fill_conesa(reverse = T)
    }

    ggp <- ggpubr::ggarrange(ggp.aux, ggp.aux2, ncol = 2, widths = c(0.5, 0.5), align = "h")
  }

  return(list(plot = ggp, lp.var = lp.new_pat_variable, norm_pat = norm_patient, pat = new_pat))

}

#' plot_cox.comparePatients.list
#'
#' @param lst_models List HDcox models
#' @param df.pat Dataframe of observations
#' @param error.bar Show error bar.
#' @param onlySig Show only significant variables.
#' @param alpha Significant value.
#' @param zero.rm Remove variables equal to 0.
#' @param auto.limits Logical. If TRUE, limits are detected for a better plot.
#' @param top Number. Show top variables.
#'
#' @export
#'
plot_cox.comparePatients.list <- function(lst_models, df.pat, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                auto.limits = T, top = NULL){

  lst_plots <- purrr::map(lst_models, ~plot_cox.comparePatients(model = ., df.pat = df.pat, error.bar = error.bar, onlySig = onlySig,
                                                           alpha = alpha, zero.rm = zero.rm,
                                                           auto.limits = auto.limits, top = top))

  return(lst_plots)
}

#' plot_cox.comparePatients
#'
#' @param model HDcox model
#' @param df.pat Dataframe of observations
#' @param error.bar Show error bar.
#' @param onlySig Show only significant variables.
#' @param alpha Significant value.
#' @param zero.rm Remove variables equal to 0.
#' @param auto.limits Logical. If TRUE, limits are detected for a better plot.
#' @param top Number. Show top variables.
#'
#' @export

plot_cox.comparePatients <- function(model, df.pat, error.bar = F, onlySig = T, alpha = 0.05, zero.rm = T,
                                auto.limits = T, top = NULL){

  #DFCALLS
  value <- patients <- NULL

  #plot
  ggp.simulated_beta <- plot_pseudobeta(model = model, error.bar = error.bar, onlySig = onlySig,
                                              alpha = alpha, zero.rm = zero.rm, auto.limits = auto.limits, top = top)

  coefficients <- ggp.simulated_beta$beta
  coefficients <- coefficients[order(coefficients, decreasing = T),,drop=F]

  if(!is.null(top)){
    if(top < nrow(coefficients)){
      aux_df <- coefficients
      aux_df[,"coef"] <- abs(aux_df[,"coef",drop=F])
      aux_df <- aux_df[order(aux_df[,"coef",drop=F], decreasing = T),,drop=F]
      aux_df <- aux_df[1:top,,drop=F]
      coefficients <- coefficients[rownames(coefficients) %in% rownames(aux_df),,drop=F]
    }
  }

  #norm patient
  if(!is.null(model$X$x.mean) & !is.null(model$X$x.sd)){
    norm_patient <- scale(df.pat, center = model$X$x.mean, scale = model$X$x.sd)
  }else if(!is.null(model$X$x.mean)){
    norm_patient <- scale(df.pat, center = model$X$x.mean, scale = F)
  }else if(!is.null(model$X$x.sd)){
    norm_patient <- scale(df.pat, center = F, scale = model$X$x.sd)
  }else{
    norm_patient <- df.pat
  }

  #lp.new_pat_manual <- norm_patient[,rownames(coefficients)] %*% coefficients #predict lp
  lp.new_pat_variable <- apply(norm_patient[,rownames(coefficients),drop=F], 1, function(x){
    x * coefficients #predict terms
  })

  #Compute LP without top variables
  #can be change for getLPforNewPatient(model = model, new_pat = patient, time = time, type = type, method = "cox")
  #for each patient on the data frame

  lp.pats <- norm_patient[,rownames(ggp.simulated_beta$beta)] %*% ggp.simulated_beta$beta
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

  lp.new_pat_variable$patients <- factor(lp.new_pat_variable$patients, levels = rownames(df.pat))

  accuracy <- 0.1
  auto.limits.flag = T
  sd.min <- ggp.simulated_beta$sd.min[rownames(coefficients),]
  sd.max <- ggp.simulated_beta$sd.max[rownames(coefficients),]
  auto.limits <- NULL
  if(auto.limits.flag){
    if(!is.null(sd.min) & !is.null(sd.max)){
      auto.limits_min <- round2any(x = max(c(abs(coefficients-sd.min),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
      auto.limits_max <- round2any(x = max(c(abs(coefficients+sd.max),abs(lp.new_pat_variable[lp.new_pat_variable$lp.flag==T,]$value))), accuracy = accuracy, f = ceiling)
      auto.limits <- max(auto.limits_min, auto.limits_max)
    }else{
      auto.limits <- round2any(max(abs(lp.new_pat_variable$value)), accuracy = accuracy, f = ceiling)
    }
  }else{
    auto.limits <- round2any(max(c(abs(sd.max), abs(sd.min), abs(lp.new_pat_variable$value))), accuracy = accuracy, f = ceiling)
  }

  ggp <- ggplot(lp.new_pat_variable[lp.new_pat_variable$lp.flag==F,], aes(x = var, y = value, fill = patients)) +
    geom_bar(stat = "identity", position = "dodge")
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
  res_lp.plot <- ggp2

  ggp <- ggp + guides(fill = "none")
  ggp2 <- ggp2 + ylab(label = "") + xlab(label = "")

  pp <- ggpubr::ggarrange(ggp, ggp2, ncol = 2, widths = c(0.8, 0.2), align = "h")

  return(list(plot = pp, var.plot = res_all.plot, lp.plot = res_lp.plot, lp = lp.pats, lp.var = lp.new_pat_variable, norm_patients = norm_patient, patients = df.pat))
}

#' patient.eventDensity
#'
#' @param patient observation
#' @param time time point of study
#' @param model HDcox model
#' @param type type between "lp", "risk, "expected", "survival"
#' @param size size
#' @param color color
#'
#' @export

plot_patient.eventDensity <- function(patient, time = NULL, model, type = "lp", size = 3, color = "red"){

  #DFCALLS
  x <- y <- event <- NULL

  pred.value <- getLPforNewPatient(model = model, new_pat = patient, time = time, type = type, method = "cox")

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
#' @param patient observation
#' @param time time point of study
#' @param model HDcox model
#' @param type type between "lp", "risk, "expected", "survival"
#' @param size size
#' @param color color
#'
#' @export

plot_patient.eventHistogram <- function(patient, time = NULL, model, type = "lp", size = 3, color = "red"){

  #DFCALLS
  x <- y <- NULL

  pred.value <- getLPforNewPatient(model = model, new_pat = patient, time = time, type = type, method = "cox")

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

#' plot_divergent.biplot
#' @description Two side plot by a qualitative and a quantitative variable and Y event matrix.
#'
#' @param X A list of HDcox objects. Each HDCox object has the attribute time measured in minutes.
#' @param Y X axis title.
#' @param NAMEVAR1 factor variable name (must be located in colnames(X)).
#' @param NAMEVAR2 numerical variable name (must be located in colnames(X)).
#' @param breaks Number of values for X numerical variable per each Y axis value. Example: breaks = 5 -> Y axis values: "0-4, 5-9, 10-14..."
#' @param x.text Text for X axis.
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
#'   plot_divergent.biplot(X, Y, NAMEVAR1, NAMEVAR2, breaks = 5, x.text = "N. of Patients")
#' }
plot_divergent.biplot <- function(X, Y, NAMEVAR1, NAMEVAR2, breaks, x.text = "N. of Samples"){
  df<-NULL

  VAR1 <- X[rownames(X), NAMEVAR1] #will be a factor
  VAR2 <- X[rownames(X), NAMEVAR2] #must be numerical

  OUTCOME <- Y[rownames(X),"event"]

  df <- data.frame(VAR1, VAR2, OUTCOME)  # merge by row names (by=0 or by="row.names")
  colnames(df) <- c(NAMEVAR1, NAMEVAR2, "event")

  #add age as category
  df.index <- NULL
  cat <- NULL
  index <- NULL

  breaks = breaks
  min <- round2any(min(VAR2), accuracy = breaks, f = floor)
  max <- round2any(max(VAR2), accuracy = breaks, f = ceiling)
  for(i in seq(min,max,breaks)){
    if(i!=max){
      new <- which(df[,NAMEVAR2]>=i & df[,NAMEVAR2]<=(i+breaks-1))
      index <- c(index, new)
      cat <- c(cat, rep(paste0(i,"-",i+breaks-1), length(new)))
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

#' plot_events
#'
#' @param Y Y matrix
#' @param roundTo Round time to which value.
#' @param categories Categories to print.
#' @param max.breaks Number maximum of breaks in X axis.
#'
#' @export

plot_events <- function(Y, roundTo = 0.25, categories = c("Censored","Death"), max.breaks = 20){

  #DFCALLS
  Category <- Time <- Values <- x.names <- breaks<- NULL

  if(class(Y[,"event"])!="logical"){
    message("Warning: Y matrix must has event column as TRUE, FALSE. as.logical() function has been used.")
    Y$event <- as.logical(Y$event)
  }

  time_aux <- round2any(as.numeric(Y$time), roundTo)

  while(length(unique(time_aux))>max.breaks){
    roundTo = roundTo + 0.05
    time_aux <- round2any(as.numeric(Y$time), roundTo)
  }

  Y$time <- round2any(as.numeric(Y$time), roundTo)
  breaks = seq(min(Y$time), max(Y$time), by=roundTo)
  breaks = round2any(breaks, 0.1)
  breaks = c(breaks,max(breaks)+roundTo)
  x.names <- cut(x = Y$time, breaks = breaks, include.lowest = T)

  Y$time_g <- x.names

  vt=NULL
  vcategory=NULL
  vvalues=NULL
  for(t in levels(x.names)){
    vt <- c(vt, t, t)
    vcategory <- c(vcategory, categories)
    vvalues<- c(vvalues, sum(Y[Y$time_g==t, "event"]==F), sum(Y[Y$time_g==t, "event"]==T))
  }

  dd <- data.frame(Time=vt, Category=vcategory, Values=vvalues)
  dd$Time <- factor(dd$Time, levels = levels(x.names))

  ggp_density <- ggplot(dd, aes(fill=Category, x=Time, y=Values)) +
    #geom_bar(position="stack", stat="identity") +
    geom_bar(stat = "identity") +
    ylab("Number of patients") +
    scale_y_continuous(n.breaks = 10) +
    guides(fill=guide_legend(title="Event type"), color = "none")

  if(requireNamespace("RColorConesa", quietly = TRUE)){
    ggp_density <- ggp_density + RColorConesa::scale_fill_conesa()
  }

  return(list(plot = ggp_density, df = dd))
}

#' plot_HDcox.PLS.model
#'
#' @param model HDcox model
#' @param comp vector of two components
#' @param mode scores, loadings o biplot
#' @param factor factor to color
#' @param legend.title Legend title
#' @param mahalanovis_limit Mahalanobis distance limit
#' @param radius Radious to plot variable names
#' @param names Logical. Show names for variables outside the radius.
#' @param allNames Logical. Show all variable names.
#' @param colorReverse Logical. Color reverse.
#' @param text.size Text size.
#'
#' @export

plot_HDcox.PLS.model <- function(model, comp = c(1,2), mode = "scores", factor = NULL, legend.title = NULL, mahalanovis_limit = 20, radius = 0.2, names = F, allNames = F, colorReverse = F, text.size = 4){

  ggp = NULL
  aux.model = model

  modes <- c("scores", "loadings", "biplot")
  if(!mode %in% modes){
    stop_quietly(paste0("mode must be one of the following: ", paste0(modes, collapse = ", ")))
  }

  if(!is.null(factor)){
    if(class(factor)!="factor" & mode %in% c("scores", "biplot")){
      stop_quietly("Factor must be a factor object.")
    }
  }

  if(!class(aux.model)==pkg.env$model_class){
    stop_quietly("aux.model must be a HDcox object.")
  }else if(!attr(aux.model, "model") %in% pkg.env$pls_methods){
    stop_quietly("aux.model must be a HDcox object pls class ('PLS-ICOX','sPLS-DRCOX','sPLS-DRCOX-MixOmics' or 'PLS-DACOX-MixOmics').")
  }

  if(mode=="scores"){
    if(ncol(aux.model$X$scores)==1){
      message("The model has only 1 component")

      aux.model$X$scores <- cbind(aux.model$X$scores, aux.model$X$scores)
      colnames(aux.model$X$scores) <- c("p1", "p2")
      df <- as.data.frame(aux.model$X$scores)

      comp = c(1, 1)
      subdata <- df[0,1:2]

      f <- as.factor(factor)

      txt.expression <- paste0("Scores (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(df) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor)) +
        # ggtitle(label = bquote("Scores (aux.model) - "~R^2 == .(max(aux.model@aux.modelDF[,2])))) +
        # xlab(label = paste0("p",as.character(comp[1]), " (", as.character(aux.model@aux.modelDF$R2X[comp[1]]*100), " %)")) +
        # ylab(label = paste0("p",as.character(comp[2]), " (", as.character(aux.model@aux.modelDF$R2X[comp[2]]*100), " %)")) +
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)")) +
        labs(color = legend.title) + theme(legend.position="bottom")

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse))
      }

      # if(allNames){
      #   ggp <- ggp + ggrepel::geom_text_repel(data = as.data.frame(df),
      #                                         aes(x = as.data.frame(df)[,comp[1]],
      #                                             y = as.data.frame(df)[,comp[2]]),
      #                                         label = rownames(as.data.frame(df)), size=text.size)
      # }else if(names){
      #   ggp <- ggp + ggrepel::geom_text_repel(data = subdata, aes(x = subdata[,comp[1]],
      #                                                             y = subdata[,comp[2]]),
      #                                         label = rownames(subdata), size=text.size)
      # }
      #
      # if(!is.null(radius) & nrow(subdata)!=0){
      #   ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
      # }

    }else{
      df <- as.data.frame(aux.model$X$scores)
      mh <- mahalanobis(x = aux.model$X$scores[,comp], center = F, cov = cov(aux.model$X$scores[,comp]))
      subdata <- as.data.frame(aux.model$X$scores)[names(mh)[mh>mahalanovis_limit],]
      for(i in 1:nrow(df)){
        if(rownames(df[i,comp,drop=F]) %in% rownames(subdata)){
          next
        }else{
          #check mahalanovis_limit
          aux <- df[i,comp,drop=F]
          y_value <- sqrt(mahalanovis_limit^2-aux[,1]^2)
          if(y_value < abs(aux[,2])){
            subdata <- rbind(subdata, df[i,,drop=F])
          }
        }
      }

      txt.expression <- paste0("Scores (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(df) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor)) +
        stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F) +
        coord_fixed(ratio=1) +
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)")) +
        labs(color = legend.title) + theme(legend.position="bottom")

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse)) +
          scale_fill_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse))
      }

      if(names){
        ggp <- ggp + ggrepel::geom_text_repel(data = subdata,
                                              aes(x = subdata[,comp[1]],
                                                  y = subdata[,comp[2]]),
                                              label = rownames(subdata), size=text.size)
      }else if(allNames){
        ggp <- ggp + ggrepel::geom_text_repel(data = as.data.frame(df),
                                              aes(x = as.data.frame(df)[,comp[1]],
                                                  y = as.data.frame(df)[,comp[2]]),
                                              label = rownames(as.data.frame(df)),
                                              size=text.size)
      }
    }
  }

  if(mode=="loadings"){

    if(ncol(aux.model$X$loadings)==1){
      message("The model has only 1 component")

      aux.model$X$loadings <- cbind(aux.model$X$loadings,aux.model$X$loadings)
      colnames(aux.model$X$loadings) <- c("p1", "p2")
      df <- as.data.frame(aux.model$X$loadings)
      #mh <- mahalanobis(x = aux.model$X$loadings[,comp], center = F, cov = cov(aux.model$X$loadings[,comp]))
      subdata <- df[apply(df[,comp],1,function(x){any(abs(x)>radius)}),]
      for(i in 1:nrow(df)){
        if(rownames(df[i,comp,drop=F]) %in% rownames(subdata)){
          next
        }else{
          #check radius
          aux <- df[i,comp,drop=F]
          y_value <- sqrt(radius^2-aux[,1]^2)
          if(y_value < abs(aux[,2])){
            subdata <- rbind(subdata, df[i,,drop=F])
          }
        }
      }

      comp = c(1, 1)

      txt.expression <- paste0("Loadings (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(df) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]])) +
        #ggtitle(label = bquote("Loadings (aux.model) - "~R^2 == .(max(aux.model@aux.modelDF[,2])))) +
        #xlab(label = paste0("p",as.character(comp[1]), " (", as.character(aux.model@aux.modelDF$R2X[comp[1]]*100), " %)")) +
        #ylab(label = paste0("p",as.character(comp[2]), " (", as.character(aux.model@aux.modelDF$R2X[comp[2]]*100), " %)"))
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)"))

      if(allNames){
        ggp <- ggp + ggrepel::geom_text_repel(data = as.data.frame(df),
                                              aes(x = as.data.frame(df)[,comp[1]],
                                                  y = as.data.frame(df)[,comp[2]]),
                                              label = rownames(as.data.frame(df)), size=text.size)
      }else if(names){
        ggp <- ggp + ggrepel::geom_text_repel(data = subdata, aes(x = subdata[,comp[1]],
                                                                  y = subdata[,comp[2]]),
                                              label = rownames(subdata), size=text.size)
      }

      if(!is.null(radius) & nrow(subdata)!=0){

        if(requireNamespace("ggforce", quietly = TRUE)){
          ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
        }

      }

    }else{
      df <- as.data.frame(aux.model$X$loadings)
      subdata <- df[apply(df[,comp],1,function(x){any(abs(x)>radius)}),] #works as square, but not as a circle
      for(i in 1:nrow(df)){
        if(rownames(df[i,comp,drop=F]) %in% rownames(subdata)){
          next
        }else{
          #check radius
          aux <- df[i,comp,drop=F]
          y_value <- sqrt(radius^2-aux[,1]^2)
          if(y_value < abs(aux[,2])){
            subdata <- rbind(subdata, df[i,,drop=F])
          }
        }
      }

      txt.expression <- paste0("Loadings (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(as.data.frame(df)) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]])) +
        coord_fixed(ratio=1) +
        # ggtitle(label = bquote("Loadings (aux.model) - "~R^2 == .(max(aux.model@aux.modelDF[,2])))) +
        # xlab(label = paste0("p",as.character(comp[1]), " (", as.character(aux.model@aux.modelDF$R2X[comp[1]]*100), " %)")) +
        # ylab(label = paste0("p",as.character(comp[2]), " (", as.character(aux.model@aux.modelDF$R2X[comp[2]]*100), " %)"))
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)"))
    }

    if(allNames){
      ggp <- ggp + ggrepel::geom_text_repel(data = as.data.frame(df),
                                            aes(x = as.data.frame(df)[,comp[1]],
                                                y = as.data.frame(df)[,comp[2]]),
                                            label = rownames(as.data.frame(df)), size=text.size)
    }else if(names){
      ggp <- ggp + ggrepel::geom_text_repel(data = subdata, aes(x = subdata[,comp[1]],
                                                                y = subdata[,comp[2]]),
                                            label = rownames(subdata), size=text.size)
    }

    if(!is.null(radius) & nrow(subdata)!=0){
      if(requireNamespace("ggforce", quietly = TRUE)){
        ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius))
      }
    }

  }

  if(mode=="biplot"){
    if(ncol(aux.model$X$loadings)==1){
      message("The model has only 1 component")
      df <- cbind(aux.model$X$scores, aux.model$X$scores)
      colnames(df) <- c("p1", "p2")
      df <- as.data.frame(df)

      comp = c(1, 1)

      txt.expression <- paste0("Biplot (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(df) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor)) +
        stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F) +
        coord_fixed(ratio=1) +
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)")) +
        labs(color = legend.title) + theme(legend.position="bottom")

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse)) +
          scale_fill_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse))
      }

      df_loading <- as.data.frame(cbind(aux.model$X$loadings, aux.model$X$loadings))
      max.loadings <- apply(abs(df_loading), 2, max)
      max.scores <- apply(abs(df), 2, max)

      ratio <- max.scores / max.loadings

      df_loading <- as.data.frame(t(apply(df_loading, 1, function(x){x * ratio})))
      subdata_loading <- df_loading[apply(df_loading[,comp],1,function(x){any(abs(x) > (radius * max(ratio[comp[1]])))}),]
      for(i in 1:nrow(df_loading)){
        if(rownames(df_loading[i,comp,drop=F]) %in% rownames(subdata_loading)){
          next
        }else{
          #check radius
          aux <- df_loading[i,comp,drop=F]
          y_value <- sqrt((radius * max(ratio[comp[1]]))^2-aux[,1]^2)
          if(y_value < abs(aux[,2])){
            subdata_loading <- rbind(subdata_loading, df_loading[i,,drop=F])
          }
        }
      }

      df.2 <- cbind(df_loading[,comp[1]], df_loading[,comp[2]])
      colnames(df.2) <- c("l1", "l2")
      df.2 <- as.data.frame(df.2)
      ggp <- ggp + geom_point(data = df.2, aes(x = df.2[,comp[1]], y = df.2[,comp[2]]), shape = 15)

    }else{
      df <- as.data.frame(aux.model$X$scores)
      mh <- mahalanobis(x = aux.model$X$scores[,comp], center = F, cov = cov(aux.model$X$scores[,comp]))
      subdata <- as.data.frame(aux.model$X$scores)[names(mh)[mh>mahalanovis_limit],]
      for(i in 1:nrow(df)){
        if(rownames(df[i,comp,drop=F]) %in% rownames(subdata)){
          next
        }else{
          #check radius
          aux <- df[i,comp,drop=F]
          y_value <- sqrt(mahalanovis_limit^2-aux[,1]^2)
          if(y_value < abs(aux[,2])){
            subdata <- rbind(subdata, df[i,,drop=F])
          }
        }
      }

      txt.expression <- paste0("Biplot (",attr(aux.model, "model"),") - ")

      ggp <- ggplot(as.data.frame(df)) +
        geom_point(aes(x = df[,comp[1]], y = df[,comp[2]], color = factor)) +
        stat_ellipse(aes(x = df[,comp[1]], y = df[,comp[2]], fill = factor), geom = "polygon", alpha = 0.1, show.legend=F) +
        coord_fixed(ratio=1) +
        ggtitle(label = bquote(.(txt.expression) ~R^2 == "not R2 yet")) +
        xlab(label = paste0("comp_",as.character(comp[1]), " (", as.character("NA"), " %)")) +
        ylab(label = paste0("comp_",as.character(comp[2]), " (", as.character("NA"), " %)")) +
        labs(color = legend.title) + theme(legend.position="bottom")

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        ggp <- ggp + scale_color_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse)) +
          scale_fill_manual(values = RColorConesa::colorConesa(length(unique(factor)), reverse = colorReverse))
      }

      df_loading <- as.data.frame(aux.model$X$loadings)
      max.loadings <- apply(abs(df_loading), 2, max)
      max.scores <- apply(abs(df), 2, max)

      ratio <- max.scores / max.loadings

      df_loading <- as.data.frame(t(apply(df_loading, 1, function(x){x * ratio})))

      subdata_loading <- df_loading[apply(df_loading[,comp],1,function(x){any(abs(x) > (radius * max(ratio[comp[1]])))}),]
      for(i in 1:nrow(df)){
        if(rownames(df[i,comp,drop=F]) %in% rownames(subdata_loading)){
          next
        }else{
          #check radius
          aux <- df[i,comp,drop=F]
          y_value <- sqrt(abs((radius * max(ratio[comp[1]]))^2-aux[,1]^2))
          if(y_value < abs(aux[,2])){
            subdata_loading <- rbind(subdata_loading, df[i,,drop=F])
          }
        }
      }

      ggp <- ggp +
        geom_point(data = df_loading, aes(x = df_loading[,comp[1]], y = df_loading[,comp[2]]), shape = 15)

      if(allNames){
        ggp <- ggp + ggrepel::geom_text_repel(data = as.data.frame(df_loading), aes(x = as.data.frame(df_loading)[,comp[1]], y = as.data.frame(df_loading)[,comp[2]]), label = rownames(as.data.frame(df_loading)), size=text.size)
      }else if(names){
        ggp <- ggp + ggrepel::geom_text_repel(data = subdata_loading, aes(x = subdata_loading[,comp[1]], y = subdata_loading[,comp[2]]), label = rownames(subdata_loading), size=text.size)
      }

      if(!is.null(radius) & nrow(subdata)!=0){
        if(requireNamespace("ggforce", quietly = TRUE)){
          ggp <- ggp + ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius * max(ratio[comp[1]])))
        }
      }
    }

  }

  return(list(plot = ggp, outliers = rownames(subdata)))
}

#' plot_cox.event.list
#'
#' @param lst_models List of HDcox models
#' @param type "lp", "risk", "expected" or "survival" methodology
#' @param h.breaks Number of time breaks.
#'
#' @export

plot_cox.event.list <- function(lst_models, type = "lp", h.breaks = 30){

  ggp_list <- purrr::map(lst_models, ~plot_cox.event(model = ., type = type, h.breaks = h.breaks))

  return(ggp_list)
}

#' plot_cox.event
#'
#' @param model HDcox model
#' @param type "lp", "risk", "expected" or "survival" methodology
#' @param h.breaks Number of time breaks.
#'
#' @export

plot_cox.event <- function(model, type = "lp", h.breaks = 30){

  #DFCALLS
  event <- NULL

  #exits
  if(all(is.null(model$survival_model$fit)) || all(is.na(model$survival_model$fit))){
    message(paste0("Survival model not found for ", attr(model, "model"), "."))
    return(NA)
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

  binwidth <- (max(df_hr[,1]) - min(df_hr[,1])) / h.breaks
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

#' plot_proportionalHazard.list
#'
#' @param lst_models List of HDcox models
#'
#' @export

plot_proportionalHazard.list <- function(lst_models){

  lst_plots <- purrr::map(lst_models, ~plot_proportionalHazard(model = .))

  return(lst_plots)
}

#' plot_proportionalHazard
#'
#' @param model HDcox model
#'
#' @export

plot_proportionalHazard <- function(model){
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

#### ### ### ###
# KAPLAN MEIER #
#### ### ### ###

#' getAutoKM.list
#'
#' @param type Kaplan Meier for linear predictors ("LP"), for PLS components ("COMP") or for original variables ("VAR").
#' @param lst_models List of HDcox models
#' @param comp vector of two components
#' @param top select top X variables
#' @param ori_data Compute the KM with the original data cutpoints or the normalize
#' @param BREAKTIME Break for time points
#' @param only_sig Return only significant log-rank test variables
#' @param alpha Significant cutoff
#' @param title Title of the plot
#' @param verbose Return messages
#'
#' @export

getAutoKM.list <- function(type = "LP", lst_models, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){
  if(!type %in% c("LP", "COMP", "VAR")){
    stop("Type parameters must be one of the following: LP, COMP or VAR")
  }

  if(type == "LP"){
    lst <- purrr::map(lst_models, ~getLPKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "COMP"){
    lst <- purrr::map(lst_models, ~getCompKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else{
    lst <- purrr::map(lst_models, ~getVarKM(model = ., comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }
  return(lst)
}

#' getAutoKM
#'
#' @param type Kaplan Meier for linear predictors ("LP"), for PLS components ("COMP") or for original variables ("VAR").
#' @param model HDcox model
#' @param comp vector of two components
#' @param top select top X variables
#' @param ori_data Compute the KM with the original data cutpoints or the normalize
#' @param BREAKTIME Break for time points
#' @param only_sig Return only significant log-rank test variables
#' @param alpha Significant cutoff
#' @param title Title of the plot
#' @param verbose Return messages
#'
#' @export

getAutoKM <- function(type = "LP", model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){
  if(!type %in% c("LP", "COMP", "VAR")){
    stop("Type parameters must be one of the following: LP, COMP or VAR")
  }

  if(type == "LP"){
    return(getLPKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else if(type == "COMP"){
    return(getCompKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }else{
    return(getVarKM(model, comp = comp, top = top, ori_data = ori_data, BREAKTIME = BREAKTIME, only_sig = only_sig, alpha = alpha, title = title, verbose = verbose))
  }
}

getLPKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

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
  vars_data <- as.data.frame(model$survival_model$lp)
  rownames(vars_data) <- rownames(model$X$data)
  colnames(vars_data) <- "LP"

  vars_num <- vars_data
  if(all(dim(vars_num)>0)){
    info_logrank_num <- getLogRank_NumVariables(data = vars_num, sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
  }else{
    info_logrank_num <- NULL
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / 20
  }

  d <- info_logrank_num$df_numASqual
  v_names <- info_logrank_num$df_nvar_lrtest[,1:2]

  LST_SPLOT <- plot_survivalplot.qual(data = d,
                                       sdata = data.frame(model$Y$data),
                                       BREAKTIME = BREAKTIME,
                                       cn_variables = v_names$Variable,
                                       name_data = NULL, title = title)

  return(list(info_logrank_num = info_logrank_num, LST_PLOTS = LST_SPLOT))

}

getCompKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  # DFCALLS
  lst_vars <- info_logrank_qual <- NULL

  if(attr(model, "model") %in% c(pkg.env$pls_methods, pkg.env$multiblock_methods)){

    if(!all(is.null(model$survival_model))){
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
    vars_data <- as.data.frame(model$X$scores[rownames(model$X$scores),unique_vars,drop=F])
  }else{
    vars_data <- list()
    for(b in names(model$X$data)){
      unique_vars <- deleteIllegalChars(unique(unlist(lst_vars[[b]])))
      vars_data[[b]] <- as.data.frame(model$X$scores[[b]][rownames(model$X$scores[[b]]),unique_vars,drop=F])
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
      vars_num[[b]] <- vars_data[[b]]

      if(all(dim(vars_num[[b]]))>0){
        info_logrank_num[[b]] <- getLogRank_NumVariables(data = vars_num[[b]], sdata = data.frame(model$Y$data), VAR_EVENT = "event", name_data = NULL, minProp = 0.1, ROUND_CP = 4)
      }else{
        info_logrank_num[[b]] <- NULL
      }
    }
  }

  if(is.null(BREAKTIME)){
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / 20
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

getVarKM <- function(model, comp = 1:2, top = 10, ori_data = T, BREAKTIME = NULL, only_sig = F, alpha = 0.05, title = NULL, verbose = FALSE){

  if(attr(model, "model") %in% pkg.env$pls_methods){

    if(all(is.null(model$survival_model))){
      if(verbose){
        message("Survival cox model not found")
      }
      return(NA)
    }

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

      if(attr(model, "model") %in% pkg.env$sb.plsicox){
        aux <- model$list_pls_models[[b]]
      }else if(attr(model, "model") %in% pkg.env$sb.splsdrcox){
        aux <- model$list_spls_models[[b]]
      }

      if(attr(model, "model") %in% c(pkg.env$sb.plsicox, pkg.env$sb.splsdrcox)){

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
    BREAKTIME <- (max(model$Y$data[,"time"]) - min(model$Y$data[,"time"])) / 20
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
    if(cn==VAR_EVENT) #skip outcome variable
      next

    variable <- data[,cn] #select the variable

    tbl <- as.data.frame(sort(table(variable)))
    if(all(dim(tbl)==c(1,1))) next #just one factor
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
    if(length(grep("-", cn, fixed = T))>0){
      cn <- gsub("-", "_", x = cn, fixed = T)
    }
    if(length(grep("+", cn, fixed = T))>0){
      cn <- gsub("+", ".", x = cn, fixed = T)
    }
    if(length(grep("*", cn, fixed = T))>0){
      cn <- gsub("*", ".star.", x = cn, fixed = T)
    }

    colnames(auxData)[3] <- cn

    # Determine the optimal cutpoint for continuous variables, using the maximally selected rank statistics from the 'maxstat' R package.
    minProp = minProp #we have to establish a minimum number of patients per group in %

    res.cut <- tryCatch(
      expr = {
        survminer::surv_cutpoint(auxData, time="time", event="event", variables = cn, minprop = minProp)
      },
      # Specifying error message
      error = function(e){
        message(paste0(cn, ": ", e))
        NA
      }
    )

    if(all(is.na(res.cut))){
      next
    }

    cutpoint_value <- round(res.cut$cutpoint[1,1], ROUND_CP)
    variable <- ifelse(variable>cutpoint_value, paste0("greater than ", cutpoint_value), paste0("lesser/equal than ", cutpoint_value))
    variable <- data.frame(factor(variable))
    colnames(variable) = cn

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
      LST_NVAR_SIG <- rbind(LST_NVAR_SIG, c(cn, NA, NA))
      next
    }else{
      pval <- surv_pvalue(kmsurvival)
      LST_NVAR_SIG <- rbind(LST_NVAR_SIG, c(cn, round(pval$pval,4), cutpoint_value))
    }

  }

  LST_NVAR_SIG <- as.data.frame(LST_NVAR_SIG)
  LST_NVAR_SIG[,2] <- as.numeric(LST_NVAR_SIG[,2])
  LST_NVAR_SIG[,3] <- as.numeric(LST_NVAR_SIG[,3])

  if(exists("VAR_DESCRIPTION")){
    colnames(LST_NVAR_SIG) <- c("Variable", "P-Val (Log Rank)", "Cutoff", "Description")
  }else{
    colnames(LST_NVAR_SIG) <- c("Variable", "P-Val (Log Rank)", "Cutoff")
  }

  LST_NVAR_SIG <- LST_NVAR_SIG[order(LST_NVAR_SIG$`P-Val (Log Rank)`),]

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

      colnames(aux)[3] <- cn

      f = as.formula(paste0("Surv(time = time, event = event) ~ ", cn))

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

      if(requireNamespace("RColorConesa", quietly = TRUE)){
        colors <- RColorConesa::colorConesa(length(levels(data[,cn])))
        names(colors) <- NULL
      }else{
        colors <- NULL
      }

      #GGSURVPLOT DOES NOT PRINT INTERVALS IF ALL DATA IS NOT SELECTED FOR RIBBON STYLE
      #IF PROBLEMS CHANGE TO STEP STYLE
      kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", palette = colors,
                                      conf.int = TRUE, ggtheme = theme_bw(), legend.labs = levels(aux[,cn]),
                                      conf.int.style = "ribbon",
                                      conf.int.alpha = 0.25,
                                      xlim = c(0, round2any(max(aux$time), 5, ceiling)),
                                      pval = T,
                                      surv.median.line = "hv", # Add medians survival
                                      risk.table = TRUE,
                                      legend.title = cn,
                                      break.time.by = BREAKTIME,
                                      font.caption = 8,
                                      font.x = 10,
                                      font.y = 10,
                                      font.tickslab = 8,
                                      font.legend = 8,
                                      title = title)

      kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
        theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))

      lst_splots[[cn]] <- kmplot
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

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", color = colors,
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
                                    font.legend = 8,
                                    title = title)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["SurvivalFunction"]] <- kmplot

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", color = colors, fun = "event",
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
                                    font.legend = 8,
                                    title = title)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["HazardCurve"]] <- kmplot

    kmplot <- survminer::ggsurvplot(fit = kmsurvival, censor.shape = "|", color = colors, fun = "cumhaz",
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
                                    font.legend = 8,
                                    title = title)

    kmplot$table <- kmplot$table + labs(title = "Patients at risk") +
      theme(axis.text = element_text(size = 8)) + theme(axis.title = element_text(size = 10))
    lst_splots[["CumulativeHazard"]] <- kmplot
  }

  return(lst_splots)
}
