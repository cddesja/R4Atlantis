#' helper functions
#'
#' Set of functions which support the package.
#'
#' @param model_path Character string of the ATLANTIS folder.
#' @param list List of dataframes read in by load_atlantis_ncdf.
#' @param n_boxes Integer giving the number of boxes (polygons). Automatically extracted from the
#' ncdf-file.
#' @param steps Integer giving the number of timesteps during the model run. Automatically extracted 
#' from the ncdf-file.
#' @param array Array.
#' @param paths Model paths defining the old and new model run in compare_atlantis.
#' @param two_dataframes Dataframes defining the old and new model run in compare_atlantis.
#' @param data Dataframe created by 'load_atlantis_ncdf'.
#' @param axis_label Character vector giving the axis labels.
#' @param which_axis Charater string specifying the x or y-axis.
#' @param plots_per_row Integer giving the number of plots per row in faceted plots.
#' 
#'
#' @details \describe{
#'
#' \item{\code{read_functionalgroups}}{
#' Read in 'functionalgroups.csv' from the model folder.}
#'
#' \item{\code{get_groups}}{
#' Extract the names of all functional groups from 'functionalgroups.csv' using the column 'Name'.}
#'
#' \item{\code{get_age_groups}}{
#' Extract the names of all age structured functional groups from 'functionalgroups.csv' using the column 'Name'.}
#'
#' \item{\code{get_acronyms}}{
#' Extract the acronyms of all functional groups from 'functionalgroups.csv' using the column 'Code'.}
#'
#' \item{\code{get_vert_acronyms}}{
#' Extract the acronyms of all age structured functional groups from 'functionalgroups.csv' using 
#' the column 'Code'. Note: The function may be renamed at some point to get_age_acronyms.}
#'
#' \item{\code{get_invert_acronyms}}{
#' Extract the acronyms of all NON age structured functional groups from 'functionalgroups.csv' using 
#' the column 'Code'. Note: The function may be renamed at some point to get_non_age_acronyms.}
#'
#' \item{\code{combine_lists_add_polygons}}{
#' Used to combine a list to a dataframe in various ploting routines.}
#'
#' \item{\code{array_to_list}}{
#' Used to convert an array to a list in various ploting routines.}
#'
#' \item{\code{add_model_run}}{
#' Add model run als linetype to plots.}
#'
#' \item{\code{combine_dataframes}}{
#' Combine multiple dataframes which were read in using 'load_atlantis_ncdf'. Used in 'compare_atlantis'.
#' Currently only 2 dataframes can be passed.}
#'
#' \item{\code{combine_output}}{
#' Combine}
#'
#' \item{\code{mean_over_ages}}{
#' Compute mean of atlantis output ('atoutput') per group, ageclass for each timestep.}
#' \item{\code{mean_over_polygons}}{
#' Compute mean of atlantis output ('atoutput') per group, ageclass and polygon for each timestep.}
#'
#' \item{\code{mean_overview}}{
#' Compute mean of atlantis output ('atoutput') per group for each timestep.}
#'
#' \item{\code{mean_over_ages_compare}}{
#' Compute mean of atlantis output ('atoutput') per group, ageclass and model run for each timestep. 
#' Only used in 'compare_atlantis'}
#'
#' \item{\code{mean_overview_compare}}{
#' Compute mean of atlantis output ('atoutput') per group and model run for each timestep. 
#' Only used in 'compare_atlantis'}
#'
#' \item{\code{plot_age}}{
#' Create ggplot for each species and ageclass.}
#'
#' \item{\code{plot_polygon}}{
#' Create ggplot for each species, ageclass and polygon.}
#'
#' \item{\code{plot_overview}}{
#' Create ggplot for each species}
#'
#' \item{\code{plot_calibrate}}{
#' Create ggplot for each species. Ageclasses are shown in different colours.}
#'
#' \item{\code{add_linetype}}{
#' Add linetype to ggplots based on model run. Used in 'compare_atlantis'.}
#'
#' \item{\code{add_label}}{
#' Add labels for ageclass and polygon in ggplots.}
#'
#' \item{\code{add_axis_label}}{
#' Add x- and y-axis labels to ggplots.}
#'
#' \item{\code{get_date_mv}}{
#' Get system date and modelversion. This is added to create output folders and names.}
#'
#' \item{\code{scale_width}}{
#' Set the width of plots depending on the number of plotted groups.}
#'
#' \item{\code{scale_height}}{
#' Set the height of plots depending on the number of plotted groups.}
#'
#' }
#'
#'
#' @name helper_functions
NULL

# Import '%>%' operator from magrittr
#' @importFrom magrittr %>%

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 1) Extract Names/Groups/Acronyms/Bboxes/Physics/Biomasspools *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions
read_functionalgroups <- function(model_path = file.path("z:", "Atlantis", "ATLANTIS NSmodel base")){
  result <- read.table(file.path(model_path, "functionalGroups.csv"), sep = ",", header = T, stringsAsFactors = F)
  return(result)
}

#' @export
#' @rdname helper_functions
get_groups <- function(){
  result <- read_functionalgroups()
  result <- result$Name
  return(result)
}

#' @export
#' @rdname helper_functions
get_age_groups <- function(){
  result <- read_functionalgroups()
  result <- subset(result, NumCohorts > 1)$Name
  return(result)
}

#' @export
#' @rdname helper_functions
get_acronyms <- function(){
  result <- read_functionalgroups()
  result <- result[, names(result) == "Code"]
  return(result)
}

#' @export
#' @rdname helper_functions
get_vert_acronyms <- function(){
  result <- read_functionalgroups()
  result <- subset(result, NumCohorts > 1, select = "Code")[,1]
  return(result)
}

#' @export
#' @rdname helper_functions
get_invert_acronyms <- function(){
  result <- read_functionalgroups()
  result <- subset(result, NumCohorts == 1, select = "Code")[,1]
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 2) NCDF-files *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions
combine_lists_add_polygons <- function(list, n_boxes, steps){
  if(steps == 1) result <- do.call(c, list)
  if(steps >= 2) result <- do.call(rbind, list)
  result <- as.data.frame(result)
  result$polygon  <- rep(c(0:(n_boxes - 1)), times = dim(result)[1] / n_boxes)
  return(result)
}

#' @export
#' @rdname helper_functions
array_to_list <- function(array){
  result <- list()
  for(i in seq(1, dim(array)[1])){
    result[[i]] <- array[i,,]
  }
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 3) Plot ATLANTIS output/input *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions
add_model_run <- function(plots){
  if(!all(sapply(plots, ggplot2::is.ggplot))) stop("Some plots are no ggplots!")
  for(i in seq_along(plots)){
    plots[[i]] <- plots[[i]] + ggplot2::aes(linetype = run)
  }
  return(plots)
}


#' @export
#' @rdname helper_functions
combine_dataframes <- function(paths, two_dataframes){
  if(!is.list(two_dataframes)) stop("Passed dataframes have to be stored as list.")
  if(length(two_dataframes) != 2) stop("Only two dataframes can be passed.")
  new_data <- two_dataframes[[which(paths == model_path)]]
  old_data <- two_dataframes[[which(paths != model_path)]]
  old_data$run <- "old"
  new_data$run <- "new"
  data <- rbind(old_data, new_data)
  return(data)
}

#' @export
#' @rdname helper_functions
# Aggregationfunctions for general output plots
mean_over_ages <- function(data){
  result <- data %>%
    dplyr::group_by(species, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  return(result)
}

#' @export
#' @rdname helper_functions
mean_over_polygons <- function(data){
  result <- data %>%
    dplyr::group_by(species, polygon, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  return(result)
}

#' @export
#' @rdname helper_functions
mean_overview <- function(data){
  result <- data %>%
    dplyr::group_by(species, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  return(result)
}

#' @export
#' @rdname helper_functions
# Aggregationfunctions for general compare plots
mean_over_ages_compare <- function(data){
  result <- data %>%
    dplyr::group_by(species, agecl, time, category_plot, run) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  return(result)
}

#' @export
#' @rdname helper_functions
mean_overview_compare <- function(data){
  result <- data %>%
    dplyr::group_by(species, time, category_plot, run) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  return(result)
}

#' @export
#' @rdname helper_functions
# Calibration functions! Divide by input year!
datatrans_calibrate <- function(data){
  min_time <- min(data$time)
  ref <- data %>%
    dplyr::filter(time == min_time)
  ref$time <- NULL
  names(ref)[names(ref) == "atoutput"] <- "atoutput_ref"
  result <- data %>%
    dplyr::left_join(ref) %>%
    dplyr::mutate(atoutput = atoutput / atoutput_ref)
  result$atoutput[result$atoutput_ref == 0] <- 0
  # Strangely ordering gets lost due to standardisation. Add it again here...
  ag_data <- result %>%
    as.data.frame() %>%
    dplyr::select(species, category_plot) %>%
    unique() %>%
    dplyr::arrange(category_plot, species)  
  levels(result$species) <- ag_data[order(ag_data$category_plot), "species"]
  
  return(result)
}

#' @export
#' @rdname helper_functions
plot_age <- function(data){
  ggplot2::ggplot(data = data, ggplot2::aes(x = time, y = atoutput, colour = category_plot)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(species ~ agecl, scales = "free") +
    theme_standard()
}

#' @export
#' @rdname helper_functions
plot_polygon <- function(data){
  ggplot2::ggplot(data = data, ggplot2::aes(x = time, y = atoutput)) +
    ggplot2::geom_line(ggplot2::aes(colour = category_plot)) +
    ggplot2::facet_grid(species ~ polygon, scales = "free") +
    theme_standard(font_tiny = 10)
}

#' @export
#' @rdname helper_functions
plot_overview <- function(data, plots_per_row){
  ggplot2::ggplot(data = data, ggplot2::aes(x = time, y = atoutput, colour = category_plot)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ species, scales = "free_y", ncol = plots_per_row) + 
    theme_standard(scale_font = 0.8)
}

#' @export
#' @rdname helper_functions
plot_calibrate <- function(data, plots_per_row){   
  anno <- c(min(data$time), max(data$time))
  ggplot2::ggplot(data = data, ggplot2::aes(x = time, y = atoutput, colour = factor(agecl))) +
    ggplot2::annotate("rect", xmin = anno[1], xmax = anno[2], ymin = 0.5, ymax = 1.5, alpha = 0.1) +
    ggplot2::annotate("rect", xmin = anno[1], xmax = anno[2], ymin = 0.8, ymax = 1.2, alpha = 0.3) +
    ggplot2::geom_line() +
#     ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::facet_wrap(~ species, scales = "free_y", ncol = plots_per_row) +
    theme_standard(scale_font = 0.8)    
}

#' @export
#' @rdname helper_functions
add_linetype <- function(plot){
  plot <- plot + ggplot2::aes(linetype = run)
  return(plot)
}

#' @export
#' @rdname helper_functions
add_label <- function(plot, plot_data, path = NULL){
  # Check if multiple columns are available as label! NOTE: If new labels are added add the correspinding column name here!
  # "path" is needed to add labels to compare-plots! Both model paths (old, new) have to be passed!
  if(sum(is.element(names(plot_data), c("polygon", "agecl"))) > 1) stop("Function 'add_label' has no unique column name!")
  plot_data <- as.data.frame(plot_data)
  
  # Aggregate labels for each species and add x-position!
  labels <- plot_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(y_pos = max(atoutput))
  labels$x_pos <- min(plot_data$time)
  
  # Add y-position for each different label column!
  if(any(is.element(names(plot_data), "polygon"))){
    if(is.null(path)){
      labels <- merge(labels, unique(subset(plot_data, select = "polygon")))
    } else {
      labels <- merge(labels, unique(subset(plot_data, select = c("polygon", "run"))))
    }
    plot <- plot + ggplot2::geom_text(data = labels, ggplot2::aes(x = x_pos, y = y_pos, label = polygon), colour = "black", size = 3, hjust = 0)     
  } 
  if(any(is.element(names(plot_data), "agecl"))){
    if(is.null(path)){
      labels <- merge(labels, unique(subset(plot_data, select = "agecl"))) 
    } else {
      labels <- merge(labels, unique(subset(plot_data, select = c("agecl", "run")))) 
    }      
    plot <- plot + ggplot2::geom_text(data = labels, ggplot2::aes(x = x_pos, y = y_pos, label = agecl), colour = "black", size = 3, hjust = 0)          
  } 
  return(plot)
}

#' @export
#' @rdname helper_functions
add_axis_label <- function(plots, axis_label, which_axis){
  if(any(!sapply(plots, ggplot2::is.ggplot))) stop("Some plots are no ggplots!")
  for(i in seq_along(plots)){
    if(which_axis == "y") plots[[i]] <- plots[[i]] + ggplot2::labs(y = axis_label[[i]])
    if(which_axis == "x") plots[[i]] <- plots[[i]] + ggplot2::labs(x = axis_label[[i]])
  }
  return(plots)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 4) General Plotting Functions *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions
get_date_mv <- function(modelversion){
  sysdate <- stringr::str_replace_all(Sys.Date(), "-", "")
  result <- paste(modelversion, sysdate, sep = "_")
  return(result)  
}

#' @export
#' @rdname helper_functions
scale_width <- function(groups, plot_scale_width = 0.16, plots_per_row = 8, ps = 14){
  ifelse(length(groups) < plots_per_row, length(groups) * plot_scale_width * ps, plots_per_row * plot_scale_width * ps)
}

#' @export
#' @rdname helper_functions
scale_height <- function(groups, plot_scale_height = 0.12, legend_height = 0.1, plots_per_row = 8, ps = 14){
  (ceiling(length(groups) / plots_per_row) * plot_scale_height + legend_height) * ps
}

