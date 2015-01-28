#' Load Atlantis outputfiles (netcdf)
#' 
#' 
#' This function loads Atlantis outputfiles (netcdf) and converts them to a dataframe.
#' @param model_path Character string of the ATLANTIS folder.
#' @param filename Character string of the general ATLANTIS output file. Usually "outputNorthSea.nc".
#' @param select_groups Vector of funtional groups which shall be plotted. Names have to match the ones used in the ncdf file. Check column "Name" in "functionalGroups.csv" as example.
#' @param select_variable Character value spefifying which variables shall be loaded. Only one variable can be loaded at a time. Currently, only "N", "Nums", "ResN", "StructN", "Eat", "Growth", "Prodn", "Grazing" can be selected.
#' @param biomasspools Character vector giving the names of biomasspools. Note this does not mean groups which are considered as biomasspools in ATLANTIS but species which are only present in the botom layer. Therefore they are found in "boxTracers.csv"
#' @param aggregate_layers Logical specifying if layers shall be aggregated or not.
#' @param load_init Logical specifying if an initial file shall be read in. Default is F wich means output files are loaded.
#' @param remove_bboxes Logical specifying if an boundary boxes  shall be excluded. Default is T.
#' @param modelstart Character value giving the start of the model run in "yyyy-mm-dd". Default is "1991-01-01".
#' @return Dataframe in long format with the following informations: Species, timestep, polygon, agecl and value (= "atoutput").
#' 
#' 
#' @details This functions converts the ATLANTIS output to a dataframe which can be processed in R.
#' @keywords gen
#' @examples 
#' load_atlantis_output(model_path = file.path("z:", "Atlantis", "ATLANTIS NSmodel base"), filename = "outputNorthSea.nc", select_groups = get_groups(), select_variable = "ResN", biomasspools = c("large_crabs", "small_epifauna", "sessile_epifauna", "epifaunal_macrobenthos"))
#' @export

load_atlantis_ncdf <- function(model_path, filename, select_groups, select_variable, biomasspools, aggregate_layers = T, load_init = F, remove_bboxes = T, modelstart = "1991-01-01"){
  # Check input!
  supported_variables <- c("N", "Nums", "ResN", "StructN", "Eat", "Growth", "Prodn", "Grazing")
  if(length(select_groups) == 0) stop("No Groups selected.")
  if(length(select_variable) == 0) stop("No Variables selected.")
  if(length(biomasspools) == 0) stop("No biomasspools selected.")
  if(length(select_variable) > 1) stop("Only one variable allowed per function call.")
  if(select_variable != "N" & all(is.element(select_groups, biomasspools))) stop("The only output for Biomasspools is N.")
  if(any(!is.element(select_variable, supported_variables))) stop(paste("Only", paste(supported_variables, collapse = ", "), "can be selected as 'select_variable'"))
  report <- 0 # check if if-statements are evaluated!
  
  # Check input structure!
  active_groups <- read.table(file.path(model_path, "functionalGroups.csv"), sep = ",", header = T)
  active_groups <- as.vector(subset(active_groups, IsTurnedOn == 1)$Name)  
  inactive_groups <- select_groups[which(!is.element(select_groups, active_groups))]
  if(length(inactive_groups) > 1){
    warning(print(paste("Some selected groups are not active in the model run. Check 'IsTurnedOn' in functionalGroups.csv in", model_path, "\n"), 
                  paste(inactive_groups, collapse = "\n")))
  }
  if(all(!is.element(select_groups, active_groups))){
    stop(paste("None of the species selected are active in the model run. Check spelling and Check 'IsTurnedOn' in functionalGroups.csv in", model_path, "\n"))
  }
    
  # Load ATLANTIS output!
  at_out <- ncdf::open.ncdf(file.path(model_path, filename))
  
  # Get info from netcdf file! (Filestructure and all variable names)
  var_names_ncdf <- names(at_out$var)
  n_timesteps <- at_out$dim[[1]]$len
  n_boxes     <- at_out$dim[[2]]$len
  n_layers    <- at_out$dim[[3]]$len
  
  # Extract data from the ncdf file! Create a vector of all potential variable names first! Only use names which
  # are available in the ncdf-file as an extraction of missing variables is not possible! Unfortunately variable
  # names for Prodn and Garzing use a "" instead of "_" as seperator... :)
  # Create vecotr of available species at the end using search_clean! This is needed to create species-names
  # lateron! This approach may seem complicate but it turns out that this approach is very robust since no
  # user input is needed as the variable names are basically extracted from the available names in the ncdf file!
  search <- unlist(lapply(select_groups, paste0, c("", 1:10)))
  if(is.element(select_variable, c("Prodn", "Grazing"))){
    search <- unlist(lapply(search, paste, select_variable, sep = ""))
  } else {
    search <- unlist(lapply(search, paste, select_variable, sep = "_"))
  }
  search_clean <- search[is.element(search, var_names_ncdf)]
  at_data <- lapply(search_clean, ncdf::get.var.ncdf, nc = at_out) 
  
  # Get final species and number of ageclasses per species
  final_species <- select_groups[sapply(lapply(select_groups, grepl, x = search_clean), any)]
  final_agecl <- read.table(file.path(model_path, "functionalGroups.csv"), sep = ",", header = T)
  final_agecl <- final_agecl$NumCohorts[is.element(final_agecl$Name, final_species)]
  
  # Exract timesetp from Parameterfile!
  toutinc <- get_timestep(model_path = model_path)
  
  # Combine List of arrays to one dataframe! Add additional informations:
  # Species-names, polygons, layers, ageclasses! This process is split into 3 categories of 
  # list-types:
  # 1) All list elements have 3 dimensions --> StructN/ResN/Nums
  # 2) All list elements have 2 dimensions --> Eat/Growth/Prodn/Grazing/Catch/Discards
  # 3) list elements have different dimensions --> N
  
  # 1) All list elements have 3 dimensions --> StructN/ResN/Nums
  if(all(sapply(lapply(at_data, dim), length) == 3) & select_variable != "N"){# datainfos: layer, boxes, time
    if(length(final_species[final_agecl == 1]) > 0) warning("Some selected groups are not age-structured and not included in the output!")
    if(aggregate_layers){
      if(select_variable == "Nums"){
        at_data <- lapply(at_data, colSums)
      } else {# select_variable == "StructN" | select_variable == "ResN"
        at_data <- lapply(at_data, colMeans)
      }    
    at_data <- combine_lists_add_polygons(at_data, n_boxes, steps = n_timesteps)
    at_data$agecl <- unlist(lapply(lapply(final_agecl, seq, from = 1, by = 1), rep, each = n_boxes))
    at_data$species <- rep(rep(final_species, times = final_agecl), each = n_boxes)
    } else {
      at_data <- lapply(at_data, array_to_list) 
      at_data <- do.call(c, at_data)
      at_data <- combine_lists_add_polygons(at_data, n_boxes, steps = n_timesteps)
      at_data$agecl <- unlist(lapply(lapply(final_agecl, seq, from = 1, by = 1), rep, each = n_layers * n_boxes))
      at_data$species <- rep(rep(final_species, times = final_agecl), each = n_boxes * n_layers) 
      # The following problem could be solved if one konws how to pass the 2nd function argument to lapply!
      if(length(unique(final_agecl)) > 1) stop("Cohort numbers differ among age based groups! Exraction not possible. Contact package development Team.")
      at_data$layer <- rep(rep(c(c((n_layers-2):0), "sed"), each = n_boxes), times = unique(final_agecl) * length(final_species))
    }
    report <- 1
  }
  
  # 2) All list elements have 2 dimensions --> Eat/Growth/Prodn/Grazing/Catch/Discards
  if(all(sapply(lapply(at_data, dim), length) == 2) & select_variable != "N"){# datainfos: boxes, time
    at_data <- combine_lists_add_polygons(at_data, n_boxes, n_timesteps)
    if(is.element(select_variable, c("Prodn", "Grazing"))){# No ageclasses!
      at_data$species <- rep(final_species, each = n_boxes) 
    } else {
      at_data$agecl <- unlist(lapply(lapply(final_agecl, seq, from = 1, by = 1), rep, each = n_boxes))
      at_data$species <- rep(rep(final_species, times = final_agecl), each = n_boxes)       
    }
    report <- 1
  }
  
  # 3) List elements have different dimensions --> N! Well if only one group es selected, e.g. "cod"
  #    this may also become FALSE... Therefore we only check for "N" as input!
  if(select_variable == "N"){
    if(aggregate_layers){
      for(i in seq_along(at_data)){
        if(length(dim(at_data[[i]])) == 3) at_data[[i]] <- colMeans(at_data[[i]])
      }
      at_data <- combine_lists_add_polygons(at_data, n_boxes, steps = n_timesteps)
      at_data$species <- rep(final_species, each = n_boxes)
    } else {# split into different lists!
      for(i in seq_along(at_data)){
        if(length(dim(at_data[[i]])) == 3){
          at_data[[i]] <- array_to_list(at_data[[i]])
          at_data[[i]] <- lapply(at_data[[i]], as.data.frame)
          at_data[[i]] <- do.call(rbind, at_data[[i]])
          at_data[[i]]$polygon <- rep(c(0:(n_boxes - 1)), times = n_layers)
          at_data[[i]]$layer <- rep(0:(n_layers - 1), each = n_boxes)
        } else {
          at_data[[i]] <- as.data.frame(at_data[[i]])
          at_data[[i]]$polygon <- c(0:(n_boxes - 1))
          at_data[[i]]$layer <- 7
        }
      }
      for(i in seq_along(at_data)){
        at_data[[i]]$species <- final_species[i]
      }
    at_data <- do.call(rbind, at_data)
    }
    report <- 1
  }
    
  if(report == 0) stop("Critical error. Contact package development team: Report = 0!")
  
  # Convert to long dataframe!
  names(at_data)[1:n_timesteps] <- c(0:(n_timesteps-1))
  at_data_long <- reshape2::melt(at_data, id.vars = names(at_data[(n_timesteps + 1):length(names(at_data))]), measure.vars = names(at_data)[1:n_timesteps], 
                       variable.name = "timestep", value.name = "atoutput")
  at_data_long$timestep <- as.numeric(levels(at_data_long$timestep))[at_data_long$timestep]
  
  # Convert timestep to actual time!
  at_data_long$time <- with(at_data_long, as.Date.numeric(timestep * toutinc, origin = modelstart))
  
  # Remove boundary boxes!
  if(remove_bboxes){
    at_data_long <- at_data_long %>%
      dplyr::filter(!is.element(polygon, get_bboxes()))
  }
  
  # Add categories and order species (factor) by category_plot!
  at_data_long <- set_categorys(at_data_long)    
  ag_data <- at_data_long %>%
    dplyr::select(species, category_plot) %>%
    unique() %>%
    dplyr::arrange(category_plot, species)
  at_data_long$species <- factor(at_data_long$species, levels = ag_data[order(ag_data$category_plot), "species"])
    
  return(at_data_long) 
}






