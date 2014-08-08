#' Load Atlantis outputfiles (netcdf)
#' 
#' 
#' This function loads Atlantis outputfiles (netcdf)
#' @param model_path Path of the ATLANTIS directory
#' @param filename Name of the ouputfile which shall be read in
#' @param select_variable ATLANTIS variables which shall be loaded! Only one variable can be loaded at a time! "Only 'N', 'Nums', 'ResN' and 'StructN' 
#' @param aggregate_layers Specifies of layers shall be aggregated or not? Currently only T is available.
#' @param select_groups Vector of groups for which ouput shall be loaded. Should be derived from "functionalGroups.csv" using the column Name.
#' @param biomasspools Vector of biomasspool names as given in functionalGroups.csv!
#' @return Dataframe in long format with the following informations: Species, timestep, polygon, agecl and value
#' 
#' @details This functions converts the ATLANTIS output to a dataframe which can be processed in R.
#' @keywords gen
#' @examples 
#' load_atlantis_output(model_path = "Z://Atlantis//ATLANTIS NSmodel", filename = "outputNorthSea.nc", select_variable = "ResN", biomasspools = c("large_crabs", "small_epifauna", "sessile_epifauna", "epifaunal_macrobenthos"))
#' @export

# Author:    Alexander Keth
# Institute: Instutute for Hydrobiology and Fisheries Science (IHF)
# Date:      12.06.2014
  
convert_atlantis_output <- function(model_path, filename, select_groups, select_variable, biomasspools, aggregate_layers = T){
  # Check input!
  if(length(select_groups) == 0) stop("No Groups selected.")
  if(length(select_variable) == 0) stop("No Variables selected.")
  if(length(biomasspools) == 0) stop("No biomasspools selected.")
  if(length(select_variable) > 1) stop("Only one variable allowed per function call.")
  if(select_variable != "N" & all(is.element(select_groups, biomasspools))) stop("The only output for Biomasspools is N.")
  if(all(select_variable != c("N", "Nums", "ResN", "StructN"))) stop("Only 'N', 'Nums', 'ResN' and 'StructN' can be selected as 'select_variable'")
  
  # Check input structure!
  group_data <- read.table(file.path(model_path, "functionalGroups.csv"), sep = ",", header = T)
  group_data <- subset(group_data, IsTurnedOn == 1)
  all_groups <- as.vector(group_data$Name)
  
  groups_bp     <- select_groups[which(is.element(select_groups, biomasspools))]  
  groups_non_bp <- select_groups[which(!is.element(select_groups, biomasspools))]
  
  if(any(!is.element(select_groups, all_groups))){
    warning(print(paste("Some selected groups are not active in the model run. Check 'IsTurnedOn' in functionalGroups.csv in", model_path)))
    if(length(groups_non_bp) >= 1) groups_non_bp <- groups_non_bp[!is.na(match(groups_non_bp, all_groups))]
    if(length(groups_bp) >= 1)     groups_bp <- groups_bp[!is.na(match(groups_bp, all_groups))]
    if(length(groups_non_bp) == 0 & length(groups_bp) == 0) stop("None of the species selected are active in the model run.")
  }
    
  # Load ATLANTIS output!
  at_out <- open.ncdf(file.path(model_path, filename))
  
  # Extract number of timesteps, polygons and layers from netcdf
  n_timesteps <- at_out$dim[[1]]$len
  n_boxes     <- at_out$dim[[2]]$len
  n_layers    <- at_out$dim[[3]]$len
  
  # Get all variable names from to netcdf file!
  var_names_ncdf <- names(at_out$var)
      
  if(select_variable == "N"){
    if(length(groups_bp) >= 1){ # Extract data for biomasspools!
      at_data_bp <- load_atlantis_output(groups_bp, select_variable, var_names_ncdf, at_out)
      if(length(groups_non_bp) == 0) at_data <- at_data_bp # Final result if only biomasspools are selected
    }
    if(length(groups_non_bp) >= 1){ # Extract data for non_biomasspools!
      at_data_non_bp <- load_atlantis_output(groups_non_bp, select_variable, var_names_ncdf, at_out)
      if(aggregate_layers) at_data_non_bp <- lapply(at_data_non_bp, colSums)
      if(length(groups_bp) == 0) at_data <- at_data_non_bp # Final result if no biomasspools are selected       
    }
    if(length(groups_non_bp) >= 1 & length(groups_bp) >= 1) at_data <- c(at_data_bp, at_data_non_bp) # Mixed scenario
    
    # Put everything together and add other columns!
    at_data <- combine_lists_add_polygons(at_data, n_boxes)
    at_data$species <- rep(c(groups_bp, groups_non_bp), each = n_boxes)
    
  } else { # select_variable not "N"
    if(length(groups_bp) >= 1) warning("Biomasspools only have N as output. They are not included in the output!")
    cohorts <- as.vector(group_data$NumCohorts[match(groups_non_bp, all_groups)])
    at_data <- load_atlantis_output(groups_non_bp, select_variable, var_names_ncdf, at_out)
    if(aggregate_layers) at_data <- lapply(at_data, colSums)
    
    # Put everything together and add other columns!
    at_data <- combine_lists_add_polygons(at_data, n_boxes)
    at_data$agecl <- unlist(lapply(lapply(cohorts, seq, from = 1, by = 1), rep, each = n_boxes))
    at_data$species <- rep(rep(groups_non_bp, times = cohorts), each = n_boxes)              
  }
  
  # Convert to long dataframe!
  names(at_data)[1:n_timesteps] <- c(0:(n_timesteps-1))
  at_data_long <- melt(at_data, id.vars = names(at_data[(n_timesteps + 1):length(names(at_data))]), measure.vars = names(at_data)[1:n_timesteps], 
                       variable.name = "timestep", value.name = "atoutput")
  return(at_data_long) 
}






