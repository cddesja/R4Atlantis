#' Load Atlantis outputfiles (netcdf)
#' 
#' 
#' This function loads Atlantis outputfiles (netcdf) and converts them to a dataframe.
#' @param model_path Character string of the ATLANTIS folder.
#' @param filename Character string of the general ATLANTIS output file. Usually "outputNorthSea.nc".
#' @param physic_variables Character value spefifying which variables shall be loaded. Only one variable can be loaded at a time. Currently, only "salt", "NO3", "NH3", "Temp", "Oxygen", "Si", "Det_Si", "DON", "Chl_a", "Denitrifiction", "Nitrification" can be selected.
#' @param aggregate_layers Logical specifying if layers shall be aggregated or not.
#' @param load_init Logical specifying if an initial file shall be read in. Default is F wich means output files are loaded.
#' @param remove_bboxes Logical specifying if an boundary boxes  shall be excluded. Default is T.
#' @param modelstart Character value giving the start of the model run in "yyyy-mm-dd". Default is "1991-01-01".
#' @return Dataframe in long format with the following informations: variable, timestep, polygon, and value (= "atoutput").

#' @details This functions converts the ATLANTIS output to a dataframe which can be processed in R.
#' @keywords gen
#' @examples 
#' load_atlantis_ncdf_physics(model_path = file.path("z:", "Atlantis", "ATLANTIS NSmodel base"), filename = "outputNorthSea.nc", physic_variables = c("salt", "Temp"))
#' @export

load_atlantis_ncdf_physics <- function(model_path, 
                                       filename, 
                                       physic_variables, 
                                       aggregate_layers = T, 
                                       load_init = F, 
                                       remove_bboxes = T, 
                                       modelstart = "1991-01-01"){  
  if(is.null(physic_variables)) stop("No physical variables selected.")
  supported_variables <- get_physics()
  
  wrong_input <- physic_variables[which(!is.element(physic_variables, supported_variables))]
  
  if(length(wrong_input) > 1){
    stop(paste(wrong_input, "not part of", paste(supported_variables, collapse = ", ")))
  } 
  
  # Load ATLANTIS output!
  at_out <- ncdf::open.ncdf(file.path(model_path, filename))
  
  # Extract number of timesteps, polygons and layers from netcdf
  n_timesteps <- at_out$dim[[1]]$len
  n_boxes     <- at_out$dim[[2]]$len
  n_layers    <- at_out$dim[[3]]$len
  
  # Exract outputlength from Parameterfile!
  toutinc <- get_timestep(model_path = model_path)
  
  physic_variables <- sort(physic_variables)
      
  physic_output <- lapply(physic_variables, ncdf::get.var.ncdf, nc = at_out)  
  if(aggregate_layers){
    physic_output <- lapply(physic_output, colMeans)
  } else {
    if(n_timesteps > 1) physic_output <- lapply(physic_output, array_to_list)
    if(n_timesteps > 1) physic_output <- do.call(c, physic_output)
  }
  if(load_init) physic_output <- lapply(physic_output, t)
  physic_output <- combine_lists_add_polygons(physic_output, n_boxes, n_timesteps)
  
  if(aggregate_layers){
    physic_output$variable <- rep(physic_variables, each = n_boxes)
  } else {
    physic_output$variable <- rep(physic_variables, each = n_boxes * n_layers)
    physic_output$layer <- rep(rep(c(0:(n_layers-1)), each = n_boxes), times = length(physic_variables))
  } 
    
  # Convert to long dataframe!
  names(physic_output)[1:n_timesteps] <- c(0:(n_timesteps-1))
  physic_output_long <- reshape2::melt(physic_output, id.vars = names(physic_output[(n_timesteps + 1):length(names(physic_output))]), 
                                       measure.vars = names(physic_output)[1:n_timesteps], 
                                       variable.name = "timestep", value.name = "atoutput")
  physic_output_long$timestep <- as.numeric(levels(physic_output_long$timestep))[physic_output_long$timestep]
  
  # Convert timestep to actual time!
  physic_output_long$time <- with(physic_output_long, as.Date.numeric(timestep * toutinc, origin = modelstart))
  
  if(remove_bboxes){
    physic_output_long <- physic_output_long %>%
      dplyr::filter(!is.element(polygon, get_bboxes()))
  }
  
  return(physic_output_long) 
}



