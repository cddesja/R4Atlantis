#' Extract fishbase IDs using the package "rfishbase"
#'
#'
#' This function extracts fishbase IDs using the package "rfishbase"
#' @param vector of fish species with genus and species
#' @param specification if subspecies (e.g. Sprattus sprattus balticus) should be excluded!
#' @return a list with species, and fishbase IDs
#'
#' @details The function depends on the package "rfishbase" which creates a local copy of the fishbase database. The IDs are needed to generate URLs to scan www.fishbase.org for detailed informations about fish growth.
#' @examples
#' get_ids_fishbase(c("Gadus morhua", "Merlangius merlangus"))
#' @export


get_ids_fishbase <- function(fish, exclude_subspecies = T){
  # Check if every fishname is composed of genus and species!
  if(any(sapply(str_split(fish, pattern = " "), length) < 2)) stop("Fishnames not complete!")
  
  fish_data <- loadCache()
  
  # Create list of all fishbase.names!
  fb_names <- sapply(fish_data, function(x)x$ScientificName)
  
  # Check if every fishname is available at fishbase!
  xx <- lapply(fish, grepl, x = fb_names)
  xx <- which(unlist(lapply(xx, any)) == F)
  if(length(xx) >= 1) warning("Not available at fishbase:\n", paste(fish[xx], collapse = "\n"))
  
  # Get list_id for every Prey.Species.Name! Thus only information for requested species is extracted!
  list_id <- unlist(lapply(fish, grep, x = fb_names))
  fb_names <- fb_names[list_id]
  
  # Extract fishbase id for every list_id
  fishbase_id <- unlist(lapply(list_id, function(x)fish_data[[x]]$id))
  
  # Exclude Subspecies if needed!
  if(exclude_subspecies){
    exclude_ids <- which(sapply(str_split(fb_names, " "), length) == 3)
    if(length(exclude_ids) >= 1){
      fb_names    <- fb_names[-exclude_ids]
      list_id     <- list_id[-exclude_ids]
      fishbase_id <- fishbase_id[-exclude_ids]     
    }
  }
  
  result <- list(fishbase_id = fishbase_id, fishbase_name = fb_names)
  return(result)
}
