#' Extract values for a and b from the length weight relationship from www.fishbase.org
#' 
#'   
#' This function extracts values for a and b from www.fishbase.org.
#' @param fish Vector of fish species with genus and species information.
#' @param exclude_subspecies Specification if subspecies (e.g. Sprattus sprattus balticus) should be excluded.
#' @return A dataframe with species, country, locality, a, b, sex, length, r2.
#'
#' @details Before the actual extraction takes place fishbase IDs for every species are extracted using the function "get_ids_fishbase". The IDs are needed to generate the URLs lateron. At the moment subspecies can only be excluded from the extraction.
#' @examples
#' extract_length_weight_fishbase(c("Gadus morhua", "Merlangius merlangus"))
#' @export

# Comment in for debugging!
# fish <- read.table(file.path(my_path, "fish_species_names_from_ibts.csv"), sep = ";", stringsAsFactors = F, header = F)[,1]
# exclude_subspecies <- T
# fish <- c("Gadus morhua", "Merlangius merlangus")

extract_length_weight_fishbase <- function(fish, exclude_subspecies = T){
  ids <- get_ids_fishbase(fish, exclude_subspecies)  
  
  # Split up Names in species and genus part to generate URLs
  ge <- sapply(str_split(ids[[2]], pattern = " "),function(x)x[1])
  sp <- sapply(str_split(ids[[2]], pattern = " "),function(x)x[2]) 
  
  urls <- paste0("http://fishbase.org/PopDyn/LWRelationshipList.php?ID=", ids[[1]], "&GenusName=", ge, "&SpeciesName=", sp, "&fc=183")
  
  # Perform extraction! May crash due to fishbase downtime...
  result <- list()
  for(i in seq_along(urls)){
    result[[i]] <- readHTMLTable(doc = urls[i], which = 3)
  }
  
  # Remove Species without Length/Growth information!
  pos_missing <- which(sapply(result, is.null) == T)
  if(length(pos_missing) >= 1){
    missing_species <- sort(ids[[2]][pos_missing])
    warning("No length growth relationship available:\n", paste(missing_species, collapse = "\n "))
    ids <- lapply(ids, function(x)x[-pos_missing]) 
    result <- result[-pos_missing]
  }
  
  # add names to dataframes
  for(i in seq_along(result)){
    result[[i]]$species <- ids[[2]][i]
  }
  
  result <- do.call(rbind, result)
  
  # Cleanup
  result_backup <- result
  result$Score <- NULL
  result$a <- as.numeric(as.character(result$a))
  result$b <- as.numeric(as.character(result$b))
  result$Country <- as.character(result$Country)
  result$Locality <- as.character(result$Locality)
  result$r2 <- as.numeric(str_replace_all(string = as.character(result$r2), pattern = "&nbsp", replacement = ""))
  
  names(result)[names(result) == "Country"] <- "country"
  names(result)[names(result) == "Locality"] <- "locality"
  
  return(result)
}



