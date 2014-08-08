#' Extract values for Linf and k from www.fishbase.org
#'
#'
#' This function extracts values for Linf and k from www.fishbase.org
#' @param vector of fish species with genus and species
#' @param specification if subspecies (e.g. Sprattus sprattus balticus) should be excluded!
#' @return a dataframe with species, country, locality, linf and k!
#'
#' @details Before the actual extraction takes place fishbaseh IDs for every species are extracted using the function "get_ids_fishbase". The IDs are needed to generate the URLs lateron. At the moment subspecies can only be excluded from the extraction.
#' @examples
#' extract_linf_k_fishbase(c("Gadus morhua", "Merlangius merlangus))
#' @export


extract_linf_k_fishbase <- function(fish, exclude_subspecies = T){
  ids <- get_ids_fishbase(fish, exclude_subspecies)  
  
  # Split up Names in species and genus part to generate URLs
  ge <- sapply(str_split(ids[[2]], pattern = " "),function(x)x[1])
  sp <- sapply(str_split(ids[[2]], pattern = " "),function(x)x[2]) 
  urls <- paste0("http://fishbase.org/PopDyn/PopGrowthList.php?ID=", ids[[1]], "&GenusName=", ge, "&SpeciesName=", sp, "&fc=183")
  
  # Extract data from fishbase!
  fishbase <- lapply(urls, readLines, warn="F")
  fishbase.backup <- fishbase
  fishbase <- fishbase.backup
  
  # First remove Species without Growth information!
  pos_missing <- which(grepl("The system found no growth information for the requested specie.", fishbase))
  if(length(pos_missing) >= 1){
    warning("No growth information available:\n", paste(ids[[2]][pos_missing], collapse = "\n "))
    ids <- lapply(ids, function(x)x[-pos_missing]) 
    fishbase <- fishbase[-pos_missing]
  }
  
  # Actual extraction is performed!
  table_start <- 142 # Based on: all(sapply(fishbase, grep, pattern = "<table cellpadding") == 142)
  table_end   <- sapply(fishbase, grep, pattern = "<table align=\"center\"") - 2
  for(i in 1:length(fishbase)){
    fishbase[[i]] <- fishbase[[i]][table_start:table_end[i]]
  }
  
  # Extract Linf and K!
  linfk_pos <- lapply(fishbase, grep, pattern = "loo")
  linfk <- list()
  for(i in 1:length(fishbase)){
    linfk[[i]] <- fishbase[[i]][linfk_pos[[i]]]
  }
  
  linf_start <- lapply(lapply(linfk, str_locate, pattern = "&loo="), function(x)x[,2] + 1)
  linf_end   <- lapply(lapply(linfk, str_locate, pattern = "&k="), function(x)x[,1] - 1)  
  linf <- list()
  for(i in 1:length(linfk)){
    linf[[i]] <- str_sub(linfk[[i]], start = linf_start[[i]], end = linf_end[[i]])
  }
  
  k_start <- lapply(lapply(linfk, str_locate, pattern = "&k="), function(x)x[,2] + 1)
  k_end   <- lapply(lapply(linfk, str_locate, pattern = "&id"), function(x)x[,1] - 1)  
  k <- list()
  for(i in 1:length(linfk)){
    k[[i]] <- str_sub(linfk[[i]], start = k_start[[i]], end = k_end[[i]])
  }
  
  # Extract Country and Locality!  
  col_pos <- lapply(fishbase, grep, pattern = "<td>")
  col_length <- sapply(col_pos, length)
  
  country_pos  <- lapply(col_length, seq, from = 11, by = 14)
  for(i in 1:length(country_pos)){
    country_pos[[i]] <- col_pos[[i]][country_pos[[i]]]
  }
  
  country <- list()
  for(i in 1:length(fishbase)){
    country[[i]] <- fishbase[[i]][country_pos[[i]]]
  }
  country <- lapply(country, str_replace_all, pattern = "\t\t\t\t<td>", replacement = "")
  country <- lapply(country, str_replace_all, pattern = "</td>", replacement = "")
  
  locality_pos <- lapply(col_length, seq, from = 12, by = 14) 
  for(i in 1:length(locality_pos)){
    locality_pos[[i]] <- col_pos[[i]][locality_pos[[i]]]
  }
  
  locality <- list()
  for(i in 1:length(fishbase)){
    locality[[i]] <- fishbase[[i]][locality_pos[[i]]]
  }
  locality <- lapply(locality, str_replace_all, pattern = "\t\t\t\t<td>", replacement = "")
  locality <- lapply(locality, str_replace_all, pattern = "</td>", replacement = "")
  
  # Check if dimensions are correct
  if(any(c(sapply(linf, length) == sapply(k, length), sapply(linf, length) == sapply(country, length), sapply(linf, length) == sapply(locality, length)) == F)){
    stop("This should not have happened. Contact package development team.")
  }
  rep_names <- sapply(linf, length)
  names <- rep(ids[[2]], times = rep_names)
  
  result <- data.frame(species = names, country = unlist(country), locality = unlist(locality), linf = unlist(linf), k = unlist(k))
  return(result)
}  
  

