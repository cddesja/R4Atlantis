load_atlantis_output <- function(select_groups, select_variable, var_names_ncdf, at_out){
  search <- unlist(lapply(select_groups, paste0, c("", 1:10)))
  search <- unlist(lapply(search, paste, select_variable, sep = "_"))
  search <- search[is.element(search, var_names_ncdf)]
  result <- lapply(search, get.var.ncdf, nc = at_out)  
  return(result)
}