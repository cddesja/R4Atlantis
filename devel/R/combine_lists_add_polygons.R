combine_lists_add_polygons <- function(list, n_boxes){
  result <- lapply(list, as.data.frame)
  result <- do.call(rbind, result)
  result$polygon  <- rep(c(0:(n_boxes - 1)), times = dim(result)[1] / n_boxes)
  return(result)
}
