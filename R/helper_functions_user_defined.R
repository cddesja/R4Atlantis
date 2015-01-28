#' helper functions user defined
#'
#' Set of functions which support the package.
#'
#' @param model_path Character string of the ATLANTIS folder.
#' @param data Dataframe created by 'load_atlantis_ncdf'.
#'
#' @details \describe{
#'
#' \item{\code{get_bps}}{
#' Set names of biomasspools. Note: These are not biomasspools as defined in ATLANTIS but all groups
#' which are only present in the sediment layer. Check 'boxTracers.csv' if you are using the initial
#' conditions generator. Use full names as in the ncdf output file.}
#' 
#' \item{\code{get_physics}}{
#' Set physical variables of the ATLANTIS model. By default all available physical parameters are set here:
#' "salt", "NO3", "NH3", "Temp", "Oxygen", "Si", "Det_Si", "DON", "Chl_a", "Denitrifiction", "Nitrification".}
#' 
#' \item{\code{get_boxes}}{
#' Set the number of polygons (boxes).}
#' 
#' \item{\code{get_bboxes}}{
#' Specify the boundary boxes.}
#' 
#' \item{\code{get_fish_acronyms}}{
#' Set the fish acronyms. By default these are extracted from 'functionalgroups.csv' using the column 'Code' 
#' filtering for 'FISH' and 'SHARK' in the InvertType column.}
#' 
#' \item{\code{get_bacteria_acronyms}}{
#' Set the bacteria acronyms. By default these are extracted from 'functionalgroups.csv' using the column 'Code' 
#' filtering for 'SED_BACT' and 'PL_BACT' in the InvertType column.}
#' 
#' \item{\code{get_timestep}}{
#' Extract timestep during model run from parameterfile.}
#' 
#' \item{\code{add_nas_layers}}{
#' WARNING: Only works for IHF members.}
#' 
#' \item{\code{remove_nas_layers}}{
#' WARNING: Only works for IHF members.}
#' 
#' \item{\code{plot_active_poylgon}}{
#' WARNING: Only works for IHF members.}
#' 
#' \item{\code{set_categorys}}{
#' Set categories for various plotting routines (plot_atlantis, compare_atlantis...). Basically these
#' define the colour settings in those plots and help ordering them.}
#' 
#' \item{\code{theme_standard}}{
#' Theme for ggplots.}
#' 
#' }
#'
#'
#' @name helper_functions_user_defined
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 1) Extract Names/Groups/Acronyms/Bboxes/Physics/Biomasspools *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions_user_defined
get_bps <- function(){
  bps <- c("large_crabs", "small_epifauna", "sessile_epifauna", "epifaunal_macrobenthos")
  return(bps)
}

#' @export
#' @rdname helper_functions_user_defined
get_physics <- function(){
  physic_var <- c("salt", "NO3", "NH3", "Temp", "Oxygen", "Si", "Det_Si", "DON", "Chl_a", "Denitrifiction", "Nitrification", "eflux", "vflux", "volume")
  return(physic_var)
}


#' @export
#' @rdname helper_functions_user_defined
get_boxes <- function(){
  c(0:25)
}
  
#' @export
#' @rdname helper_functions_user_defined
get_bboxes <- function(){
  c(0, 16, 17, 23, 24, 25)  
}

#' @export
#' @rdname helper_functions_user_defined
get_fish_acronyms <- function(){
  result <- read_functionalgroups()
  result <- subset(result, InvertType %in% c("FISH", "SHARK"), select = "Code")[,1]
  return(result)  
}

#' @export
#' @rdname helper_functions_user_defined
get_bacteria_acronyms <- function(){
  result <- read_functionalgroups()
  result <- subset(result, InvertType %in% c("SED_BACT", "PL_BACT"), select = "Code")[,1]
  return(result)    
}

#' @export
#' @rdname helper_functions_user_defined
get_timestep <- function(model_path){
  toutinc <- readLines(file.path(model_path, "NorthSea_run_fishing_F.prm"))
  toutinc <- toutinc[grep("toutinc", toutinc)]
  toutinc <- stringr::str_split(toutinc, pattern = "\t")[[1]][2]
  toutinc <- as.numeric(stringr::str_split(toutinc, pattern = " ")[[1]][1])
  return(toutinc)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 3) Plot ATLANTIS output/input *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions_user_defined
# Add NAs for missing layers!
add_nas_layers <- function(data){
#   layers <- load("boxTracers.RData")
  layers <- readLines(file.path("z:", "Atlantis/Initial_conditions_generator/atlantisNCGen/NorthSea/boxTracers.csv"), n = 1)
  layers <- unlist(stringr::str_split(layers, pattern = ","))
  layers <- as.numeric(layers[2:length(layers)])
  layers <- lapply(layers-1, seq, from = 0)
  layers <- lapply(layers, function(x)c(x,7))
  layers[[17]] <- 0
  layers[[18]] <- 0
  for(i in seq_along(layers)){
    layers[[i]] <- data.frame(polygon = i-1, layer = layers[[i]])
  }
  layers <- do.call(rbind, layers)
  layers$keep <- 1
  data <- merge(data, layers, all = T)
  data$atoutput[is.na(data$keep)] <- NA
  return(data)
}

#' @export
#' @rdname helper_functions_user_defined
# Add NAs for missing layers!
remove_nas_layers <- function(data){
  layers <- readLines(file.path("z:", "Atlantis/Initial_conditions_generator/atlantisNCGen/NorthSea/boxTracers.csv"), n = 1)
  layers <- unlist(stringr::str_split(layers, pattern = ","))
  layers <- as.numeric(layers[2:length(layers)])
  layers <- lapply(layers-1, seq, from = 0)
  layers[[17]] <- 0
  layers[[18]] <- 0
  for(i in seq_along(layers)){
    layers[[i]] <- data.frame(polygon = i-1, layer = layers[[i]])
  }
  layers <- do.call(rbind, layers)
  layers <- subset(layers, !is.element(polygon, get_bboxes()))
  layers$keep <- 1
  data <- merge(data, layers, all = T)
  data <- subset(data, keep == 1)
  return(data)
}

#' @export
#' @rdname helper_functions_user_defined
plot_active_poylgon <- function(){
  wuwu <- data.frame(polygon = get_boxes()[!is.element(get_boxes(), get_bboxes())], active = 1)
  poly_plot <- read.table(file.path("z:", "Atlantis", "Polygons_and_bgm_new", "NSpolygons_new_channel_raw_R_Contour_test.txt"), sep = "", header = F)
  names(poly_plot) <- c("region", "long", "lat")
  plot <- ggplot2::ggplot() +
    ggplot2::geom_map(data = wuwu, ggplot2::aes(map_id = polygon, fill = active), map = poly_plot) +
    ggplot2::geom_path(data = poly_plot, ggplot2::aes(x = long, y = lat), colour = "black", lineend = "round") +
    ggplot2::facet_wrap(~ polygon, nrow = 1) +
    theme_standard() +
    ggplot2::theme(legend.position = "none", 
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(), 
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank())
  return(plot)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# *** 4) General Plotting Functions *** ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname helper_functions_user_defined
theme_standard <- function(font = 22, font_small = 18, font_tiny = 14, scale_font = 1, rot_xaxis_text = T, rot_strips_y = T){
  ggplot2::theme(
    text                = ggplot2::element_text(family="sans", size = font * scale_font),
    title               = ggplot2::element_text(hjust=.5), 
    axis.title.x        = ggplot2::element_text(hjust=.5),
    axis.title.y        = ggplot2::element_text(hjust=.5, vjust=0.3),
    #    axis.text  =  element_text(),  #	inherits from text   
    axis.text.x         = ggplot2::element_text(angle = ifelse(rot_xaxis_text, 45, 0), hjust = 1, size = scale_font * font_tiny, colour = "black"),
    axis.text.y         = ggplot2::element_text(size = font_tiny * scale_font, colour = "black"),
    axis.line           = ggplot2::element_blank(),
    #    axis.line.x  =	element_line(),	#	inherits from axis.line
    #    axis.line.y	=	element_line(),	#	inherits from axis.line
    #    axis.ticks   =	element_line(),	#	inherits from line
    #    axis.ticks.x	=	element_line(),	#	inherits from axis.ticks
    #    axis.ticks.y	=	element_line(),	#	inherits from axis.ticks
    #    axis.ticks.length	=	unit(),		
    #    axis.ticks.margin	=	unit(),		
    plot.margin         = grid::unit(c(1,1,1,1), "mm"),
    #    plot.background     =	element_rect(),	#	inherits from rect
    #    plot.title	         =	element_text(),	#	 inherits from title
    panel.grid          = ggplot2::element_blank(),
    panel.border        = ggplot2::element_rect(fill = NA, colour="black"),
    panel.background    = ggplot2::element_blank(),
    #    panel.margin  =	unit	,		
    #    panel.grid.major	=	element_line(),	#	inherits from panel.grid
    #    panel.grid.minor	=	element_line(),	#	inherits from panel.grid
    #    panel.grid.major.x	=	element_line(),	#	inherits from panel.grid.major
    #    panel.grid.major.y	=	element_line(),	#	inherits from panel.grid.major
    #    panel.grid.minor.x	=	element_line(),	#	inherits from panel.grid.minor
    #    panel.grid.minor.y	=	element_line(),	#	inherits from panel.grid.minor
    legend.position     = "bottom",
    legend.text         = ggplot2::element_text(size = font_tiny * scale_font),
    legend.key.width    = grid::unit(0.75, "cm"),
    legend.title        = ggplot2::element_text(size = font_small * scale_font),
    #    legend.background   =	element_rect(),	#	inherits from rect
    #    legend.margin       =	unit(),	
    #    legend.key          =	element_rect(fill = NULL, colour = NULL, size = NULL, linetype = NULL, color = NULL)	,	#	inherits from rect
    #    legend.key.size	   =	unit,	#	inherits from legend.key.size
    #    legend.key.height	 =	unit,	#	inherits from legend.key.size
    #    legend.text.align	 =	,	#	number from 0 (left) to 1 (right)
    #    legend.title.align	 =	,	#	number from 0 (left) to 1 (right)
    #     legend.direction	   =	"horizontal",
    #    legend.justification	=	,	#	center or two-element numeric vector
    legend.box	         =	"horizontal",
    #    legend.box.just	   =	,	#	top, "bottom", "left", or "right"
    strip.background    = ggplot2::element_blank(),
    strip.text          = ggplot2::element_text(size = font_small),
    strip.text.x        = ggplot2::element_text(size = scale_font * font_small),	#	inherits from strip.text
    strip.text.y	      = ggplot2::element_text(size = scale_font * font_small, angle = ifelse(rot_strips_y, 0, 90))	#	inherits from strip.text
  )
}

#' @export
#' @rdname helper_functions_user_defined
set_categorys <- function(data){
  cats <- list(data.frame(species = "baleen_whales", category_plot= "mammal"),
               data.frame(species = "toothed_whales", category_plot= "mammal"),
               data.frame(species = "seals", category_plot= "mammal"),
               data.frame(species = "seabirds", category_plot= "bird"),
               data.frame(species = "piscivorous_sharks", category_plot= "shark & ray"),
               data.frame(species = "other_sharks", category_plot= "shark & ray"),
               data.frame(species = "skate_ray", category_plot= "shark & ray"),
               data.frame(species = "cod", category_plot= "fish"),
               data.frame(species = "whiting", category_plot= "fish"),
               data.frame(species = "haddock", category_plot= "fish"),
               data.frame(species = "saithe", category_plot= "fish"),
               data.frame(species = "hake", category_plot= "fish"),
               data.frame(species = "blue_whiting", category_plot= "fish"),
               data.frame(species = "norway_pout", category_plot= "fish"),
               data.frame(species = "other_large_demersals", category_plot= "fish"),
               data.frame(species = "other_small_demersals", category_plot= "fish"),
               data.frame(species = "monkfish", category_plot= "fish"),
               data.frame(species = "gurnards", category_plot= "fish"),
               data.frame(species = "herring", category_plot= "fish"),
               data.frame(species = "sprat", category_plot= "fish"),
               data.frame(species = "mackerel", category_plot= "fish"),
               data.frame(species = "horse_mackerel", category_plot= "fish"),
               data.frame(species = "sandeel", category_plot= "fish"),
               data.frame(species = "plaice", category_plot= "fish"),
               data.frame(species = "dab", category_plot= "fish"),
               data.frame(species = "witch", category_plot= "fish"),
               data.frame(species = "sole", category_plot= "fish"),
               data.frame(species = "turbot", category_plot= "fish"),
               data.frame(species = "bass", category_plot= "fish"),
               data.frame(species = "red_mullet", category_plot= "fish"),
               data.frame(species = "small_pelagic_filterfeeders", category_plot= "fish"),
               data.frame(species = "squid", category_plot= "shrimp & squid"),
               data.frame(species = "pandalus", category_plot= "shrimp & squid"),
               data.frame(species = "crangon", category_plot= "shrimp & squid"),
               data.frame(species = "nephrops", category_plot= "shrimp & squid"),
               data.frame(species = "large_crabs", category_plot= "invertebrate"),
               data.frame(species = "small_epifauna", category_plot= "invertebrate"),
               data.frame(species = "sessile_epifauna", category_plot= "invertebrate"),
               data.frame(species = "epifaunal_macrobenthos", category_plot= "invertebrate"),
               data.frame(species = "infaunal_macrobenthos", category_plot= "invertebrate"),
               data.frame(species = "small_infauna", category_plot= "invertebrate"),
               data.frame(species = "gelantinous_zoo", category_plot= "invertebrate"),
               data.frame(species = "benthic_microflora", category_plot= "small stuff"),
               data.frame(species = "meiofauna", category_plot= "small stuff"),
               data.frame(species = "micro_zoo", category_plot= "small stuff"),
               data.frame(species = "meso_zoo", category_plot= "small stuff"),
               data.frame(species = "planktonic_microflora", category_plot= "small stuff"),
               data.frame(species = "diatoms", category_plot= "primary producer"),
               data.frame(species = "other_phytoplankton", category_plot= "primary producer"),
               data.frame(species = "labile_detritus", category_plot= "detritus"),
               data.frame(species = "refractory_detritus", category_plot= "detritus"),
               data.frame(species = "carrion", category_plot= "detritus"))
  cats <- do.call(rbind, cats)
  if(!is.null(intersect(names(data), names(cats)))) data <- dplyr::left_join(data, cats)
  return(data)
}


