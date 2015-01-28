#' Plot Atlantis output (netcdf)
#' 
#' 
#' This function plots Atlantis output (netcdf)
#' @param model_path Character string of the ATLANTIS folder.
#' @param filename_output Character string of the general ATLANTIS output file. Usually "outputNorthSea.nc".
#' @param filename_prod Character string of the productivity ATLANTIS output file. Usually "outputNorthSeaPROD.nc".
#' @param filename_init Character string of the initial ATLANTIS file. Usually "initNorthSea.nc".
#' @param select_groups Vector of funtional groups which shall be plotted. Names have to match the ones used in the ncdf file. Check column "Name" in "functionalGroups.csv" as example.
#' @param output_path Character string of the destination the plots shall be saved into.
#' @param mv Character value giving the modelversion. E.g. "v.1.0.0".
#' @param ps Numeric value giving the size of the plots. Default = 14.
#' @param plots_per_row Integer value giving the number of plots in each row in faceted plots. Default is 8.
#' @param report Logical If TRUE status of plotting procedure is commented in detail.
#' @return 20 Plots as pdf in output_path! 

#' @details This functions plots nums, n and physical variables (Oxygen, Salinity...) for each functional group!
#' @keywords gen
#' @examples 
#' plot_atlantis <- function(model_path = file.path("z:", "Atlantis", "ATLANTIS NSmodel base"), filename_output = "outputNorthSea.nc", filename_prod = "outputNorthSeaPROD.nc", select_groups = get_groups(), output_path =  file.path("z:", "Atlantis", "ATLANTIS NSmodel Plots"))
#' @export

plot_atlantis_new <- function(model_path, 
                              filename_output = "outputNorthSea.nc", 
                              filename_prod = "outputNorthSeaPROD.nc", 
                              filename_init = "init_NorthSea.nc", 
                              select_groups, 
                              output_path, 
                              mv, 
                              ps = 14, 
                              plots_per_row = 8,
                              report = T){

  age_groups <- get_age_groups()
  select_age_groups <- select_groups[is.element(select_groups, age_groups)] 
  if(length(select_age_groups) == 0) stop("At least one age-structured group has to be selected!")
  
  select_other_groups <- select_groups[!is.element(select_groups, age_groups)]
  bps <- get_bps()
  physic_var <- get_physics()[-which(is.element(get_physics(), c("vflux", "eflux", "volume")))]
  bio_conv <- 5.7 * 20 / 1000000000
  
  if(report) print("*** Start: reading in data! ***")
  # "load_atlantis_ncdf" and "load_atlantis_ncdf_physics" are independent functions in seperate R-files!
  # NOTE: Data for new plots has to be added here if not already available!
  at_n       <- load_atlantis_ncdf(model_path, filename_output, select_groups,     "N",       bps)
  
  at_structn_l <- load_atlantis_ncdf(model_path, filename_output, select_age_groups, "StructN", bps, aggregate_layers = F)
  at_resn_l    <- load_atlantis_ncdf(model_path, filename_output, select_age_groups, "ResN"   , bps, aggregate_layers = F)
  at_nums_l    <- load_atlantis_ncdf(model_path, filename_output, select_age_groups, "Nums"   , bps, aggregate_layers = F)
  if(length(select_other_groups) >= 1){
    at_n_pools   <- load_atlantis_ncdf(model_path, filename_output, select_other_groups, "N"    , bps, aggregate_layers = F)
  }
  at_eat     <- load_atlantis_ncdf(model_path, filename_prod,   select_age_groups, "Eat",     bps)
  at_growth  <- load_atlantis_ncdf(model_path, filename_prod,   select_age_groups, "Growth",  bps)
  
  flux       <- load_atlantis_ncdf_physics(model_path, filename_output, c("eflux", "vflux"), aggregate_layers = F)
  physics    <- load_atlantis_ncdf_physics(model_path, filename_output, physic_var) 
  
  biomass_timeseries <- read.csv(file.path(model_path, "biomass_timeseries.csv"))  
  biomass_timeseries <- set_categorys(data = biomass_timeseries)
  
  if(report) print("*** Start: data transformations! ***")
  # Aggregate Layers for Nums, ResN, StructN
  at_resn <- at_resn_l %>%
    dplyr::group_by(species, polygon, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = mean(atoutput))
  at_structn <- at_structn_l %>%
    dplyr::group_by(species, polygon, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = mean(atoutput))
  at_nums <- at_nums_l %>%
    dplyr::group_by(species, polygon, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  
  # Aggregate Numbers! This is done seperately since numbers need to be summed!
  at_nums_age <- at_nums %>%
    dplyr::group_by(species, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput)) 
  at_nums_polygon <- at_nums %>%
    dplyr::group_by(species, polygon, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  at_nums_overview <- at_nums %>%
    dplyr::group_by(species, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(atoutput))
  
  # Calculate biomass for age-groups
  at_structn_l$biomass_ind <- (at_resn_l$atoutput + at_structn_l$atoutput) * at_nums_l$atoutput * bio_conv
  biomass <- at_structn_l %>%
    dplyr::group_by(species, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(biomass_ind)) %>%
    dplyr::mutate(model = "atlantis")
  # Biomass per ageclass
  biomass_ages <- at_structn_l %>%
    dplyr::group_by(species, agecl, time, category_plot) %>%
    dplyr::summarise(atoutput = sum(biomass_ind))
  # Calculate relative change to allow plotting on same y-axis!
#   biomass_ages <- datatrans_calibrate(data = biomass_ages)
  # Calculate biomass for non-age-groups
  if(length(select_other_groups) >= 1){
    vol <- load_atlantis_ncdf_physics(model_path, filename = filename_init, physic_variables = c("volume", "dz"), 
                                      aggregate_layers = F, load_init = T, remove_bboxes = T)
    vol <- reshape2::dcast(vol, polygon + layer ~ variable, value.var = "atoutput")
    at_n_pools <- dplyr::left_join(at_n_pools, vol)
    at_n_pools$biomass_ind <- with(at_n_pools, ifelse(species %in% bps, atoutput * volume / dz * bio_conv, atoutput * volume * bio_conv))
    biomass_pools <- at_n_pools %>%
      dplyr::group_by(species, time, category_plot) %>%
      dplyr::summarise(atoutput = sum(biomass_ind)) %>%
      dplyr::mutate(model = "atlantis")
    
    # Put everyhing together!
    biomass <- rbind(biomass, biomass_pools, biomass_timeseries)
  } else {
    biomass <- rbind(biomass, biomass_timeseries) 
  }
  
  # NOTE: New dataframes also have to be added here depending on the calculations needed!
  agg_age      <- lapply(list(at_structn, at_resn, at_eat, at_growth), mean_over_ages)
  agg_polygon  <- lapply(list(at_n), mean_over_polygons)
  agg_overview <- lapply(list(at_n, at_eat, at_growth), mean_overview)
  
  # Divide by initial value
  at_calibrate <- lapply(list(at_nums_age, agg_age[[1]], agg_age[[2]]), datatrans_calibrate)
    
  # Calculate agestructure per group and timestep!
  at_agestructure <- at_nums_age %>%
    dplyr::group_by(species, time) %>%
    dplyr::mutate(atoutput = atoutput / sum(atoutput))
  
  flux <- remove_nas_layers(data = flux)
  
  # Combine Numbers-data with other dataframes (age, overview, polygon) because the plotting routines work with all dataframes!
  output <- list("overview" = c(list(at_nums_overview, subset(biomass, model == "atlantis")), agg_overview),
                 "age"      = c(list(at_nums_age, biomass_ages), agg_age), 
                 "polygon"  = c(list(at_nums_polygon), agg_polygon)) 
    
  if(report) print("*** Start: Plotting! ***")  
  plots_overview  <- lapply(output$overview, plot_overview, plots_per_row = plots_per_row)  
  plots_age       <- lapply(output$age,      plot_age)
  plots_polygon   <- lapply(output$polygon,  plot_polygon)
  plots_calibrate <- lapply(at_calibrate,    plot_calibrate, plots_per_row = plots_per_row)  
  
  plot_agestructure <- ggplot2::ggplot(at_agestructure, ggplot2::aes(x = time, y = atoutput, fill = factor(agecl))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~ species, ncol = plots_per_row) +  
    theme_standard(scale_font = 0.8)
  
  plot_benchmark <- ggplot2::ggplot(data = biomass, ggplot2::aes(x = time, y = atoutput, colour = model)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ species, scales = "free_y", ncol = plots_per_row) + 
    theme_standard(scale_font = 0.8)
  
  plot_flux <- ggplot2::ggplot(flux, ggplot2::aes(x = time, y = atoutput, colour = variable, linetype = variable)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(layer ~ polygon, scales = "free") +
    theme_standard(font_tiny = 10)
        
  plot_physics <- ggplot2::ggplot(physics, ggplot2::aes(x = time, y = atoutput)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(variable ~ polygon, scales = "free") +
    theme_standard(font_tiny = 10)
  
  # Add Labels for each row in age/polygon plots!
  for(i in seq_along(plots_age)){
    plots_age[[i]] <- add_label(plot = plots_age[[i]], plot_data = output$age[[i]])
  }
  for(i in seq_along(plots_polygon)){
    plots_polygon[[i]] <- add_label(plot = plots_polygon[[i]], plot_data = output$polygon[[i]])
  }
      
  # Combine all plots to add labels and save them!
  # NOTE: The order of the plots is defined here! Listplots are prioritized to singleplots
  # NOTE: Data for new plots has to be added here --> y_labels, plotnames, plotwidths, plotheights
  # 1) Overview  --> Numbers, Nitrogen, Eat, Growth
  # 2) Age       --> Numbers, StructN, ResN, Eat, Growth
  # 3) Polygon   --> Numbers, Nitrogen
  # 4) Calibrate --> Numbers, StructN, ResN
  # 5) Single    --> Agestructure, Benchmark, Physics, Flux
  plots <- c(plots_overview, 
             plots_age, 
             plots_polygon, 
             plots_calibrate,
             list(plot_agestructure, plot_benchmark, plot_physics, plot_flux))
  
  y_labels <- list("Numbers", "Relative Biomass", "Nitrogen", "Eat", "Growth",                                # overview
                   "Numbers", "Relative Biomass", "Structural nitrogen", "Reserve nitrogen", "Eat", "Growth", # ageclasses
                   "Numbers", "Nitrogen",                                                 # polygon
                   expression("Numbers/Numbers"[initial]),                                # calibrate
                   expression("StructN/StructN"[initial]),
                   expression("ResN/ResN"[initial]),
                   "Numbers [%]", "Biomass [t]", "Physical parameters", "Flux")           # Single
  
  plots <- add_axis_label(plots, axis_label = y_labels, which_axis = "y")
  
  plotnames <- paste0(mv,
                      c(paste0("_at_output_", "overview_", c("nums", "biomass", "nitrogen", "eat", "growth")),
                        paste0("_at_output_", "ages_",     c("nums", "biomass", "structn", "resn", "eat", "growth")),
                        paste0("_at_output_", "polygons_", c("nums", "nitrogen")),
                        paste0("_at_calibrate_", "reference_", c("nums", "structn", "resn")),
                        paste0("_at_calibrate_", "agestructure"),
                        paste0("_at_benchmark_", "biomass"), 
                        paste0("_at_output_", c("physics", "flux"))),
                      ".pdf")

  plotwidths <- c(scale_width(select_age_groups),                                                      # overview nums
                  rep(scale_width(select_groups), times = 2),                                          # overview_biomass and n
                  rep(scale_width(select_age_groups), times = 2),                                      # overview_eat/growth
                  rep(1.6 * ps, times = length(plots_age)),                                            # ageclasses
                  rep(2.4 * ps, times = length(plots_polygon)),                                        # polygon
                  rep(scale_width(select_age_groups), times = length(plots_calibrate) + 1),            # calibrate + agestructure
                  ifelse(length(select_groups) >= 28, scale_width(select_groups), scale_width(1:28)),  # benchmark
                  rep(1.6 * ps, times = 2))                                                            # physics + flux
    
  plotheights <- c(scale_height(select_age_groups),                                                       # overview_nums
                   rep(scale_height(select_groups), times = 2),                                           # overview_biomass and n
                   rep(scale_height(select_age_groups), times = 2),                                       # overview_eat/growth
                   rep((0.08 * length(select_age_groups) + 2 * 0.1) * ps, times = length(plots_age)),     # ageplots
                   (0.08 * length(select_age_groups) + 0.1) * ps,                                         # polygon_nums
                   (0.08 * length(select_groups) + 0.1) * ps,                                             # polygon_n 
                   rep(scale_height(select_age_groups), times = length(plots_calibrate) + 1),             # calibrate + agestructure
                   ifelse(length(select_groups) >= 28, scale_height(select_groups), scale_height(1:28)),  # benchmark
                   0.84 * ps,                                                                             # physics
                   0.68 * ps)                                                                             # flux

  if(report) print("*** Start: Saving to pdf! ***")
    
  for(i in seq_along(plots)){
    pdf(file = file.path(output_path, plotnames[i]), width = plotwidths[i], height = plotheights[i])
    print(plots[[i]])   
    dev.off()
    if(report) print(paste("*** plot", i ,"of", length(plots), "done!***"))
  }
}
                  







