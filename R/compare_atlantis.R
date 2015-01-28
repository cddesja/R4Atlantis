#' Plot Atlantis output (netcdf)
#' 
#' 
#' This function plots Atlantis output (netcdf)
#' @param model_paths Vector of character strings of the ATLANTIS folders.
#' @param runs Vector if character strings giving the ID of each model run.
#' @param filename_output Character string of the general ATLANTIS output file. Usually "outputNorthSea.nc".
#' @param filename_prod Character string of the productivity ATLANTIS output file. Usually "outputNorthSeaPROD.nc".
#' @param filename_init Character string of the initial ATLANTIS file. Usually "initNorthSea.nc".
#' @param select_groups Vector of funtional groups which shall be plotted. Names have to match the ones used in the ncdf file. Check column "Name" in "functionalGroups.csv" as example.
#' @param output_path Character string of the destination the plots shall be saved into.
#' @param mv Character value giving the modelversion. E.g. "v.1.0.0".
#' @param ps Numeric value giving the size of the plots. Default = 14.
#' @param plots_per_row Integer value giving the number of plots in each row in faceted plots. Default is 8.
#' @param report Logical If TRUE status of plotting procedure is commented in detail.
#' @return 14 Plots as pdf in output_path! 

#' @details This functions plots nums, n and physical variables (Oxygen, Salinity...) for each functional group!
#' @keywords gen
#' @examples 
#' plot_atlantis <- function(model_path = file.path("z:", "Atlantis", "ATLANTIS NSmodel base"), filename_output = "outputNorthSea.nc", filename_prod = "outputNorthSeaPROD.nc", select_groups = get_groups(), output_path =  file.path("z:", "Atlantis", "ATLANTIS NSmodel Plots"))
#' @export
#' 
compare_atlantis_new <- function(model_paths,
                                 runs,
                                 filename_output = "outputNorthSea.nc", 
                                 filename_prod = "outputNorthSeaPROD.nc", 
                                 filename_init = "init_NorthSea.nc", 
                                 select_groups, 
                                 output_path, 
                                 mv, 
                                 ps = 14, 
                                 plots_per_row = 8,
                                 report = T){
  if(length(model_paths) == 1) stop("At least 2 model runs have to be selected.")
  if(length(model_paths) != length(runs)) stop("Number of model_paths and runs do not match!")
  
  age_groups <- get_age_groups()
  select_age_groups <- select_groups[is.element(select_groups, age_groups)] 
  if(length(select_age_groups) == 0) stop("At least one age-structured group has to be selected!")
  
  select_other_groups <- select_groups[!is.element(select_groups, age_groups)]
  bps <- get_bps()
  physic_var <- get_physics()[-which(is.element(get_physics(), c("vflux", "eflux", "volume")))]
  bio_conv <- 5.7 * 20 / 1000000000  
  
  biomass_timeseries <- read.csv(file.path(model_paths[[1]], "biomass_timeseries.csv"))  
  biomass_timeseries <- set_categorys(data = biomass_timeseries)
  biomass_timeseries$run <- "base"
  
  result_runs <- list()
  for(i in seq_along(model_paths)){
    if(report) print(paste("*** Start run", i, "of", length(model_paths)))                             
    if(report) print("*** Start: reading in data! ***")
    
    # "load_atlantis_ncdf" and "load_atlantis_ncdf_physics" are independent functions in seperate R-files!
    # NOTE: Data for new plots has to be added here if not already available!
    at_n       <- load_atlantis_ncdf(model_paths[i], filename_output, select_groups, "N", bps)
    
    at_structn_l <- load_atlantis_ncdf(model_paths[i], filename_output, select_age_groups, "StructN", bps, aggregate_layers = F)
    at_resn_l    <- load_atlantis_ncdf(model_paths[i], filename_output, select_age_groups, "ResN"   , bps, aggregate_layers = F)
    at_nums_l    <- load_atlantis_ncdf(model_paths[i], filename_output, select_age_groups, "Nums"   , bps, aggregate_layers = F)
    if(length(select_other_groups) >= 1){
      at_n_pools   <- load_atlantis_ncdf(model_paths[i], filename_output, select_other_groups, "N", bps, aggregate_layers = F)
    }
    at_eat     <- load_atlantis_ncdf(model_paths[i], filename_prod,   select_age_groups, "Eat",     bps)
    at_growth  <- load_atlantis_ncdf(model_paths[i], filename_prod,   select_age_groups, "Growth",  bps)
    
    physics    <- load_atlantis_ncdf_physics(model_paths[i], filename_output, physic_var) 
      
    if(report) print("*** Start run: data transformations! ***")
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
      vol <- load_atlantis_ncdf_physics(model_paths[i], filename = filename_init, physic_variables = c("volume", "dz"), 
                                        aggregate_layers = F, load_init = T, remove_bboxes = T)
      vol <- reshape2::dcast(vol, polygon + layer ~ variable, value.var = "atoutput")
      at_n_pools <- dplyr::left_join(at_n_pools, vol)
      at_n_pools$biomass_ind <- with(at_n_pools, ifelse(species %in% bps, atoutput * volume / dz * bio_conv, atoutput * volume * bio_conv))
      biomass_pools <- at_n_pools %>%
        dplyr::group_by(species, time, category_plot) %>%
        dplyr::summarise(atoutput = sum(biomass_ind)) %>%
        dplyr::mutate(model = "atlantis")
      # Put everyhing together!
      biomass <- rbind(biomass, biomass_pools)
    }
     
  
    # NOTE: New dataframes also have to be added here depending on the calculations needed!
    agg_age      <- lapply(list(at_structn, at_resn, at_eat, at_growth), mean_over_ages)
    agg_polygon  <- lapply(list(at_n), mean_over_polygons)
    agg_overview <- lapply(list(at_n, at_eat, at_growth), mean_overview) 
  
    result_runs[[i]] <- c(list(at_nums_overview, subset(biomass, model == "atlantis")), agg_overview,  # overview
                          list(at_nums_age, biomass_ages), agg_age,                                    # age
                          list(at_nums_polygon), agg_polygon,                                          # polygon
                          list(biomass),                                                               # benchmark
                          list(physics))                                                               # physics
  
    if(i == length(model_paths)) print(" *** All runs read in! ***")
  }
  
  if(any(sapply(result_runs, length) != 15)) stop("Output from model runs not equal!")
  
  # rbind output from different runs add column runs to dataframe!
  combine_results <- list()
  for(i in seq_along(result_runs[[1]])){
    dummy <- list()
    for(j in seq_along(runs)){
      result_runs[[j]][[i]]$run <- runs[j]
      dummy[[j]] <- result_runs[[j]][[i]]
    }
    combine_results[[i]] <- do.call(rbind, dummy)
  }
    
  # Combine Numbers-data with other dataframes (age, overview, polygon) because the plotting routines work with all dataframes!
  output <- list("overview"  = list(combine_results[[1]], combine_results[[2]], combine_results[[3]], combine_results[[4]], combine_results[[5]]),
                 "age"       = list(combine_results[[6]], combine_results[[7]], combine_results[[8]], combine_results[[9]], combine_results[[10]], combine_results[[11]]), 
                 "polygon"   = list(combine_results[[12]], combine_results[[13]]))
  
  if(report) print("*** Start: Plotting! ***")  
  plots_overview  <- lapply(output$overview, plot_overview, plots_per_row = plots_per_row)  
  plots_age       <- lapply(output$age,      plot_age)
  plots_polygon   <- lapply(output$polygon,  plot_polygon)
  
  biomass <- rbind(combine_results[[14]], biomass_timeseries)
  plot_benchmark <- ggplot2::ggplot(data = biomass, ggplot2::aes(x = time, y = atoutput, colour = model, linetype = run)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ species, scales = "free_y", ncol = plots_per_row) + 
    theme_standard(scale_font = 0.8)
  
  plot_physics <- ggplot2::ggplot(combine_results[[15]], ggplot2::aes(x = time, y = atoutput, linetype = run)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(variable ~ polygon, scales = "free") +
    theme_standard(font_tiny = 10)
  
  # Add model run as linetype!
  plots_overview <- add_model_run(plots = plots_overview)
  plots_age      <- add_model_run(plots = plots_age)
  plots_polygon  <- add_model_run(plots = plots_polygon)

  # Add Labels for each row in age/polygon plots!
  for(i in seq_along(plots_age)){
    plots_age[[i]] <- add_label(plot = plots_age[[i]], plot_data = output$age[[i]], path = model_paths)
  }
  for(i in seq_along(plots_polygon)){
    plots_polygon[[i]] <- add_label(plot = plots_polygon[[i]], plot_data = output$polygon[[i]], path = model_paths)
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
             list(plot_benchmark, plot_physics))
    
  y_labels <- list("Numbers", "Relative Biomass", "Nitrogen", "Eat", "Growth",                                # overview
                   "Numbers", "Relative Biomass", "Structural nitrogen", "Reserve nitrogen", "Eat", "Growth", # ageclasses
                   "Numbers", "Nitrogen",                                                                     # polygon
                   "Numbers [%]", "Biomass [t]", "Physical parameters")                                       # Single
  
  plots <- add_axis_label(plots, axis_label = y_labels, which_axis = "y")
  
  plotnames <- paste0(mv,
                      c(paste0("_at_compare_", "overview_", c("nums", "biomass", "nitrogen", "eat", "growth")),
                        paste0("_at_compare_", "ages_",     c("nums", "biomass", "structn", "resn", "eat", "growth")),
                        paste0("_at_compare_", "polygons_", c("nums", "nitrogen")),
                        paste0("_at_compare_", "biomass"), 
                        "_at_compare_physics"),
                      ".pdf")

  plotwidths <- c(scale_width(select_age_groups),                                                      # overview nums
                  rep(scale_width(select_groups), times = 2),                                          # overview_biomass and n
                  rep(scale_width(select_age_groups), times = 2),                                      # overview_eat/growth
                  rep(1.6 * ps, times = length(plots_age)),                                            # ageclasses
                  rep(2.4 * ps, times = length(plots_polygon)),                                        # polygon
                  ifelse(length(select_groups) >= 28, scale_width(select_groups), scale_width(1:28)),  # benchmark
                  rep(1.6 * ps, times = 1))                                                            # physics + flux
    
  plotheights <- c(scale_height(select_age_groups),                                                       # overview_nums
                   rep(scale_height(select_groups), times = 2),                                           # overview_biomass and n
                   rep(scale_height(select_age_groups), times = 2),                                       # overview_eat/growth
                   rep((0.08 * length(select_age_groups) + 2 * 0.1) * ps, times = length(plots_age)),     # ageplots
                   (0.08 * length(select_age_groups) + 0.1) * ps,                                         # polygon_nums
                   (0.08 * length(select_groups) + 0.1) * ps,                                             # polygon_n 
                   ifelse(length(select_groups) >= 28, scale_height(select_groups), scale_height(1:28)),  # benchmark
                   0.84 * ps)                                                                             # physics

  if(report) print("*** Start: Saving to pdf! ***")
    
  for(i in seq_along(plots)){
    pdf(file = file.path(output_path, plotnames[i]), width = plotwidths[i], height = plotheights[i])
    print(plots[[i]])   
    dev.off()
    if(report) print(paste("*** plot", i ,"of", length(plots), "done!***"))
  }
}
                  







