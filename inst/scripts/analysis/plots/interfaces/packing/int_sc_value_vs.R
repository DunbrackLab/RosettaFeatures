# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)
library(plyr)
library(grid)

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "int_packing-sc_value_vs",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs Interface metrics such as packstat vs other metrics",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
    SELECT
      sc_value,
      packstat,
      dSASA,
      dG,
      dG_cross,
      delta_unsatHbonds,
      interface
    FROM
      interfaces
  "
  
  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())
  
  
  data = query_sample_sources(sample_sources, sele)
  
  #Scatterplots
  #sc_value vs packstat
  parts = list(
    geom_point(size=1.5, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  p <- ggplot(data=data, aes(x =sc_value, y=packstat)) + parts + 
    ggtitle("sc_value vs packstat") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("packstat", limit = c(0, 1.0))
  plot_field(p, "sc_value_vs_packstat_by_all", grid = sample_source ~ .)
  plot_field(p, "sc_value_vs_packstat_by_interface", grid=interface ~ sample_source)
  
  #sc_value vs dSASA
  p <- ggplot(data=data, aes(x =sc_value, y=dSASA)) + parts + 
    ggtitle("sc_value vs dSASA") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("Buried SASA")
  plot_field(p, "sc_value_vs dSASA_all", grid = sample_source ~ .)
  plot_field(p, "sc_value_vs_dSASA_by_interface", grid=interface ~ sample_source)
  
  data_rm_out <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$dG <= quantile(d2$dG, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$dG <= quantile(d2$dG, .10))) #Top 10 percent
  })
  
  #sc_value vs dG
  p <- ggplot(data = data_rm_out, aes(x=sc_value, y=dG)) + parts +
    ggtitle("sc_value_vs_dG") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "sc_value_vs_dG-(top_90_percent)-by_all", grid=sample_source ~ .)
  plot_field(p, "sc_value_vs_dG_(top_90_percent)-by_interface", grid=interface ~ sample_source)
  
  p <- ggplot(data = data_top, aes(x=sc_value, y=dG)) + parts +
    ggtitle("sc_value_vs_dG") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "top_10_percent_dG-sc_value_vs_dG_by_all", grid=sample_source ~ .)
  plot_field(p, "top_10_percent_dG-sc_value_vs_dG_by_interface", grid=interface ~ sample_source)
  
  #deltaUnsatHbonds vs sc_value
  p <- ggplot(data = data, aes(x=sc_value, y=delta_unsatHbonds)) + parts +
    ggtitle("sc_value vs Interface unsatisfied polar atoms") +
    scale_y_continuous("n") +
    scale_x_continuous("sc_value", limit = c(0, 1.0))
  plot_field(p, "sc_value_vs_delta_unsat_polars_by_all", grid=sample_source ~ .)
  plot_field(p, "sc_value_vs_delta_unsat_polars_by_interface", grid=interface ~ sample_source)
  

  #By Native
  
  sele <- "
    SELECT
  sc_value,
  packstat,
  dSASA,
  dG,
  dG_cross,
  delta_unsatHbonds,
  interface,
  natives.native as native
  FROM
  interfaces,
  natives,
  WHERE
  interfaces.struct_id = natives.struct_id
  "
  
  data_rm_out <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$dG <= quantile(d2$dG, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$dG <= quantile(d2$dG, .10))) #Top 10 percent
  })
  
  #sc_value vs dG
  p <- ggplot(data = data_rm_out, aes(x=sc_value, y=dG)) + parts +
    ggtitle("sc_value_vs_dG") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "sc_value_vs_dG(top_90_percent)-by_native_by_all", grid=sample_source ~ .)
  plot_field(p, "sc_value_vs_dG(top_90_percent)-by_native_interface", grid=interface ~ sample_source)
  
  p <- ggplot(data = data_top, aes(x=sc_value, y=dG)) + parts +
    ggtitle("sc_value_vs_dG") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "top_10_percent_dG-sc_value_vs_dG_by_native_by_all", grid=sample_source ~ .)
  plot_field(p, "top_10_percent_dG-sc_value_vs_dG_by_native_by_interface", grid=interface ~ sample_source)
  
  
  
})) # end FeaturesAnalysis