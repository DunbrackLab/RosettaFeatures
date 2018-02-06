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

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "int_SASA-dSASA_vs",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic dSASA and SASA information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    dSASA,
    dSASA_hphobic,
    dSASA_polar,
    dG,
    interface
  FROM
    interfaces"
  
  int_data = query_sample_sources(sample_sources, sele)
 

  
  #dSASA sides
  sele = "
  SELECT
    dG,
    dSASA,
    dSASA_sc,
    dSASA - dSASA_sc as dSASA_bb,
    dhSASA,
    dhSASA_sc,
    dhSASA - dhSASA_sc as dhSASA_bb,
    dSASA-dhSASA as dpSASA,
    dhSASA_rel_by_charge,
    aromatic_dSASA_fraction,
    interface_nres,
    interface,
    side
  FROM
    interface_sides
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
  
  data = query_sample_sources(sample_sources, sele)
  #print(data)
  
  
  
  #ScatterPlots
  
  
  parts = list(
    geom_point(size=.75),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(size=.5),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  #JAB - commenting out the hydrophobic dSASA.  Unclear if this is useful or not.  I don't think it is very much.
  #fields = c("dSASA", "dSASA_bb", "dSASA_sc", "dhSASA", "dhSASA_bb", "dhSASA_sc", "dhSASA_rel_by_charge")
  
  #fields = c("dSASA", "dSASA_bb", "dSASA_sc")
  fields = c("dSASA")
  
  data_rm_out = subset(int_data, subset=(int_data$dG <= quantile(int_data$dG, .90))) #Remove high energy outliers
  data_top = subset(int_data, subset=(int_data$dG <= quantile(int_data$dG, .10))) #Top 10 percent
  
  for (f in fields){
    
    #dSASA vs dG
    

    
    p <- ggplot(data=data_rm_out, aes(x = f, y = dG, colour=sample_source)) + parts +
      ggtitle(paste(f,"vs dG"))
    
    plot_field(p, paste(f, "vs_dG_by_all", sep="_"))
    plot_field(p, paste(f, "vs_dG_by_interface", sep="_"), grid=interface ~ .)
  
    data_rm_out$e_density = data_rm_out$dG/data_rm_out$dSASA
    field = c("dSASA")
    p <- ggplot(data = data_rm_out, aes(x = dSASA, y = e_density, colour=sample_source)) + parts +
      ggtitle(paste(field, "vs Interface energy density")) +
    plot_field(p, paste("control", field, "vs_energy_density", sep="_"), grid=side ~ .)
  
    #Top 10 Percent
    p <- ggplot(data=data_top, aes(x = f, y = dG, colour=sample_source)) + parts +
      ggtitle(paste(f,"vs dG"))
    
    plot_field(p, paste(f, "vs_dG_top_10_percent_by_all", sep="_"))
    plot_field(p, paste(f, "vs_dG_top_10_percent_by_interface", sep="_"), grid=interface ~ .)
    
    data_top$e_density = data_top$dG/data_top$dSASA
    field = c("dSASA")
    p <- ggplot(data = data_top, aes(x = dSASA, y = e_density, colour=sample_source)) + parts +
      ggtitle(paste(field, "vs Interface energy density")) +
      plot_field(p, paste("control", field, "vs_energy_density_top_10_percent", sep="_"), grid=side ~ .)
    
  }
  
})) # end FeaturesAnalysis