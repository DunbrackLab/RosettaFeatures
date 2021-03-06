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
id = "int_SASA_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic dSASA and SASA information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

 plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  sele = "
  SELECT
    dSASA,
    dSASA_sc,
    dSASA - dSASA_sc as dSASA_bb,
    dhSASA,
    dhSASA_sc,
    dhSASA - dhSASA_sc as dhSASA_bb,
    dhSASA_rel_by_charge,
    aromatic_dSASA_fraction,
    interface,
    side
  FROM
    interface_sides
  ORDER BY dSASA DESC
  "
  
  #Polar fraction - from Ben Strange's Paper
  

  parts = list(plot_parts, xlab("SASA"))
  
  data = query_sample_sources(sample_sources, sele)
  
  data$polar_fraction = (data$dSASA - data$dhSASA)/data$dSASA
  field = "polar_fraction"
  
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("dSASA(Polar)/dSASA") +
    ggtitle("Polar dSASA Fraction")
  plot_field(p, "dSASA_polar_fraction_den_by_all", grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("dSASA(Polar)/dSASA") +
    ggtitle("Polar dSASA Fraction")
  plot_field(p, "dSASA_polar_fraction_den_by_interface", grid=side~interface)
  
  #print(data)
  
  #Backbone SASA may not be interesting, but I want I still want to know for now.
  #JAB - Commenting out hydrophibic sasa.  Not very useful from my experience and it makes too many plots.
  fields = c("dSASA", "dSASA_bb", "dSASA_sc", "dhSASA", "dhSASA_bb", "dhSASA_sc", "dhSASA_rel_by_charge")
  fields = c("dSASA")
  
  data_rm_out <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$dSASA <= quantile(d2$dSASA, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$dSASA <= quantile(d2$dSASA, .10))) #Top 10 percent
  })
  
  for (field in fields){

    parts = list(plot_parts, scale_x_continuous("SASA"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_90_percent_den_sides_by_all", sep="_"), grid=side ~ .)
  
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_90_percent_den_sides","by_interface", sep="_"), grid=side~interface)
    
    parts = list(plot_parts, scale_x_continuous("SASA"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data_top,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_10_percent_den_sides_by_all", sep="_"), grid=side ~ .)
    
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data_top,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_10_percent_den_sides","by_interface", sep="_"), grid=side~interface)
  
  }
  
  
  ####  Means  #########
  fields = c("dSASA")
  for (field in fields){

  avgs <- ddply(data, .(sample_source, side, field), function(d2){
    data.frame(m = mean(d2[,field]), std_dev = sd(d2[,field]), m_top10 = mean(d2[1:10,field]), std_dev_top_10 = sd(d2[1:10,field]), top = d2[1,field])
  })
    
  p <- ggplot(data=avgs ) +
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field,"Average", sep=" "))+
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
  plot_field(p, paste("avg_sides_by_all", field, sep = "_"), grid=side ~ .)
  
  #Average Top 10
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "Average Best 10",sep=" ")) +
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
  plot_field(p, paste("avg_sides_top_10_by_all", field, sep="_"), grid=side ~ .)
  
  #Best
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "top", sep=" ")) +
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
  plot_field(p, paste("sides_top_by_all", field, sep = "_"), grid=side ~ .)
  
  avgs <- ddply(data, .(sample_source, side, field, interface), function(d2){
    data.frame(m = mean(d2[,field]), std_dev = sd(d2[,field]), m_top10 = mean(d2[1:10,field]), std_dev_top_10 = sd(d2[1:10,field]), top = d2[1,field])
  })
      
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field,"Average", sep=" ")) +
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
  plot_field(p, paste("avg_sides_by_interface", field, sep="_"), grid=side ~ interface)
  
  #Average Top 10
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "Average Best 10",sep=" ")) +
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17)) +
    ylab(field)
  plot_field(p, paste("avg_sides_top_10_by_interface", field, sep="_"), grid=side ~ interface)
  
  #Best
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top , fill=sample_source)) +
    theme_bw() +
    ggtitle(paste("Buried", field, "top", sep=" ")) +
    scale_x_discrete(labels=function(x) abbreviate(x, minlength=17)) +
    ylab(field)
  plot_field(p, paste("sides_top_by_interface", field, sep = "_"), grid=side ~ interface)
    
  } #End each side
  
  
  #Fractions
  field = "aromatic_dSASA_fraction"
  parts = list(plot_parts, scale_x_continuous("fraction", limit=c(0, 1.0)))
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction") +
  plot_field(p, "dSASA_aromatic_fraction_den_by_all", grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic dSASA Fraction")
  plot_field(p, "dSASA_aromatic_fraction_den_by_interface", grid=side~interface)
  
  #Natives
  
  sele = "
  SELECT
  dSASA,
  dSASA_sc,
  dSASA - dSASA_sc as dSASA_bb,
  dhSASA,
  dhSASA_sc,
  dhSASA - dhSASA_sc as dhSASA_bb,
  dhSASA_rel_by_charge,
  aromatic_dSASA_fraction,
  interface,
  side,
  natives.native as native
  FROM
  interface_sides,
  natives
  WHERE
  interface_sides.struct_id = natives.struct_id
  ORDER BY dSASA DESC
  "
  data = query_sample_sources(sample_sources, sele)
  
  fields = c("dSASA")
  
  data_rm_out <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$dSASA <= quantile(d2$dSASA, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$dSASA <= quantile(d2$dSASA, .10))) #Top 10 percent
  })
  
  for (field in fields){
    
    parts = list(plot_parts, scale_x_continuous("SASA"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_90_percent_den_sides_by_native_by_all", sep="_"), grid=side ~ .)
    
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data_rm_out,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_90_percent_den_sides_by_native","by_interface", sep="_"), grid=side~interface)
    
    parts = list(plot_parts, scale_x_continuous("SASA"))
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data_top,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_10_percent_den_sides_by_native_by_all", sep="_"), grid=side ~ .)
    
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data_top,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste("Buried", field, sep=" "))
    plot_field(p, paste(field, "top_10_percent_den_sides_by_native","by_interface", sep="_"), grid=side~interface)
    
  }
  
  
})) # end FeaturesAnalysis