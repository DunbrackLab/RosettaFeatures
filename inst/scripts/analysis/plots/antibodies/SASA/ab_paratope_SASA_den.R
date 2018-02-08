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
id = "ab_SASA-paratope_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Various statistics on the H3 Kink",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    paratope_SASA,
    paratope_hSASA,
    paratope_SASA - paratope_hSASA as paratope_pSASA
  FROM
    ab_metrics
    "
  
  data = query_sample_sources(sample_sources, sele)
 
  parts <- list(
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
  
  
  #Paratope SASA
  data_rm_out <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$paratope_SASA <= quantile(d2$paratope_SASA, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source), function(d2){
    subset(d2, subset=(d2$paratope_SASA <= quantile(d2$paratope_SASA, .10))) #Top 10 percent
  })
  
  group = c("sample_source")
  dens <- estimate_density_1d(data, group, c("paratope_SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR Paratope SASA")
  plot_field(p, "paratope_sasa_den")
  
  group = c("sample_source")
  dens <- estimate_density_1d(data_rm_out, group, c("paratope_SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR Paratope SASA")
  plot_field(p, "top_90_percent_paratope_sasa_den")

  group = c("sample_source")
  dens <- estimate_density_1d(data_top, group, c("paratope_SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR Paratope SASA")
  plot_field(p, "paratope_sasa_top_10_percent_den")
  
  #Natives
  
  sele = "
  SELECT
  paratope_SASA,
  paratope_hSASA,
  paratope_SASA - paratope_hSASA as paratope_pSASA,
  natives.native as native
  FROM
  ab_metrics,
  natives
  WHERE
  ab_metrics.struct_id = natives.struct_id
  "
  
  data = query_sample_sources(sample_sources, sele)
  
  data_rm_out <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$paratope_SASA <= quantile(d2$paratope_SASA, .90))) #Remove high energy outliers
  })
  
  data_top <- ddply(data, .(sample_source, native), function(d2){
    subset(d2, subset=(d2$paratope_SASA <= quantile(d2$paratope_SASA, .10))) #Top 10 percent
  })
  
  group = c("sample_source")
  dens <- estimate_density_1d(data_rm_out, group, c("paratope_SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR Paratope SASA")
  plot_field(p, "top_90_percent_paratope_sasa_den_by_native")
  
  group = c("sample_source")
  dens <- estimate_density_1d(data_top, group, c("paratope_SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR Paratope SASA")
  plot_field(p, "paratope_sasa_top_10_percent_den_by_native")
  
})) # end FeaturesAnalysis