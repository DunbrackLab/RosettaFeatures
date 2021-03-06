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
id = "ab_dSASA-ab_cdr_den",
author = "Jared Adolf-Bryfogle",
brief_description = "CDR Sasas",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  
  #if ("FALSE" %in% opt$options$include_cdr4 & "FALSE" %in% opt$options$cdr4_only){
  sele = "
  SELECT
    ag_ab_dSASA as dSASA,
    ag_ab_dSASA_sc as dSASA_sc,
    ag_ab_dhSASA as dhSASA,
    ag_ab_dhSASA_sc as dhSASA_sc,
    ag_ab_dhSASA_rel_by_charge as dhSASA_rel_by_charge,
    struct_id,
    CDR,
    length
  FROM
    cdr_metrics
  WHERE
    dSASA > 0 and
    CDR NOT LIKE '%Proto%'
    "
  
  data = query_sample_sources(sample_sources, sele)
 
  parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_wrap(facets=grid, ncol=3)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  
  
  data$polar_fraction = (data$dSASA - data$dhSASA)/data$dSASA

  #Avg CDR dSASA
  means_save <- ddply(data, .(sample_source, CDR), function(data){
    data.frame(sample_source = data$sample_source, CDR = data$CDR, m=mean(data$dSASA))
  })
  p <- ggplot(data=means_save, na.rm=T) +
    geom_bar(position="dodge", stat='identity', aes(x=CDR, y=m, fill=sample_source)) +
    ggtitle("Average Antigen Buried SASA") +
    xlab("CDR") +
    ylab("dSASA") +
    theme_bw()
  plot_field(p, "avg_cdr_dSASA_hist")
  
  
  #Avg CDR dSASA by length
  means_save <- ddply(data, .(sample_source, CDR, length), function(data){
    data.frame(sample_source = data$sample_source, CDR = data$CDR, m=mean(data$dSASA))
  })
  p <- ggplot(data=means_save, na.rm=T) +
    geom_bar(position="dodge", stat='identity', aes(x=length, y=m, fill=sample_source)) +
    ggtitle("Average Antigen Buried SASA") +
    xlab("CDR Length") +
    ylab("dSASA") +
    theme_bw()
  plot_field(p, "avg_cdr_dSASA_hist_by_length", ~CDR)
  
  
  

  
  #CDR dSASA
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(data, group, c("dSASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("dSASA") +
    ggtitle("Antigen Buried Solvent Accessible Surface Area")
  plot_field(p, "cdr_dSASA_den", ~CDR)
  
})) # end FeaturesAnalysis