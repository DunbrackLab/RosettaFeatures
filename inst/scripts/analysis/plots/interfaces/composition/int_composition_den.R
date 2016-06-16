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
id = "int_composition_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic composition of the interfaces, restypes, etc",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures","ResidueFeatures", "ResidueTypesFeatures", "PdbDataFeatures"),

#Because R is stupid and can't find this function.



run=function(self, sample_sources, output_dir, output_formats){

  #Aromatic Composition
  sele <- "
  SELECT
    aromatic_fraction,
    interface_nres,
    interface,
    side
  FROM
    interface_sides
  "

  capwords <- function(s, strict = FALSE)
  {
    cap <- function(s) paste(toupper(substring(s, 1, 1)), {
      s <- substring(s, 2)
      if (strict)
        tolower(s)
      else s
    }, sep = "", collapse = " ")
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }
  
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

  plot_field = function(p, plot_id, grid = NULL, ssLegend=T){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(ssLegend){
      if(nrow(sample_sources) <= 3){
        p <- p + theme(legend.position="bottom", legend.direction="horizontal")
      }
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  fields = c("aromatic_fraction")
  data = query_sample_sources(sample_sources, sele)
  
  for(field in fields){
    fieldSP = unlist(strsplit(field, split="_"))
    parts = list(plot_parts, xlab("fraction"))
    
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides_by_all", sep="_"), grid=side ~ .)
    
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
  }
  
  parts = list(plot_parts, scale_x_continuous("number of interface residues"))
  
  field = "interface_nres"
  group = c("sample_source", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface nres")
  plot_field(p, paste(field, "den_sides_by_all", sep="_"), grid=side ~ .)
  
  group = c("sample_source", "interface", "side")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Interface nres")
  plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
  
  #Restype composition - Overall interface 
  sele <-"
  SELECT
    interface_residues.interface as interface,
    residues.name3 as restype,
    residue_type.name1 as restype1,
    interface_residues.SASA_int as SASA_int,
    interface_residues.dSASA as dSASA,
    interface_residues.dSASA - interface_residues.dSASA_sc as dSASA_bb,
    interface_residues.dSASA_sc as dSASA_sc,
    interface_residues.dhSASA as dhSASA,
    interface_residues.dG as dG,
    
    interface_residues.relative_dSASA_fraction as dSASA_fraction,
    interface_residues.struct_id as struct_id
  FROM
    residues,
    interface_residues,
    residue_type
  WHERE
    interface_residues.struct_id == residues.struct_id and
    interface_residues.resNum == residues.resNum and
    residues.name3==residue_type.name3
  "
  res_data = query_sample_sources(sample_sources, sele)
  
  

  ##Histogram - only plot residues that have a dSASA fraction > 5 % -  change this for sidechains once we have that data

  ##### Typical way is not working, so we will have to do it manually. #####
  
  #Restype Composition - Classical, not working!
#  p <- ggplot(data=res_data, aes(x=restype1)) + 
#    geom_bar(position="dodge", aes(y = ..density.., fill=sample_source), binwidth=1)+ 
#    theme_bw() +
#    ggtitle("Interface ResType Composition") +
#    scale_y_continuous(label=percent)
#  plot_field(p, "restype_composition_by_all_test") 
#  plot_field(p, "restype_composition_by_interface_test", grid=interface ~ .)
    
  get_percent <- function(d) {
    d_per <- ddply(d, .(sample_source, interface, struct_id), function(per_struct_id){
      d_per_restype <- ddply(per_struct_id, .(restype1), function(per_restype){
        #print(head(per_restype))
        perc = length(per_restype$restype1)/length(per_struct_id$struct_id) 
        df = data.frame(perc = perc)
      })
    })
    d_per
  }
    
        
  #Restype Composition
  
  p <- ggplot(data=get_percent(res_data), aes(x=restype1)) + 
    geom_bar(position="dodge", stat="identity", aes(y=perc, fill=sample_source))+ 
    theme_bw() +
    ggtitle("Interface ResType Composition") +
    scale_y_continuous(label="percent") +
    ylab("% of Sample Source")
  plot_field(p, "restype_composition_by_all") 
  plot_field(p, "restype_composition_by_interface", grid=interface ~ .)
      
      
    
  
})) # end FeaturesAnalysis