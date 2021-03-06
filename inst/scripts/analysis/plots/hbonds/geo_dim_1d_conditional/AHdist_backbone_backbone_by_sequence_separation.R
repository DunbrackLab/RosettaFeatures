# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "AHdist_backbone_backbone_by_squence_separation",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	geom.AHdist,
	CASE don_site.resNum - acc_site.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4' WHEN 5 THEN '5'
		ELSE 'long' END AS seq_sep
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	hbond.struct_id = geom.struct_id AND
	hbond.hbond_id =  geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND
	hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND
	hbond.acc_id = acc_site.site_id AND
	(acc_site.HBChemType == 'hbacc_PBA' OR don_site.HBChemType == 'hbdon_PBA') AND
	don_pdb.struct_id = hbond.struct_id AND don_pdb.site_id = hbond.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hbond.struct_id AND acc_pdb.site_id = hbond.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;"
f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	seq_sep = factor(f$seq_sep,
		levels = c("-4", "-3", "-2", "-1", "2", "3", "4", "5", "long"),
		labels = c("-4", "-3", "-2", "-1", "2", "3", "4", "5", "long")))

f <- na.omit(f, method="r")

plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=x, y=y)),
	geom_indicator(aes(indicator=counts), group=1),
	scale_y_continuous("FeatureDensity", breaks=c(1,3,5,7)),
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), breaks=c(1.6, 1.9, 2.2, 2.5)),
	coord_trans(limx=c(1.4,2.7), limy=c(0,7.5)))

dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep"),
	"AHdist", weight_fun = radial_3d_normalization)

plot_id <- "hbond_AHdist_backbone_backbone_by_sequence_separation"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ seq_sep ) +
	ggtitle("H-Bond Backbone-Backbone A-H Distance by SeqSep: DonRes - AccRes") +
	scale_y_continuous("FeatureDensity", limits=c(0,6), breaks=c(1,3,5)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.5))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
