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
#library(grid)

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "total_score",
author = "Jared Adolf-Bryfogle", 
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  structure_scores.struct_id as struct_id, 
  structure_scores.score_value as total_score, 
  structure_scores.score_type_id as score_type 
FROM 
	structure_scores, 
  score_types 
  
WHERE 
	score_types.score_type_name='total_score' AND 
  structure_scores.score_type_id = score_types.score_type_id 
ORDER BY score_value;"


data <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = data,
  ids = c("sample_source"),
  variable = "total_score")

plot_id <- "total_score"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x, y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Total Rosetta Energy") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

f <- ddply(data, .(sample_source), function(d2){
  data.frame(total_score = d2[1:20,]$total_score)
})

data_rm_out <- ddply(data, .(sample_source), function(d2){
  subset(d2, subset=(d2$total_score <= quantile(d2$total_score, .90))) #Remove high energy outliers
})

data_top <- ddply(data, .(sample_source), function(d2){
  subset(d2, subset=(d2$total_score <= quantile(d2$total_score, .10))) #Top 10 percent
})

dens <- estimate_density_1d(f, ids = c("sample_source"), variable = "total_score")

plot_id <- "total_score_top_20"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Rosetta Structure Score - Top 20") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(
  data = data_rm_out,
  ids = c("sample_source"),
  variable = "total_score")

plot_id <- "total_score_top_90_percent"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Rosetta Structure Score") +
  labs(x="Rosetta Energy Units") +
  scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(
  data = data_top,
  ids = c("sample_source"),
  variable = "total_score")

plot_id <- "total_score_top_10_percent"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Rosetta Structure Score") +
  labs(x="Rosetta Energy Units") +
  scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#Averages Scoring
avgs <- ddply(data, .(sample_source), function(d2){
  data.frame(m = mean(d2$total_score), std_dev = sd(d2$total_score), m_top10 = mean(d2[1:10,]$total_score), std_dev_top_10 = sd(d2[1:10,]$total_score), top = d2[1,]$total_score)
})

p <- ggplot(data=avgs ) + 
  geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
  #geom_errorbar(aes(ymin = m-std_dev, ymax=m+std_dev) +
  theme_bw() +
  ggtitle("Average Total Score") +
  ylab("REU") +
  scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
save_plots(self, "avg_total_score", sample_sources, output_dir, output_formats)
  
#Avg Top 10 Scoring
p <- ggplot(data=avgs ) + 
  geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m_top10 , fill=sample_source)) +
  theme_bw() +
  ggtitle("Average Best 10 Score") +
  ylab("REU") +
  scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
save_plots(self, "avg_top_10_total_score", sample_sources, output_dir, output_formats)

#Top Scoring
p <- ggplot(data=avgs ) + 
  geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= top, fill=sample_source)) +
  theme_bw() +
  ggtitle("Best Score") +
  ylab("REU") +
  scale_x_discrete(labels=function(x) abbreviate(x, minlength=17))
save_plots(self, "best_total_score", sample_sources, output_dir, output_formats)


#By Native

sele <-"
SELECT
structure_scores.struct_id as struct_id, 
natives.native as native,
structure_scores.score_value as total_score, 
structure_scores.score_type_id as score_type 
FROM 
structure_scores, 
score_types,
natives

WHERE 
score_types.score_type_name='total_score' AND 
structure_scores.score_type_id = score_types.score_type_id AND
natives.struct_id = structure_scores.struct_id
ORDER BY score_value;"


data <- query_sample_sources(sample_sources, sele)

data_rm_out <- ddply(data, .(sample_source, native), function(d2){
  
  subset(d2, subset=(d2$total_score <= quantile(d2$total_score, .90))) #Remove high energy outliers
})

data_top <- ddply(data, .(sample_source, native), function(d2){
  subset(d2, subset=(d2$total_score <= quantile(d2$total_score, .10))) #Top 10 percent
})

f <- ddply(data, .(sample_source, native), function(d2){
  data.frame(total_score = d2[1:5,]$total_score)
})



dens <- estimate_density_1d(f, ids = c("sample_source"), variable = "total_score")

plot_id <- "total_score_top_5_by_native"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Rosetta Structure Score - Top 20") +
  labs(x="Rosetta Energy Units") +
  scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(
  data = data_rm_out,
  ids = c("sample_source"),
  variable = "total_score")

plot_id <- "total_score_top_90_percent_by_native"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Rosetta Structure Score") +
  labs(x="Rosetta Energy Units") +
  scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(
  data = data_top,
  ids = c("sample_source"),
  variable = "total_score")

plot_id <- "total_score_top_10_percent_by_native"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x, y, colour=sample_source), size=1.4) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Rosetta Structure Score") +
  labs(x="Rosetta Energy Units") +
  scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
