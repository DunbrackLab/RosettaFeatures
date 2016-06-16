#!/usr/bin/env Rscript

#NOTE: You must run this script with SUDO!!


required_packages <- c("devtools", "roxygen2")

if (length(setdiff(required_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(required_packages, rownames(installed.packages())), repos="http://cran.rstudio.com/")  
}


devtools::install_github("momeara/RosettaFeatures")