#!/usr/bin/env Rscript


#NOTE: You must run this script with SUDO!!

#install.packages("devtools", repos="http://cran.rstudio.com/")

required_packages <- c("devtools", "roxygen2")

if (length(setdiff(required_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(required_packages, rownames(installed.packages())), repos="http://cran.rstudio.com/")  
}

library(methods)

#devtools::document() # if you changed function signatures
devtools::build()

devtools::install_local(".")
