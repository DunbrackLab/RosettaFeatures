#!/usr/bin/env Rscript

#NOTE: You must run this script with SUDO!!


#install.packages("devtools", repos="http://cran.rstudio.com/")
devtools::document() # if you changed function signatures
devtools:build()

devtools::install_local(".")