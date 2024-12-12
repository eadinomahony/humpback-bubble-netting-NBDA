#########################################################
#########################################################
# The diffusion of cooperative and solo bubble net feeding in Canadian Pacific humpback whales
# v. Éadin O'Mahony, Oct 2024

#set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load custom functions:
source("00 Source functions.R")

#install.packages("assocInd")
library(dplyr)
library(tidyverse)
library(assocInd)
library(lubridate)
library(devtools)
library(rootSolve)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#download NBDA package from github:
#devtools::install_github("whoppitt/NBDA")
library(NBDA)

################################################################################
################################################################################

### THE WHOLE STUDY HINGES UPON THREE DATASETS: 1) events.csv; 2) ILV.csv and 3) assoc-final.csv
### These are all the sightings and each individual's ILVs 
### (sex information has been gleaned from blow sampling data)
### 1) Each row is one sighting
### 2) Each row is one individual
### 3) assoc-final.csv is a matrix with dyadic SRI values

### NBDA models using NBDA package: https://github.com/whoppitt/NBDA

################################################################################
################################################################################

