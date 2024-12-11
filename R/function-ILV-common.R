#########################################################
#########################################################
# Functions
# Eric Keen, v. June 2020

#########################################################
#########################################################
# Test data for function development

testdata <- FALSE
if(testdata){

  # Set working directory to folder that this R file is in
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  source("00 Source functions.R")

  mr <- read.csv("../data/hw-master.csv",stringsAsFactors=FALSE) ; nrow(mr)
  head(mr)
  nrow(mr)

  summ <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  nrow(summ)

  sits <- summ

  id <- "BCX0711"

}

#########################################################
#########################################################

filter2id <- function(id,mr,summ){
  mri <- mr[mr$id==id,] ; mri
  sitids <- unique(mri$groupid) ; sitids
  sits <- summ[summ$sit%in%sitids,] ; sits
  years <- unique(sits$year) ; years
  return(list(id=id,mri=mri,sitids=sitids,sits=sits,years=years))
}

#filter2id(id=id,mr=mr,summ=summ)


#########################################################
#########################################################

sampling.occasions <- function(sits){
  sitdates <- sits$ymd
  #sitdates <- substr(sits$sit,1,8) ; sitdates
  udates <- unique(sitdates) ; udates
  return(list(n=length(udates),dates=udates))
}


#########################################################
#########################################################

get.ids <- function(siti){
  j=1
  ids <- c()
  for(j in 1:nrow(siti)){
    sitj <- siti[j,] ; sitj
    idj <- sitj$id ; idj
    idj <- strsplit(as.character(idj)," ")[[1]] ; idj
    ids <- c(ids,idj)
  }
  ids
  return(ids)
}
