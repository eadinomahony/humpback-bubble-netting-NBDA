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

  load('../data/events.rds')
  mr <- events

  SUMM <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  nrow(SUMM)

  id <- "BCX0711"

  idata <- filter2id(id,mr,summ) ; idata

}

#########################################################
#########################################################

ILV.context <- function(id){

  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  load('../data/events.rds')
  mr <- events


  # old code from 2020/2021
  # idkey <- read.csv("../data/id-key.csv",stringsAsFactors=FALSE) ; idkey
  # keyi <- idkey[idkey$id==id,] ; keyi
  # rakes <- lone.bnfe <- damage <- white <- stage <-  NA
  # if(nrow(keyi)>0){
  #   rakes <- keyi$kw.rakes ; rakes
  #   lone.bnfe <- keyi$lone.bnfe ; lone.bnfe
  #   damage <- keyi$damage ; damage
  #   white <- keyi$white ; white
  #   stage <- as.character(keyi$first.stage)
  # }

  context <- data.frame(stage,white,damage,rakes,lone.bnfe) ; context
  return(context)
}

#ILV.context(id)

#########################################################
#########################################################
