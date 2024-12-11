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
  head(events)
  #mr <- read.csv("../data/events.csv",stringsAsFactors=FALSE) ; nrow(mr)
  #head(mr)

}

#########################################################
#########################################################

known.mothers <- function(){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  #mr <- read.csv("../data/events.csv",stringsAsFactors=FALSE) ; nrow(mr)
  load('../data/events.rds')
  head(mr)

  ms <- grep("M",mr$stage) ; ms
  cs <- grep("C",mr$stage) ; cs
  mcs <- unique(c(ms,cs)) ; mcs
  mids <- mr$id[mcs] ; mids
  moms <- sort(unique(mids)); moms
  moms <- moms[moms != 'CALF']
  moms

  return(moms)
}



#########################################################
#########################################################
