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

known.bnfers <- function(){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  load('../data/events.rds')
  #mr <- read.csv("../data/events.csv",stringsAsFactors=FALSE) ; nrow(mr)
  head(mr)
  #mr$b <- paste(mr$bhvr,mr$bhvr2,mr$bhvr3,sep="-")
  bnfs <- grep("BN",mr$bhvr) ; bnfs
  bnfers <- mr$id[bnfs] ; bnfers
  bnfers <- sort(unique(bnfers)); bnfers
  bnfers <- bnfers[bnfers != 'CALF']
  bnfers

  return(bnfers)
}



#########################################################
#########################################################
