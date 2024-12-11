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
  
  mr <- read.csv("../data/hw-master.csv",stringsAsFactors=FALSE) ; nrow(mr)
  head(mr)
  nrow(mr)
  
  SUMM <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  nrow(SUMM)

  id <- "BCX0711"
  
  idata <- filter2id(id,mr,summ) ; idata
  
}

#########################################################
#########################################################

ILV.geography <- function(idata){
  
  #idata <- filter2id(id,mr,summ) ; idata
  id <- idata$id
  mri <- idata$mri
  years <- idata$years
  sits <- idata$sits
  
  # Spatial patterns
  summi <- sits
  head(summi)
  
  x <- mean(summi$x,na.rm=TRUE) ; x
  if(!is.finite(x)){x <- NA} ; x
  xsd <- sd(summi$x,na.rm=TRUE) ; xsd
  
  y <- mean(summi$y,na.rm=TRUE) ; y
  if(!is.finite(y)){y <- NA} ; y
  ysd <- sd(summi$y,na.rm=TRUE) ; ysd
  
  km.out <- mean(summi$km.out,na.rm=TRUE)
  if(!is.finite(km.out)){km.out <- NA} ; km.out
  
  km.sd <- sd(summi$km.out,na.rm=TRUE)
  if(!is.finite(km.sd)){km.sd <- NA} ; km.sd
  
  km.min <- min(summi$km.out,na.rm=TRUE)
  if(!is.finite(km.min)){km.min <- NA} ; km.min
  
  km.max <- max(summi$km.out,na.rm=TRUE)
  if(!is.finite(km.max)){km.max <- NA} ; km.max
  
  ilvs <- data.frame(x,
                     xsd,
                     y,
                     ysd,
                     km.out,
                     km.sd,
                     km.min,
                     km.max
                     )
    
  ilvs
  
  return(ilvs)
}



#########################################################
#########################################################
