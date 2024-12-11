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

  sits <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  nrow(sits)

  id <- "BCX0711"

  idata <- filter2id(id,mr,summ) ; idata

}

#########################################################
#########################################################
# Site fidelity metrics
#########################################################
# Individual Occurrence (IO)

occurrence <- function(id,sits){
  siti <- sits[grep(id,sits$id),] ; siti
  ioi <- NA

  if(nrow(siti)>0){
  Ti <- sampling.occasions(sits)$n ; Ti # define sampling occasions as the number of days of effort
  ci <- sampling.occasions(siti)$n ; ci
  ioi <- (ci-1) / (Ti - 1) ; ioi
  }

  return(ioi)
}

#occurrence(id=id,sits=sits)

#########################################################
#########################################################
# Permanence (IT)

permanence <- function(id,sits){
  (siti <- sits[grep(id,sits$id),])
  iti <- NA

  if(nrow(siti)==1){
    iti <- 0
  }

  if(nrow(siti)>1){
    siti
    Fi <- diff(range(siti$doy, na.rm=TRUE)) ; Fi
    Ftot <- diff(range(sits$doy, na.rm=TRUE)) ; Ftot
    iti <- Fi/Ftot
  }
  iti

  return(iti)
}

#permanence(id=id,sits=sits)

#########################################################
#########################################################
# Periodicity

periodicity <- function(id,sits){
  siti <- sits[grep(id,sits$id),]
  iti <- NA

  if(nrow(siti)==1){
    iti <- 0
  }

  if(nrow(siti)>1){
    Fi <- diff(range(siti$doy, na.rm=TRUE)) ; Fi
    ci <- sampling.occasions(siti)$n ; ci
    iti <- (Fi/(ci-1))^(-1)
  }
  iti

  return(iti)
}

#periodicity(id=id,sits=sits)

#########################################################
#########################################################
# Standardised Site Fidelity Index (SSFI)

ssfi <- function(id,sits){
  ssfii <- NA

  IT <- permanence(id,sits) ; IT
  It <- periodicity(id,sits) ; It

  ssfii <- 2/((1/IT) + (1/It)) ; ssfii

  return(ssfii)
}

#ssfi(id,sits)

#########################################################
#########################################################
# Annual capture rate (overall)

loyalty.overall <- function(id,sits){
  all.years <- unique(sits$year) ; all.years
  siti <- sits[grep(id,sits$id),] ; siti
  id.years <- unique(siti$year) ; id.years
  ao <- length(id.years) / length(all.years) ; ao
  return(ao)
}

#loyalty.overall(id,sits)

#########################################################
#########################################################
# Loyalty since first capture

loyalty.sincefirst <- function(id,sits){
  all.years <- unique(sits$year) ; all.years
  siti <- sits[grep(id,sits$id),] ; siti
  id.years <- unique(siti$year) ; id.years
  first.year <- min(id.years) ; first.year
  possible.years  <- all.years[which(all.years >= first.year)] ; possible.years
  ao <- length(id.years) / length(possible.years) ; ao
  return(ao)
}

#loyalty.sincefirst(id,sits)

#########################################################
#########################################################
