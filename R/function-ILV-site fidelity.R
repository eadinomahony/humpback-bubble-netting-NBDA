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
  head(mr)

  summ <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  nrow(summ)

  id <- "BCX0711"
  idata <- filter2id(id,mr,summ) ; idata
}

#########################################################
#########################################################

ILV.fidelity <- function(idata,summ){

  #idata <- filter2id(id,mr,summ) ; idata
  id <- idata$id
  mri <- idata$mri ; mri
  years <- sort(idata$years) ; years
  sits <- idata$sits ; sits

  # Overall patterns

  earliest <- min(sits$doy,na.rm=TRUE) ; earliest
  latest <- max(sits$doy,na.rm=TRUE) ; latest
  SF.total <- loyalty.overall(id,sits=summ) ; SF.total
  SF.since <- loyalty.sincefirst(id,sits=summ) ; SF.since

  # Year-averaged patterns
  overview <- data.frame()
  y=1
  for(y in 1:length(years)){
    yi <- years[y] ;  yi
    sity <- sits[sits$year==yi,] ; sity
    summy <- summ[summ$year==yi,] ; nrow(summy)

    arrival <- min(sity$doy,na.rm=TRUE) ; arrival
    departure <- max(sity$doy,na.rm=TRUE) ; departure
    min.stay <- departure - arrival ; min.stay
    if(min.stay==0){min.stay <- 1} ; min.stay

    #sits <- sity
    IO <- occurrence(id=id,sits=summy) ; IO
    IT <- permanence(id=id,sits=summy) ; IT
    It <- periodicity(id=id,sits=summy) ; It
    SSFI <- ssfi(id,sity) ; SSFI

    # Make dataframe row for year
    dfy <- data.frame(arrival,
                      departure,
                      min.stay,
                      IO,
                      IT,
                      It,
                      SSFI
                      )
    dfy
    #print(dfy)
    # Add to annual dataframe
    overview <- rbind(overview,dfy)
  }
  overview
  oversumm <- apply(overview,2,function(x){mean(x,na.rm=TRUE)}) ; oversumm
  oversumm <- t(data.frame(oversumm))

  ilvs <- data.frame(first=years[1],
                     SF.total,
                     SF.since,
                     earliest,
                     latest,
                     oversumm
                     )
  row.names(ilvs) <- NULL
  ilvs

  return(list(ilv=ilvs,overview=overview))
}

#ILV.fidelity(idata,summ)

#########################################################
#########################################################
