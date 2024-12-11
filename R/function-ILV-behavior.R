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

  idata <- filter2id(id,mr,SUMM) ; idata

}

#########################################################
#########################################################

ILV.behavior <- function(idata){

  #idata <- filter2id(id,mr,summ) ; idata
  id <- idata$id
  mri <- idata$mri
  years <- idata$years
  summi <- idata$sits
  sits <- idata$sitids ; sits

  # Behavior variables

  # Motherhood
  moms <- known.mothers() ; moms
  if(id %in% moms){mom <- 1}else{mom <- 0}
  mom

  momyears <- grep("M",mri$stage) ; momyears
  momyears <- which(mri$stage[momyears] != "FEMALE")
  momyears <- unique(mri$year[momyears]) ; momyears

  # Group size
  grps <- as.numeric(as.character(summi$n)) ; grps
  grp_mean <- mean(grps,na.rm=TRUE) ; grp_mean

  # BNFE
  bnfers <- known.bnfers() ; bnfers
  if(id %in% bnfers){bnfe <- 1}else{bnfe <- 0}
  bnfe

  (summi <- summi %>% arrange(ymd))
  (summbnf <- summi %>% filter(grepl('BN', bhv)))
  if(nrow(summbnf)>0){
    bnfyears <- unique(summbnf$year) ; bnfyears
    bnfsits <- summbnf$sit ; bnfsits
  }else{
    bnfyears <- c()
    bnfsits <- c()
  }
  bnfyears
  bnfsits

  # First year seen BNFE
  bnf.year1 <- bnf.sit1 <- NA
  if(length(bnfyears)>0){
    bnf.year1 <- min(bnfyears)
    bnf.sit1 <- bnfsits[1]
  }
  bnf.year1
  bnf.sit1

  # BNF demonstrator? (doing this in 2004?)
  bnf.demon <- 0
  if(!is.na(bnf.year1)){
    if(as.character(bnf.year1)=="2004"){bnf.demon <- 1}
  }
  bnf.demon

  # BNFE rate
  (bnf.rate <- length(bnfsits) / length(sits))

  # BNF Solo
  (bnfsolo <- summbnf %>% filter(n == 1))
  if(nrow(bnfsolo)>0){
    soloyears <- unique(bnfsolo$year) ; soloyears
    solosits <- bnfsolo$sit ; solosits
  }else{
    soloyears <- c()
    solosits <- c()
  }
  if(nrow(bnfsolo)>0){bnfsolo <- 1}else{bnfsolo <- 0}
  bnfsolo
  soloyears
  solosits

  solo.year1 <- NA
  if(length(soloyears)>0){solo.year1 <- min(soloyears)}

  solo.sit1 <- NA
  if(length(solosits)>0){solo.sit1 <- solosits[1]}

  # Feeding rate (other than BNF)
  (bnfsitid <- summbnf$sit)
  (summfe <- summi %>% filter(grepl('F', bhv)))
  (summnobnf <- summfe[! summfe$sit %in% bnfsitid, ])
  (fe.rate <- nrow(summnobnf)/ length(sits))

  # Old variables not used in 2023
  # Sea lion
  #sl.rate <- length(which(summi$sealion != 0)) / length(sits) ; sl.rate

  # Robust
  #robust.rate <- length(which(summi$robust != 0)) / length(sits) ; robust.rate

  # Kelp roll
  #kelp.rate <- length(which(summi$kelp != 0)) / length(sits) ; kelp.rate

  # Orca
  #orca.rate <- length(which(summi$orca != 0)) / length(sits) ; orca.rate

  # Rest
  #rest.rate <- length(which(summi$rest != 0)) / length(sits) ; rest.rate

  # Social
  #social.rate <- length(which(summi$social != 0)) / length(sits) ; social.rate

  ilvs <- data.frame(n.years=length(years),
                     n.sits=length(sits),
                     grp_mean,
                     mom,
                     n.calves = length(momyears),
                     bnfe,
                     bnf.year1,
                     bnf.sit1,
                     bnf.demon,
                     bnf.rate,
                     bnfesolo=bnfsolo,
                     solo.year1,
                     solo.sit1,
                     fe.rate,
                     #social.rate,
                     #rest.rate,
                     #sl.rate,
                     #kelp.rate,
                     #orca.rate,
                     years=paste(years,collapse=" "),
                     momyears=paste(momyears,collapse=" "),
                     bnfyears=paste(bnfyears,collapse=" "),
                     bnfsits=paste(bnfsits,collapse=" "),
                     soloyears=paste(soloyears,collapse=" "),
                     solosits=paste(solosits,collapse=" ")
                     )
  ilvs

  return(ilvs)
}



#########################################################
#########################################################
