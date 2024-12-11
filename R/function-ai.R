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

  dyads <- read.csv("../data/dyads-1.csv",stringsAsFactor=FALSE) ; dyads

  SUMM <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  ntot <- nrow(SUMM) ; ntot

  #head(dyads)
  #di <- dyads[2,] ; di
  #xbar <- di$X / ntot ; xbar
  #yabar <- di$ya / ntot ; yabar
  #ybbar <- di$yb / ntot ; ybbar
  #yabbar <- 0

  # Calculate c/k (the proportion of compartments surveyed)
  # Compartment is equivalent to the group rule
  #ck <- xbar + yabbar + 0.5*(yabar + ybbar) ; ck

  # If we assume that c (the compartments surveyed) = T, our total number of encounters
  #C <- ntot
  #k <- C / ck ; k

  # w = uniform probability of identification
  # Eq. 5 in Weko 2018
  #numerator <- xbar + ((k-1)/(C-1))*yabbar ; numerator
  #denominator <- xbar + yabbar + ((yabar + ybbar)/2) ; denominator
  #w <- numerator / denominator ; w

  #GLECI(x=di$X,Ya=di$ya,Yb=di$yb,Yab=0,Ynull=di$ynull,w=w)
  #HWI(x=di$X,Ya=di$ya,Yb=di$yb,Yab=0)

}

#########################################################
#########################################################

ai <- function(dyads){
  # Association indices
  # half weight, simple ratio, and simple gleci (assuming w=0)

  library(assocInd)
  library(rootSolve)

  mr <-
    dyads %>%
    rowwise() %>%
    mutate(hwi = HWI(x=X, Ya=ya, Yb=yb, Yab=0)[1]) %>%
    mutate(sri = SRI(x=X, Ya=ya, Yb=yb, Yab=0)[1])

  # old way from 2020 -- slower but works
  #mr <- dyads
  #mr$hwi <- NA
  #mr$gleci <- NA
  #mr$sri <- NA
  #for(i in 1:nrow(mr)){
  #  dati <- mr[i,] ; dati
    #haii <- hai(x=dati$X, ya=dati$ya, yb=dati$yb)
    #mr$hai[i] <- haii
    #dati

  #  hwi <- HWI(x=dati$X,Ya=dati$ya,Yb=dati$yb,Yab=0) ; hwi
  #  mr$hwi[i] <- hwi[1]

    #totsits <- 3959
    #ynull <- totsits - (dati$ya + dati$yb + dati$X) ; ynull
    #gleci <- GLECI(x=dati$X,Ya=dati$ya,Yb=dati$yb,Yab=0,Ynull=ynull,w=0) ; gleci
    #mr$gleci[i] <- gleci[1]

  #  sri <- SRI(x=dati$X,Ya=dati$ya,Yb=dati$yb,Yab=0) ; sri
  #  mr$sri[i] <- sri[1]
  #}
  #head(mr)

  #par(mfrow=c(2,2)) ; breaks <- seq(0,1,by=.05)
  #hist(mr$hwi,breaks=breaks)
  #hist(mr$sri,breaks=breaks)
  #hist(mr$gleci,breaks=breaks)
  #par(mfrow=c(1,1))

  return(mr)
}


#########################################################
#########################################################
# Simulation summary function

gai.summ <- function(gai.obs){
  names(gai.obs)
  head(gai.obs)

  # Randomization columns
  rcol <- which(substr(names(gai.obs),1,1)=="b") ; rcol

  # Dyadic medians and 95% CIs
  med <- apply(gai.obs[,rcol],1,function(x){quantile(x,.5,na.rm=TRUE)}) ; med
  lci <- apply(gai.obs[,rcol],1,function(x){quantile(x,.025,na.rm=TRUE)}) ; lci #; hist(lci)
  uci <- apply(gai.obs[,rcol],1,function(x){quantile(x,.975,na.rm=TRUE)}) ; uci #; hist(uci)

  # Dyadic significance
  pu <- apply(gai.obs,1,function(x){length(which(x[rcol] >= x[19])) / B}) ; head(pu) #; hist(pu)
  pl <- apply(gai.obs,1,function(x){length(which(x[rcol] <= x[19])) / B}) ; head(pl) #; hist(pl)

  # Create test statistics
  gmean <- apply(gai.obs[,rcol],1,function(x){mean(x,na.rm=TRUE)}) #; hist(gai.obs$mean)
  gsd <- apply(gai.obs[,rcol],1,function(x){sd(x,na.rm=TRUE)}) #; hist(gai.obs$sd)
  gcv <- gsd / gmean
  names(gai.obs)

  # Add to dataframe
  gai.obs$med <- med
  gai.obs$lci <- lci
  gai.obs$uci <- uci
  gai.obs$pu <- pu
  gai.obs$pl <- pl
  gai.obs$mean <- gmean
  gai.obs$sd <- gsd
  gai.obs$cv <- gcv
  head(gai.obs)

  return(gai.obs)
}
