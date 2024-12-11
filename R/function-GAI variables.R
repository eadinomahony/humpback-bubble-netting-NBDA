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

  mr <- read.csv("../data/ILV.csv",stringsAsFactor=FALSE) ; mr
  sits <- read.csv("../data/Sightings.csv",stringsAsFactor=FALSE) ; sits
  dyads <- read.csv("../data/dyads-0.csv",stringsAsFactor=FALSE) ; dyads

  sit.cutoff <- 10
  subd <- dyads[dyads$nA >= sit.cutoff & dyads$nB >= sit.cutoff,] ; nrow(subd)
  subd <- ai(subd)

}

#########################################################
#########################################################

#########################################################
#########################################################
# Spatial overlap =======================================
# Proportion of those years in which both were identified,
# that both were identified within 20km of swimming distance from the inner fjord.

geoverlap <- function(di,km.cutoff=15){
  #i = which.max(dyads$X)
  #km.cutoff = 20 # km
  #di <- dyads[i,]  ; di

  a <- di[2] ; a
  b <- di[3] ; b

  #a <- di$A ; a
  #b <- di$B ; b
  sita <- sits[grep(a,sits$id),] ; nrow(sita)
  sitb <- sits[grep(b,sits$id),] ; nrow(sitb)

  yrs <- unique(sita$year[which(sita$year %in% sitb$year)]) ; yrs
  geoverlap <- NA
  if(length(yrs)>0){
    kms <- c()
    y <- 3
    for(y in 1:length(yrs)){
      yi <- yrs[y] ; yi
      ya <- sita[sita$year==yi,] ; ya
      yb <- sitb[sitb$year==yi,] ; yb
      kmdf <- data.frame()
      for(a in 1:nrow(ya)){
        kmdf <- rbind(kmdf, data.frame(a=ya$km.out[a],b=yb$km.out))
      }
      kmdf$diff <- abs(kmdf$a - kmdf$b)  ; kmdf
      kms[y] <- min(kmdf$diff,na.rm=TRUE)
    }
    kms
    geoverlap <- length(which(kms <= km.cutoff)) / length(yrs)
  }
  return(geoverlap)
}

#########################################################
#########################################################

geo.wrapper <- function(dyads,sits){
  geos <- rep(NA,times=nrow(dyads))
  progbar <- 0
  for(i in 1:nrow(dyads)){
    geos[i] <- geoverlap(dyads[i,])
    #print(geos[i])
    newprog <- round((i / nrow(dyads))*100)
    if(newprog != progbar){
      progbar <- newprog
      print(paste("Geographic overlap :: ", progbar," % complete . . . "))
    }
  }
  return(geos)
}

#########################################################
#########################################################
# Temporal overlap ======================================
# custom SRI calculated on whether
# pairs were observed in the study area within 14 days of each
# other, within sampling periods of 60 days

# Simple ratio:
# N periods i and j associated /
# N periods i or j or both are seen

toverlap <- function(di,sits){
  sits$yrdec <- sits$year + (sits$doy/365) ; sits$yrdec
  dint <- 14/365 ; dint
  starts <- seq(2004,2020,by=dint) ; length(starts)
  ends <- c((starts[2:length(starts)]),(starts[length(starts)] + dint)) ; length(ends)
  bins <- data.frame(starts,ends) ; head(bins)

  #i=911
  #di <- dyads[i,]  ; di
  a <- di$A ; a
  b <- di$B ; b

  sita <- sits[grep(a,sits$id),] ; nrow(sita)
  sitb <- sits[grep(b,sits$id),] ; nrow(sitb)

  bins$a <- 0
  bins$b <- 0
  b <- 10
  for(b in 1:nrow(bins)){
    bini <- bins[b,] ; bini

    ain <- which(sita$yrdec >= bini$starts & sita$yrdec < bini$ends) ; ain
    if(length(ain)>0){bins$a[b] <- 1}

    bin <- which(sitb$yrdec >= bini$starts & sitb$yrdec <= bini$ends) ; bin
    if(length(bin)>0){bins$b[b] <- 1}
  }
  bins

  bins$sum <- bins$a + bins$b ; bins$sum
  sri <- length(which(bins$sum==2)) / length(which(bins$sum>0)) ; sri
  return(sri)
}

#toverlap(di)

#########################################################
#########################################################

toverlap.wrapper <- function(dyads,sits){
  tovers <- rep(NA,times=nrow(dyads))
  progbar <- 0
  for(i in 1:nrow(dyads)){
    tovers[i] <- toverlap(dyads[i,],sits)
    #print(geos[i])
    newprog <- round((i / nrow(dyads))*100)
    if(newprog != progbar){
      progbar <- newprog
      print(paste("Temporal overlap :: ", progbar," % complete . . . "))
    }
  }
  return(tovers)
}

#########################################################
#########################################################
# Gregariousness ========================================
# Sum of SRIs for A
# Sum of SRIs for B
# Multiply them
# Take log of product
# Whitehead and James 2015 recommend removing the SRI of the dyad of interest to avoid circularity

gregariousness <- function(dyads){
  dyads %>% head

  dyads %>%
    group_by(A) %>%
    # Get sum of SRI for A and B
    mutate(Asum = sum(sri)) %>%
    ungroup() %>%
    group_by(B) %>%
    mutate(Bsum = sum(sri)) %>%
    ungroup() %>%
    # Revise by subtracting each dyad's sri
    mutate(Asum = Asum - sri,
           Bsum = Bsum - sri) %>%
    mutate(GI = log10(Asum*Bsum)) %>%
    select(-Asum, -Bsum)
}

# Old way (slow)
gregariousness_old <- function(di,dyads){
  #dyads <- ai(dyads)
  #i=977
  #di <- dyads[i,]  ; di

  a <- di$A ; a
  b <- di$B ; b
  sriab <- di$sri ; sriab

  sria <- dyads$sri[grep(a,dyads$dyad)] ; sria
  srib <- dyads$sri[grep(b,dyads$dyad)] ; srib

  suma <- sum(sria) - sriab ; suma
  sumb <- sum(srib) - sriab ; sumb
  gi <- log10(suma*sumb) ; gi

  return(gi)
}

#gregariousness(di,dyads=dyads)

#########################################################
#########################################################

# Old way (slow)
greg.wrapper <- function(dyads){
  gregs <- rep(NA,times=nrow(dyads))
  progbar <- 0
  for(i in 1:nrow(dyads)){
    gregs[i] <- gregariousness(di=dyads[i,],dyads=dyads)
    newprog <- round((i / nrow(dyads))*100)
    if(newprog != progbar){
      progbar <- newprog
      print(paste("Gregariousness :: ", progbar," % complete . . . "))
    }
  }
  return(gregs)
}

#########################################################
#########################################################
# Not used  =============================================
# Gender / class similarity
# Social unit membership
# Kinship


