#########################################################
#########################################################
# Generalized Affiliation Index Function

# Eric Keen, v. June 2020

#########################################################
#########################################################
# Test data for function development

testdata <- FALSE
if(testdata){

  # Set working directory to folder that this R file is in
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

  dsv <- read.csv("../data/GAI/dyads-10-structural.csv",stringsAsFactors=FALSE)  ; head(dsv)
  obs <- dsv

  hw <- read.csv("../data/hw-master.csv",stringsAsFactor=FALSE) ; hw
  hwi <- hw
  hwi$id <- hwi$id[sample(1:nrow(hwi),nrow(hwi),replace=FALSE)]

  gai.obs <- read.csv("../data/GAI/dyads-10-GAI-siml.csv", stringsAsFactors=FALSE)

}

#########################################################
#########################################################
# Matrix Formatting Function

mat2edgl <- function(M){
  #M <- My
  test <- lower.tri(M) # Extract matrix lower triangle
  # Extract matrix cells
  weight <- M[test]
  # Extract cells ids
  tmp <- which(test, arr.ind = TRUE)
  # Create an edgelist of actor, receiver and interactions weights
  DF <- cbind("from" = colnames(M)[tmp[, 1]], "to" = colnames(M)[tmp[, 2]])
  DF <- data.frame(DF, weight)
  names(DF) <- c("A","B","var")
  head(DF)
  return(DF)
}

#########################################################
#########################################################
# Base GAI function

gai.calc2 <- function(hwi,obs,sit.floor=5){
  nrow(hwi)
  dyads <- fast.dyads2(hwi,sit.floor=sit.floor)
  dyads$dyad <- as.character(dyads$dyad)
  dyads <- ai(dyads)
  nrow(dyads)

  # Remove navel-gazers and duplicated rows
  dyads <- dyads[dyads$A != dyads$B,]
  nrow(dyads)
  dsort <- apply(dyads,1,function(x){paste(sort(c(as.character(x[2]),as.character(x[3]))),collapse="-")})
  head(dsort)
  dyads$dsort <- dsort
  dyads <- dyads[duplicated(dyads$dsort),]
  nrow(dyads)

  # Match the new dyads to the observed dsv
  dyads <- dyads[order(dyads$dyad),]
  obs <- obs[order(obs$dyad),]
  dyads %>% nrow
  obs %>% nrow
  head(dyads$dyad)
  head(obs$dyad)
  #test <- data.frame(new=dyads$dyad,old=obs$dyad) ; head(test)
  #which(test$new != test$old)

  # Add same old structural variables to new dyads df
  #dyads$geoverlap <- obs$geoverlap
  #dyads$temporal <- obs$temporal

  # Recalculate gregariousness
  dyads <- gregariousness(dyads)
  #gregs <- greg.wrapper(dyads)
  #dyads$GI <- gregs ; head(dyads)
  head(dyads)

  infs <- which(!is.finite(dyads$GI)) ; infs
  if(length(infs)>0){dyads$GI[infs] <- -5}
  range(dyads$GI)

  # Remove incomplete cases
  dyads <- dyads[complete.cases(dyads),]
  nrow(dyads)

  # GLM
  res <- glm(sri ~ GI + geoverlap + temporal,
             data=dyads,
             family = binomial(link = "logit"))$residuals
  head(res,50)
  length(res)
  dyads$gnum <- res
  head(dyads)

  # Calculate GAI
  dyads$gai <- dyads$gnum / dyads$yabx
  head(dyads)
  #plot(gai~sri,data=dyads,xlim=c(0,.3))

  #hist(dyads$gai)
  range(dyads$gai,na.rm=TRUE)
  #length(which(dyads$gai>0)) / nrow(dyads)
  return(dyads)
}

#gai.calc2(hwi,obs=dsv)

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

#########################################################
#########################################################
# Signifiance tests and summaries

gai.test <- function(gai.obs){
  gai <- gai.obs[!is.na(gai.obs$gai),] ; head(gai)

  # Randomization columns
  rcol <- which(substr(names(gai.obs),1,1)=="b") ; rcol

  par(mfrow=c(2,2))
  hist(gai$gai,main="Observed GAI",breaks=20,col="grey")

  print("Observed: N dyads = ")
  print(nrow(gai))
  print(".")

  print("Observed: N individuals = ")
  print(length(unique(c(gai$A,gai$B))))
  print(".")

  print("Observed: Median GAI = ")
  print(median(gai$gai))
  print(".")

  print("Observed: Mean GAI = ")
  print(mean(gai$gai))
  print(".")

  print("Observed: SD GAI = ")
  print(sd(gai$gai))
  print(".")

  print("Observed: Proportion of GAIs that are negative = ")
  print(length(which(gai$gai < 0)) / nrow(gai))
  print(".")

  print("Observed: Max GAI:")
  print(gai[which.max(gai$gai),1:22])
  print(".")

  # Significance tests
  means <- apply(gai.obs[,rcol],2,function(x){mean(x,na.rm=TRUE)})
  hist(means,main="Mean of randomizations",breaks=20,col="grey")

  sds <- apply(gai.obs[,rcol],2,function(x){sd(x,na.rm=TRUE)})
  hist(sds,main="SD of randomizations",breaks=20,col="grey")

  cvs <- sds / means
  hist(cvs,main="CV of randomizations",breaks=20,col="grey")

  # Test for social preferences
  print("P value: Are there social preferences?")
  cvobs <- sd(gai$gai) / mean(gai$gai)  ; cvobs
  print(
    length(which(cvs >= cvobs)) / nrow(gai)
  )
  print(".")

  # Test for shorter-term than expected preferences
  print("P value: Are there preferences weaker than expected, on average?")
  print(
    length(which(means <= mean(gai$gai))) / nrow(gai)
  )
  print(".")

  # Test for presence of strong social preference relationships
  print("P value: Are there strong social preferences in the case of some dyads?")
  print(
    length(which(sds >= sd(gai$gai))) / nrow(gai)
  )
  print(".")

  # Percentage of preferred dyadic associations
  print("What percentage of dyadic preferences are significant?")
  #sigs <- which(gai$pu >= 0.975) ; sigs
  sigs <- which(gai$pl >= 0.95) ; sigs
  print(
    length(sigs) / nrow(gai)
  )
  print(".")

}

#########################################################
#########################################################
# SRI test

sri.test <- function(gai.obs){
  gai <- gai.obs[!is.na(gai.obs$sri),] ; head(gai)

  # Randomization columns
  rcol <- which(substr(names(gai.obs),1,1)=="b") ; rcol

  par(mfrow=c(2,2))
  hist(gai$sri,main="Observed SRI",breaks=20,col="grey")

  print("Observed: N dyads = ")
  print(nrow(gai))
  print(".")

  print("Observed: N individuals = ")
  print(length(unique(c(gai$A,gai$B))))
  print(".")

  print("Observed: Median SRI = ")
  print(median(gai$sri))
  print(".")

  print("Observed: Mean SRI = ")
  print(mean(gai$sri))
  print(".")

  print("Observed: SD SRI = ")
  print(sd(gai$sri))
  print(".")

  print("Observed: Median SRI (non-zero)= ")
  print(median(gai$sri[which(gai$sri>0)]))
  print(".")

  print("Observed: Mean SRI (non-zero)= ")
  print(mean(gai$sri[which(gai$sri>0)]))
  print(".")

  print("Observed: SD SRI (non-zero)= ")
  print(sd(gai$sri[which(gai$sri>0)]))
  print(".")

  print("Observed: Percentage > 0")
  print(length(which(gai$sri > 0)) / nrow(gai))
  print(".")

  nonzeroes <- gai[gai$sri > 0,]
  print("Observed: Percentage <= 0.01")
  print(length(which(nonzeroes$sri <= 0.01)) / nrow(nonzeroes))
  print(".")

  print("Observed: Percentage <= 0.02")
  print(length(which(nonzeroes$sri <= 0.02)) / nrow(nonzeroes))
  print(".")

  print("Observed: Percentage <= 0.05")
  print(length(which(nonzeroes$sri <= 0.05)) / nrow(nonzeroes))
  print(".")

  print("Observed: Max SRI:")
  print(gai[which.max(gai$sri),1:20])
  print(".")

  # Significance tests
  means <- apply(gai.obs[,rcol],2,function(x){mean(x,na.rm=TRUE)})
  hist(means,main="Mean of randomizations",breaks=20,col="grey")

  sds <- apply(gai.obs[,rcol],2,function(x){sd(x,na.rm=TRUE)})
  hist(sds,main="SD of randomizations",breaks=20,col="grey")

  cvs <- sds / means
  hist(cvs,main="CV of randomizations",breaks=20,col="grey")

  # Test for social preferences
  print("P value: Are strong social associations more common than expected?")
  cvobs <- sd(gai$sri) / mean(gai$sri)
  print(
    length(which(cvs >= cvobs)) / nrow(gai)
  )
  print(".")

  # Percentage of signficant dyadic associations
  print("What percentage of dyadic assocations are significant?")
  nonzeroes <- gai[gai$sri > 0,]
  #sigs <- which(gai$pu >= 0.975) ; sigs
  sigs <- which(nonzeroes$pl >= 0.95) ; sigs
  print(
    length(sigs) / nrow(nonzeroes)
  )
  print(".")

}



