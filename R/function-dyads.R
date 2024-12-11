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

  # Needs the ILV dataset
  mr <- read.csv("../data/ILV.csv",stringsAsFactors=FALSE)

  SUMM <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
  ntot <- nrow(SUMM) ; ntot

  hw <- read.csv("../data/hw-master.csv",stringsAsFactors=FALSE)
}

#########################################################
#########################################################
# Functions

fast.dyads2 <- function(hw,sit.floor=10){
  #sit.floor <- 5

  head(hw)
  uuid <- names(which(table(hw$id)>=sit.floor)) ; length(uuid)
  uuid

  # Put all necessary info in a string for each species
  uid <-
    hw %>%
    filter(id %in% uuid) %>%
    group_by(id) %>%
    summarize(uid = paste0(id[1],'---',
                           length(unique(groupid)),'---',
                           paste(unique(groupid), collapse=' '))) %>%
    #ungroup() %>%
    pull(uid)

  uid %>% length
  uid[1]

  # Create a set of all potential dyads
  #uid
  df <- expand.grid(A=uid,B=uid) ; nrow(df)
  #head(df)

  # Extra info from strings
  dfex <-
    df %>%
    mutate(Aid =    lapply(str_split(A,'---'), '[[',1) %>% unlist %>% as.character) %>%
    mutate(Bid =    lapply(str_split(B,'---'), '[[',1) %>% unlist %>% as.character) %>%
    #mutate(nay =    lapply(str_split(A,'---'), '[[',2) %>% unlist %>% as.numeric) %>%
    #mutate(nby =    lapply(str_split(B,'---'), '[[',2) %>% unlist %>% as.numeric) %>%
    mutate(na =    lapply(str_split(A,'---'), '[[',2) %>% unlist %>% as.numeric) %>%
    mutate(nb =    lapply(str_split(B,'---'), '[[',2) %>% unlist %>% as.numeric) %>%
    #mutate(Ayears = lapply(str_split(A,'---'), '[[',3) %>% unlist) %>%
    #mutate(Byears = lapply(str_split(B,'---'), '[[',3) %>% unlist) %>%
    mutate(Asits =  lapply(str_split(A,'---'), '[[',3) %>% unlist) %>%
    mutate(Bsits =  lapply(str_split(B,'---'), '[[',3) %>% unlist)

  #dfex %>% head

  # Calculate X
  dfex1 <-
    dfex %>%
    rowwise() %>%
    mutate(X = length(which(unique(unlist(str_split(Asits,' '))) %in% unique(unlist(str_split(Bsits,' ')))))) %>%
    select(A=Aid, B=Bid, na:nb, X) # get rid of hefty columns

  names(dfex1)

  # Calculate final metrics & format dataframe
  (nsits <- length(unique(hw$groupid)))
  dyads <-
    dfex1 %>%
    mutate(dyad = paste0(A,'-',B),
           ya = na - X,
           yb = nb - X,
           yabx = ya + yb + X,
           ynull = nsits - (na + nb - X)) %>%
    select(dyad, A, B, na, nb, X, ya, yb, yabx, ynull)

  head(dyads)
  return(dyads)
}


