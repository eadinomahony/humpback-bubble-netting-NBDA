#########################################################
#########################################################
# The diffusion of cooperative and solo bubble net feeding in Canadian Pacific humpback whales
# v. Éadin O'Mahony, Oct 2024

#set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load custom functions:
source("00 Source functions.R")

#install.packages("assocInd")
library(dplyr)
library(tidyverse)
library(assocInd)
library(lubridate)
library(devtools)
library(rootSolve)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#download NBDA package from github:
#devtools::install_github("whoppitt/NBDA")
library(NBDA)

################################################################################
################################################################################

### THE WHOLE STUDY HINGES UPON THREE DATASETS: 1) events.csv; 2) ILV.csv and 3) assoc-final.csv
### These are all the sightings and each individual's ILVs 
### (sex information has been gleaned from blow sampling data)
### 1) Each row is one sighting
### 2) Each row is one individual
### 3) assoc-final.csv is a matrix with dyadic SRI values

### NBDA models using NBDA package: https://github.com/whoppitt/NBDA

################################################################################
################################################################################

##### Load events.csv and clean it up a little:

hw <- read.csv("../data/events.csv",stringsAsFactor=FALSE) ; hw

hw.unfiltered <- hw # all sightings retained

# this sightings csv has sightings of calves in it which need to be removed 
# because they don't have IDs:
nrow(hw) # 7697

# Remove calves:
hw <- hw %>% filter(id != 'CALF')
nrow(hw) # 7490

hw.nocalves <- hw; nrow(hw.nocalves)
# 7697 - 7490 = 207 calf sightings removed

# check for other discrepancies:
hw[nchar(hw$id) < 7,] 

# there are a few erroneous IDs, remove them
hw <- hw[nchar(hw$id) >= 7, ]
nrow(hw) # 7485 photo IDs

################################################################################
################################################################################
# Methods 
# (a) Data collection

### earliest and latest dates for each platform:
date_summary <- hw.unfiltered %>%
  group_by(platform) %>%
  summarize(
    earliest_date = min(ymd, na.rm = TRUE),  
    latest_date = max(ymd, na.rm = TRUE) 
  ); date_summary


################################################################################
################################################################################
# Results
# Paragraph 1

#Of the 526 identified individual whales, 250 were encountered on 5 occasions. 
#58% of whales were encountered across multiple years (32% in 5 or more years, 
#15% in 10 or more years). Bubble netting was observed in 254 individuals (48% 
#of the identified population, hereafter referred to as ‘bubble netters’) on 635 occasions.

# Load ILV.csv and clean it up a little:
mr <- read.csv("../data/ILV-update.csv", stringsAsFactors = F); head(mr[,1:13])

#extract the earliest sighting ID for each individual (makes 'sit1' column)
mr <- mr %>%
  mutate(sit1 = sapply(strsplit(sits, " "), 
                       function(x) x[which.min(as.Date(substr(x, 5, 12), 
                                                       format = "%Y%m%d"))]))
head(mr)

# what % of whales were seen in multiple years?
mr %>% filter(n.years > 1) %>% select(id) # 306
306 / nrow(mr) # 0.5806

# what % of whales were seen in 5 or more years?
mr %>% filter(n.years > 4) %>% select(id)
169 / nrow(mr) #0.3207

# what % of whales were seen in 10 or more years?
mr %>% filter(n.years > 9) %>% select(id)
81 / nrow(mr) #0.1537

# Bring in sightings data
sits_raw <- read.csv("../data/Sightings.csv",stringsAsFactors=FALSE)
sits_raw %>% head

# 1. Total # of photo identifications and total number of unique IDs?
nrow(hw) # 7485 photo IDs, when filtered for 'CALF' ID
length(unique(hw$id)) # 526 IDs

# 2. Total number of encounters?
head(sits_raw)
nrow(unique(sits_raw)) # 4235
sits_raw[nchar(sits_raw$sit) < 14,] # wonky sit id because the date info is missing in raw dataset
length(unique(sits_raw$ymd)) # 1357
length(unique(sits_raw$sit)) # 4062

length(unique(hw$groupid)) # 4056 ---
# there should be 14 characters in groupids, clean this up in hw:
hw[nchar(hw$groupid) < 14,] 
hw <- hw[nchar(hw$groupid) >= 14, ]
length(unique(hw$groupid)) # encounters 

# 3. Number of days of sampling?
head(hw)
length(unique(hw$ymd)) # 1356 days of sampling

################################################################################
################################################################################
# Results
# Paragraph 2
# How many bubble netters are there in total?

#how many bubble netters are there in total in mr?
nrow(mr[mr$bnfe == 1,]) # 254
(254 / 526)*100 # 48.29% of the total identified population

################################################################################
################################################################################
# Results
# Paragraph 2
# How many bubble netters were bubble netting on our first encounter of them?

#the number of rows where 'bnf.sit1' equals 'sit1'
num_equal_rows <- sum(mr$bnf.sit1 == mr$sit1, na.rm = TRUE); num_equal_rows 
# 82 whales bubble net fed on 1st sighting

num_unequal_rows <- sum(mr$bnf.sit1 != mr$sit1, na.rm = TRUE); num_unequal_rows 
# 172 bubble netters don't, subset out these 172 whales:

subset <- mr[mr$bnf.sit1 != mr$sit1 & !is.na(mr$bnf.sit1) & !is.na(mr$sit1), ]

#now find these sit IDs from bnf.sit1 in hw, and check if they were solo or cooperative on their first bubble net:
bnf_explore <- subset %>%
  dplyr::select(bnf.sit1) %>%                  
  dplyr::distinct() %>%                       
  dplyr::left_join(hw %>%                     
                     dplyr::select(groupid, n),
                   by = c("bnf.sit1" = "groupid")) %>%
  dplyr::arrange(bnf.sit1)                   

print(bnf_explore)

# how many of these are solo bnf events?
bnf_explore %>%
  dplyr::filter(n == 1) %>%
  dplyr::count()
# 24

# manually checked these 24 s=bnf.sit1 IDs in 'subset', 21 of the 172 whales solo bubble net fed on their first bnf sit
172 - 21 # 151 group bubble net fed on their first bnf sit

# % of whales seen any number of times who bubble net fed IN A GROUP on their first sighting (learning by doing?)
# i.e. how many of the 82 whales that BNFed on their first sit, did so in a group?
homophilic_whales <- subset(mr, bnf.sit1 == sit1 & !is.na(bnf.sit1)); nrow(homophilic_whales) # 82

#explore
head(homophilic_whales[,1:15])

# 65 individuals were observed cooperatively bubble netting on 1st encounter
# 17 were solo bubble netting on 1st encounter
17+65 # = 82
(65/254)*100 # 25.59% of all bubble netters were cooperatively bubble netting on 1st encounter
(17/254)*100 # 6.69% of all bubble netters were solo bubble netting on 1st encounter

################################################################################
################################################################################
# Results
# Paragraph 2
# How many bubble netting events did we see in total? 
# How many of these were cooperative?
# How much of all behaviour was feeding behaviour?

unique(hw$bhvr) 
count_bnf <- sum(hw$bhvr == "BNF", na.rm = TRUE); count_bnf
count_fe <- sum(hw$bhvr %in% c("BNF", "FE", "TRAPFEED", "L-FE", "SKIM FE", "LU-FE", 
                               "TR-FE", "BN", "LF", "SG-FE", "TF", "MI-FE", 
                               "MI-FE-SG", "SG-MI-FE", "MIFE", "TR-FE", "MI-FE-SOC", 
                               "TR-LF", "SG-FE-SL", "MI-LFF", "MI-FE-RE"), na.rm = TRUE)
count_fe

count_fe / nrow(hw) # 0.639; proportion of sightings assumed to be feeding associated

################################################################################
################################################################################

# Restrict sightings to individuals seen at least 5 times

idtab<-table(hw$id); idtab
hw <- subset(hw, id %in% names(idtab[idtab >= 5]))
table(hw$id) # check that it worked

nrow(hw) # 6984 sightings left when filtered to whales seen 5 or more times

################################################################################
################################################################################


# the # of bubble netters cooperating on first sighting but who eventually solo BNF also:
first_feed_in_group <- sum(homophilic_whales$bnf.sit1 != homophilic_whales$solo.sit1, na.rm = TRUE); first_feed_in_group 
first_feed_in_group_subset <- subset(homophilic_whales, bnf.sit1 != solo.sit1)

first_sit_soloBNF <- sum(homophilic_whales$bnf.sit1 == homophilic_whales$solo.sit1, na.rm = TRUE); first_sit_soloBNF # 17

sum(!is.na(homophilic_whales$solo.sit1)) # 26 soloists who bubble net fed on their first sighting (either in a group or as solo)

### explore the whales that BOTH group and solo bnf:

BHVmorphing <- mr[mr$bnfsits != "" & mr$solosits!="", c("id", "bnfsits", "solosits", "sit1", "solo.sit1", "bnf.sit1")] # subset to include ALL whales known to solo bnf
nrow(BHVmorphing)

soloIDs.2 <- BHVmorphing$id

PUREsoloists <- BHVmorphing[BHVmorphing$bnfsits == BHVmorphing$solosits,]; nrow(PUREsoloists) # whales that ONLY solo bnf

(nrow(BHVmorphing)-nrow(PUREsoloists))/nrow(BHVmorphing)*100 # 74.1573% of solo bubble netters perform both solo and group bnf


################################################################################
################################################################################
### Now filter dataset for minimum of 5 sightings per individual: 

### % of whales seen 5 or more times who BNF on first sit (sit1)?

# reload data to subset properly:

mr <- read.csv("../data/ILV-update.csv", stringsAsFactors = F); head(mr[,1:13])

# Restrict to individuals seen at least 5 times
# this removes indivs seen less than 5 times, so ROA will no longer have consecutive numbers!
mr <- mr[mr$n.sits >= 5,] 

# Extract the earliest sighting ID for each row
mr <- mr %>%
  mutate(sit1 = sapply(strsplit(sits, " "), 
                       function(x) x[which.min(as.Date(substr(x, 5, 12), 
                                                       format = "%Y%m%d"))]))

# Count the number of rows where 'bnf.sit1' equals 'sit1'
num_equal_rows_min5sits <- sum(mr$bnf.sit1 == mr$sit1, na.rm = TRUE); num_equal_rows_min5sits # 28 whales BNFed on their first sighting
num_unequal_rows_min5sits <- sum(mr$bnf.sit1 != mr$sit1, na.rm = TRUE); num_unequal_rows_min5sits # 151

num_equal_rows_min5sits / sum(!is.na(mr$bnf.sit1)) * 100 # 15.64% of bubble netting whales seen 5 or more times were bubble netting on their first sighting

# How many bubble netters are there that were seen at least 5 times?
sum(!is.na(mr$bnf.sit1)) # 179
# How many of these were encountered BEFORE we saw them bubble netting
sum(mr$bnf.sit1 != mr$sit1, na.rm = T)

#########################################################
#########################################################

# all of BNF sightings, how many were group and how many were solo?

bnf.events <- hw.unfiltered %>%
  filter(grepl("^BN", bhvr))

nrow(bnf.events) # 2201
max(bnf.events$n, na.rm = T) # 16

length(unique(bnf.events$groupid)) # 635 unique encounters of BNF

# of whales sighted 5 or more times, 
bnf.events2 <- hw %>%
  filter(grepl("^BN", bhvr))

length(unique(bnf.events2$groupid)) # 612 unique encounters when filters to whales seen 5 or more times

# now from this, how many have group size of 1?
solo.bnf.events <- bnf.events[bnf.events$n == 1 & !is.na(bnf.events$n),]
nrow(solo.bnf.events) # 169, although two of these are not true solo BNF, so = 167

(nrow(solo.bnf.events)-2) / nrow(bnf.events) * 100 # 7.58% of bnf events are solo events

group.bnf <- nrow(bnf.events) - (nrow(solo.bnf.events)-2) # = 2034

group.bnf / nrow(bnf.events) * 100 # 92.41 %

# how many bubble netters do we have in the study of those whales seen 5 or more times?
bnfers <- mr[!is.na(mr$bnf.sit1),]; nrow(bnfers)

length(unique(solo.bnf.events$id)) # 95 soloists, 93 true soloists (2 events not solo bnf)

soloIDs <- unique(solo.bnf.events$id)

#########################################################
#########################################################
# Prepare final ILV info we want to use in later NBDA models:

# summarise data in a tidy dataframe called roa:
roa <- data.frame(id=mr$id,
                  ROA=mr$ROA,
                  first.sit=mr$sit1, # first sighting
                  first.year=mr$first,
                  first.bnf.sit=mr$bnf.sit1, # first BNF sighting
                  first.bnf.year=mr$bnf.year1, # first year seen BNFing
                  n.sits=mr$n.sits, # Total number of encounters
                  bnfrate=mr$bnf.rate, # Proportion of sightings that involved BNF
                  momyears=mr$momyears, # include mom years to check mom binary variable
                  ssfi=mr$SSFI, # Site Fidelity index
                  nyears=mr$n.years,
                  genetic.sex=mr$Sex) # genetic sex is female = -1, male = 1 and NA = 0

# Standardize n.sits to have mean of 0 & sd of 1:
roa <- roa %>% mutate_at(c('n.sits'), ~(scale(.) %>% as.vector))
# check
sd(roa$n.sits); mean(roa$n.sits); sd(roa$ssfi); mean(roa$ssfi) # mean is basically zero, sd = 1


### Combine momyears and genetic.sex into final sex ILV:
# First a for loop that assigns a -1 to each whale that HAS momyears (i.e. momyears > 0)
roa$female <- 0
for(i in 1:nrow(roa)){
  # if momyears has info in it, then put a -1 in the corresponding female column
  if (nchar(as.character(roa$momyears[i])) > 0) {
    roa$female[i] <- -1
  }
}


# Now need to bring columns 'female' and 'genetic.sex' together
roa <- roa %>%
  mutate(final.sex = case_when(
    genetic.sex==1 & female==0 ~ 1,
    genetic.sex==-1 & female==0 ~ -1,
    genetic.sex==-1 & female==-1 ~ -1,
    is.na(genetic.sex) & female==-1 ~ -1,
    genetic.sex==0 & female==-1 ~ -1,
    genetic.sex==0 & female==0 ~ 0,
    is.na(genetic.sex) & female==0 ~ 0,
    TRUE~genetic.sex
  ))


################################################################################
################################################################################
#mr has 250 individuals but hw has 254 -- which ones are missing ILVs?

missing_individuals2 <- setdiff(hw$id, mr$id); print(missing_individuals2) # same 4 individuals as in missing_individuals

# and remove from hw:
hw <- hw[!(hw$id %in% missing_individuals2), ]
length(unique(hw$id)) # 250

################################################################################
################################################################################

