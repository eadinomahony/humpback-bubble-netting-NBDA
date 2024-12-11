#########################################################
#########################################################
# Format raw datasets

#########################################################
#########################################################

library(tidyverse)
library(readxl)
library(lubridate)

# If working in R Studio, shortcut to set wd to folder that this R file is in
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################################################
#########################################################
# Bring in datasets

land_file <- '../data/raw/HW-LAND-2014-2023-v20231023.xlsx'
land <- read_xlsx(land_file)

boat_file <- '../data/raw/HW-BOAT-2004-2023-v20231023.xlsx'
boat <- read_xlsx(boat_file)

#########################################################
#########################################################
# Format land data

land %>% head

# Format column names
names(land) <- tolower(names(land)) #lower case
names(land) <- gsub(' ','_', names(land)) #replace space with underscore
land %>% names

# Remove letters from group_#
groups <- land$`group_#`
groups <- sapply(groups, function(x){substr(x, 1, nchar(x)-1)})
groups %>% unique
land$group <- groups

# Remove, add, rename, rearrange/select
land <-
  land %>%
  mutate(dt = paste(year, month, day, sep='-')) %>%
  mutate(ymd = ymd(dt)) %>%
  mutate(doy = yday(ymd)) %>%
  mutate(yfrac = year + (doy / 365)) %>%
  mutate(time = NA, lat = NA, lon = NA) %>%
  dplyr::rename(n = true_grp_sz,
         n_bnf = total_bnf,
         bhvr = `final-bhv`,
         stage = `life_stage_or_sex`) %>%
  mutate(groupid = paste0(substr(platform, 1, 4), gsub('-','',ymd),str_pad(group, width=2, side='left', pad='0'))) %>%
  mutate(n = as.numeric(n)) %>%
  mutate(n = as.numeric(n_bnf)) %>%
  mutate(bhvr = as.character(bhvr)) %>%
  mutate(stage = as.character(stage)) %>%
  mutate(group = as.numeric(group)) %>%
  select(groupid, id, ymd, year, month, day, doy, yfrac, time, platform, group, n, n_bnf, lat, lon, bhvr, stage)

land %>% head
land %>% tail

#########################################################
#########################################################

boat %>% head

# Format column names
boat <- boat[,1:17]
names(boat) <- tolower(names(boat))
names(boat) <- gsub(' ','_', names(boat))
boat %>% names


# Deal with Verney group in 2020 -- treat as one large group
(groups <- boat$`grp_#`)
i <- grep('-', groups)
groups[i]
groups[i] <- '1A'

# Deal with long group names e.g. 1A/B
i <- grep('/', groups)
groups[i]
groups[i] <- substr(groups[i], 1, 3)
groups[i]

# Remove letters
groups <- sapply(groups, function(x){substr(x, 1, nchar(x)-1)})
groups <- gsub('A', '', groups)
groups <- gsub('B', '', groups)
groups <- gsub(' ?','', groups) # Remove question marks
groups <- str_pad(groups, width=2, side='left', pad='0')
groups %>% unique
boat$grp <- groups

boat %>% names

# Remove extraneous columns & add date
boat <-
  boat %>%
  mutate(dt = paste(year, month, day, sep='-')) %>%
  mutate(ymd = ymd(dt)) %>%
  mutate(platform = 'Elemiah') %>%
  mutate(doy = yday(ymd)) %>%
  mutate(yfrac = year + (doy / 365)) %>%
  mutate(time = NA, lat = NA, lon = NA) %>%
  dplyr::rename(n='size',
         n_bnf='total_bnf',
         bhvr = 'bhv-final',
         stage = 'type') %>%
  mutate(groupid = paste0(substr(platform, 1, 4), gsub('-','', ymd), grp)) %>%
  mutate(n = as.numeric(n)) %>%
  mutate(n = as.numeric(n_bnf)) %>%
  mutate(bhvr = as.character(bhvr)) %>%
  mutate(stage = as.character(stage)) %>%
  mutate(group = as.numeric(grp)) %>%
  select(groupid, id, ymd, year, month, day, doy, yfrac, time, platform, group, n, n_bnf, lat, lon, bhvr, stage)

boat %>% head
boat %>% tail

#########################################################
#########################################################
# Combine datasets

mr <- rbind(land, boat)
nrow(mr)

#########################################################
#########################################################
# More cleaning

# Deal with calf names
(i <- grep('CALF', toupper(mr$id)))
mr$id[i] %>% unique()
mr$id[i] <- 'CALF'
(i <- grep('CLAF', toupper(mr$id)))
mr$id[i] <- 'CALF'
(i <- grep('BCX0239-C', toupper(mr$id)))
mr$id[i] <- 'CALF'

# Deal with problematic IDs
mr$id %>% unique
mr$id <- gsub(' ','', mr$id)
mr <- mr[mr$id != 'BCY', ]
bads <- c('NOID', 'NEWID', 'NOISPS', 'POORPIC', 'POORPC', 'RE-CHK', 'RECHK',
          'CHAMP', 'CHANP', 'NA', '\\?', 'JUVXX', 'BELOWSM', 'BXC?')
badi <- sapply(bads, function(x){grep(x, mr$id)}) %>% unlist %>% unique
(goodi <- which(! 1:nrow(mr) %in% badi))
mr <- mr[sort(goodi), ]
mr$id %>% unique
nchar(mr$id) %>% table
nrow(mr)

# Remove invalid groups
mr$group %>% unique
mr$group %>% table(useNA='always')
mr <- mr %>% filter(!is.na(group))
nrow(mr)
head(mr)

# Remove FL behavior
(fl <- grep('FL', mr$bhvr))
mr$bhvr[fl] <- NA


#########################################################
#########################################################
# Save dataset

save(mr, file='../data/events.rds')
write.csv(mr, '../data/events.csv', quote=FALSE, row.names=FALSE)

#########################################################
#########################################################
# Review / explore / misc

mr %>% filter(groupid == 'Fin-2018051901')

