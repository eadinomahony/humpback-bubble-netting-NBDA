#### NBDA Figure S3

#set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# reload data 

mr <- read.csv("../data/ILV-update.csv", stringsAsFactors = F); head(mr[,1:13])

### distribution of observation frequencies to justify threshold of 5:
library(dplyr)

mr %>%
  count(n.sits) %>%
  arrange(n.sits)

sum(mr$n.sits >= 5)
sum(mr$n.sits < 5)

library(ggplot2)

sit.freq <- ggplot(mr, aes(x = n.sits)) +
  geom_histogram(binwidth = 1, fill = "cornflowerblue", color = "black", alpha=0.5) +
  geom_vline(xintercept = 4, color = "red", size = 0.6) +
  labs(
    x = "Number of sightings",
    y = "Number of individuals") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust= 0.5),
        axis.text = element_text(color = "black", size = 18),  
        axis.title = element_text(color = "black", size = 20))

ggsave("../figures/resubmission-plots/sightings-frequency-distribution.pdf", sit.freq, width = 12, height = 6, dpi= 300)
ggsave("../figures/resubmission-plots/sightings-frequency-distribution.png", sit.freq, width = 12, height = 6, dpi= 300)
ggsave("../figures/resubmission-plots/sightings-frequency-distribution.svg", sit.freq, width = 12, height = 6, dpi= 300)
