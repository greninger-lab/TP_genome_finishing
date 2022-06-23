#By: Nicole Lieberman, PhD
#naliebe@uw.edu

library(tidyverse)
library(reshape2)

setwd("~/Desktop/")

sample <- "samplename" #strain of interest for plotting

#ARP

arp_repeats <- read_csv("~/Desktop/arp_full_testing.csv") #path to arp_full csv 

arp_alignments <- filter(arp_repeats, Sample == sample) 
n <- as.numeric(gsub("ARP", "", arp_alignments$Ref))
gg_arp <- data.frame(Reads = arp_alignments$total_reads, Repeats = n)

arp <- ggplot(gg_arp, aes(x=Repeats, y=Reads, label=Reads)) + 
  geom_bar(stat="identity") +
  labs(title=sample)
arp

png(paste0("~/Desktop", sample, "_arp.png"), width=3, height=3, units="in", res=300)
arp
dev.off()


#tp0470

tp0470_repeats <- read_csv("~/Desktop/tp0470_full_testing.csv") #path to tp0470_full csv 

tp0470_alignments <- filter(tp0470_repeats, Sample == sample) 
n <- as.numeric(gsub("TP0470_", "", tp0470_alignments$Ref))
gg_tp0470 <- data.frame(Reads = tp0470_alignments$total_reads, Repeats = n)

tp0470 <- ggplot(gg_tp0470, aes(x=Repeats, y=Reads, label=Reads)) + 
  geom_bar(stat="identity") +
  labs(title=sample)
tp0470

png(paste0("~/Desktop", sample, "_tp0470.png"), width=3, height=3, units="in", res=300)
tp0470
dev.off()
