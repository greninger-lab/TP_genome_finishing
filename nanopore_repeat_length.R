#By: Nicole Lieberman, PhD
#naliebe@uw.edu
#this script is used to parse all tp0470 and arp mapping statistics

library(tidyverse)

##### SETUP #################

args <- (commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
  run_name="testing"
}else{
  run_name=args[[1]]
}

path <- "/PATH/TO/fastq/DI/bam/" #containing directories to arp and tp0470 cov stats. NEEDS TO END IN /
key <- read_csv("/PATH/TO/fastq/sample_sheet.csv")

######################
cov_path <- paste0(path,"stats/") #folder containing coverage stats
arp_cov_ext <- ".arp.covstats-all.txt" #arp covstats file name extension
arp_glob_cov <- paste0("*",arp_cov_ext)
tp0470_cov_ext <- ".tp0470.covstats.txt" #tp0470 covstats file name extension 
tp0470_glob_cov <- paste0("*",tp0470_cov_ext)

########################### 
#ARP:

#read in coverage stats and merge:
arp_cov_files <- fs::dir_ls(cov_path, glob=arp_glob_cov)   
arp_cov_stats <- arp_cov_files %>% purrr::map_dfr(read_tsv, .id = "Barcode") 

#tidy Barcode names, add sample names, etc:
arp_cov_stats$Barcode <- gsub(cov_path, "",arp_cov_stats$Barcode)
arp_cov_stats$Barcode <- gsub(arp_cov_ext, "",arp_cov_stats$Barcode)

arp_cov_stats <- inner_join(key, arp_cov_stats)
colnames(arp_cov_stats)[3] <- "Ref"

arp_cov_stats$total_reads <- arp_cov_stats$Plus_reads + arp_cov_stats$Minus_reads
arp_reads <- arp_cov_stats %>% group_by(Barcode) %>% summarize(sum(total_reads)) 
arp_cov_stats <- inner_join(arp_cov_stats, arp_reads)
arp_cov_stats$percent_of_reads <- (arp_cov_stats$total_reads / arp_cov_stats$`sum(total_reads)`)*100

write.csv(arp_cov_stats, paste0("~/Desktop/arp_full_", run_name, ".csv"))

arp_top <- arp_cov_stats %>% group_by(Barcode) %>% arrange(Barcode, desc(percent_of_reads)) %>% filter(row_number()==1)
arp_top$`Manual Fix` <- ""
write.csv(arp_top, paste0("~/Desktop/arp_top_", run_name, ".csv"))

##############################
#TP0470

#read in coverage stats and merge:
tp0470_cov_files <- fs::dir_ls(cov_path, glob=tp0470_glob_cov)   
tp0470_cov_stats <- tp0470_cov_files %>% purrr::map_dfr(read_tsv, .id = "Barcode") #Barcode is sample name for now

#tidy Barcode names, add sample names, etc:
tp0470_cov_stats$Barcode <- gsub(cov_path, "",tp0470_cov_stats$Barcode)
tp0470_cov_stats$Barcode <- gsub(tp0470_cov_ext, "",tp0470_cov_stats$Barcode)

tp0470_cov_stats <- inner_join(key, tp0470_cov_stats)
colnames(tp0470_cov_stats)[3] <- "Ref"

tp0470_cov_stats$total_reads <- tp0470_cov_stats$Plus_reads + tp0470_cov_stats$Minus_reads
tp0470_reads <- tp0470_cov_stats %>% group_by(Barcode) %>% summarize(sum(total_reads)) #Export this after adding name - total tp0470 reads mapped
tp0470_cov_stats <- inner_join(tp0470_cov_stats, tp0470_reads)
tp0470_cov_stats$percent_of_reads <- (tp0470_cov_stats$total_reads / tp0470_cov_stats$`sum(total_reads)`)*100

write.csv(tp0470_cov_stats, paste0("~/Desktop/tp0470_full_", run_name, ".csv"))

tp0470_top <- tp0470_cov_stats %>% group_by(Barcode) %>% arrange(Barcode, desc(percent_of_reads)) %>% filter(row_number()==1)
tp0470_top$`Manual Fix` <- ""
write.csv(tp0470_top, paste0("~/Desktop/tp0470_top_", run_name, ".csv"))

