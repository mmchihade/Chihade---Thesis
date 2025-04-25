library(tidyverse)
library(janitor)
library(here)
library(skimr)
library(vegan)
library(phyloseq)
library(patchwork)
library(agricolae)
library(ggpubr)
library(ggplot2)
library(wesanderson)
library(scales)
library(data.table)

data_otu <- read.csv("data/OTU_tables/kraken2_READS_sg_2024_decontam_G.csv", header = TRUE)


#I want to add the number of zeros present for each genus: 

otu_w_zero <- data_otu %>% 
  sjmisc::rotate_df() %>% 
  row_to_names(1) %>% 
  mutate(sum_na = rowSums(is.na(.))) %>% 
  rownames_to_column(var="genus") %>% 
  mutate(rownumber=row_number())

#while we are at it let's add a variable that counts the number of low values: 

otu_w_info <- otu_w_zero %>%
  pivot_longer(names_to = "lab_id", values_to = "count", CABX004:TZHZ469) %>% 
  mutate(count = as.numeric(count)) %>%
  filter(count <= 50) %>% 
  group_by(genus) %>% 
  mutate(low_value_count = n()) %>% 
  select(genus, low_value_count) %>% 
  distinct() 

otu <- merge(otu_w_zero, otu_w_info, by="genus", all=TRUE) %>% 
  arrange(rownumber) %>% 
  select(-rownumber) %>%
  replace(is.na(.), 0)

write_csv(otu, "data/OTU_tables/zerocount_kraken2_READS_sg_2024_decontam_G.csv")

#I also want to make a new otu table that have high and low instead of values: 

#for this I need the percent version of the otu table: 

data_otu_percent <- read.csv("data/OTU_tables/kraken2_PERCENT_sg_2024_decontam_G.csv", header = TRUE)

options(scipen = 999)

high_low_otu_1 <- data_otu_percent %>% 
  column_to_rownames(var = "subject") %>% 
  mutate_if(is.character, as.numeric) %>%
  mutate(across(everything(), ~replace(., . >= 0.1, "high"))) %>%
  mutate(across(everything(), ~replace(., . < 9e-7, "super_low"))) %>% 
  mutate(across(everything(), ~replace(., . < 1e-1, "low"))) %>% 
  rownames_to_column(var="subject") 

write_csv(high_low_otu_1, "data/OTU_tables/high_low_values_0.1_threshold.csv")


