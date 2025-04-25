#load packages ----
library(tidyverse)
library(patchwork)
library(agricolae)
library(ggpubr)
library(ggplot2)
library(janitor)
library(wesanderson)
library(scales)
library(data.table)
library(ggspatial)
library(patchwork)
library(scico)


#load data ----

data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)

sumstats <- data_grp %>% 
  group_by(country) %>%
  mutate(samples_per_country = n()) %>% 
  ungroup() %>% 
  group_by(ethnicity) %>% 
  mutate(samples_per_ethnicity = n()) %>% 
  ungroup() %>% 
  group_by(is_adult) %>% 
  mutate(samples_per_age = n()) %>% 
  group_by(year) %>% 
  mutate(samples_per_year = n()) %>% 
  group_by(study_group) %>% 
  mutate(samples_per_group = n()) %>%
  select(country, samples_per_country, ethnicity, samples_per_ethnicity, is_adult, samples_per_age, year, samples_per_year, study_group, samples_per_group)

write_csv(sumstats, "output/sumstats.csv")

numofeth <- data_grp %>%
 select(ethnicity) %>% 
 mutate(num = n_distinct(ethnicity))

data_grp %>% 
  summarise(max = max(age), 
            min = min(age)) 
  
