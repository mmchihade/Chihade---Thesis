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

metadata <- read_csv("data/Subject_Data.TishkoffLab.nonPHI.2024_freeze_3_MOD_SHOTGUN_CURATE.csv")


#clean up data ---- 

cleanmetadata <- metadata %>% #we are creating a new dataframe based on our beaches raw data  
  clean_names()

write_csv(cleanmetadata, "clean_metadata_tishkofflab_2024.csv")

clean_expanded_metadata <- cleanmetadata %>%
  select(-subject_cm8, -subject_cm7, -abbpop, -contains('admixture'), -contains('plain'), -contains('639'), -contains('genetic'), -contains('plate'), -contains('dna'), -contains('library'), -udi_coord, -contains('barcode'), -contains('reads'), -contains('flow'), -contains('nova'), -run_start_date, -index_id, -vol_ep_motion_ul, -contains('bracken')) %>% 
  mutate(across(.cols = 11:27, .fns = tolower)) %>%
  mutate(across(.cols = 29:226, .fns = tolower))

write_csv(clean_expanded_metadata, "data/clean_expanded_metadata_tishkofflab_2024.csv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

clean_simple_metadata <- clean_expanded_metadata %>%
  select(record_id, lab_id, field_id, topmed_id, collection_name, country, region, sample_site, longitude, latitude, date_collected_full, sample_year, sex, age, is_adult, study_group, country_group, lifestyle_subsistence_major,  subject_ethnicity, waist_circumference, body_fat_percent, bmi) %>%
  rename(date = date_collected_full) %>% 
  rename(year = sample_year) %>% 
  rename(ethnicity = subject_ethnicity) %>% 
  rename(subsistance = lifestyle_subsistence_major) %>% 
  mutate(ethnicity = recode(ethnicity,
                          "ju||hoan" = "khoesan",
                          "kaukau" = "khoesan",
                          "ju|'hoan/yeyi" = "khoesan",
                          "naro" = "khoesan",
                          "kaukau/naro" = "khoesan",
                          "!xoo" = "khoesan",
                          "ju|'hoan" = "khoesan",
                          "kaukau/!xoo" = "khoesan", 
                          "naro/kaukau" = "khoesan", 
                          "bangbelle - see subject notes" = "bantu General",
                          "(bantu) ngoumba" = "bantu General",
                          "(bantu)" = "bantu General",
                          "bakoko - see subject notes" = "bantu General",
                          "bantu - see subject notes" = "bantu General",
                          "bantu" = "bantu General",
                          "bantu ngoumba" = "bantu General",
                          "bulu" = "bantu General",
                          "mabea" = "bantu General",
                          "mpompong" = "bantu General",
                          "ndjem" = "bantu General", 
                          "bangbele" = "bantu General",
                          "tswana" =  "bantu General",
                          "mophadima" = "bantu General",
                          "mpompong" = "bantu General",
                          "maka" = "bantu General",
                          "tswana (mophadima)/tswana" = "bantu General", 
                          "aku mbororo fulani - see subject notes" = "mbororo Fulani", 
                          "bahurutshe (tswana)" = "bantu General", 
                          "mabea - see subject notes" = "bantu General",
                          "baka - see subject notes" = "baka", 
                          "bantu fang - see subject notes" = "fang", 
                          "bantu fang" = "fang",
                          "(bantu) fang" = "fang",
                          "bikele maka" = "bantu General", 
                          "ewondo - see subject notes" = "bantu General", 
                          "fang - subject notes" = "fang", 
                          "fulani" = "mbororo Fulani",
                          "mbororo" = "mbororo Fulani",
                          "kgalagadi (mopebana)" = "kgalagadi", 
                          "ngoumba - see subject notes" = "ngoumba", 
                          "Bandem" = "bantu General", 
                          "mbororo fulani - see subject notes" = "mbororo Fulani",
                          "mbororo fulani" = "mbororo Fulani", 
                          "zime" = "nzime",
                          "bandem" = "bantu General", 
                          "ewondo" = "bantu General")) %>% 
  mutate_if(is.character, str_replace_all, '_', ' ') %>% 
  mutate_if(is.character, str_replace_all, ';', ' ') %>% 
  mutate(lab_id=as.factor(lab_id)) %>% 
  mutate_if(is.character, firstup) %>%
  select(-country_group) %>%
  unite(col = "country_group", country, study_group, sep = " ", remove = FALSE) %>%
  arrange(sample_site) 

write_csv(clean_simple_metadata, "data/clean_simple_metadata_tishkofflab_2024.csv")

#let's add some alpha diversity values to this table ----

#to do this I actually need to load the values: 

alpha <- read_csv("data/alpha_diversity.csv") %>% 
  select(record_id, S.obs, Shannon)

clean_simple_metadata_alpha <- merge(alpha, clean_simple_metadata, by = "record_id")

clean_simple_metadata_alpha <- clean_simple_metadata_alpha %>%
  filter(!row_number() %in% c(413, 414)) %>%
  select(record_id, lab_id, field_id, topmed_id, collection_name, country, region, sample_site, longitude, latitude, date, year, sex, age, is_adult, study_group, country_group, subsistance,  ethnicity, waist_circumference, body_fat_percent, bmi, S.obs, Shannon)  %>%
  arrange(sample_site)

write_csv(clean_simple_metadata_alpha, "data/clean_simple_metadata_tishkofflab_alpha_2024.csv")

