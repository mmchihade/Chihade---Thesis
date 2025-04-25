library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(janitor)
library(wesanderson)
library(scales)
library(data.table)
library(purrr)
library(forcats)
library(ggthemes)
library(TSP)
library(Polychrome)
library(rcartocolor)
library(stringr)


#load data ----

data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE) %>%
  select(-record_id, -lab_id, -topmed_id, -collection_name, -sample_site, -longitude, -latitude, -date, -year, -is_adult, -age, -waist_circumference, -body_fat_percent, -bmi) %>% 
  rename(lab_id = field_id) %>%
  mutate(lab_id = toupper(lab_id))

#create a dataframe with all of the cog tally files, with the file name as a column so that we can match column to tally

merged_df <- list.files(path = "data/cog_tally_overall",
                        pattern = "\\.csv$", 
                        full.names = TRUE) %>% 
  set_names() %>%
  map_dfr(read_csv, .id = "file_name") 

#we need the cog categories too! 

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#now let's clean the data up 

#first this guy is HUGE so let's cut her down: 

cog_gen_trim <- clean_names(merged_df) %>% 
  rename(lab_id = file_name) %>%
  mutate(across(everything(), gsub, pattern = "data/cog_tally_overall/", replacement = "")) %>%
  mutate(across(everything(), gsub, pattern = ".csv", replacement = "")) %>%
  mutate(n = as.numeric(n)) %>% 
  filter(genus %in% c("WG-1", "UBA3375","Anaerovibrio", "Treponema_D", "UBA2804", "Campylobacter", "RUG410", "RF16", "CAG-568", "Campylobacter_D", "UBA636", "Succinivibrio", "UBA1436", "UBA71", "RUG572", "Butyricicoccus_A", "CAG-45", "Lachnospira", "Rs-D84", "Brachyspira", "Slackia_A", "CAG-533", "UBA4644"))

#I also have to replace all of the TZ codes with B codes so that I can divide by country and ethnicity 

tally <- cog_gen_trim %>% 
  mutate(lab_id = recode(lab_id,
                         "TZBG067"="B607",
                         "TZBG069"="B609",
                         "TZBG070"="B610",
                         "TZBG072"="B612",
                         "TZBG073"="B613",
                         "TZBG074"="B614",
                         "TZBG076"="B616",
                         "TZBG0771"="B617.1",
                         "TZBG0772"="B617.2",
                         "TZBG079"="B619",
                         "TZBG080"="B620",
                         "TZHZ336"="B082",
                         "TZHZ337"="B083",
                         "TZHZ338"="B084",
                         "TZHZ339"="B085",
                         "TZHZ340"="B086",
                         "TZHZ343"="B089",
                         "TZHZ344"="B090",
                         "TZHZ345"="B091",
                         "TZHZ346"="B092",
                         "TZHZ453"="B202",
                         "TZHZ457"="B206",
                         "TZHZ463"="B212",
                         "TZHZ469"="B218",
                         "TZHZ472"="B221",
                         "TZHZ474"="B223",
                         "TZHZ488"="B238",
                         "TZHZ559"="B312",
                         "TZHZ586"="B339",
                         "TZHZ587"="B340",
                         "TZMS164"="B464",
                         "TZMS177"="B477",
                         "TZMS178"="B478",
                         "TZMS183"="B483",
                         "TZMS190"="B490",
                         "TZMS192"="B492",
                         "TZMS195"="B495",
                         "TZMS196"="B496",
                         "TZMS197"="B497",
                         "TZMS198"="B498",
                         "TZMS199"="B499")) %>% 
  mutate(cog_category = recode(cog_category,
                               "?" = "S",
                               "!" = "S")) %>% 
  group_by(lab_id, cog_category, genus) %>% 
  mutate(n = sum(n, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#I also want to separate all of the multiple letter columns. They will probably be useful later on, but it is hella overwhelming now

full_tally <- tally %>% 
  separate(cog_category, into = c("cog_1", "cog_2", "cog_3", "cog_4", "cog_5", "cog_6", "cog_7"), sep = "") %>%
  select(-cog_1) %>% 
  pivot_longer(names_to="cog_name", values_to = "cog_category", cog_2:cog_7)%>%
  na.omit() %>%
  select(-cog_name) %>%
  group_by(lab_id, cog_category) %>% 
  mutate(tally = sum(n, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#umm, I am going to try to make a combined dataset, let's see if R can handle it ----

merge <- merge(full_tally, data_grp, by = "lab_id")

merge_labels <- merge(merge, cog, by = "cog_category")

#yay!!!!!

#now I want to divide this up by genus: 

WG1 <- merge_labels %>% 
  filter(genus %in% c("WG-1"))

UBA3375 <- merge_labels %>% 
  filter(genus %in% c("UBA3375"))
  
  
Anaerovibrio <- merge_labels %>% 
  filter(genus %in% c("Anaerovibrio"))
  
Treponema_D <- merge_labels %>% 
  filter(genus %in% c("Treponema_D"))

UBA2804 <- merge_labels %>% 
  filter(genus %in% c("WG-1"))

Campylobacter <- merge_labels %>% 
  filter(genus %in% c("Campylobacter"))
  
RUG410 <-  merge_labels %>% 
  filter(genus %in% c("RUG410"))

RF16 <- merge_labels %>% 
  filter(genus %in% c("RF16"))

CAG568 <- merge_labels %>% 
  filter(genus %in% c("CAG-568"))

Campylobacter_D <- merge_labels %>% 
  filter(genus %in% c("Campylobacter_D"))

UBA636 <-  merge_labels %>% 
  filter(genus %in% c("UBA636"))

Succinivibrio <- merge_labels %>% 
  filter(genus %in% c("Succinivibrio"))
  
UBA1436 <- merge_labels %>% 
  filter(genus %in% c("UBA1436"))
  
UBA71 <- merge_labels %>% 
  filter(genus %in% c("UBA71"))

RUG572 <- merge_labels %>% 
  filter(genus %in% c("RUG572"))

Butyricicoccus_A <- merge_labels %>% 
  filter(genus %in% c("Butyricicoccus_A"))

CAG45 <- merge_labels %>% 
  filter(genus %in% c("CAG-45"))
  
Lachnospira <- merge_labels %>% 
  filter(genus %in% c("Lachnospira"))

RS_D84 <- merge_labels %>% 
  filter(genus %in% c("Rs-D84"))

Brachyspira <- merge_labels %>% 
  filter(genus %in% c("Brachyspira"))

Slackia_A <- merge_labels %>% 
  filter(genus %in% c("Slackia_A"))

CAG533 <-  merge_labels %>% 
  filter(genus %in% c("CAG-533"))

UBA4644 <- merge_labels %>% 
  filter(genus %in% c("UBA4644"))


#let's group by a bunch of factors and average within them ----

#okay this is gonna be so long

country_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA4644 <- merge(country_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA4644 <- merge(sex_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA4644 <- merge(ethnicity_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA4644 <- merge(country_group_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA4644 <- merge(region_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA4644 <- UBA4644 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA4644 <- merge(study_group_mean_UBA4644, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_CAG533 <- merge(country_mean_CAG533, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_CAG533 <- merge(sex_mean_CAG533, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_CAG533 <- merge(ethnicity_mean_CAG533, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_CAG533 <- merge(country_group_mean_CAG533, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_CAG533 <- merge(region_mean_CAG533, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_CAG533 <- CAG533 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_CAG533 <- merge(study_group_mean_CAG533, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Slackia_A <- merge(country_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Slackia_A <- merge(sex_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Slackia_A <- merge(ethnicity_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Slackia_A <- merge(country_group_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Slackia_A <- merge(region_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Slackia_A <- Slackia_A %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Slackia_A <- merge(study_group_mean_Slackia_A, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Brachyspira <- merge(country_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Brachyspira <- merge(sex_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Brachyspira <- merge(ethnicity_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Brachyspira <- merge(country_group_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 


region_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Brachyspira <- merge(region_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Brachyspira <- Brachyspira %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Brachyspira <- merge(study_group_mean_Brachyspira, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_RS_D84 <- merge(country_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_RS_D84 <- merge(sex_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_RS_D84 <- merge(ethnicity_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_RS_D84 <- merge(country_group_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_RS_D84 <- merge(region_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_RS_D84 <- RS_D84 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_RS_D84 <- merge(study_group_mean_RS_D84, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 



country_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Lachnospira <- merge(country_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Lachnospira <- merge(sex_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Lachnospira <- merge(ethnicity_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Lachnospira <- merge(country_group_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Lachnospira <- merge(region_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Lachnospira <- Lachnospira %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Lachnospira <- merge(study_group_mean_Lachnospira, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 



country_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_CAG45 <- merge(country_mean_CAG45, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_CAG45 <- merge(sex_mean_CAG45, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_CAG45 <- merge(ethnicity_mean_CAG45, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_CAG45 <- merge(country_group_mean_CAG45, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_CAG45 <- merge(region_mean_CAG45, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_CAG45 <- CAG45 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_CAG45 <- merge(study_group_mean_CAG45, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 



country_mean_Butyricicoccus_A <- Butyricicoccus_A %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Butyricicoccus_A <- merge(country_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Butyricicoccus_A <- Butyricicoccus_A %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Butyricicoccus_A <- merge(sex_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Butyricicoccus_A <- Butyricicoccus_A %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Butyricicoccus_A <- merge(ethnicity_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Butyricicoccus_A <- Butyricicoccus_A %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Butyricicoccus_A <- merge(country_group_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Butyricicoccus_A <- Butyricicoccus_A %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Butyricicoccus_A <- merge(region_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Butyricicoccus_A <- UBA4644 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Butyricicoccus_A <- merge(study_group_mean_Butyricicoccus_A, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_WG1 <- WG1 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_WG1 <- merge(country_mean_WG1, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_WG1 <- WG1 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_WG1 <- merge(sex_mean_WG1, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_WG1 <- WG1 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_WG1 <- merge(ethnicity_mean_WG1, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_WG1 <- WG1 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_WG1 <- merge(country_group_mean_WG1, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_WG1 <- WG1 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_WG1 <- merge(region_mean_WG1, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_WG1 <- WG1 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_WG1 <- merge(study_group_mean_WG1, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA3375 <- merge(country_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA3375 <- merge(sex_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA3375 <- merge(ethnicity_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA3375 <- merge(country_group_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA3375 <- merge(region_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA3375 <- UBA3375 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA3375 <- merge(study_group_mean_UBA3375, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Anaerovibrio <- merge(country_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Anaerovibrio <- merge(sex_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Anaerovibrio <- merge(ethnicity_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Anaerovibrio <- merge(country_group_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Anaerovibrio <- merge(region_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Anaerovibrio <- Anaerovibrio %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Anaerovibrio <- merge(study_group_mean_Anaerovibrio, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Treponema_D <- merge(country_mean_Treponema_D, cog, by = "cog_category") %>%
  select(country, cog_category, mean_per_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Treponema_D <- merge(sex_mean_Treponema_D, cog, by = "cog_category") %>%
  select(sex, cog_category, mean_per_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Treponema_D <- merge(ethnicity_mean_Treponema_D, cog, by = "cog_category") %>%
  select(ethnicity, cog_category, mean_per_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Treponema_D <- merge(country_group_mean_Treponema_D, cog, by = "cog_category") %>%
  select(country_group, cog_category, mean_per_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Treponema_D <- merge(region_mean_Treponema_D, cog, by = "cog_category") %>%
  select(region, cog_category, mean_per_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Treponema_D <- Treponema_D %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Treponema_D <- merge(study_group_mean_Treponema_D, cog, by = "cog_category") %>%
  select(study_group, cog_category, mean_per_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA2804 <- merge(country_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA2804 <- merge(sex_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA2804 <- merge(ethnicity_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA2804 <- merge(country_group_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA2804 <- merge(region_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA2804 <- UBA2804 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA2804 <- merge(study_group_mean_UBA2804, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Campylobacter <- merge(country_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Campylobacter <- merge(sex_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Campylobacter <- merge(ethnicity_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Campylobacter <- merge(country_group_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Campylobacter <- merge(region_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Campylobacter <- Campylobacter %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Campylobacter <- merge(study_group_mean_Campylobacter, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_RUG410 <- merge(country_mean_RUG410, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_RUG410 <- merge(sex_mean_RUG410, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_RUG410 <- merge(ethnicity_mean_RUG410, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_RUG410 <- merge(country_group_mean_RUG410, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_RUG410 <- merge(region_mean_RUG410, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_RUG410 <- RUG410 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_RUG410 <- merge(study_group_mean_RUG410, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_RF16 <- RF16 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_RF16 <- merge(country_mean_RF16, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_RF16 <- RF16 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_RF16 <- merge(sex_mean_RF16, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_RF16 <- RF16 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_RF16 <- merge(ethnicity_mean_RF16, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_RF16 <- RF16 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_RF16 <- merge(country_group_mean_RF16, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_RF16 <- RF16 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_RF16 <- merge(region_mean_RF16, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_RF16 <- RF16 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_RF16 <- merge(study_group_mean_RF16, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_CAG568 <- merge(country_mean_CAG568, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_CAG568 <- merge(sex_mean_CAG568, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_CAG568 <- merge(ethnicity_mean_CAG568, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_CAG568 <- merge(country_group_mean_CAG568, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_CAG568 <- merge(region_mean_CAG568, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_CAG568 <- CAG568 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_CAG568 <- merge(study_group_mean_CAG568, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Campylobacter_D <- merge(country_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Campylobacter_D <- merge(sex_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Campylobacter_D <- merge(ethnicity_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Campylobacter_D <- merge(country_group_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Campylobacter_D <- merge(region_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Campylobacter_D <- Campylobacter_D %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Campylobacter_D <- merge(study_group_mean_Campylobacter_D, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA636 <- merge(country_mean_UBA636, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA636 <- merge(sex_mean_UBA636, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA636 <- merge(ethnicity_mean_UBA636, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA636 <- merge(country_group_mean_UBA636, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA636 <- merge(region_mean_UBA636, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA636 <- UBA636 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA636 <- merge(study_group_mean_UBA636, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_Succinivibrio <- merge(country_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_Succinivibrio <- merge(sex_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_Succinivibrio <- merge(ethnicity_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_Succinivibrio <- merge(country_group_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_Succinivibrio <- merge(region_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_Succinivibrio <- Succinivibrio %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_Succinivibrio <- merge(study_group_mean_Succinivibrio, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA1436 <- merge(country_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA1436 <- merge(sex_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA1436 <- merge(ethnicity_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA1436 <- merge(country_group_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA1436 <- merge(region_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA1436 <- UBA1436 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA1436 <- merge(study_group_mean_UBA1436, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_UBA71 <- merge(country_mean_UBA71, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_UBA71 <- merge(sex_mean_UBA71, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_UBA71 <- merge(ethnicity_mean_UBA71, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_UBA71 <- merge(country_group_mean_UBA71, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_UBA71 <- merge(region_mean_UBA71, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_UBA71 <- UBA71 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_UBA71 <- merge(study_group_mean_UBA71, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 

country_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_RUG572 <- merge(country_mean_RUG572, cog, by = "cog_category") %>%
  group_by(country) %>%
  mutate(total_per_country = sum(mean_per_country, group_by = TRUE)) %>%
  mutate(percent_by_country = mean_per_country/total_per_country) %>% 
  group_by(country) %>% 
  select(country, cog_category, mean_per_country, total_per_country, percent_by_country, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 

sex_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

sex_RUG572 <- merge(sex_mean_RUG572, cog, by = "cog_category") %>%
  group_by(sex) %>%
  mutate(total_per_sex = sum(mean_per_sex, group_by = TRUE)) %>%
  mutate(percent_by_sex = mean_per_sex/total_per_sex) %>% 
  group_by(sex) %>% 
  select(sex, cog_category, mean_per_sex, total_per_sex, percent_by_sex, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 

ethnicity_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

ethnicity_RUG572 <- merge(ethnicity_mean_RUG572, cog, by = "cog_category") %>%
  group_by(ethnicity) %>%
  mutate(total_per_ethnicity = sum(mean_per_ethnicity, group_by = TRUE)) %>%
  mutate(percent_by_ethnicity = mean_per_ethnicity/total_per_ethnicity) %>% 
  group_by(ethnicity) %>% 
  select(ethnicity, cog_category, mean_per_ethnicity, total_per_ethnicity, percent_by_ethnicity, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

country_group_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

country_group_RUG572 <- merge(country_group_mean_RUG572, cog, by = "cog_category") %>%
  group_by(country_group) %>%
  mutate(total_per_country_group = sum(mean_per_country_group, group_by = TRUE)) %>%
  mutate(percent_by_country_group = mean_per_country_group/total_per_country_group) %>% 
  group_by(country_group) %>% 
  select(country_group, cog_category, mean_per_country_group, total_per_country_group, percent_by_country_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 

region_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

region_RUG572 <- merge(region_mean_RUG572, cog, by = "cog_category") %>%
  group_by(region) %>%
  mutate(total_per_region = sum(mean_per_region, group_by = TRUE)) %>%
  mutate(percent_by_region = mean_per_region/total_per_region) %>% 
  group_by(region) %>% 
  select(region, cog_category, mean_per_region, total_per_region, percent_by_region, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 

study_group_mean_RUG572 <- RUG572 %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(tally)) %>% #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 
  mutate_if(is.numeric, signif, digits = 3)

study_group_RUG572 <- merge(study_group_mean_RUG572, cog, by = "cog_category") %>%
  group_by(study_group) %>%
  mutate(total_per_study_group = sum(mean_per_study_group, group_by = TRUE)) %>%
  mutate(percent_by_study_group = mean_per_study_group/total_per_study_group) %>% 
  group_by(study_group) %>% 
  select(study_group, cog_category, mean_per_study_group, total_per_study_group, percent_by_study_group, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


#okay let's plot ----

site_color <- (c(wes_palette("Darjeeling1", 5, type = c("discrete")), wes_palette("Darjeeling2", 4, type = c("discrete")), wes_palette("Royal1", 4, type = c("discrete")), wes_palette("Royal2", 5, type = c("discrete")), wes_palette("Rushmore1", 3, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete")), wes_palette("Moonrise3", 5, type = c("discrete"))))

show_col(site_color)

#"WG-1", "UBA3375","Anaerovibrio", "Treponema_D", "UBA2804", "Campylobacter", "RUG410", "RF16", "CAG-568", "Campylobacter_D", "UBA636", "Succinivibrio", "UBA1436", "UBA71", "RUG572", "Butyricicoccus_A", "CAG-45", "Lachnospira", "RS-D84", "Brachyspira", "Slackia_A", "CAG-533", "UBA4644"

#UBA4644

#plot by country

country_UBA4644 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA4644.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA4644 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA4644.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA4644 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA4644.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA4644 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA4644.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA4644 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA4644.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA4644 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA4644.pdf", width = 40, height = 20, units = "cm")



#CAG533

#plot by country

country_CAG533 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_CAG533.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_CAG533 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_CAG533.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_CAG533 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_CAG533.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_CAG533 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_CAG533.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_CAG533 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_CAG533.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_CAG533 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_CAG533.pdf", width = 40, height = 20, units = "cm")



#Slackia_A

#plot by country

country_Slackia_A %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Slackia_A.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Slackia_A %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Slackia_A.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Slackia_A %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Slackia_A.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Slackia_A %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Slackia_A.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Slackia_A %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Slackia_A.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Slackia_A %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Slackia_A.pdf", width = 40, height = 20, units = "cm")



#Brachyspira

#plot by country

country_Brachyspira %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Brachyspira.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Brachyspira %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Brachyspira.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Brachyspira %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Brachyspira.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Brachyspira %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Brachyspira.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Brachyspira %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Brachyspira.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Brachyspira %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Brachyspira.pdf", width = 40, height = 20, units = "cm")



#RS_D84

#plot by country

country_RS_D84 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_RS_D84.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_RS_D84 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_RS_D84.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_RS_D84 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_RS_D84.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_RS_D84 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_RS_D84.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_RS_D84 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_RS_D84.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_RS_D84 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_RS_D84.pdf", width = 40, height = 20, units = "cm")



#Lachnospira

#plot by country

country_Lachnospira %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Lachnospira.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Lachnospira %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Lachnospira.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Lachnospira %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Lachnospira.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Lachnospira %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Lachnospira.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Lachnospira %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Lachnospira.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Lachnospira %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Lachnospira.pdf", width = 40, height = 20, units = "cm")



#CAG45

#plot by country

country_CAG45 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_CAG45.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_CAG45 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_CAG45.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_CAG45 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_CAG45.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_CAG45 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_CAG45.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_CAG45 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_CAG45.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_CAG45 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_CAG45.pdf", width = 40, height = 20, units = "cm")



#Butyricicoccus_A

#plot by country

country_Butyricicoccus_A %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Butyricicoccus_A %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Butyricicoccus_A %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Butyricicoccus_A %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Butyricicoccus_A %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Butyricicoccus_A %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Butyricicoccus_A.pdf", width = 40, height = 20, units = "cm")



#RUG572

#plot by country

country_RUG572 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_RUG572.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_RUG572 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_RUG572.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_RUG572 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_RUG572.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_RUG572 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_RUG572.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_RUG572 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_RUG572.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_RUG572 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_RUG572.pdf", width = 40, height = 20, units = "cm")



#UBA71

#plot by country

country_UBA71 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA71.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA71 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA71.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA71 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA71.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA71 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA71.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA71 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA71.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA71 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA71.pdf", width = 40, height = 20, units = "cm")



#UBA1436

#plot by country

country_UBA1436 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA1436.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA1436 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA1436.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA1436 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA1436.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA1436 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA1436.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA1436 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA1436.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA1436 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA1436.pdf", width = 40, height = 20, units = "cm")



#Succinivibrio

#plot by country

country_Succinivibrio %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Succinivibrio.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Succinivibrio %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Succinivibrio.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Succinivibrio %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Succinivibrio.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Succinivibrio %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Succinivibrio.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Succinivibrio %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Succinivibrio.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Succinivibrio %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Succinivibrio.pdf", width = 40, height = 20, units = "cm")



#UBA636

#plot by country

country_UBA636 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA636.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA636 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA636.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA636 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA636.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA636 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA636.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA636 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA636.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA636 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA636.pdf", width = 40, height = 20, units = "cm")



#Campylobacter_D

#plot by country

country_Campylobacter_D %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Campylobacter_D %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Campylobacter_D %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Campylobacter_D %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Campylobacter_D %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Campylobacter_D %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Campylobacter_D.pdf", width = 40, height = 20, units = "cm")



#CAG568

#plot by country

country_CAG568 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_CAG568.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_CAG568 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_CAG568.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_CAG568 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_CAG568.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_CAG568 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_CAG568.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_CAG568 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_CAG568.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_CAG568 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_CAG568.pdf", width = 40, height = 20, units = "cm")



#RF16

#plot by country

country_RF16 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_RF16.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_RF16 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_RF16.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_RF16 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_RF16.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_RF16 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_RF16.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_RF16 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_RF16.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_RF16 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_RF16.pdf", width = 40, height = 20, units = "cm")




#RUG410

#plot by country

country_RUG410 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_RUG410.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_RUG410 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_RUG410.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_RUG410 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_RUG410.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_RUG410 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_RUG410.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_RUG410 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_RUG410.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_RUG410 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_RUG410.pdf", width = 40, height = 20, units = "cm")



#Campylobacter

#plot by country

country_Campylobacter %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Campylobacter.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Campylobacter %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Campylobacter.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Campylobacter %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Campylobacter.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Campylobacter %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Campylobacter.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Campylobacter %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Campylobacter.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Campylobacter %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Campylobacter.pdf", width = 40, height = 20, units = "cm")



#UBA2804

#plot by country

country_UBA2804 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA2804.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA2804 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA2804.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA2804 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA2804.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA2804 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA2804.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA2804 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA2804.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA2804 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA2804.pdf", width = 40, height = 20, units = "cm")



#Treponema_D

#plot by country

country_Treponema_D %>% ggplot(aes(x=factor(country, levels = unique(country)), y=mean_per_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=mean_per_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Treponema_D.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Treponema_D %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Treponema_D.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Treponema_D %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=mean_per_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=mean_per_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Treponema_D.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Treponema_D %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Treponema_D.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Treponema_D %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=mean_per_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=mean_per_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Treponema_D.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Treponema_D %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_country_group_Treponema_D.pdf", width = 40, height = 20, units = "cm")



#UBA3375

#plot by country

country_UBA3375 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_UBA3375.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_UBA3375 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_UBA3375.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_UBA3375 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_UBA3375.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_UBA3375 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_UBA3375.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_UBA3375 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_UBA3375.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_UBA3375 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_UBA3375.pdf", width = 40, height = 20, units = "cm")



#WG1

#plot by country

country_WG1 %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_WG1.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_WG1 %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_WG1.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_WG1 %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_WG1.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_WG1 %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_WG1.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_WG1 %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_WG1.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_WG1 %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_WG1.pdf", width = 40, height = 20, units = "cm")



#Anaerovibrio

#plot by country

country_Anaerovibrio %>% ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_country_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")


#plot by sex

sex_Anaerovibrio %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=percent_by_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_sex_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")

#plot by ethnicity general

ethnicity_Anaerovibrio %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=percent_by_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) + 
  geom_text(aes(label=percent_by_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_ethnicity_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")

#plot by region

region_Anaerovibrio %>% ggplot(aes(x=factor(region, levels = unique(region)), y=percent_by_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))


ggsave("output/function_DA/cog_function_by_region_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")

#plot by study_group

study_group_Anaerovibrio %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=percent_by_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  geom_text(aes(label=percent_by_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_DA/cog_function_by_study_group_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")

#plot by country_group

country_group_Anaerovibrio %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=percent_by_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 15)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=percent_by_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 1, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 

ggsave("output/function_DA/cog_function_by_country_group_Anaerovibrio.pdf", width = 40, height = 20, units = "cm")

