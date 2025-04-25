# load packages ---- 
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
  rename(id_lab = field_id) %>%
  mutate(id_lab = toupper(id_lab))

#create a dataframe with all of the cog tally files, with the file name as a column so that we can match column to tally

merged_df <- list.files(path = "data/cog_tally_overall",
                        pattern = "\\.csv$", 
                        full.names = TRUE) %>% 
  set_names() %>%
  map_dfr(read_csv, .id = "file_name") 

#we need the cog categories too! 

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#now let's clean the data up  

cog_tally <- clean_names(merged_df) %>% 
  rename(lab_id = file_name) %>% 
  select(-genus) %>% 
  mutate(across(everything(), gsub, pattern = "data/cog_tally_overall/", replacement = "")) %>%
  mutate(across(everything(), gsub, pattern = ".csv", replacement = "")) %>%
  mutate(n = as.numeric(n)) %>%
  group_by(lab_id, cog_category) %>% 
  mutate(tally = sum(n, na.rm =TRUE, group_by=TRUE)) %>% 
  select(-n) %>%
  distinct()

#I also have to replace all of the TZ codes with B codes so that I can divide by country and ethnicity 

tally <- cog_tally %>% 
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
                         "TZHZ340"="B091",
                         "TZHZ346"="B092",
                         "TZHZ403"="B202",
                         "TZHZ407"="B206",
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
  rename(id_lab = lab_id) %>% 
  mutate(cog_category = recode(cog_category,
                                    "?" = "S",
                                    "!" = "S")) %>% 
  group_by(id_lab, cog_category) %>% 
  mutate(tally = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#I also want to separate all of the multiple letter columns. They will probably be useful later on, but it is hella overwhelming now
full_tally <- tally %>% 
  separate(cog_category, into = c("cog_1", "cog_2", "cog_3", "cog_4", "cog_5", "cog_6", "cog_7"), sep = "") %>%
  select(-cog_1) %>% 
  pivot_longer(names_to="cog_name", values_to = "cog_category", cog_2:cog_7)%>%
  na.omit() %>%
  select(-cog_name) %>%
  group_by(id_lab, cog_category) %>% 
  mutate(tally = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()
  
#umm, I am going to try to make a combined dataset, let's see if R can handle it ----

merge <- merge(full_tally, data_grp, by = "id_lab")

merge_labels <- merge(merge, cog, by = "cog_category") %>% 
  filter(!cog_category %in% c("Y", "B", "Z"))

#yay!!!!!

#let's load the p values for all our factors of interest  ----

pv_all <-  read.csv("data/cog_variance_values.csv", header=TRUE)

#and merge them with the labels

pv_merge <- merge(merge_labels, pv_all, by = "cog_category") %>% 
  mutate_if(is.numeric, signif, digits = 3) %>% 
  mutate(label = stringr::str_wrap(label, width = 30)) 
  

#let's group by a bunch of factors and average within them ----

individual_label <- merge_labels %>% #and adding back the important info here
  select(cog_category, id_lab, tally, label) %>%
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3)

individual <- individual_label %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(id_lab) 

country_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, country, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, country) %>% 
  summarise(mean_per_country = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

country <- merge(country_mean, pv_merge, by=c("cog_category", "country")) %>% #and adding back the important info here
  select(cog_category, country, mean_per_country, label, pvcountry) %>% 
  mutate(label = paste(label, "(p =", pvcountry, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country) 


sex_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, sex, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, sex) %>% 
  summarise(mean_per_sex = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

sex <- merge(sex_mean, pv_merge, by=c("cog_category", "sex")) %>% #and adding back the important info here
  select(cog_category, sex, mean_per_sex, label, pvsex) %>% 
  mutate(label = paste(label, "(p =", pvsex, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(sex) 


ethnicity_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, ethnicity, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, ethnicity) %>% 
  summarise(mean_per_ethnicity = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

ethnicity <- merge(ethnicity_mean, pv_merge, by=c("cog_category", "ethnicity")) %>% #and adding back the important info here
  select(cog_category, ethnicity, mean_per_ethnicity, label, pvethnicity) %>% 
  mutate(label = paste(label, "(p =", pvethnicity, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(ethnicity) 

  
region_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, region, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, region) %>% 
  summarise(mean_per_region = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

region <- merge(region_mean, pv_merge, by=c("cog_category", "region")) %>% #and adding back the important info here
  select(cog_category, region, mean_per_region, label, pvregion) %>% 
  mutate(label = paste(label, "(p =", pvregion, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(region) 


study_group_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, study_group, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, study_group) %>% 
  summarise(mean_per_study_group = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

study_group <- merge(study_group_mean, pv_merge, by=c("cog_category", "study_group")) %>% #and adding back the important info here
  select(cog_category, study_group, mean_per_study_group, label, pvstudy_group) %>% 
  mutate(label = paste(label, "(p =", pvstudy_group, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(study_group) 


country_group_mean <- pv_merge %>%
  group_by(id_lab) %>%
  mutate(total_per_ind = sum(tally, group_by = TRUE)) %>%
  mutate(percent_by_ind = tally/total_per_ind) %>% 
  group_by(id_lab) %>% 
  select(id_lab, country_group, cog_category, tally, total_per_ind, percent_by_ind, label) %>% 
  distinct() %>%
  group_by(cog_category, country_group) %>% 
  summarise(mean_per_country_group = mean(percent_by_ind)) #for some reason, if I use the mutate function to find the mean, I cannot get an accurate total. For this reason, I am using summarise here. 

country_group <- merge(country_group_mean, pv_merge, by=c("cog_category", "country_group")) %>% #and adding back the important info here
  select(cog_category, country_group, mean_per_country_group, label, pvcountry_group) %>% 
  mutate(label = paste(label, "(p =", pvcountry_group, ")", sep = " ")) %>% 
  distinct() %>%
  mutate_if(is.numeric, signif, digits = 3) %>% 
  arrange(country_group) 



#okay let's plot ----

site_color <- (c(wes_palette("Darjeeling1", 4, type = c("discrete")), wes_palette("GrandBudapest2", 2, type = c("discrete")),  wes_palette("Royal2", 1, type = c("discrete")), wes_palette("Darjeeling2", 3, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete"))))

show_col(site_color)

#plot by individual 

individual %>% ggplot(aes(x=factor(id_lab, levels = unique(id_lab)), y=percent_by_ind, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Individual",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.text.x=element_blank(), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null"))

ggsave("output/function_overall/cog_function_by_ind.pdf", width = 40, height = 20, units = "cm")

#plot by country

country %>% ggplot(aes(x=factor(country, levels = unique(country)), y=mean_per_country, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Country",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 30)) + 
  geom_text(aes(label=mean_per_country), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.8, label.theme = element_text(size = 10))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13)))


ggsave("output/function_overall/cog_function_by_country.pdf", width = 40, height = 30, units = "cm")

#plot by sex

sex %>% ggplot(aes(x=factor(sex, levels = unique(sex)), y=mean_per_sex, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Sex",y="Protein Function", fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  geom_text(aes(label=mean_per_sex), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13)))


ggsave("output/function_overall/cog_function_by_sex.pdf", width = 40, height = 30, units = "cm")

#plot by ethnicity general

ethnicity %>% ggplot(aes(x=factor(ethnicity, levels = unique(ethnicity)), y=mean_per_ethnicity, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Ethnicity",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  geom_text(aes(label=mean_per_ethnicity), size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 10)))


ggsave("output/function_overall/cog_function_by_ethnicity.pdf", width = 40, height = 25, units = "cm")

#plot by region

region %>% ggplot(aes(x=factor(region, levels = unique(region)), y=mean_per_region, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Region",y="Protein Function", fill = "COG Category and p value") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  geom_text(aes(label=mean_per_region), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20),, axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 10)))


ggsave("output/function_overall/cog_function_by_region.pdf", width = 40, height = 25, units = "cm")

#plot by study_group

study_group %>% ggplot(aes(x=factor(study_group, levels = unique(study_group)), y=mean_per_study_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Subsistance Strategy",y="Protein Function", fill = "COG Category and p value") + 
  scale_x_discrete(labels = label_wrap(5)) + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  geom_text(aes(label=mean_per_study_group), size = 7, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13))) +    guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 13)))

ggsave("output/function_overall/cog_function_by_study_group.pdf", width = 40, height = 25, units = "cm")

#plot by country_group

country_group %>% ggplot(aes(x=factor(country_group, levels = unique(country_group)), y=mean_per_country_group, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Country and Subsistance Strategy",y="Protein Function", , fill = "COG Category and p value") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.text.x=element_text(size=7)) +
  geom_text(aes(label=mean_per_country_group), size = 5, position = position_stack(vjust = 0.5), check_overlap = TRUE) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +    theme(legend.key.height = unit(1, "null"))

ggsave("output/function_overall/cog_function_by_country_group.pdf", width = 40, height = 25, units = "cm")


test <- data_otu %>% 
 rotate_df() %>% 
  rownames_to_column("taxon") %>% 
  select(taxon) %>%
  mutate(across(everything(), gsub, pattern = "g_", replacement = ""))
  

test2 <- genus_infotest2 <- merge_labels %>% 
  select(genus) %>% 
  mutate(taxon = tolower(genus)) %>% 
  distinct()


 test <- merge(test2, test, by="taxon", all = TRUE)

