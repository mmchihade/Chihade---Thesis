#load packages 
library(ggplot2)
library(reshape2)
library(ggrepel)
library(grid)
library(gridExtra)
library(pals)
library(ggcorrplot)
library(ggpmisc)
suppressMessages(library(Hmisc)) #not sure if needed, but there is some ggplot interaction
suppressMessages(library(tidyverse)) #note that installing tidyverse takes a LONG time!
library(dplyr)
library(purrr)
library(forcats)
library(ggthemes)
library(data.table)
library(TSP)
library(Polychrome)
library(rcartocolor)

#make color vector 

safe_pal <- carto_pal(12, "Safe")

thecolors = createPalette(65,  safe_pal)

show_col(thecolors)



#load data 

data_otu <- read.csv("data/OTU_tables/bracken_PERCENT_sg_2024_decontam_G.csv", header = TRUE)
data_meta <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)

#create a version of the otu table that only shows genera with low percentages and then create a new colunm that totals these low values. finally, select only that column and the lab_id 

data_otu_low <- data_otu %>% 
  select(-total_root) %>% 
  clean_names() %>% 
  rename(lab_id = subject) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  column_to_rownames("lab_id") %>% 
  mutate(across(everything(), ~replace(., . > 0.05 , 0))) %>%
  mutate(all_low = rowSums(pick(where(is.numeric)))) %>%
  select(all_low) %>%
  rownames_to_column("lab_id")

#now create a version of the otu table that only shows genera with high percentages

data_otu_high <- data_otu %>% 
  select(-total_root) %>% 
  clean_names() %>% 
  rename(lab_id = subject) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  column_to_rownames("lab_id") %>% 
  mutate(across(everything(), ~replace(., . < 0.05 , 0))) %>%
  rownames_to_column("lab_id")

#add those two together 

otu_percent <- merge(data_otu_high, data_otu_low, by = "lab_id")

#add this with metadata

otu_meta <- merge(otu_percent, data_meta, by = "lab_id") %>%
  select(-record_id, -field_id, -topmed_id, -collection_name, -sample_site, -longitude, -latitude, -date, -year, -sex, -age, -is_adult, -subsistance, -waist_circumference, -bmi, -ethnicity, -body_fat_percent)

#use pivot longer to make a plotable data frame 

stacked <- otu_meta %>% 
  arrange(all_low) %>% #arranging the study participants by how much low variablity they have
  pivot_longer(names_to="species", values_to = "percent", g_prevotella:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


#plot! 

stacked %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) +
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors), breaks=c("all_low", "g_prevotella", "g_escherichia", "g_bacteroides", "g_wg_1", "g_treponema_d", "g_phocaeicola", "g_agathobacter", "g_akkermansia", "g_alistipes", "g_bifidobacterium", "g_blautia_a", "g_butyrivibrio_a", "g_cag_103", "g_cag_110", "g_cag_127", "g_cag_170", "g_cag_177", "g_cag_180", "g_cag_279", "g_cag_303", "g_cag_349", "g_cag_353", "g_cag_475", "g_cag_510", "g_cag_568", "g_cag_632", "g_cag_83", "g_cag_95", "g_coe1", "g_coprococcus", "g_dialister", "g_eisenbergiella", "g_enterocloster", "g_er4", "g_f082", "g_faecalibacterium", "g_gemmiger", "g_klebsiella", "g_lachnospira", "g_megamonas", "g_mitsuokella", "g_parabacteroides", "g_paraprevotella", "g_phocaeicola", "g_prevotellamassilia", "g_rc9", "g_rf16", "g_roseburia", "g_rug115", "g_rug410", "g_ruminiclostridium_e", "g_ruminococcus_d", "g_ruminococcus_e", "g_sfdb01", "g_sodaliphilus", "g_succinivibrio", "g_treponema_d", "g_uba10281", "g_uba11524", "g_uba1777", "g_uba4248", "g_uba4372", "g_uba7173", "g_uba9732", "g_vsob01")) +
   theme(axis.text.x=element_blank()) +    theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"))


ggsave("output/g_stackedbar0.05.pdf", width = 40, height = 20, units = "cm")

#by country 

stacked_country <- otu_meta %>% 
  group_by(country) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", g_prevotella:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_country %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) +
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors), breaks=c("all_low", "g_prevotella", "g_escherichia", "g_bacteroides", "g_wg_1", "g_treponema_d", "g_phocaeicola", "g_agathobacter", "g_akkermansia", "g_alistipes", "g_bifidobacterium", "g_blautia_a", "g_butyrivibrio_a", "g_cag_103", "g_cag_110", "g_cag_127", "g_cag_170", "g_cag_177", "g_cag_180", "g_cag_279", "g_cag_303", "g_cag_349", "g_cag_353", "g_cag_475", "g_cag_510", "g_cag_568", "g_cag_632", "g_cag_83", "g_cag_95", "g_coe1", "g_coprococcus", "g_dialister", "g_eisenbergiella", "g_enterocloster", "g_er4", "g_f082", "g_faecalibacterium", "g_gemmiger", "g_klebsiella", "g_lachnospira", "g_megamonas", "g_mitsuokella", "g_parabacteroides", "g_paraprevotella", "g_phocaeicola", "g_prevotellamassilia", "g_rc9", "g_rf16", "g_roseburia", "g_rug115", "g_rug410", "g_ruminiclostridium_e", "g_ruminococcus_d", "g_ruminococcus_e", "g_sfdb01", "g_sodaliphilus", "g_succinivibrio", "g_treponema_d", "g_uba10281", "g_uba11524", "g_uba1777", "g_uba4248", "g_uba4372", "g_uba7173", "g_uba9732", "g_vsob01")) +
   theme(axis.text.x=element_blank()) +    theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_grid(~fct_inorder(country),switch="x",scales="free",space="free")

ggsave("output/g_stackedbar0.05country.pdf", width = 40, height = 20, units = "cm")

#by region

stacked_region <- otu_meta %>% 
  group_by(region) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", g_prevotella:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_region %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) +
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors), breaks=c("all_low", "g_prevotella", "g_escherichia", "g_bacteroides", "g_wg_1", "g_treponema_d", "g_phocaeicola", "g_agathobacter", "g_akkermansia", "g_alistipes", "g_bifidobacterium", "g_blautia_a", "g_butyrivibrio_a", "g_cag_103", "g_cag_110", "g_cag_127", "g_cag_170", "g_cag_177", "g_cag_180", "g_cag_279", "g_cag_303", "g_cag_349", "g_cag_353", "g_cag_475", "g_cag_510", "g_cag_568", "g_cag_632", "g_cag_83", "g_cag_95", "g_coe1", "g_coprococcus", "g_dialister", "g_eisenbergiella", "g_enterocloster", "g_er4", "g_f082", "g_faecalibacterium", "g_gemmiger", "g_klebsiella", "g_lachnospira", "g_megamonas", "g_mitsuokella", "g_parabacteroides", "g_paraprevotella", "g_phocaeicola", "g_prevotellamassilia", "g_rc9", "g_rf16", "g_roseburia", "g_rug115", "g_rug410", "g_ruminiclostridium_e", "g_ruminococcus_d", "g_ruminococcus_e", "g_sfdb01", "g_sodaliphilus", "g_succinivibrio", "g_treponema_d", "g_uba10281", "g_uba11524", "g_uba1777", "g_uba4248", "g_uba4372", "g_uba7173", "g_uba9732", "g_vsob01")) +
   theme(axis.text.x=element_blank()) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_wrap(~fct_inorder(region),strip.position ="bottom",scales="free")

ggsave("output/g_stackedbar0.05region.pdf", width = 40, height = 20, units = "cm")

#by study group

stacked_group <- otu_meta %>% 
  group_by(study_group) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", g_prevotella:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_group %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) +
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors), breaks=c("all_low", "g_prevotella", "g_escherichia", "g_bacteroides", "g_wg_1", "g_treponema_d", "g_phocaeicola", "g_agathobacter", "g_akkermansia", "g_alistipes", "g_bifidobacterium", "g_blautia_a", "g_butyrivibrio_a", "g_cag_103", "g_cag_110", "g_cag_127", "g_cag_170", "g_cag_177", "g_cag_180", "g_cag_279", "g_cag_303", "g_cag_349", "g_cag_353", "g_cag_475", "g_cag_510", "g_cag_568", "g_cag_632", "g_cag_83", "g_cag_95", "g_coe1", "g_coprococcus", "g_dialister", "g_eisenbergiella", "g_enterocloster", "g_er4", "g_f082", "g_faecalibacterium", "g_gemmiger", "g_klebsiella", "g_lachnospira", "g_megamonas", "g_mitsuokella", "g_parabacteroides", "g_paraprevotella", "g_phocaeicola", "g_prevotellamassilia", "g_rc9", "g_rf16", "g_roseburia", "g_rug115", "g_rug410", "g_ruminiclostridium_e", "g_ruminococcus_d", "g_ruminococcus_e", "g_sfdb01", "g_sodaliphilus", "g_succinivibrio", "g_treponema_d", "g_uba10281", "g_uba11524", "g_uba1777", "g_uba4248", "g_uba4372", "g_uba7173", "g_uba9732", "g_vsob01")) +
  theme(axis.text.x=element_blank()) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_grid(~fct_inorder(study_group),switch="x",scales="free",space="free")

ggsave("output/g_stackedbar0.05group.pdf", width = 40, height = 20, units = "cm")

stacked_prevotella <- otu_meta %>% 
  arrange(all_low) %>%
  pivot_longer(names_to="species", values_to = "percent", g_prevotella:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_prevotella %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) +
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors), breaks=c("all_low", "g_prevotella")) +
  theme(axis.text.x=element_blank()) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_grid(~fct_inorder(study_group),switch="x",scales="free",space="free")

ggsave("output/g_stackedbar0.05prev.pdf", width = 40, height = 20, units = "cm")