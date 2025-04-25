#load packages 
library(ggplot2)
library(reshape2)
library(pals)
library(ggcorrplot)
library(ggpmisc)
suppressMessages(library(Hmisc)) #not sure if needed, but there is some ggplot interaction
suppressMessages(library(tidyverse)) #note that installing tidyverse takes a LONG time!
library(dplyr)
library(purrr)
library(forcats)
library(data.table)
library(TSP)
library(rcartocolor)

#make color vector 

safe_pal <- carto_pal(12, "Safe")

thecolors = as.vector(createPalette(65,  safe_pal))

show_col(thecolors)

#load data 

data_otu <- read.csv("data/OTU_tables/bracken_PERCENT_sg_2024_decontam_G.csv", header = TRUE) %>% 
  clean_names() %>% 
  select(-total_root)

data_meta <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)

#create a version of the otu table that only shows genera with low percentages and then create a new colunm that totals these low values. finally, select only that column and the lab_id 

data_otu_low <- data_otu %>% 
  column_to_rownames("subject") %>% 
  mutate(across(everything(), ~replace(., . > 0.1 , 0))) %>%
  mutate(all_low = rowSums(across(where(is.numeric)), na.rm = T)) %>%
  select(all_low) %>%
  rownames_to_column("lab_id")

#now create a version of the otu table that only shows genera with high percentages

data_otu_high <- data_otu %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  column_to_rownames("subject") %>% 
  mutate(across(everything(), ~replace(., . < 0.1 , 0))) %>%
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
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD", 
                               g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481", 
                               g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9" , g_acetatifactor="#B30063"),
                    breaks=c("all_low", "g_prevotella", "g_escherichia", "g_bacteroides", "g_rc9", "g_treponema_d")) +
  theme(axis.text.x=element_blank()) +  
  theme(axis.text.y=element_text(size=20), axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(keywidth = 1.5, label.theme = element_text(size = 15))) + 
  theme(legend.position="bottom")


ggsave("output/g_stackedbar0.1.pdf", width = 40, height = 20, units = "cm")

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
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD", g_faecalibacterium="#664F26", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481", 
                               g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9", , g_acetatifactor="#B30063"), breaks=c("all_low", "g_prevotella","g_bacteroides", "g_treponema_d")) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_grid(~fct_inorder(country),switch="x",scales="free",space="free") + 
  theme(legend.position="bottom")

ggsave("output/g_stackedbar0.1country.pdf", width = 40, height = 20, units = "cm")

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
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD", g_faecalibacterium="#664F26", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481", 
                               g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9", , g_acetatifactor="#B30063"), breaks=c("all_low", "g_prevotella","g_bacteroides", "g_treponema_d")) +
  theme(axis.text.x=element_blank()) +   
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_wrap(~fct_inorder(region),strip.position ="bottom",scales="free") + 
  theme(legend.position="bottom")

ggsave("output/g_stackedbar0.1region.pdf", width = 40, height = 20, units = "cm")

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
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD", g_faecalibacterium="#664F26", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481", 
                              g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9", , g_acetatifactor="#B30063"), breaks=c("all_low", "g_prevotella","g_bacteroides", "g_treponema_d")) +
  theme(axis.text.x=element_blank()) +   
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) + 
  facet_grid(~fct_inorder(study_group),switch="x",scales="free",space="free") + 
  guides(fill = guide_legend(keywidth = 1.5, label.theme = element_text(size = 15))) + 
  theme(strip.text = element_text(size = 15))  + 
  theme(legend.position="bottom")

ggsave("output/g_stackedbar0.1group.pdf", width = 40, height = 20, units = "cm")

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

ggsave("output/g_stackedbar0.1prev.pdf", width = 40, height = 20, units = "cm")
