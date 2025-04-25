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

thecolors = createPalette(130,  safe_pal)

show_col(thecolors)



#load data 

data_otu <- read.csv("data/OTU_tables/bracken_PERCENT_sg_2024_decontam_S.csv", header = TRUE)
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
  pivot_longer(names_to="species", values_to = "percent", s_prevotella_sp900557255:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


#plot! 

stacked %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() +  
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors)) +
  theme(axis.text.x=element_blank()) + 
  theme(legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 0.5))) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"))

ggsave("output/s_stackedbar0.05.pdf", width = 40, height = 20, units = "cm")

#by country 

stacked_country <- otu_meta %>% 
  group_by(country) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", s_prevotella_sp900557255:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_country %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) + 
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors)) +
  theme(axis.text.x=element_blank()) + 
  facet_grid(~fct_inorder(country),switch="x",scales="free",space="free") + 
  theme(legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 0.5))) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"))

ggsave("output/s_stackedbar0.05country.pdf", width = 40, height = 20, units = "cm")

#by region

stacked_region <- otu_meta %>% 
  group_by(region) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", s_prevotella_sp900557255:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_region %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) + 
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors)) +
  theme(axis.text.x=element_blank()) + 
  facet_wrap(~fct_inorder(region),strip.position ="bottom",scales="free") + 
  theme(legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 0.5))) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"))

ggsave("output/s_stackedbar0.05region.pdf", width = 40, height = 20, units = "cm")

#by study group

stacked_group <- otu_meta %>% 
  group_by(study_group) %>%
  arrange(all_low, .by_group = TRUE) %>%
  pivot_longer(names_to="species", values_to = "percent", s_prevotella_sp900557255:all_low) %>% 
  filter_if(is.numeric, all_vars((.) != 0))


stacked_group %>% ggplot(aes(x=factor(lab_id, levels = unique(lab_id)),y=percent, fill=species)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(text = element_text(size = 15)) + 
  labs(x="Study participants",y="Microbiome composition", fill="Microbial taxa") + 
  scale_fill_manual(values=as.vector(thecolors)) +
  theme(axis.text.x=element_blank()) + 
  facet_grid(~fct_inorder(study_group),switch="x",scales="free",space="free") + 
  theme(legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 0.5))) +    
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"))

ggsave("output/s_stackedbar0.05group.pdf", width = 40, height = 20, units = "cm")

