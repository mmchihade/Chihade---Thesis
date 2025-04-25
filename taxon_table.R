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
library(stringr)

#load data ----

metadata <- read.csv("data/genomes-all_metadata.csv", header=TRUE)

#let's split data up by the taxonomic categories ----

tax_table <- metadata %>% 
  select(Lineage) %>%
  separate(Lineage, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  clean_names()

tax_table_big <- metadata %>% 
  select(Lineage, Genome) %>%
  separate(Lineage, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  clean_names()

#the big error that we care about here is that prevotella is in the wrong family, let's try to fix that

test <- tax_table %>% 
  filter(family %in% c("f__Bacteroidaceae")) %>% 
  filter(genus %in% c("g__Alloprevotella", "g__Hallella", "g__Hoylesella", "g__Ihuprevotella", "g__Leyella", "g__Marseilla", "g__Massiliprevotella", "g__Metaprevotella", "g__Palleniella", "g__Paraprevotella", "g__Prevotella", "g__Prevotellamassilia", "g__Pseudoprevotella", "g__Segatella", "g__Xylanibacter")) %>% 
  mutate(across('family', str_replace, 'f__Bacteroidaceae', 'f__Prevotellaceae'))

test_big <- tax_table_big %>% 
  filter(family %in% c("f__Bacteroidaceae")) %>% 
  filter(genus %in% c("g__Alloprevotella", "g__Hallella", "g__Hoylesella", "g__Ihuprevotella", "g__Leyella", "g__Marseilla", "g__Massiliprevotella", "g__Metaprevotella", "g__Palleniella", "g__Paraprevotella", "g__Prevotella", "g__Prevotellamassilia", "g__Pseudoprevotella", "g__Segatella", "g__Xylanibacter")) %>% 
  mutate(across('family', str_replace, 'f__Bacteroidaceae', 'f__Prevotellaceae'))

#okay now let's put these fixed taxa into the big table 

taxon_wo <- tax_table %>% 
  filter(!genus %in% c("g__Alloprevotella", "g__Hallella", "g__Hoylesella", "g__Ihuprevotella", "g__Leyella", "g__Marseilla", "g__Massiliprevotella", "g__Metaprevotella", "g__Palleniella", "g__Paraprevotella", "g__Prevotella", "g__Prevotellamassilia", "g__Pseudoprevotella", "g__Segatella", "g__Xylanibacter"))

taxon_wo_big <- tax_table_big %>% 
  filter(!genus %in% c("g__Alloprevotella", "g__Hallella", "g__Hoylesella", "g__Ihuprevotella", "g__Leyella", "g__Marseilla", "g__Massiliprevotella", "g__Metaprevotella", "g__Palleniella", "g__Paraprevotella", "g__Prevotella", "g__Prevotellamassilia", "g__Pseudoprevotella", "g__Segatella", "g__Xylanibacter"))

taxon <- rbind(taxon_wo, test)

write_csv(taxon, "data/taxon_table.csv")
  
taxon <- rbind(taxon_wo_big, test_big)

write_csv(taxon, "data/taxon_table_big.csv")