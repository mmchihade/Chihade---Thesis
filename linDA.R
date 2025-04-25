#load packages ----
library(wesanderson)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggcorrplot)
library(ggpmisc)
library(umap)
suppressMessages(library(Hmisc)) #not sure if needed, but there is some ggplot interaction
suppressMessages(library(tidyverse)) #note that installing tidyverse takes a LONG time!
library(dplyr)
suppressMessages(library(factoextra)) #don't know how much of these packages are necessary; this is drawn from lipid_micro_EXTRA's packages
# library(FactoMineR) #do we have this?
library(igraph)
library(stringi)
library(pals)
library(phyloseq)
library(MicrobiomeStat)
library(sjmisc)
library(janitor)
library(phyloseq)
library(devtools)
library(remotes)
library(microbiome)
library(microViz)
library(eatTools)


#load data ---- 

#to keep ourselves organized, let's create a phyloseq 

data_otu <- read.csv("data/OTU_tables/kraken2_READS_sg_2024_decontam_G.csv", header = TRUE)
data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_alpha_2024.csv", header=TRUE)
data_taxo <- read.csv("data/taxon_table.csv", header=TRUE)

data_grp <- data_grp %>%
  tibble::column_to_rownames("lab_id") 

data_taxo <- data_taxo %>% 
  dplyr::select(-species) %>%
  distinct() %>% 
  filter(!genus %in% c("g__")) %>% 
  distinct() %>%
  tibble::column_to_rownames("genus") %>% 
  rotate_df() %>%
  clean_names() %>% 
  rotate_df()


in_otu_core <- data_otu %>% 
  mutate_all(~ifelse(is.na(.), "0", .)) %>% 
  column_to_rownames("lab_id") %>%
  mutate_at(vars(g_prevotella:g_methanobacterium), as.numeric) %>% 
  rotate_df()

sum <- in_otu_core %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  arrange(-sum) %>% 
  select(sum) %>% 
  tibble::rownames_to_column("genera")


# now put this data in a phyloseq 
OTU <- otu_table(as.matrix(in_otu_core), taxa_are_rows = TRUE)
SAM <- sample_data(data_grp) 
TAX <- tax_table(as.matrix(data_taxo))
data_phylo <- phyloseq(OTU, SAM, TAX)

phyloseq::taxa_are_rows(data_phylo) #testing to see if taxa are rows

#okay, now we are supposed to transform the phyloseq into a MicrobiomeStat data set 

data.obj <- mStat_convert_phyloseq_to_data_obj(data_phylo) #yay! 

#let's make sure all of our factors are the type of variable they should be: 

data.obj$meta.dat$sex <- factor(as.factor(data.obj$meta.dat$sex), levels = c("M", "F"))

data.obj$meta.dat$country <- factor(as.factor(data.obj$meta.dat$country), levels = c("Cameroon", "Botswana", "Tanzania"))

data.obj$meta.dat$country_group <- factor(as.factor(data.obj$meta.dat$country_group), levels = c("Cameroon Agropastoralist", "Cameroon Hunter-gatherer","Cameroon Pastoralist", 
                                                                                                 "Botswana Agropastoralist", "Botswana Hunter-gatherer",
                                                                                                 "Botswana Pastoralist", "Tanzania Agropastoralist", 
                                                                                                 "Tanzania Hunter-gatherer", "Tanzania Pastoralist"))

data.obj$meta.dat$ethnicity <- factor(as.factor(data.obj$meta.dat$ethnicity), levels = c("Bantu General","Baka", "Khoesan", "Kgalagadi","Bagyeli","Fang", "Ngoumba", "Ewondo", "Nzime", "Hadzabe", "Burunge", "Herero",
 "Maasai", "Mbororo Fulani", "Tikari"))

data.obj$meta.dat$study_group <- factor(as.factor(data.obj$meta.dat$study_group), levels = c("Agropastoralist", "Hunter-gatherer", "Pastoralist"))

data.obj$meta.dat$region <- factor(as.factor(data.obj$meta.dat$region), levels = c("CM-East", "CM-Southwest", "CM-Northwest", "BW-Central West",
"BW-West", "BW-Northwest", "TZ-North Central", "TZ-Central", "TZ-North"))

data.obj$meta.dat$waist_circumference <- as.numeric(data.obj$meta.dat$waist_circumference)

data.obj$meta.dat$bmi <- as.numeric(data.obj$meta.dat$bmi)

data.obj$meta.dat$S.obs <- as.numeric(data.obj$meta.dat$S.obs)

data.obj$meta.dat$Shannon <- as.numeric(data.obj$meta.dat$Shannon)

#now let's try to use LinDA to find some DA genera and families ----

test.list.country <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "country", 
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.country_group <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "country_group",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.ethnicity <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "ethnicity",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.group <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "study_group",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.region <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "region",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

#now let's try to use LinDA to find some DA genera and families ----

test.list.country <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "country", 
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.country_group <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "country_group",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.ethnicity <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "ethnicity",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.group <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "study_group",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

test.list.region <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "region",
  adj.vars = c("sex", "age", "bmi", "S.obs"),
  feature.dat.type = "count",
  feature.level = c("original", "family"),
  prev.filter = 0.1,
  abund.filter = 0.0001)

#now let's generate a volcano plots ----

volcano_plots.country <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "country",
  test.list = test.list.country,
  feature.sig.level = 0.05,
  feature.mt.method = "none")

volcano_plots.country[["original"]][["Botswana vs Cameroon (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))    + 
  labs(y="-log(p-value)", x="Fold Change") + 
  scale_size_binned(
    range = c(1, 9),
    n.breaks = 10) + 
  theme(legend.key.width = unit(0.8, "cm"))
  

ggsave("output/LinDA/g_botswana.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country[["original"]][["Tanzania vs Cameroon (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))    

ggsave("output/LinDA/g_tanzania.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country[["family"]][["Botswana vs Cameroon (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_botswana.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country[["family"]][["Tanzania vs Cameroon (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tanzania.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "country_group",
  test.list = test.list.country_group,
  feature.sig.level = 0.05,
  feature.mt.method = "none")


volcano_plots.country_group[["original"]][["Cameroon Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_chg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Cameroon Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_cp.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Botswana Agropastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bag.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Botswana Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bhg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Botswana Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bp.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Tanzania Agropastoralist vs Cameroon Agropastoralist (Reference)"]] + labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm")) 

ggsave("output/LinDA/g_tag.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Tanzania Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    
  labs(y="-log(p-value)", x="Fold Change") +    
  scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    
  theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_thg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["original"]][["Tanzania Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_tp.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Cameroon Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_chg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Cameroon Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_cp.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Botswana Agropastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bag.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Botswana Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bhg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Botswana Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bp.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Tanzania Agropastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tag.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Tanzania Hunter-gatherer vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_thg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.country_group[["family"]][["Tanzania Pastoralist vs Cameroon Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tp.pdf", width = 20, height = 20, units = "cm")


volcano_plots.ethnicity <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "ethnicity",
  test.list = test.list.ethnicity,
  feature.sig.level = 0.05,
  feature.mt.method = "none")

volcano_plots.ethnicity[["original"]][["Baka vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_baka.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Khoesan vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_khoesan.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Kgalagadi vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_kgalagadi.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Bagyeli vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bagyeli.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Fang vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_fang.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Ngoumba vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_ngoumba.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Ewondo vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_ewondo.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Nzime vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_nzime.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Hadzabe vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_hadza.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Burunge vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_burunge.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Herero vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_herero.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Maasai vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_maasai.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Mbororo Fulani vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_fulani.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["original"]][["Tikari vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_tikari.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Baka vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_baka.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Khoesan vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_khoesan.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Kgalagadi vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_kgalagadi.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Bagyeli vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bagyeli.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Fang vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_fang.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Ngoumba vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_ngoumba.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Ewondo vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_ewondo.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Nzime vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_nzime.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Hadzabe vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_hadza.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Burunge vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_burunge.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Herero vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_herero.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Maasai vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_maasai.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Mbororo Fulani vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_fulani.pdf", width = 20, height = 20, units = "cm")

volcano_plots.ethnicity[["family"]][["Tikari vs Bantu General (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tikari.pdf", width = 20, height = 20, units = "cm")

volcano_plots.group <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "study_group",
  test.list = test.list.group,
  feature.sig.level = 0.05,
  feature.mt.method = "none")

volcano_plots.group[["original"]][["Hunter-gatherer vs Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_hg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.group[["original"]][["Pastoralist vs Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_p.pdf", width = 20, height = 20, units = "cm")

volcano_plots.group[["family"]][["Hunter-gatherer vs Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_hg.pdf", width = 20, height = 20, units = "cm")

volcano_plots.group[["family"]][["Pastoralist vs Agropastoralist (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_p.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "region",
  test.list = test.list.region,
  feature.sig.level = 0.05,
  feature.mt.method = "none")

volcano_plots.region[["original"]][["CM-Southwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_cme.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["CM-Northwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_cmn.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["BW-Central West vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bcw.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["BW-West vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bww.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["BW-Northwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_bwn.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["TZ-North Central vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_tznc.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["TZ-Central vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/g_tzc.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["original"]][["TZ-North vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   
 
ggsave("output/LinDA/g_tzn.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["CM-Southwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_cme.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["CM-Northwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_cmn.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["BW-Central West vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bcw.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["BW-West vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bww.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["BW-Northwest vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_bwn.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["TZ-North Central vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tznc.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["TZ-Central vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tzc.pdf", width = 20, height = 20, units = "cm")

volcano_plots.region[["family"]][["TZ-North vs CM-East (Reference)"]] +    labs(y="-log(p-value)", x="Fold Change") +    scale_size_binned(     range = c(1, 9),     n.breaks = 10) +    theme(legend.key.width = unit(0.8, "cm"))   

ggsave("output/LinDA/f_tzn.pdf", width = 20, height = 20, units = "cm")

#finally, I want to be able to make a list of all of the significant DA results ----

country_all <- unlist(test.list.country, recursive = FALSE) %>% 
  map2_df(., names(.), ~ mutate(.x, ID = .y)) %>%
  separate(ID, into = c("taxon", "country", "vs", "reference")) %>% 
  select(-vs) %>% 
  mutate(taxon = recode(taxon, "original" = "genus")) %>% 
  filter(Adjusted.P.Value <= 0.05)
  

country_group_all <- unlist(test.list.country_group, recursive = FALSE) %>% 
  map2_df(., names(.), ~ mutate(.x, ID = .y)) %>%
  separate(ID, into = c("taxon", "country", "group", "group2", "vs", "reference", "reference2"))

hg <- country_group_all %>% 
  filter(vs %in% c("vs")) %>% 
  unite(group, c("group", "group2")) %>% 
  unite(reference, c("reference", "reference2")) %>% 
  select(-vs)

ag_p <- country_group_all %>% 
  filter(group2 %in% c("vs")) %>%
  unite(reference, c("vs", "reference")) %>% 
  select(-group2, -reference2)

country_group_all <- rbind(hg, ag_p) %>%
  mutate(taxon = recode(taxon, "original" = "genus")) %>% 
  filter(Adjusted.P.Value <= 0.05)


ethnicity_all <- unlist(test.list.ethnicity, recursive = FALSE) %>% 
  map2_df(., names(.), ~ mutate(.x, ID = .y)) %>%
  separate(ID, into = c("taxon", "ethnicity", "vs", "reference")) %>% 
  select(-vs) %>% 
  mutate(taxon = recode(taxon, "original" = "genus"))  %>% 
  filter(Adjusted.P.Value <= 0.05)
  

group_all <- unlist(test.list.group, recursive = FALSE) %>% 
  map2_df(., names(.), ~ mutate(.x, ID = .y)) %>%
  separate(ID, into = c("taxon", "group", "group2", "vs", "reference")) 


hg2 <- group_all %>% 
  filter(vs %in% c("vs")) %>% 
  unite(group, c("group", "group2")) %>%
  select(-vs)

ag_p2 <- group_all %>% 
  filter(group2 %in% c("vs")) %>%
  mutate(reference = vs) %>% 
  select(-group2, -vs)

group_all <- rbind(hg2, ag_p2) %>%
  mutate(taxon = recode(taxon, "original" = "genus"))  %>% 
  filter(Adjusted.P.Value <= 0.05)

region_all <- unlist(test.list.region, recursive = FALSE) %>% 
  map2_df(., names(.), ~ mutate(.x, ID = .y)) %>%
  separate(ID, into = c("taxon", "region", "region2", "vs", "reference", "reference2", "reference3"))

region1 <- region_all %>% 
  filter(vs %in% c("vs")) %>% 
  unite(region, c("region", "region2")) %>% 
  unite(reference, c("reference", "reference2")) %>%
  select(-vs, -reference3) 

region2 <- region_all %>% 
  filter(reference %in% c("vs")) %>%
  unite(region, c("region", "region2", "vs")) %>% 
  unite(reference, c("reference2", "reference3")) 

region_all <- rbind(region1, region2) %>%
  mutate(taxon = recode(taxon, "original" = "genus"))  %>% 
  filter(Adjusted.P.Value <= 0.05)


#okay, I also want to be able to make a list of the top 20 most interesting DA results ----

country_20 <- country_all %>%
  group_by(country, taxon) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=20)

write_csv(country_20, "output/DAranks/country.csv")

country_group_20 <- country_group_all %>%
  group_by(country, group, taxon) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=20)

write_csv(country_group_20, "output/DAranks/country_group.csv")

ethnicity_20 <- ethnicity_all %>% 
  group_by(ethnicity, taxon) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=20)

write_csv(ethnicity_20, "output/DAranks/ethnicity.csv")

group_20 <- group_all %>% 
  group_by(group, taxon) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=20)

write_csv(group_20, "output/DAranks/group.csv")

region_20 <- region_all %>%
  group_by(region, taxon) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=20)

write_csv(region_20, "output/DAranks/region.csv")

#to make sure that I can bind everything together, I need to bind together the comparison columns

country_for_all <- country_20 %>% 
  unite(comparison, c("taxon", "country", "reference"))

country_group_for_all <- country_group_20 %>% 
  unite(comparison, c("taxon", "country", "group", "reference"))

group_for_all <- group_20 %>% 
  unite(comparison, c("taxon", "group", "reference"))

ethnicity_for_all <- ethnicity_20 %>% 
  unite(comparison, c("taxon", "ethnicity", "reference"))

region_for_all <- region_20 %>% 
  unite(comparison, c("taxon", "region", "reference"))

region_everything <- region_all %>% 
  unite(comparison, c("taxon", "region", "reference"))

country_everything <- country_all %>% 
  unite(comparison, c("taxon", "country", "reference"))

country_group_everything <- country_group_all %>% 
  unite(comparison, c("taxon", "country", "group", "reference"))

group_everything <- group_all %>% 
  unite(comparison, c("taxon", "group", "reference"))

ethnicity_everything <- ethnicity_all %>% 
  unite(comparison, c("taxon", "ethnicity", "reference"))

region_everything <- region_all %>% 
  unite(comparison, c("taxon", "region", "reference"))

#now I want to make two data frames, one with all the significant DAs ordered by which ones present across the most factors, and then one that is that but just the ones with the largest fold changes  

every_DA <- rbind(country_everything, country_group_everything, group_everything, ethnicity_everything, region_everything) %>% 
  ungroup() %>%
  unite(foldChange_Pvalue, c("Coefficient", "Adjusted.P.Value")) %>% 
  select(Variable, foldChange_Pvalue, comparison) %>% 
  pivot_wider(names_from = comparison,
              values_from = foldChange_Pvalue) %>%
  rotate_df() %>% 
  row_to_names(1)

sum_NA_2 <- data.frame(colSums(is.na(every_DA))) %>% 
  rename(na_sum = colSums.is.na.every_DA..) %>% 
  rotate_df()

every_DA_ordered <- rbind(every_DA, sum_NA_2) %>% 
  rotate_df() %>% 
  arrange(na_sum) %>% 
  rownames_to_column("taxon")

write_csv(every_DA_ordered, "output/DAranks/everyone.csv")

all_DAs <- rbind(country_for_all, country_group_for_all, group_for_all, ethnicity_for_all, region_for_all) %>% 
  ungroup() %>%
  unite(foldChange_Pvalue, c("Coefficient", "Adjusted.P.Value")) %>% 
  select(Variable, foldChange_Pvalue, comparison) %>%
  distinct() %>% 
  pivot_wider(names_from = comparison,
              values_from = foldChange_Pvalue) %>% 
  rotate_df() %>% 
  row_to_names(1)
  

sum_NA <- data.frame(colSums(is.na(all_DAs))) %>% 
  rename(na_sum = colSums.is.na.all_DAs..) %>% 
  rotate_df()

all_DAs_ordered <- rbind(all_DAs, sum_NA) %>% 
  rotate_df() %>% 
  arrange(na_sum) %>% 
  rownames_to_column("taxon")

write_csv(all_DAs_ordered, "output/DAranks/all.csv")

#okay, finally, let's make these prettier 

#for these I want the number of reads for each taxon, so let's make something I can add to these for that: 

sum <- in_otu_core %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  arrange(-sum) %>% 
  select(sum) %>% 
  tibble::rownames_to_column("taxon")

data_phylo_family <- tax_agg(ps = data_phylo, "family")

sum_fam <- otu_table(data_phylo_family) %>%
  as.data.frame() %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  arrange(-sum) %>% 
  select(sum) %>% 
  tibble::rownames_to_column("taxon")

sum_pretty <- rbind(sum, sum_fam) %>% 
  rotate_df() %>% 
  row_to_names(1)

all_pretty <- all_DAs_ordered %>% 
  rotate_df() %>% 
  row_to_names(1) %>% 
  rownames_to_column("comparison") %>%
  pivot_longer(f__Treponemataceae:g_umgs1795, names_to = "taxon", values_to = "foldChange_Pvalue") %>% 
  na.omit() %>% 
  unite(all_info, c("comparison", "foldChange_Pvalue")) %>% 
  group_by(taxon) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = all_info) %>%
  select(-row)

all_pretty <- rbind_common(all_pretty, sum_pretty) 
  
write_csv(all_pretty, "output/DAranks/all_pretty.csv") 

everyone_pretty <- every_DA_ordered %>% 
  rotate_df() %>% 
  row_to_names(1) %>% 
  rownames_to_column("comparison") %>%
  pivot_longer(g_rf16:f__Porphyromonadaceae, names_to = "taxon", values_to = "foldChange_Pvalue") %>% 
  na.omit() %>% 
  unite(all_info, c("comparison", "foldChange_Pvalue")) %>% 
  group_by(taxon) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = all_info) %>%
  select(-row)

everyone_pretty <- rbind_common(everyone_pretty, sum_pretty) 

write_csv(everyone_pretty, "output/DAranks/everyone_pretty.csv")

#okay, this is great, but I think I want my genera of interest to be those with the highest p-values, not the most significant values. 

#for this I want to be able to filter by the sum, number of zeros, and number of ind with low read counts for each genus

sum <- in_otu_core %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  arrange(-sum) %>% 
  select(sum) %>% 
  tibble::rownames_to_column("taxon")

low_zero <- read.csv("data/OTU_tables/zerocount_kraken2_READS_sg_2024_decontam_G.csv", header=TRUE) %>% 
  clean_names() %>%
  rename(taxon = genus) %>%
  select(taxon, sum_na, low_value_count, low_value_count)


genus_info <- merge(sum, low_zero, by = "taxon")

#now let's untie all our data frames: 

everyone_genus <- rbind(country_everything, country_group_everything, group_everything, ethnicity_everything, region_everything) %>% 
  rename("taxon" = "Variable")
  
everyone_genus <- merge(everyone_genus, genus_info, by = "taxon") %>% 
  select(-contains('family')) %>% 
  filter(sum_na <= 1 & low_value_count < 15) %>% 
  group_by(comparison) %>% 
  arrange(desc(abs(Coefficient)), Adjusted.P.Value, .by_group = TRUE) %>% 
  slice_head(n=3) %>% 
  unite(all_info, c("comparison", "Coefficient", "Adjusted.P.Value")) %>% 
  select(taxon, all_info) %>% 
  group_by(taxon) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = all_info) %>%
  select(-row)

write_csv(everyone_genus, "output/DAranks/everyone_pretty_genus_filtered.csv")


everyone_pretty_genus <- merge(every_DA_ordered, genus_info, by = "taxon")%>% 
  select(-contains('family')) %>% 
  arrange(na_sum) %>% 
  filter(sum_na <= 1 & low_value_count < 15) %>%
  rotate_df() %>% 
  row_to_names(1) %>% 
  rownames_to_column("comparison") %>%
  pivot_longer(g_butyricicoccus_a:g_cag_873, names_to = "taxon", values_to = "foldChange_Pvalue") %>% 
  unite(all_info, c("comparison", "foldChange_Pvalue")) %>% 
  group_by(taxon) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = all_info) %>%
  select(-row)

write_csv(everyone_pretty_genus, "output/DAranks/everyone_pretty_genus_filtered.csv")


all_pretty_genus <- merge(all_DAs_ordered, genus_info, by = "taxon")%>% 
  select(-contains('family')) %>% 
  arrange(na_sum) %>% 
  filter(sum_na <= 1 & low_value_count < 15) %>%
  rotate_df() %>% 
  row_to_names(1) %>% 
  rownames_to_column("comparison") %>%
  pivot_longer(g_treponema_d:g_zag111, names_to = "taxon", values_to = "foldChange_Pvalue") %>% 
  unite(all_info, c("comparison", "foldChange_Pvalue")) %>% 
  group_by(taxon) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = taxon, values_from = all_info) %>%
  select(-row)

write_csv(all_pretty_genus, "output/DAranks/all_pretty_genus_filtered.csv")
