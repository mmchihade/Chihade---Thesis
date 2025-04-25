#load packages ----
library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(ggpubr)
library(ggplot2)
library(janitor)
library(wesanderson)
library(scales)
library(data.table)
library(sjmisc)
library(purrr)
library(MicrobiomeStat)
library(microViz)
library(ggbreak)
library(ggh4x)


#load data ----
data_otu <- read.csv("data/OTU_tables/kraken2_READS_sg_2024_decontam_G.csv", header = TRUE)
data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)
data_taxo <- read.csv("data/taxon_table.csv", header=TRUE)

data_grp <- data_grp %>%
  column_to_rownames("lab_id") 

data_taxo <- data_taxo %>% 
  select(-species) %>%
  distinct() %>% 
  filter(!genus %in% c("g__")) %>% 
  column_to_rownames("genus") %>% 
  rotate_df() %>%
  clean_names() %>% 
  rotate_df()

data_otu_clean <- clean_names(data_otu) %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

in_otu_core <- data_otu %>% 
  mutate_all(~ifelse(is.na(.), "0", .)) %>% 
  column_to_rownames("lab_id") %>%
  mutate_at(vars(g_prevotella:g_methanobacterium), as.numeric) %>% 
  rotate_df()

sum <- in_otu_core %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  arrange(-sum)

#we can also put this data in a phyloseq ----
OTU <- otu_table(as.matrix(in_otu_core), taxa_are_rows = TRUE)
SAM <- sample_data(data_grp) 
TAX <- tax_table(as.matrix(data_taxo))
data_phylo <- phyloseq(OTU, SAM, TAX)

#I want to find the relative abundance: 

ps_rel_abund <- phyloseq::transform_sample_counts(data_phylo, function(x){x / sum(x)})

#for my plots here, I want to highlight certain taxons. To do this, I need to make new phyloseqs with just certain taxon

#first let's do top 10 genera:

top_ten = filter_taxa(data_phylo, function(x) sum(x) > 148000000, TRUE)

ps_topten_rel_abund <- phyloseq::transform_sample_counts(top_ten, function(x){x / sum(x)}) 


#how about some differently abundant taxa 

top = (prune_taxa(c("g_butyricicoccus_a",
                    "g_treponema_d",
                    "g_eubacterium_j",
                    "g_umgs1484",
                    "g_uba2882",
                    "g_cag_45",
                    "g_lachnospira",
                    "g_uba3282",
                    "g_uba4644",
                    "g_cag_110",
                    "g_cag_238",
                    "g_oribacterium",
                    "g_slackia_a",
                    "g_uba1740",
                    "g_uba738"), data_phylo))

top_rel_abund <- phyloseq::transform_sample_counts(top, function(x){x / sum(x)}) 

top_big =  (prune_taxa(c("g_treponema_d",
                         "g_uba4372",
                         "g_lachnospira",
                         "g_sodaliphilus",
                         "g_succinivibrio",
                         "g_klebsiella",
                         "g_v9d3004",
                         "g_butyricicoccus_a",
                         "g_odoribacter",
                         "g_rc9",
                         "g_sutterella",
                         "g_alistipes",
                         "g_bacteroides",
                         "g_cag_791",
                         "g_oribacterium"), data_phylo))

top_big_rel_abund <- phyloseq::transform_sample_counts(top_big, function(x){x / sum(x)}) 

top_both <- (prune_taxa(c("g_alistipes",
                          "g_bacteroides",
                          "g_butyricicoccus_a",
                          "g_cag_110",
                          "g_cag_238",
                          "g_cag_45",
                          "g_cag_791",
                          "g_eubacterium_j",
                          "g_klebsiella",
                          "g_lachnospira",
                          "g_odoribacter",
                          "g_oribacterium",
                          "g_rc9",
                          "g_slackia_a",
                          "g_sodaliphilus",
                          "g_succinivibrio",
                          "g_sutterella",
                          "g_treponema_d",
                          "g_uba1740",
                          "g_uba2882",
                          "g_uba3282",
                          "g_uba4372",
                          "g_uba4644",
                          "g_uba738",
                          "g_umgs1484",
                          "g_v9d3004"), data_phylo))

top_both_rel_abund <- phyloseq::transform_sample_counts(top_both, function(x){x / sum(x)}) 

#let's plot! ----

#by top ten 

#relative abundance

phyloseq::plot_bar(ps_topten_rel_abund, fill = "OTU") +
  geom_bar(aes(fill = OTU), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_grid(~ country, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(.1,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",
                               g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481", 
                               g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",
                               g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),
  )

ggsave("output/rel_abun/phylo_rel_abun_top10_country.pdf", width = 40, height = 20, units = "cm")


#melt 

phyloseq::psmelt(top_ten) %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(0.1,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))
 
ggsave("output/rel_abun/phylo_melt_topten_country.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top_ten) %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))
 

ggsave("output/rel_abun/phylo_melt_topten_studygroup.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top_ten) %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_cag_83="#B1EEAF", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_succinivibrio="#385640", g_acetatifactor="#B30063"))  +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/phylo_melt_topten_ethnicity.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top_ten) %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))  +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/phylo_melt_topten_countrygroup.pdf", width = 40, height = 40, units = "cm")


#by top five genera based on all my DA values

#relative abundance

phyloseq::plot_bar(top_rel_abund, fill = "OTU") +
  geom_bar(aes(fill = OTU), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ country, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))

ggsave("output/rel_abun/top_rel_abun_country.pdf", width = 40, height = 40, units = "cm")


#melt 

phyloseq::psmelt(top) %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))
  
ggsave("output/rel_abun/phylo_melt_top_country.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top) %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))
 
ggsave("output/rel_abun/phylo_melt_top_studygroup.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top) %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" )) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_ethnicity.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top) %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" )) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/phylo_melt_top_countrygroup.pdf", width = 40, height = 40, units = "cm")


#by top five genera based on my DA values with the largest fold changes

#relative abundance

phyloseq::plot_bar(top_big_rel_abund, fill = "OTU") +
  geom_bar(aes(fill = OTU), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ country, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))

ggsave("output/rel_abun/top_big_rel_abun_country.pdf", width = 40, height = 40, units = "cm")


#melt 

phyloseq::psmelt(top_big) %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) +  
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))      
 

ggsave("output/rel_abun/phylo_melt_top_big_country.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top_big) %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))
  
ggsave("output/rel_abun/phylo_melt_top_big_studygroup.pdf", width = 40, height = 40, units = "cm")

phyloseq::psmelt(top_big) %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) + 
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" )) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_big_ethnicity.pdf", width = 60, height = 40, units = "cm")

phyloseq::psmelt(top_big) %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_big_countrygroup.pdf", width = 60, height = 40, units = "cm")

#By the two combined: 

#relative abundance

phyloseq::plot_bar(top_both_rel_abund, fill = "OTU") +
  geom_bar(aes(fill = OTU), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ country, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))

ggsave("output/rel_abun/top_both_rel_abun_country.pdf", width = 40, height = 40, units = "cm")


#melt 

phyloseq::psmelt(top_both) %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=4,label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) +  
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_cag_791 = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))      


ggsave("output/rel_abun/phylo_melt_top_both_country.pdf", width = 50, height = 60, units = "cm")

phyloseq::psmelt(top_both) %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=4,label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_cag_791 = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_both_studygroup.pdf", width = 50, height = 60, units = "cm")

phyloseq::psmelt(top_both) %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=4,label.x.npc = 0.2, label.y.npc = 0.9) + 
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" )) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_both_ethnicity.pdf", width = 60, height = 40, units = "cm")

phyloseq::psmelt(top_both) %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +   stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=4,label.x.npc = 0.2, label.y.npc = 0.9) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/phylo_melt_top_both_countrygroup.pdf", width = 60, height = 40, units = "cm")


#wait I am sorry I want to see my means more clearly so I am gonna do all these plots without the super high values ----

topten_plots <- phyloseq::psmelt(top_ten) %>% 
  group_by(OTU) %>%
  group_split()

P1 <-topten_plots[[3]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n",
       tag = "A") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none")

P2 <-topten_plots[[3]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "C") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none")

P3 <-topten_plots[[3]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "B") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), angle=90, hjust=-0.1, vjust=0.5, size = 4, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

P4 <- topten_plots[[3]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "D") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), fontface = "bold", angle=90, hjust=-0.1, vjust=0.5, size = 4) +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

(P1|P3)/(P2|P4)

ggsave("output/rel_abun/topten_bac_supertrim.pdf", width = 40, height = 20, units = "cm")

P1 <-topten_plots[[7]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "D") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,20000000)) + theme(legend.position="none")

P2 <-topten_plots[[7]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "C") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,20000000)) + theme(legend.position="none")

P3 <-topten_plots[[7]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "B") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), angle=90, hjust=-0.1, vjust=0.5, size = 4, fontface = "bold") +  
  coord_cartesian(ylim = c(0,20000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

P4 <-topten_plots[[7]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n", 
       tag = "D") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), fontface = "bold", angle=90, hjust=-0.1, vjust=0.5, size = 4) +  
  coord_cartesian(ylim = c(0,20000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

(P1|P3)/(P2|P4)

ggsave("output/rel_abun/topten_prev_supertrim.pdf", width = 40, height = 20, units = "cm")

#by top five genera based on all my DA values

#melt 

top_plots <- phyloseq::psmelt(top) %>% 
  group_by(OTU) %>%
  group_split()

P1 <- top_plots[[9]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n", 
       tag = "A") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none")

P2 <- top_plots[[9]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n", 
       tag = "C") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,800000)) + theme(legend.position="none")

P3 <- top_plots[[9]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n", 
       tag = "B") + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), fontface = "bold", angle=90, hjust=-0.1, vjust=0.5, size = 4) +   
  coord_cartesian(ylim = c(0,1000000)) + theme(legend.position="none")  + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

P4 <- top_plots[[9]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n", 
       tag = "D") +
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), fontface = "bold", angle=90, hjust=-0.1, vjust=0.5, size = 4) +  
  coord_cartesian(ylim = c(0,800000)) + theme(legend.position="none")  + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

(P1|P3)/(P2|P4)

ggsave("output/rel_abun/top_trep_supertrim.pdf", width = 40, height = 20, units = "cm")



top_plots[[1]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_bac_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[2]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_buty_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[3]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_cag238_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[4]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,13000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_cag45_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[5]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_eubf_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[6]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,25000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_eubj_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[7]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,400000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_lach_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[8]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,15000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_slack_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[9]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_trep_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[10]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,25000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba1740_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[11]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba2821_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[12]] %>%
  ggplot(aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5,  label.x.npc = 0.2, label.y.npc = 0.9, label.x = 1, label.y = 40000) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 50000)) + theme(legend.position="none") 

ggsave("output/rel_abun/top_uba2882_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[13]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,30000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba3282_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[14]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba4644_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[15]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_umgs1484_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_plots[[1]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_bac_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[2]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,30000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_buty_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[3]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_cag238_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[4]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,13000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_cag45_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[5]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,70000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_eubf_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[6]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,25000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_eubj_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[7]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,400000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_lach_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[8]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,10000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_slack_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[9]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,800000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_trep_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[10]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,25000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba1740_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[11]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,4000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba2821_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[12]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 30000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba2882_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[13]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,10000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba3282_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[14]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_uba4644_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[15]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,7000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_umgs1484_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[1]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_bac_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[2]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_buty_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[3]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_cag238_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[4]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,13000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_cag45_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[5]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_eubf_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[6]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_eubj_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[7]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,400000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_lach_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[8]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,15000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_slack_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[9]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,1100000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 


ggsave("output/rel_abun/top_trep_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[10]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba1740_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[11]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba2821_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[12]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 60000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba2882_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[13]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,30000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba3282_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[14]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba4644_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[15]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,8000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_umgs1484_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_plots[[1]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,5000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_bac_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[2]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_buty_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[3]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_cag238_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[4]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,13000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_cag45_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[5]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_eubf_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[6]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_eubj_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[7]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,400000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_lach_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[8]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,15000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_slack_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[9]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,1000000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_trep_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[10]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,25000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_uba1740_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[11]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,4000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba2821_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[12]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 60000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba2882_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[13]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,30000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggsave("output/rel_abun/top_uba3282_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[14]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_uba4644_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_plots[[15]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,7000)) + theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_umgs1484_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")


#by top five genera based on my DA values with the largest fold changes

#melt 

top_big_plots <- phyloseq::psmelt(top_big) %>% 
  group_by(OTU) %>%
  group_split()

top_big_plots[[1]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 500000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_ali_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[2]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,2000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_bac_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[3]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,20000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_buty_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[4]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_cag791_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[5]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,23000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_kleb_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[6]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_lach_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[7]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_ord_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[8]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_rc9_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[9]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,75000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_rug115_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[10]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_sod_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[11]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,750000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_suc_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[12]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,300000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_sutt_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[13]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_trep_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[14]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_uba2282_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[15]] %>%
  ggplot(data = ., aes(x = country, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,300000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_uba4372_supertrim_country.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[1]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 900000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_ali_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[2]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_bac_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[3]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,20000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_buty_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[4]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,50000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_cag791_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[5]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,130000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_kleb_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[6]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_lach_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[7]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_ord_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[8]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_rc9_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[9]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,75000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_rug115_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[10]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_sod_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[11]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,750000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_suc_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[12]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_sutt_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[13]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,1000000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_trep_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[14]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_uba2282_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[15]] %>%
  ggplot(data = ., aes(x = study_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,30000)) + theme(legend.position="none")

ggsave("output/rel_abun/top_big_uba4372_supertrim_study_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[1]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 1500000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_ali_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[2]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_bac_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[3]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_buty_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[4]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,90000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_cag791_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[5]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,150000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_kleb_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[6]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,300000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_lach_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[7]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_ord_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[8]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_rc9_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[9]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,75000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_rug115_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[10]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_sod_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[11]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,900000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_suc_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[12]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_sutt_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[13]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,1000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_trep_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[14]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,80000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_uba2282_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[15]] %>%
  ggplot(data = ., aes(x = ethnicity, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,700000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_uba4372_supertrim_ethnicity.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[1]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0, 2000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_ali_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[2]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_bac_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[3]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,40000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_buty_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[4]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,100000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_cag791_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[5]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,130000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_kleb_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[6]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_lach_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[7]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,200000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_ord_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[8]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_fill_manual(values = c(all_low = "#8BD0F0", g_prevotella = "#D26E7A", g_escherichia = "#DDCC75", g_bacteroides = "#0D763B", g_wg_1 = "#423B94", g_treponema_d = "#AA459B", g_phocaeicola = "#40AA99", g_agathobacter = "#9E9D40", g_akkermansia = "#993563", g_alistipes = "#843835", g_bifidobacterium = "#6599CC", g_blautia_a="#868486", g_butyrivibrio_a="#FC1626", g_cag_103="#16FC22", g_cag_110="#CB00FF", g_cag_127="#FD009B", g_cag_170="#264DFC", g_cag_177="#FC7E16", g_cag_180="#FB0DDD", g_cag_279="#22FEAF", g_cag_303="#BFF01C", g_cag_349="#F6C8FD", g_cag_353="#CD8AFD", g_cag_475="#D42E49", g_cag_510="#F9B000", g_cag_568="#1C9BFD", g_cag_632="#FBC9BF", g_cag_83="#B1EEAF", g_cag_95="#2EAD0D", g_coe1="#FD68C0", g_coprococcus="#00F9FC", g_dialister="#A9731C", g_eisenbergiella="#68406D", g_enterocloster="#FDE416", g_er4="#B600BF", g_f082="#732ACD",                                g_faecalibacterium="#3DC273", g_gemmiger="#CDE8E0", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_megamonas= "#F79BCA", g_mitsuokella="#BA0D62", g_parabacteroides="#FF7EEE", g_paraprevotella="#FCA481",                                 g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_rf16="#F64F22", g_roseburia="#A8717F", g_rug115="#C8C5A0", g_rug410="#9BCB56", g_ruminiclostridium_e="#FC0D6A", g_ruminococcus_d="#9100FF", g_ruminococcus_e="#0DC7FC", g_sfdb01="#0D51C5", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba9732="#FD7A78", g_uba10281="#F435FE", g_uba11524="#ADB2D1", g_uba1777="#0DB97C", g_uba4248="#C600A5", g_uba4372="#852EA9", g_uba7173="#F3DBE9",                                g_vsob01="#C4CD1C", g_acetatifactor="#B30063"),   )+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,3000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_rc9_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[9]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,75000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_rug115_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[10]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,500000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_sod_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[11]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,750000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_suc_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[12]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,300000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_sutt_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[13]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,1000000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_trep_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[14]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,60000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_uba2282_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")

top_big_plots[[15]] %>%
  ggplot(data = ., aes(x = country_group, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2, show.legend = FALSE) +  
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +   theme(strip.placement = "outside",         strip.background = element_rect(fill = NA, color = "white"),         panel.spacing = unit(-.01,"cm")) +    theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  scale_color_manual(values = c(g_prevotella = "#D26E7A", g_bacteroides = "#0D763B", g_treponema_d = "#AA459B", g_butyricicoccus_a = "#40AA99", g_agathobacter = "#9E9D40", g_odoribacter = "#993563", g_alistipes = "#843835", g_sutterella = "#6599CC", g_eubacterium_j ="#868486", g_eubacterium_f ="#FC1626", g_cag_110="#16FC22", g_cag_238="#CB00FF", g_cag_45="#FD009B", g_slackia_a ="#264DFC", g_cag_83="#B1EEAF", g_faecalibacterium="#3DC273", g_klebsiella="#00626E", g_lachnospira="#AF4900", g_phocaeicola="#9CA1FE", g_prevotellamassilia="#53F6D2", g_rc9="#9F73BE", g_oribacterium="#C8C5A0", g_umgs1484 ="#0DC7FC", g_sodaliphilus="#88F958", g_succinivibrio="#385640", g_uba2882="#FD7A78", g_uba1740="#F435FE", g_uba2821="#ADB2D1", g_uba3282="#0DB97C", g_uba4644="#C600A5", g_uba738="#852EA9", g_v9d3004 ="#F3DBE9",  g_vsob01="#C4CD1C", g_acetatifactor="#B30063" ))+
  stat_summary(fun = mean, geom="point",colour="black", size=1) +      
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, aes(label=signif(after_stat(y), digits = 3)), size = 4, vjust = -0.7, fontface = "bold") +  
  coord_cartesian(ylim = c(0,800000)) + theme(legend.position="none") +    theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/rel_abun/top_big_uba4372_supertrim_country_group.pdf", width = 40, height = 20, units = "cm")
