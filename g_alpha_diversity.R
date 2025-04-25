#load packages ----
library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(janitor)
library(wesanderson)
library(scales)
library(data.table)



#load data ----
data_otu <- read.csv("data/OTU_tables/kraken2_READS_sg_2024_decontam_G.csv", header = TRUE)
data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)


data_grp <- data_grp %>%
  column_to_rownames("lab_id")


in_otu_core <- data_otu %>%
  mutate_all(~ifelse(is.na(.), "0", .)) %>% 
  column_to_rownames("sample") %>%
  mutate_at(vars(g_prevotella:g_methanobacterium), as.numeric)


#we can also put this data in a phyloseq. I am not sure how this is useful yet but i'll learn! ----
OTU <- otu_table(as.matrix(in_otu_core), taxa_are_rows = FALSE)
SAM <- sample_data(data_grp)                
data_phylo <- phyloseq(OTU, SAM)


#Okay, now let's actually find the alpha diversity! ----

data_richness <- estimateR(in_otu_core)                                            # calculate richness and Chao1 using vegan package

total_alpha <- estimate_richness(OTU, measures=c("Shannon", "Simpson", "InvSimpson", "Fisher"))

data_alphadiv <- cbind(data_grp, t(data_richness), total_alpha)   # combine all indices in one data table

#and make it pretty 

data_alphadiv_tidy <- data_alphadiv %>%
  select(-S.chao1, -se.chao1) %>% 
  pivot_longer(names_to = "alphadiv_index",
               values_to = "obs_values",
               S.obs:Fisher)

#now let's plot! ----

site_color <- (c(wes_palette("Darjeeling1", 5, type = c("discrete")), wes_palette("Darjeeling2", 4, type = c("discrete")), wes_palette("Royal1", 4, type = c("discrete")), wes_palette("Royal2", 5, type = c("discrete")), wes_palette("Rushmore1", 5, type = c("discrete")), wes_palette("BottleRocket1", 6, type = c("discrete")), wes_palette("BottleRocket2", 5, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete")), wes_palette("Moonrise3", 5, type = c("discrete"))))

show_col(site_color)

#By one thing ----

P1 <- ggplot(data_alphadiv, aes(x=sample_site, y=S.obs, fill = sample_site)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
            x = "site",
            y = "richness", 
            tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 750) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P2 <- ggplot(data_alphadiv, aes(x=sample_site, y=Shannon, fill = sample_site )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'site', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



P3 <- ggplot(data_alphadiv, aes(x=sample_site, y=Fisher, fill = sample_site)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'site', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



P4 <- ggplot(data_alphadiv, aes(x=sample_site, y=InvSimpson, fill = sample_site )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'site', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 35) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


 (P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_site.pdf", width = 40, height = 20, units = "cm")


P1 <- ggplot(data_alphadiv, aes(x=country, y=S.obs, fill = country)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "country",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 750) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=country, y=Shannon, fill = country )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'country', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=country, y=Fisher, fill = country)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'country', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=country, y=InvSimpson, fill = country )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'country', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 35) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_country.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=sex, y=S.obs, fill = sex)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "sex",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 750) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

P2 <- ggplot(data_alphadiv, aes(x=sex, y=Shannon, fill = sex )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +   
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

P3 <- ggplot(data_alphadiv, aes(x=sex, y=Fisher, fill = sex)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'sex', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +   
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P4 <- ggplot(data_alphadiv, aes(x=sex, y=InvSimpson, fill = sex )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'sex', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 35) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_sex.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs, fill = ethnicity)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "ethnicity",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 4, label.y = 850) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon, fill = ethnicity )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 4, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Fisher, fill = ethnicity)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'ethnicity', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 4, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=ethnicity, y=InvSimpson, fill = ethnicity )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'ethnicity', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 38) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_ethnicity.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=study_group, y=S.obs, fill = study_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "study group",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 750) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=study_group, y=Shannon, fill = study_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'study group', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=study_group, y=Fisher, fill = study_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'study group', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=study_group, y=InvSimpson, fill = study_group )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'study group', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_study_group.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=country_group, y=S.obs, fill = country_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "country group",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1.3, label.y = 800) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P2 <- ggplot(data_alphadiv, aes(x=country_group, y=Shannon, fill = country_group )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'country group', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 0.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P3 <- ggplot(data_alphadiv, aes(x=country_group, y=Fisher, fill = country_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'country group', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1.3, label.y = 45) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



P4 <- ggplot(data_alphadiv, aes(x=country_group, y=InvSimpson, fill = country_group )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'country_group', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1.3, label.y = 35) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_country_group.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=region, y=S.obs, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "region",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 800) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P2 <- ggplot(data_alphadiv, aes(x=region, y=Shannon, fill = region )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'region', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P3 <- ggplot(data_alphadiv, aes(x=region, y=Fisher, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Fisher', x= 'region', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +  
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 45) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



P4 <- ggplot(data_alphadiv, aes(x=region, y=InvSimpson, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'region', y= 'Inverse Simpson', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 2.5, label.y = 35) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_region.pdf", width = 40, height = 20, units = "cm")

#By two things ---- 

pv_kw_country <- data_alphadiv %>% group_by(country) %>%
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
            shannon = kruskal.test(Shannon ~ ethnicity)$p.value,
            fisher = kruskal.test(Fisher ~ ethnicity)$p.value,
            inverse_simpson = kruskal.test(InvSimpson ~ ethnicity)$p.value) %>%
  column_to_rownames("country") %>%
  signif(digits = 3) %>% 
  rownames_to_column("country")
  

P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "ethnicity",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(.~country, scales = "free") + 
  geom_text(data=pv_kw_country, aes(x=0, y=950, label=paste0("p=",richness)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~country, scales = "free") + 
  geom_text(data=pv_kw_country, aes(x=0, y=1.5, label=paste0("p=",shannon)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


P3 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Fisher)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Fisher', x= 'ethnicity', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~country, scales = "free") + 
  geom_text(data=pv_kw_country, aes(x=0, y=70, label=paste0("p=",fisher)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P4 <- ggplot(data_alphadiv, aes(x=ethnicity, y=InvSimpson)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Inverse Simpson', x= 'ethnicity', y= 'Inverse Simpson', tag = "D") +
  geom_point(show.legend = FALSE) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  facet_grid(~country, scales = "free")  + 
  geom_text(data=pv_kw_country, aes(x=0, y=37, label=paste0("p=",inverse_simpson)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



(P1 | P2) / (P3 | P4) # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_ethnicity_country.pdf", width = 40, height = 20, units = "cm")

pv_kw_region <- data_alphadiv %>% 
  filter(region %in% c("BW-Central West","BW-West", "CM-East", "CM-Southwest")) %>% 
  group_by(region) %>% 
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
                   shannon = kruskal.test(Shannon ~ ethnicity)$p.value,
                   fisher = kruskal.test(Fisher ~ ethnicity)$p.value,
                   inverse_simpson = kruskal.test(InvSimpson ~ ethnicity)$p.value)

anova_bw <- data_alphadiv %>%
  filter(region %in% c("BW-Northwest"))

anova_bw <- data.table(region = "BW-Northwest",
                       richness = summary(aov(S.obs ~ ethnicity, data = anova_bw))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ ethnicity, data = anova_bw))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ ethnicity, data = anova_bw))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ ethnicity, data = anova_bw))[[1]][["Pr(>F)"]][[1]])

anova_cm <- data_alphadiv %>%
  filter(region %in% c("CM-Northwest"))

anova_cm <- data.table(region = "CM-Northwest",
                       richness = summary(aov(S.obs ~ ethnicity, data = anova_cm))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ ethnicity, data = anova_cm))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ ethnicity, data = anova_cm))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ ethnicity, data = anova_cm))[[1]][["Pr(>F)"]][[1]])


pv_region <-rbind(pv_kw_region, anova_bw, anova_cm) %>% 
  column_to_rownames("region") %>%
  signif(digits = 3) %>% 
  rownames_to_column("region") 

rm(anova_cm, anova_bw)

P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "ethnicity",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(.~region, scales = "free") + 
  geom_text(data=pv_region, aes(x=0, y=800, label=paste0("p=",richness)),  size = 7/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +  
  facet_grid(~region, scales = "free") + 
  geom_text(data=pv_region, aes(x=0, y=6, label=paste0("p=",shannon)),  size = 7/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 





P3 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Fisher)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Fisher', x= 'ethnicity', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +   
  facet_grid(~region, scales = "free") + 
  geom_text(data=pv_region, aes(x=0, y=75, label=paste0("p=",fisher)),  size = 7/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




P4 <- ggplot(data_alphadiv, aes(x=ethnicity, y=InvSimpson)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE, position = position_dodge(0.8)) +
  labs(title= 'Inversion Simpson', x= 'ethnicity', y= 'Inverse Simpson', tag = "C") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +   
  facet_grid(.~region, scales = "free") + 
  geom_text(data=pv_region, aes(x=0, y=45, label=paste0("p=",inverse_simpson)),  size = 7/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



(P1 | P2) / (P3 | P4) # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_ethnicity_region.pdf", width = 40, height = 20, units = "cm")

by <- data_alphadiv %>%
  filter(ethnicity %in% c("Bagyeli"))


anova_by <- data.table(ethnicity = "Bagyeli",
                       richness = summary(aov(S.obs ~ sex, data = by))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = by))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = by))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = by))[[1]][["Pr(>F)"]][[1]])

bk <- data_alphadiv %>%
  filter(ethnicity %in% c("Baka"))

anova_bk <- data.table(ethnicity = "Baka",
                       richness = summary(aov(S.obs ~ sex, data = bk))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = bk))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = bk))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = bk))[[1]][["Pr(>F)"]][[1]])

br <- data_alphadiv %>%
  filter(ethnicity %in% c("Burunge"))

anova_br <- data.table(ethnicity = "Burunge",
                       richness = summary(aov(S.obs ~ sex, data = br))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = br))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = br))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = br))[[1]][["Pr(>F)"]][[1]])

fg <- data_alphadiv %>%
  filter(ethnicity %in% c("Fang"))

anova_fg <- data.table(ethnicity = "Fang",
                       richness = summary(aov(S.obs ~ sex, data = fg))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = fg))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = fg))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = fg))[[1]][["Pr(>F)"]][[1]])

fl <- data_alphadiv %>%
  filter(ethnicity %in% c("Mbororo Fulani"))

anova_fl <- data.table(ethnicity = "Mbororo Fulani",
                       richness = summary(aov(S.obs ~ sex, data = fl))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = fl))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = fl))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = fl))[[1]][["Pr(>F)"]][[1]])

hd <- data_alphadiv %>%
  filter(ethnicity %in% c("Hadzabe"))

anova_hd <- data.table(ethnicity = "Hadzabe",
                       richness = summary(aov(S.obs ~ sex, data = hd))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = hd))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = hd))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = hd))[[1]][["Pr(>F)"]][[1]])

kg <- data_alphadiv %>%
  filter(ethnicity %in% c("Kgalagadi"))

anova_kg <- data.table(ethnicity = "Kgalagadi",
                       richness = summary(aov(S.obs ~ sex, data = kg))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = kg))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = kg))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = kg))[[1]][["Pr(>F)"]][[1]])

ms <- data_alphadiv %>%
  na.omit() %>%
  filter(ethnicity %in% c("Maasai"))

anova_ms <- data.table(ethnicity = "Maasai",
                       richness = summary(aov(S.obs ~ sex, data = ms))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = ms))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = ms))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = ms))[[1]][["Pr(>F)"]][[1]])

ng <- data_alphadiv %>%
  filter(ethnicity %in% c("Ngoumba"))

anova_ng <- data.table(ethnicity = "Ngoumba",
                       richness = summary(aov(S.obs ~ sex, data = ng))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = ng))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = ng))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = ng))[[1]][["Pr(>F)"]][[1]])

nz <- data_alphadiv %>%
  filter(ethnicity %in% c("Nzime"))

anova_nz <- data.table(ethnicity = "Nzime",
                       richness = summary(aov(S.obs ~ sex, data = nz))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = nz))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = nz))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = nz))[[1]][["Pr(>F)"]][[1]])

tk <- data_alphadiv %>%
  filter(ethnicity %in% c("Tikari"))

anova_tk <- data.table(ethnicity = "Tikari",
                       richness = summary(aov(S.obs ~ sex, data = tk))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = tk))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = tk))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = tk))[[1]][["Pr(>F)"]][[1]])

bt <- data_alphadiv %>%
  filter(ethnicity %in% c("Bantu General"))

anova_bt <- data.table(ethnicity = "Bantu General",
                       richness = summary(aov(S.obs ~ sex, data = bt))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = bt))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = bt))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = bt))[[1]][["Pr(>F)"]][[1]])

ks <- data_alphadiv %>%
  filter(ethnicity %in% c("Khoesan"))

anova_ks <- data.table(ethnicity = "Khoesan",
                       richness = summary(aov(S.obs ~ sex, data = ks))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ sex, data = ks))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ sex, data = ks))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ sex, data = ks))[[1]][["Pr(>F)"]][[1]])



pv_sex <-rbind(anova_by, anova_bk, anova_br, anova_fg, anova_fl, anova_hd, anova_kg, anova_ms, anova_ng, anova_nz, anova_tk, anova_bt, anova_ks) %>% 
  column_to_rownames("ethnicity") %>%
  signif(digits = 3) %>% 
  rownames_to_column("ethnicity") %>%
  arrange(ethnicity)

rm(anova_by, anova_bk, anova_br, anova_fg, anova_fl, anova_hd, anova_mk, anova_kg, anova_ms, anova_ng, anova_nz, anova_tk, anova_bt, anova_ks, by, bk, br, fg, fl, hd, mk, kg, ms, ng, nz, tk, bt, ks)

data_sex <- data_alphadiv %>% 
  filter(ethnicity != "tswana") %>%
  group_by(ethnicity) %>% 
  filter(n() >= 7) %>% 
  na.omit()
 
P1 <- data_sex %>% 
  ggplot(aes(x=sex, y=S.obs)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title = "Richness ", 
       x = "sex",
       y = "richness", 
       tag = "A") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) + 
  facet_grid(~ethnicity, scales = "free", labeller = label_wrap_gen()) + 
  geom_text(data=pv_sex, aes(x=0, y=800, label=paste0("p=",richness)),  size = 7/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P2 <- data_sex %>% 
  ggplot(aes(x=sex, y=Shannon)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Shannon', x= 'sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +  
  facet_grid(~ethnicity, scales = "free")+ 
  geom_text(data=pv_sex, aes(x=0, y=4.5, label=paste0("p=",shannon)),  size = 7/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P3 <- data_sex %>% 
  ggplot(aes(x=sex, y=Fisher)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Fisher', x= 'sex', y= 'Fisher', tag = "C") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +   
  facet_grid(~ethnicity, scales = "free") + 
  geom_text(data=pv_sex, aes(x=0, y=70, label=paste0("p=",fisher)),  size = 7/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




P4 <- data_sex %>% 
  ggplot(aes(x=sex, y=InvSimpson)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Inverse Simpson', x= 'sex', y= 'Inverse Simpson', tag = "C") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +   
  facet_grid(~ethnicity, scales = "free") + 
  geom_text(data=pv_sex, aes(x=0, y=38, label=paste0("p=",inverse_simpson)),  size = 7/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



(P1 | P2) / (P3 | P4) # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_sex_ethnicity.pdf", width = 40, height = 20, units = "cm")

#by numeric values ----

linreg_age <- data.table(richness = summary(lm(S.obs~age, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~age, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~age, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~age, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=age, y=S.obs)) +
  labs(title = "Richness ", 
       x = "age",
       y = "richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=890, label=paste0("R^2=",richness)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P2 <- ggplot(data_alphadiv, aes(x=age, y=Shannon)) +
  labs(title= 'Shannon', x= 'age', y= 'Shannon', tag = "B") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=5, label=paste0("R^2=",shannon)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=age, y=Fisher)) +
  labs(title= 'Fisher', x= 'age', y= 'Fisher', tag = "C") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15))  +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=70, label=paste0("R^2=",fisher)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=age, y=InvSimpson)) +
  labs(title= 'Inverse Simpson', x= 'age', y= 'Inverse Simpson', tag = "D") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=38, label=paste0("R^2=",inverse_simpson)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_age.pdf", width = 40, height = 20, units = "cm")

linreg_bmi <- data.table(richness = summary(lm(S.obs~waist_circumference, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~waist_circumference, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~waist_circumference, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~waist_circumference, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=bmi, y=S.obs)) +
  labs(title = "Richness ", 
       x = "bmi",
       y = "richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_bmi, aes(x=50, y=890, label=paste0("R^2=",richness)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P2 <- ggplot(data_alphadiv, aes(x=bmi, y=Shannon)) +
  labs(title= 'Shannon', x= 'bmi', y= 'Shannon', tag = "B") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_bmi, aes(x=50, y=5, label=paste0("R^2=",shannon)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=bmi, y=Fisher)) +
  labs(title= 'Fisher', x= 'bmi', y= 'Fisher', tag = "C") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15))  +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_bmi, aes(x=50, y=70, label=paste0("R^2=",fisher)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=bmi, y=InvSimpson)) +
  labs(title= 'Inverse Simpson', x= 'bmi', y= 'Inverse Simpson', tag = "D") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_bmi, aes(x=50, y=38, label=paste0("R^2=",inverse_simpson)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_bmi.pdf", width = 40, height = 20, units = "cm")

linreg_wc <- data.table(richness = summary(lm(S.obs~waist_circumference, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~waist_circumference, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~waist_circumference, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~waist_circumference, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=S.obs)) +
  labs(title = "Richness ", 
       x = "waist circumference",
       y = "richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_wc, aes(x=50, y=890, label=paste0("R^2=",richness)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P2 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=Shannon)) +
  labs(title= 'Shannon', x= 'waist circumference', y= 'Shannon', tag = "B") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_wc, aes(x=50, y=5, label=paste0("R^2=",shannon)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P3 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=Fisher)) +
  labs(title= 'Fisher', x= 'waist circumference', y= 'Fisher', tag = "C") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15))  +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_wc, aes(x=50, y=70, label=paste0("R^2=",fisher)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P4 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=InvSimpson)) +
  labs(title= 'Inverse Simpson', x= 'waist circumference', y= 'Inverse Simpson', tag = "D") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18)) +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth(colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_wc, aes(x=50, y=38, label=paste0("R^2=",inverse_simpson)), size = 10/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2) / (P3 | P4)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_4/alpha_diversity_by_waist_circumference.pdf", width = 40, height = 20, units = "cm")



