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



#load data ----
data_otu <- read.csv("data/OTU_tables/kraken2_READS_sg_2024_decontam_G.csv", header = TRUE)
data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE)
  

data_grp <- data_grp %>%
  column_to_rownames("lab_id")


in_otu_core <- data_otu %>%
  mutate_all(~ifelse(is.na(.), "0", .)) %>% 
  column_to_rownames("lab_id") %>%
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
            y = "Richness", 
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

 (P1) / (P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_site.pdf", width = 40, height = 20, units = "cm")

country_comp <- list( c("Cameroon", "Tanzania"), c("Cameroon", "Botswana"), c("Botswana", "Tanzania"))
  
P1 <- ggplot(data_alphadiv, aes(x=country, y=S.obs, fill = country)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Country",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=13, label.x = 1.1, label.y = 750) + 
  stat_compare_means(comparisons = country_comp) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=country, y=Shannon, fill = country )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'country', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=13, label.x = 1, label.y = 0.5) + 
  stat_compare_means(comparisons = country_comp)+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_country.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=sex, y=S.obs, fill = sex)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Sex",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) + 
  stat_compare_means(method = "anova", label.y = 700, size = 10) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

P2 <- ggplot(data_alphadiv, aes(x=sex, y=Shannon, fill = sex )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) + 
  stat_compare_means(method = "anova", label.y = 5, size = 10) +  
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_sex.pdf", width = 40, height = 20, units = "cm")

P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs, fill = ethnicity)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 4, label.y = 850) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon, fill = ethnicity )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 4, label.y = 0.5) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))




(P1) / (P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity.pdf", width = 40, height = 20, units = "cm")

sg_comp <- list( c("Agropastoralist", "Pastoralist"), c("Agropastoralist", "Hunter-gatherer"), c("Hunter-gatherer", "Pastoralist"))

P1 <- ggplot(data_alphadiv, aes(x=study_group, y=S.obs, fill = study_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Subsistance Strategy",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +    
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=13, label.x = 1, label.y = 750) + 
  stat_compare_means(comparisons = sg_comp) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=study_group, y=Shannon, fill = study_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Subsistance Strategy', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  scale_x_discrete(labels = label_wrap(10)) +  
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=13, label.x = 1, label.y = 0.5) + 
  stat_compare_means(comparisons = sg_comp) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_study_group.pdf", width = 40, height = 20, units = "cm")

cg_comp <- list(c("Tanzania Agropastoralist", "Tanzania Hunter-gatherer"), c("Tanzania Agropastoralist", "Tanzania Pastoralist"), c("Tanzania Pastoralist", "Tanzania Hunter-gatherer")) 

P1 <- ggplot(data_alphadiv, aes(x=country_group, y=S.obs, fill = country_group)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Country group",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 1.3, label.y = 800) + 
  stat_compare_means(comparisons = cg_comp) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


P2 <- ggplot(data_alphadiv, aes(x=country_group, y=Shannon, fill = country_group )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Country group', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_compare_means(comparisons = cg_comp) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 0.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))




(P1) / (P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_country_group.pdf", width = 40, height = 20, units = "cm")

bw_comp <- list(c("BW-Central West", "BW-Northwest"), c("BW-Central West", "BW-West"), c("BW-Northwest", "BW-West")) 
                    
cm_comp <- list(c("CM-East", "CM-Northwest"), c("CM-East", "CM-Southwest"), c("CM-Southwest", "CM-Northwest")) 
                    
tz_comp <- list(c("TZ-Central", "TZ-North"), c("TZ-Central", "TZ-North Central"), c("TZ-North", "TZ-North Central"))


ggplot(data_alphadiv, aes(x=region, y=S.obs, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Region",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 950) +
  stat_compare_means(comparisons = bw_comp) +
  stat_compare_means(comparisons = cm_comp) +
  stat_compare_means(comparisons = tz_comp) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  coord_cartesian(ylim = c(900, 1270)) 


P2 <- ggplot(data_alphadiv, aes(x=region, y=Shannon, fill = region )) +
  geom_boxplot(show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Region', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 7)) +
  scale_fill_manual(values=as.vector(site_color)) +
  stat_kruskal_test(p.adjust.method="BH", na.rm=T, size=5, label.x = 3, label.y = 0.5) + 
  stat_compare_means(comparisons = bw_comp) +
  stat_compare_means(comparisons = cm_comp) +
  stat_compare_means(comparisons = tz_comp) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  coord_cartesian(ylim = c(0, 6.5))

(P1) / (P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_region.pdf", width = 40, height = 20, units = "cm")

#By two things ---- 

pv_kw_country <- data_alphadiv %>% group_by(country) %>%
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
            shannon = kruskal.test(Shannon ~ ethnicity)$p.value,
            fisher = kruskal.test(Fisher ~ ethnicity)$p.value,
            inverse_simpson = kruskal.test(InvSimpson ~ ethnicity)$p.value) %>%
  column_to_rownames("country") %>%
  signif(digits = 3) %>% 
  rownames_to_column("country")

pv_kw_country$richness <- p.adjust(pv_kw_country$richness, method = "BH")

pv_kw_country$shannon <- p.adjust(pv_kw_country$shannon, method = "BH")
  
P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(.~country, scales = "free") + 
  geom_text(data=pv_kw_country, aes(x=0, y=950, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 15/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~country, scales = "free") + 
  geom_text(data=pv_kw_country, aes(x=-2, y=2, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity_country.pdf", width = 40, height = 20, units = "cm")

pv_kw_studygroup <- data_alphadiv %>% group_by(study_group) %>%
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
                   shannon = kruskal.test(Shannon ~ ethnicity)$p.value,
                   fisher = kruskal.test(Fisher ~ ethnicity)$p.value,
                   inverse_simpson = kruskal.test(InvSimpson ~ ethnicity)$p.value) %>%
  column_to_rownames("study_group") %>%
  signif(digits = 3) %>% 
  rownames_to_column("study_group")

pv_kw_studygroup$richness <- p.adjust(pv_kw_studygroup$richness, method = "BH")

pv_kw_studygroup$shannon <- p.adjust(pv_kw_studygroup$shannon, method = "BH")

P1 <- ggplot(data_alphadiv, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~study_group, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 12, face = "bold")) +  
  geom_text(data=pv_kw_studygroup, aes(x=0, y=1000, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 15/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  coord_cartesian(ylim = c(940, 1150)) 


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~study_group, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 12, face = "bold")) +  
  geom_text(data=pv_kw_studygroup, aes(x=0, y=1, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 14/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  coord_cartesian(ylim = c(-1, 5)) 

(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity_studygroup.pdf", width = 40, height = 20, units = "cm")


pv_kw_region <- data_alphadiv %>% 
  filter(region %in% c("BW-Central West","BW-West", "CM-East", "CM-Southwest")) %>% 
  group_by(region) %>% 
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
            shannon = kruskal.test(Shannon ~ ethnicity)$p.value,
            fisher = kruskal.test(Fisher ~ ethnicity)$p.value,
            inverse_simpson = kruskal.test(InvSimpson ~ ethnicity)$p.value)

pv_kw_region$richness <- p.adjust(pv_kw_region$richness, method = "BH")

pv_kw_region$shannon <- p.adjust(pv_kw_region$shannon, method = "BH")

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
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~ region, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.001,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  geom_text(data=pv_region, aes(x=0, y=1015, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 13/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  coord_cartesian(ylim = c(950, 1130)) 


P2 <- ggplot(data_alphadiv, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +  
  facet_grid(~ region, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  geom_text(data=pv_region, aes(x=0, y=6, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 13/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity_region.pdf", width = 45, height = 20, units = "cm")

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

pv_sex$richness <- p.adjust(pv_sex$richness, method = "BH")

pv_sex$shannon <- p.adjust(pv_sex$shannon, method = "BH")

rm(anova_by, anova_bk, anova_br, anova_fg, anova_fl, anova_hd, anova_mk, anova_kg, anova_ms, anova_ng, anova_nz, anova_tk, anova_bt, anova_ks, by, bk, br, fg, fl, hd, mk, kg, ms, ng, nz, tk, bt, ks)

data_sex <- data_alphadiv %>% 
  filter(ethnicity != "Tswana") %>%
  group_by(ethnicity) %>% 
  filter(n() >= 7) %>% 
  na.omit()
 
P1 <- data_sex %>% 
  ggplot(aes(x=sex, y=S.obs)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title = "Richness ", 
       x = "Sex",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) + 
  facet_grid(~ethnicity, switch = "y", scales = "free", labeller = label_wrap_gen()) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 12, face = "bold")) + 
  geom_text(data=pv_sex, aes(x=0, y=900, label=paste0("p=",richness)),  size = 13/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P2 <- data_sex %>% 
  ggplot(aes(x=sex, y=Shannon)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Shannon', x= 'Sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +  
  facet_grid(~ethnicity, labeller = label_wrap_gen(), switch = "y", scales = "free") + 
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  geom_text(data=pv_sex, aes(x=0, y=5, label=paste0("p=",shannon)),  size = 13/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_sex_ethnicity.pdf", width = 40, height = 20, units = "cm")

pv_kw_region_sg <- data_alphadiv %>% 
  filter(region %in% c("BW-Central West","BW-West", "CM-East", "CM-Southwest")) %>% 
  group_by(region) %>% 
  dplyr::summarize(richness = kruskal.test(S.obs ~ study_group)$p.value,
                   shannon = kruskal.test(Shannon ~ study_group)$p.value,
                   fisher = kruskal.test(Fisher ~ study_group)$p.value,
                   inverse_simpson = kruskal.test(InvSimpson ~ study_group)$p.value)

anova_bw_sg <- data_alphadiv %>%
  filter(region %in% c("BW-Northwest"))

anova_bw_sg <- data.table(region = "BW-Northwest",
                       richness = summary(aov(S.obs ~ study_group, data = anova_bw_sg))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ study_group, data = anova_bw_sg))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ study_group, data = anova_bw_sg))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ study_group, data = anova_bw_sg))[[1]][["Pr(>F)"]][[1]])

anova_cm_sg <- data_alphadiv %>%
  filter(region %in% c("CM-Northwest"))

anova_cm_sg <- data.table(region = "CM-Northwest",
                       richness = summary(aov(S.obs ~ study_group, data = anova_cm_sg))[[1]][["Pr(>F)"]][[1]],
                       shannon=summary(aov(Shannon ~ study_group, data = anova_cm_sg))[[1]][["Pr(>F)"]][[1]],
                       fisher=summary(aov(Fisher ~ study_group, data = anova_cm_sg))[[1]][["Pr(>F)"]][[1]], 
                       inverse_simpson=summary(aov(InvSimpson ~ study_group, data = anova_cm_sg))[[1]][["Pr(>F)"]][[1]])


pv_region_sg <-rbind(pv_kw_region_sg, anova_bw_sg, anova_cm_sg) %>% 
  column_to_rownames("region") %>%
  signif(digits = 3) %>% 
  rownames_to_column("region") 

pv_region_sg$richness <- p.adjust(pv_region_sg$richness, method = "BH")

pv_region_sg$shannon <- p.adjust(pv_region_sg$shannon, method = "BH")



P1 <- ggplot(data_alphadiv, aes(x=study_group, y=S.obs)) +
  geom_boxplot(aes(fill = study_group), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~ region, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.001,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  geom_text(data=pv_region_sg, aes(x=0, y=1015, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 13/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) + 
  coord_cartesian(ylim = c(950, 1130)) 


P2 <- ggplot(data_alphadiv, aes(x=study_group, y=Shannon)) +
  geom_boxplot(aes(fill = study_group ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +  
  facet_grid(~ region, switch = "y", scales = "free") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  geom_text(data=pv_region_sg, aes(x=0, y=6, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 13/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_study_group_region.pdf", width = 45, height = 20, units = "cm")



#by numeric values ----

linreg_age <- data.table(richness = summary(lm(S.obs~age, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~age, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~age, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~age, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=age, y=S.obs)) +
  labs(title = "Richness ", 
       x = "Age",
       y = "Richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = FALSE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=890, label=paste0("R^2=",richness)), size = 17/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P2 <- ggplot(data_alphadiv, aes(x=age, y=Shannon)) +
  labs(title= 'Shannon', x= 'Age', y= 'Shannon', tag = "B") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  geom_text(data=linreg_age, aes(x=50, y=5, label=paste0("R^2=",shannon)), size = 17/.pt) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_age.pdf", width = 40, height = 20, units = "cm")

linreg_bmi <- data.table(richness = summary(lm(S.obs~waist_circumference, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~waist_circumference, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~waist_circumference, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~waist_circumference, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=bmi, y=S.obs)) +
  labs(title = "Richness ", 
       x = "BMI",
       y = "Richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = FALSE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) + 
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  stat_poly_eq(use_label(c("eq", "R2", "P")), size = 6, label.y = "bottom") +  
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P2 <- ggplot(data_alphadiv, aes(x=bmi, y=Shannon)) +
  labs(title= 'Shannon', x= 'BMI', y= 'Shannon', tag = "B", 
       color = "Country and 
Subsistance Strategy", 
       shape = "Country and 
Subsistance Strategy") +
  geom_point(aes(color = country_group, shape = country_group)) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth(colour="black", method='lm', se=FALSE, size=0.7) + 
  stat_poly_eq(use_label(c("eq", "R2", "P")), size = 6, label.y = "bottom") +  
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_bmi.pdf", width = 40, height = 20, units = "cm")

linreg_wc <- data.table(richness = summary(lm(S.obs~waist_circumference, data=data_alphadiv))$r.squared,
                         shannon = summary(lm(Shannon~waist_circumference, data=data_alphadiv))$r.squared,
                         fisher = summary(lm(Fisher~waist_circumference, data=data_alphadiv))$r.squared,
                         inverse_simpson = summary(lm(InvSimpson~waist_circumference, data=data_alphadiv))$r.squared) %>%
  signif(digits = 3)

P1 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=S.obs)) +
  labs(title = "Richness ", 
       x = "Waist circumference",
       y = "Richness", 
       tag = "A") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = FALSE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  stat_poly_eq(use_label(c("eq", "R2", "P")), size = 5, label.y = "bottom") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



P2 <- ggplot(data_alphadiv, aes(x=waist_circumference, y=Shannon)) +
  labs(title= 'Shannon', x= 'Waist circumference', y= 'Shannon', tag = "B", 
       color = "Country and 
Subsistance Strategy", 
       shape = "Country and 
Subsistance Strategy") +
  geom_point(aes(color = country_group, shape = country_group), show.legend = TRUE) +
  scale_shape_manual(values=c(1, 2, 3, 4, 9, 15, 16, 17, 18))+
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = label_wrap(8)) +    
  scale_color_manual(values=as.vector(site_color)) +
  geom_smooth( colour="black", method='lm', se=FALSE, size=0.7) + 
  stat_poly_eq(use_label(c("eq", "R2", "P")), size = 5, label.y = "bottom") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))



(P1 | P2)  # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_waist_circumference.pdf", width = 40, height = 20, units = "cm")



#by just one region or group ----

data_alphadiv_cm_nw <- data_alphadiv %>% 
  filter(region %in% c("CM-Northwest"))

pv_region_cm_nw <- pv_region %>% 
  filter(region %in% c("CM-Northwest"))

P1 <- ggplot(data_alphadiv_cm_nw, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#046C9A","#F8A")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +
  geom_text(data=pv_region_cm_nw, aes(x=0, y=1000, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




P2 <- ggplot(data_alphadiv_cm_nw, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#046C9A","#F8A")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +  
  geom_text(data=pv_region_cm_nw, aes(x=0, y=4, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 10/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity_region_cm_nw.pdf", width = 20, height = 20, units = "cm")

data_alphadiv_cm_nw<- data_alphadiv %>% 
  filter(region %in% c("CM-Northwest")) %>% 
  rownames_to_column("lab_id")


pv_region_cm_nw <- pv_region %>% 
  filter(region %in% c("CM-Northwest"))

P1 <- ggplot(data_alphadiv_cm_nw, aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#C93312", "#9A8822")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +
  geom_text(data=pv_region_cm_nw, aes(x=0, y=1000, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 20/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




P2 <- ggplot(data_alphadiv_cm_nw, aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity ), show.legend = FALSE) +
  labs(title= 'Shannon', x= 'Ethnicity', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#C93312", "#9A8822")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  +  
  geom_text(data=pv_region_cm_nw, aes(x=0, y=4, label=paste0("Kruskal-Wallis\n p=",shannon)),  size = 20/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_ethnicity_region_cm_nw.pdf", width = 20, height = 20, units = "cm")

data_sex_hd <- data_sex %>% 
  filter(ethnicity %in% c("Hadzabe"))

pv_sex_hd <- pv_sex %>% 
  filter(ethnicity %in% c("Hadzabe"))

P1 <- data_sex_hd %>% 
  ggplot(aes(x=sex, y=S.obs)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title = "Richness ", 
       x = "Sex",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) + 
  facet_grid(~ethnicity, scales = "free", labeller = label_wrap_gen()) + 
  geom_text(data=pv_sex_hd, aes(x=0, y=1070, label=paste0("p=",richness)),  size = 15/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 



P2 <- data_sex_hd %>% 
  ggplot(aes(x=sex, y=Shannon)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Shannon', x= 'Sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +  
  facet_grid(~ethnicity, scales = "free")+ 
  geom_text(data=pv_sex_hd, aes(x=0, y=4, label=paste0("p=",shannon)),  size = 15/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 




(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_sex_ethnicity_hd.pdf", width = 20, height = 20, units = "cm")

data_sex_ms <- data_sex %>% 
  filter(ethnicity %in% c("Maasai"))

pv_sex_ms <- pv_sex %>% 
  filter(ethnicity %in% c("Maasai"))

P1 <- data_sex_ms %>% 
  ggplot(aes(x=sex, y=S.obs)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title = "Richness ", 
       x = "Sex",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) + 
  facet_grid(~ethnicity, scales = "free", labeller = label_wrap_gen()) + 
  geom_text(data=pv_sex_ms, aes(x=0, y=1110, label=paste0("p=",richness)),  size = 15/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


P2 <- data_sex_ms %>% 
  ggplot(aes(x=sex, y=Shannon)) +
  geom_boxplot(show.legend = FALSE, aes(fill = sex)) +
  labs(title= 'Shannon', x= 'Sex', y= 'Shannon', tag = "B") +
  geom_point(show.legend = FALSE, aes(fill = sex)) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  scale_x_discrete(labels = label_wrap(10)) +  
  facet_grid(~ethnicity, scales = "free")+ 
  geom_text(data=pv_sex_ms, aes(x=0, y=4, label=paste0("p=",shannon)),  size = 15/.pt, vjust = 0.9, hjust = -0.1) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20)) 


(P1) / (P2) # all plots together using the patchwork package

ggsave("output/alpha_diversity_2/alpha_diversity_by_sex_ethnicity_ms.pdf", width = 20, height = 20, units = "cm")


pv_kw_cam <- data_alphadiv %>% filter(country == "Cameroon") %>%
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
                   shannon = kruskal.test(Shannon ~ ethnicity)$p.value, 
                   country = "Cameroon") %>%
  column_to_rownames("country") %>%
  signif(digits = 3) %>% 
  rownames_to_column("country")

cam_comp1 <- list(c("Baka", "Bantu General"),  
                 c("Fang", "Mbororo Fulani"), 
                 c("Mbororo Fulani", "Ngoumba")) 

cam_comp2 <- list(c("Baka", "Fang"),   
                 c("Mbororo Fulani", "Nzime")) 

cam_comp3 <- list(c("Ngoumba", "Tikari"), 
                 c("Baka", "Mbororo Fulani")) 

cam_comp4 <- list(c("Bantu General", "Mbororo Fulani"),
                 c("Mbororo Fulani", "Tikari")) 

cam_comp5 <- list(c("Bagyeli", "Mbororo Fulani"), 
                 c("Baka", "Ngoumba"), 
                 c("Fang", "Tikari"), 
                 c("Bantu General", "Tikari"),  
                 c("Baka", "Tikari"), 
                 c("Bagyeli", "Tikari"))

data_alphadiv %>% 
  filter(country == "Cameroon") %>%
  ggplot(aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#FF0000","#00A08A","#F2AD00","#5BBCD6","#C93312"
                             ,"#FAEFD1","#DC863B", "#9A8822")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  facet_grid(~ country, switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) + 
  stat_compare_means(comparisons = cam_comp1, aes(label=..p.signif..), step.increase = 0, size =5) +
  stat_compare_means(comparisons = cam_comp2, aes(label=..p.signif..), step.increase = 0, label.y = 1150) +
  stat_compare_means(comparisons = cam_comp3, aes(label=..p.signif..), step.increase = 0, label.y = 1190) +
  stat_compare_means(comparisons = cam_comp4, aes(label=..p.signif..), step.increase = 0, label.y = 1220) +
  stat_compare_means(comparisons = cam_comp5, aes(label=..p.signif..), label.y = c(1250, 1270, 1300, 1320, 1340)) +
  geom_text(data=pv_kw_cam, aes(x=0, y=900, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 18/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/alpha_diversity_2/alpha_diversity_by_country_ethnicity_cam.pdf", width = 30, height = 30, units = "cm")


pv_kw_tan <- data_alphadiv %>% filter(country == "Tanzania") %>%
  dplyr::summarize(richness = kruskal.test(S.obs ~ ethnicity)$p.value,
                   shannon = kruskal.test(Shannon ~ ethnicity)$p.value, 
                   country = "Tanzania") %>%
  column_to_rownames("country") %>%
  signif(digits = 3) %>% 
  rownames_to_column("country")


tan_comp <- list( c("Burunge", "Hadzabe"), c("Hadzabe", "Maasai"), c("Burunge", "Maasai"))

P1 <- data_alphadiv %>% 
  filter(country == "Tanzania") %>%
  ggplot(aes(x=ethnicity, y=S.obs)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Richness ", 
       x = "Ethnicity",
       y = "Richness", 
       tag = "A") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#F98400", "#ECCBAE", "#899DA4")) + 
  facet_grid(~ country, switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) +  
  stat_compare_means(comparisons = tan_comp, label.y = c(1130, 1130, 1140)) +
  geom_text(data=pv_kw_tan, aes(x=0,y=1090, label=paste0("Kruskal-Wallis\n p=",richness)),  size = 18/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

P2 <- data_alphadiv %>% 
  filter(country == "Tanzania") %>%
  ggplot(aes(x=ethnicity, y=Shannon)) +
  geom_boxplot(aes(fill = ethnicity), show.legend = FALSE) +
  labs(title = "Shannon", 
       x = "Ethnicity",
       y = "Shannon", 
       tag = "B") +
  geom_point(show.legend = FALSE) +
  theme(text = element_text(size = 10)) +
  scale_fill_manual(values=c("#F98400", "#ECCBAE", "#899DA4")) + 
  facet_grid(~ country, switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 15, face = "bold")) +  
  stat_compare_means(comparisons = tan_comp, label.y = c(4.5, 4.5, 5)) +
  geom_text(data=pv_kw_tan, aes(x=0, y=5, label=paste0("Kruskal-Wallis\n p=", shannon)),  size = 18/.pt, vjust = 0.9, hjust = -0.2) + 
  theme(axis.text.y=element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

(P1)/(P2)

ggsave("output/alpha_diversity_2/alpha_diversity_by_country_ethnicity_tan.pdf", width = 40, height = 30, units = "cm")

