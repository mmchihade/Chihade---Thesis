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
library(sjmisc)
library(RColorBrewer)
library(vegan)
library(ggrepel)
library(ggpmisc)

#first I need to load all my values ---

cog_variance <- read.csv("data/cog_variance_values.csv", header = TRUE) %>% 
  pivot_longer(names_to = "grouping", values_to = "pvalue", pvcountry:stdev) %>% 
  mutate(grouping = recode(grouping,
                         "stdev"="sd", 
                         "pvcountry_group" = "pvcountrygroup", 
                         "pvstudy_group"="pvstudygroup", 
                         "Fstudy_group"="Fstudygroup", 
                         "Fcountry_group"="Fcountrygroup"))

cog_genus_variance <- read.csv("data/cog_genus_variance_means.csv", header = TRUE)

cog_count <- read.csv("data/cog_count.csv", header = TRUE)

#we need the cog categories too! 

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#now let's bind these together ----

variance <- merge(cog_genus_variance, cog_variance, by = c("cog_category", "grouping")) 

variance <- merge(variance, cog, by= "cog_category")

variance <- merge(variance, cog_count, by = "cog_category") %>% 
  filter(!cog_category %in% c("Y", "B", "Z"))

#and now we plot! 

site_color <- (c(wes_palette("Darjeeling1", 4, type = c("discrete")), wes_palette("GrandBudapest2", 2, type = c("discrete")),  wes_palette("Royal2", 1, type = c("discrete")), wes_palette("Darjeeling2", 3, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete"))))
show_col(site_color)

variance %>% 
  filter(grouping == "sd") %>%
  ggplot(aes(x=pvalue, y=mean)) +
  geom_point(aes(size=count, colour = label)) + 
  geom_smooth(colour="black", method='lm', se=FALSE, size=0.7) +  
  stat_poly_eq(use_label(c("eq", "R2")), size = 6) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Individual", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 8,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc")), 
                   nudge_x = 1) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +    
  guides(color = guide_legend(override.aes = list(size = 6))) + 
  theme(plot.title = element_text(face = "bold", size = 15)) + 
  theme(text = element_text(size = 15))
  


ggsave("output/q3/ind.pdf", width = 40, height = 20, units = "cm")


variance %>% 
  filter(grouping == "pvcountry") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Country", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 2,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") + 
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave("output/q3/country.pdf", width = 40, height = 20, units = "cm")

variance %>% 
  filter(grouping == "pvsex") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Sex", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 2,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +    guides(color = guide_legend(override.aes = list(size = 6)))

ggsave("output/q3/sex.pdf", width = 40, height = 20, units = "cm")

variance %>% 
  filter(grouping == "pvstudygroup") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Study Group", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 3,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +    guides(color = guide_legend(override.aes = list(size = 6)))

ggsave("output/q3/studygroup.pdf", width = 40, height = 20, units = "cm")

variance %>% 
  filter(grouping == "pvregion") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Region", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 3,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 6)))

ggsave("output/q3/region.pdf", width = 40, height = 20, units = "cm")

variance %>% 
  filter(grouping == "pvethnicity") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Ethnicity", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 3,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +    guides(color = guide_legend(override.aes = list(size = 6))) 

ggsave("output/q3/ethnicity.pdf", width = 40, height = 20, units = "cm")



variance %>% 
  filter(grouping == "pvcountrygroup") %>% 
  ggplot(aes(x=-log(pvalue), y=-log(mean), colour = label)) +
  geom_point(aes(size=count)) + 
  scale_color_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) + 
  labs(title = "Country Group", x= 'Variance  within a function category', y= 'Average variance in functional fraction of each taxons proteome', color = "COG Category") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=15,face="bold")) + 
  guides(color = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  geom_label_repel(aes(label = cog_category),
                   box.padding   = 0.35, 
                   point.padding = 0.2,
                   force = 3,
                   segment.color = 'black', 
                   color = "black",
                   max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.015, "npc"))) + 
  scale_size_binned(range = c(1, 15), n.breaks = 14, guide = "none") +    guides(color = guide_legend(override.aes = list(size = 6))) 

ggsave("output/q3/countrygroup.pdf", width = 40, height = 20, units = "cm")