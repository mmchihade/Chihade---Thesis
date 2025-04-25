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
library(sp)
library(sf)
library(rnaturalearth)
library(ggsflabel)
library(colorspace)


#load data ----
data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs("EPSG:4326")) 

data_grp <- data_grp %>%
  mutate(study_group = recode(study_group,
                            "Agropastoralist" = "Agriculturalists - raise livestock and grow their own crops; no regular dairy consumption (n = 224)",
                            "Hunter-gatherer" = "Hunter-gatherers - rely on the environment for their dietary needs; high pathogen exposure (n = 136)",
                            "Pastoralist" = "Pastoralists - raise livestock, diet high in dairy (n = 110)")) %>% 
  mutate(study_group = str_wrap(study_group, width = 23)) %>%
  mutate(ethnicity = recode(ethnicity,
                              "Mbororo Fulani" = "Mbororo Fulani (n=97)",
                              "Herero" = "Herero (n=2)",
                              "Maasai" = "Maasai (n=11) ",
                              "Bantu General" = "Bantu Speakers (n=54)",
                              "Fang" = "Fang (n=34)",
                              "Ngoumba"="Ngoumba (n=26)",
                              "Nzime" = "Nzime (n=9)",
                              "Tikari" = "Tikari (n=80)",
                              "Kgalagadi" = "Kgalagadi (n=10)",
                              "Burunge" = "Burunge (n=11)",
                              "Baka" = "Baka (n=69)",
                              "Bagyeli" = "Bagyeli (n=31)",
                              "Khoesan" = "Khoisan (n=17)",
                              "Hadzabe" = "Hadzabe (n=19)")) %>%
  select("lab_id", "country", "region", "sample_site", "year", "study_group", "ethnicity","ethnicity" ,"geometry")

world <- read_sf("data/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp") 

#subset data and turn it into spatial objects ----

africa <- subset(world, ADMIN %in% c("Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cabo Verde","Cameroon","Central African Republic","Chad","Comoros",
                                      "Democratic Republic of the Congo","Republic of the Congo","Ivory Coast","Djibouti","Egypt","Equatorial Guinea","Eritrea","Swaziland","Ethiopia","Gabon","Gambia","Ghana","Guinea","Guinea-Bissau","Kenya","Lesotho","Liberia","Libya","Madagascar","Malawi","Mali","Mauritania","Mauritius","Morocco","Mozambique","Namibia","Niger","Nigeria","Rwanda","Sao Tome and Principe","Senegal","Seychelles","Sierra Leone","Somaliland","Somalia","South Africa","South Sudan","Sudan","United Republic of Tanzania","Togo","Tunisia","Uganda","Zambia","Zimbabwe"))

africa <- africa %>% 
  clean_names() %>%
  select(-contains('name'), -contains('formal'), -contains('fclass'), -featurecla, -adm0_dif, -level, -tlc, -geou_dif, -su_dif, -contains('brk'), -contains("note"), -tlc_diff, -contains("adm0"))

study_countries <- africa %>% 
  filter(admin %in% c("Cameroon", "Botswana", "United Republic of Tanzania"))

cameroon <- africa %>% 
  filter(admin %in% c("Cameroon"))

tanzania <- africa %>% 
  filter(admin %in% c("United Republic of Tanzania"))

botswana <- africa %>% 
  filter(admin %in% c("Botswana"))


site_cam <- data_grp %>% 
  filter(country %in% c("Cameroon")) %>% 
  select(sample_site, geometry) 

cam <- data_grp %>% 
  filter(country %in% c("Cameroon")) %>% 
  select(ethnicity, study_group, geometry) %>% 
  distinct(ethnicity, geometry, .keep_all = TRUE) %>% 
  filter(!row_number() %in% c(9, 18, 11, 13, 10, 14, 17, 19, 16, 23, 15, 12, 22)) %>% 
  mutate(study_group = str_wrap(study_group, width = 20))
  

site_bot <- data_grp %>% 
  filter(country %in% c("Botswana")) %>% 
  select(sample_site, geometry) 


bot <- data_grp %>% 
  filter(country %in% c("Botswana")) %>% 
  select(ethnicity, study_group, geometry) %>% 
  distinct(ethnicity, geometry, .keep_all = TRUE) %>% 
  filter(!row_number() %in% c(2, 4, 9, 10))

site_tan <- data_grp %>% 
  filter(country %in% c("Tanzania")) %>% 
  select(sample_site, geometry)


tan <- data_grp %>% 
  filter(country %in% c("Tanzania")) %>% 
  select(ethnicity, study_group, geometry) %>% 
  distinct(ethnicity, geometry, .keep_all = TRUE)%>% 
  filter(!row_number() %in% c(1, 3, 6,7))


#load colors ---- 
site_color <- (c(wes_palette("Darjeeling1", 5, type = c("discrete")), wes_palette("Darjeeling2", 4, type = c("discrete")), wes_palette("Royal1", 4, type = c("discrete")), wes_palette("Royal2", 5, type = c("discrete")), wes_palette("Rushmore1", 5, type = c("discrete")), wes_palette("BottleRocket1", 6, type = c("discrete")), wes_palette("BottleRocket2", 5, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete")), wes_palette("Moonrise3", 5, type = c("discrete"))))

site_color <- lighten(
  site_color,
  amount = 0.5,
  method = c("relative", "absolute"),
  space = c("HCL", "HLS", "combined"),
  fixup = TRUE
)

show_col(site_color)

#let's plot! 

ggplot() + 
  geom_sf(data = africa, fill = "white", color = "black") +
  geom_sf(data = study_countries, fill = "#A5D46A") + 
  geom_sf_label_repel(data = study_countries, aes(label = geounit, fontface = "bold"), force = 100, nudge_x = -15, seed = 10, size = 10) + 
  theme_void() + 
  theme(rect = element_rect(fill = "transparent")) 

ggsave("output/largemap.png", bg = "transparent")


ggplot() + 
  geom_sf(data = cameroon, fill = "#A5D46A") + 
  geom_sf(data = site_cam, aes(shape = sample_site), size = 3, show.legend = FALSE) + 
  ggrepel::geom_label_repel(data = cam, aes(label = ethnicity, geometry = geometry, fill = study_group), stat = "sf_coordinates", box.padding = 1.2, max.overlaps = Inf, point.padding = 0.6, force = 3, arrow = arrow(length = unit(0.015, "npc")), nudge_x = -0.2, size = 6) + 
  guides(fill = guide_legend(title = "Subsistance Strategy", override.aes = list(label = "", size=6))) + 
  scale_fill_manual(values = as.vector(site_color)) +
  scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  ggtitle("Cameroon") + 
  theme(rect = element_rect(fill = "transparent")) + 
  theme_void() + 
  theme(plot.title = element_text(size = 30, face = "bold"))  + 
  theme(legend.text = element_text(size=20), legend.position=c(1.6, 0.5), legend.title=element_text(size=25), legend.key.width = unit(0.5,"cm")) + 
  theme(legend.key.height = unit(6, "cm"))

ggsave("output/cameroon.png", bg = "transparent", width = 20, height = 30, units = "cm")  


 
ggplot() + 
  geom_sf(data = botswana, fill = "#A5D46A") + 
  geom_sf(data = site_bot, aes(shape = sample_site), size = 3, show.legend = FALSE) + 
  ggrepel::geom_label_repel(data = bot, aes(label = ethnicity, geometry = geometry, fill = study_group ), stat = "sf_coordinates", box.padding = 1.2, max.overlaps = Inf, point.padding = 0.6, force = 3, arrow = arrow(length = unit(0.015, "npc")), nudge_x = -0.2, size = 5) + 
  guides(fill = guide_legend(title = "Subsistance Strategy", override.aes = aes(label = ""), size = 7)) + 
  scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  ggtitle("Botswana") + 
  theme(rect = element_rect(fill = "transparent")) + 
  theme_void() + 
  theme(plot.title = element_text(size = 30, face = "bold"), legend.key.width = unit(0.5,"cm")) + 
  theme(legend.key.height = unit(6, "cm"))


ggsave("output/botswana.png", bg = "transparent") 


ggplot() + 
  geom_sf(data = tanzania, fill = "#A5D46A") + 
  geom_sf(data = site_tan, aes(shape = sample_site), size = 3) + 
  ggrepel::geom_label_repel(data = tan, aes(label = ethnicity, geometry = geometry, fill = study_group ), stat = "sf_coordinates", box.padding = 1.2, max.overlaps = Inf, point.padding = 0.6, force = 3, arrow = arrow(length = unit(0.015, "npc")), nudge_x = -0.2, size = 7) + 
  guides(fill = guide_legend(title = "Subsistance Strategy", override.aes = aes(label = ""), size = 7)) + 
  scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_fill_manual(values=as.vector(site_color)) +
  ggtitle("Tanzania") + 
  theme(rect = element_rect(fill = "transparent")) + 
  theme_void()


ggsave("output/tanzania.png", bg = "transparent") 






