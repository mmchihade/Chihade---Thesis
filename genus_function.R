# load packages ---- 
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

#load data that I will to make various groupings ----

data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_alpha_2024.csv", header=TRUE) %>%
  select(-record_id, -lab_id, -topmed_id, -collection_name, -sample_site, -longitude, -latitude, -date, -year, -is_adult, -age, -waist_circumference, -body_fat_percent, -bmi, -subsistance) %>% 
  rename(lab_id = field_id) %>%
  mutate(lab_id = toupper(lab_id)) 

high_low <- read.csv("data/OTU_tables/high_low_values_0.1_threshold.csv", header = TRUE) %>% 
  select(-true_total, -total_root, -unclassified) %>%
  pivot_longer(names_to = "genus", values_to = "abundance", g__Prevotella:g__Methanobacterium) %>% 
  mutate(across(everything(), gsub, pattern = "g__", replacement = "")) %>% 
  mutate(across(everything(), gsub, pattern = "\\.", replacement = "-"))

id_merge <- read.csv("data/clean_simple_metadata_tishkofflab_alpha_2024.csv", header=TRUE) %>% 
  select(field_id, lab_id) %>% 
  rename(subject = lab_id) %>% 
  rename(lab_id = field_id)

high_low <- merge(id_merge, high_low, by = "subject") %>% 
  select(-subject) 

#create a dataframe with all of the cog tally files, with the file name as a column so that we can match column to tally

merged_df <- list.files(path = "data/cog_tally_overall",
                        pattern = "\\.csv$", 
                        full.names = TRUE) %>% 
  set_names() %>%
  map_dfr(read_csv, .id = "file_name") 

#we need the cog categories too! 

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#now let's clean the data up  

#I also have to replace all of the TZ codes with B codes so that I can divide by country and ethnicity 

tally <- clean_names(merged_df) %>% 
  rename(lab_id = file_name) %>%
  mutate(across(everything(), gsub, pattern = "data/cog_tally_overall/", replacement = "")) %>%
  mutate(across(everything(), gsub, pattern = ".csv", replacement = "")) %>%
  mutate(n = as.numeric(n)) %>%
  mutate(lab_id = recode(lab_id,
                         "TZBG067"="B607",
                         "TZBG069"="B609",
                         "TZBG070"="B610",
                         "TZBG072"="B612",
                         "TZBG073"="B610",
                         "TZBG074"="B614",
                         "TZBG076"="B616",
                         "TZBG0771"="B617.1",
                         "TZBG0772"="B617.2",
                         "TZBG079"="B619",
                         "TZBG080"="B620",
                         "TZHZ336"="B082",
                         "TZHZ337"="B083",
                         "TZHZ338"="B084",
                         "TZHZ339"="B085",
                         "TZHZ340"="B086",
                         "TZHZ343"="B089",
                         "TZHZ344"="B090",
                         "TZHZ345"="B091",
                         "TZHZ346"="B092",
                         "TZHZ453"="B202",
                         "TZHZ457"="B206",
                         "TZHZ463"="B212",
                         "TZHZ469"="B218",
                         "TZHZ472"="B221",
                         "TZHZ474"="B223",
                         "TZHZ488"="B238",
                         "TZHZ559"="B312",
                         "TZHZ586"="B339",
                         "TZHZ587"="B340",
                         "TZMS164"="B464",
                         "TZMS177"="B477",
                         "TZMS178"="B478",
                         "TZMS183"="B483",
                         "TZMS190"="B490",
                         "TZMS192"="B492",
                         "TZMS195"="B495",
                         "TZMS196"="B496",
                         "TZMS197"="B497",
                         "TZMS198"="B498",
                         "TZMS199"="B499")) %>% 
  mutate(cog_category = recode(cog_category,
                               "?" = "S",
                               "!" = "S")) %>% 
  group_by(lab_id, cog_category, genus) %>% 
  mutate(tally = sum(n, na.rm =TRUE, group_by=TRUE)) %>%
  select(-n) %>%
  distinct()

#I also want to separate all of the multiple letter columns. They will probably be useful later on, but it is hella overwhelming now

full_tally <- tally %>% 
  separate(cog_category, into = c("cog_1", "cog_2", "cog_3", "cog_4", "cog_5", "cog_6", "cog_7"), sep = "") %>%
  select(-cog_1) %>% 
  pivot_longer(names_to="cog_name", values_to = "cog_category", cog_2:cog_7)%>%
  na.omit() %>%
  select(-cog_name) %>%
  group_by(lab_id, cog_category, genus) %>% 
  mutate(tally = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#umm, I am going to try to make a combined dataset, let's see if R can handle it ----

merge <- merge(full_tally, data_grp, by = "lab_id")

merge <-merge(merge, high_low, by = c("lab_id", "genus"))

merge_labels <- merge(merge, cog, by = "cog_category") %>% 
  filter(!cog_category %in% c("Y", "B", "Z"))

#yay!!!!!

#now let's split this up by a bunch of variables 

cameroon_label <- merge_labels %>% filter(country == "Cameroon")

botswana_label <- merge_labels %>% filter(country == "Botswana")

tanzania_label <- merge_labels %>% filter(country == "Tanzania") 

ag_label <- merge_labels %>% filter(study_group == "Agropastoralist")

pas_label <- merge_labels %>% filter(study_group == "Pastoralist")

hg_label <- merge_labels %>% filter(study_group == "Hunter-gatherer")

high_label <- merge_labels %>% filter(abundance == "high")

low_label <- merge_labels %>% filter(abundance == "low")

super_low_label <- merge_labels %>% filter(abundance == "super_low")

alpha_high <- merge_labels %>% filter(S.obs > 1100) 

alpha_low <- merge_labels %>% filter(S.obs <= 1100)

all <- merge_labels %>% filter(genus %in% c("Alistipes",
                                            "Bacteroides",
                                            "CAG-791",
                                            "Klebsiella",
                                            "Lachnospira",
                                            "Odoribacter",
                                            "RC9",
                                            "Oribacterium",
                                            "Sodaliphilus",
                                            "Succinivibrio",
                                            "Sutterella",
                                            "UBA4372",
                                            "V9D3004",
                                                  "Eubacterium_J",
                                                  "Lachnospira",
                                                  "Slackia_A",
                                                  "Treponema_D",
                                                  "UBA1740",
                                                  "UBA2882",
                                                  "UBA3282",
                                                  "UBA4644",
                                                  "UMGS1484",
                                                  "CAG-45",
                                                  "Butyricicoccus_A",
                                                  "CAG-110",
                                                  "CAG-238"))


#now let's get rid of the individual groupings and just look at function per genus ---- 

#I am feeling lazy so let's use lapply here: 

merge_list <- list(merge_labels, cameroon_label, botswana_label, tanzania_label, ag_label, pas_label, hg_label, high_label, low_label, super_low_label, alpha_high, alpha_low, all)

#here is my function 

delete_id <- function(i) { 
  options(digits=15)
  sum <- i %>% 
    group_by(genus, cog_category) %>% 
    summarise(sum = sum(tally))
  
  merge(sum, cog, by=c("cog_category")) %>% #and adding back the important info here
    select(genus, cog_category, sum, label)
}

merge_list_final <- lapply(merge_list, delete_id)

#wait we can do the same thing to get the percent of each function for each genus

findPerecent <- function(i) { 
   i %>% 
    group_by(genus) %>%
    mutate(total_per_genus = sum(sum, group_by = TRUE)) %>%
    mutate(percent_by_genus = sum/total_per_genus) %>% 
    select(genus, cog_category, sum, total_per_genus, percent_by_genus, label) %>% 
    distinct() %>%
    mutate(id = cur_group_id()) %>% 
    mutate(test = sum(percent_by_genus))
}

merge_list_percent <- lapply(merge_list_final, findPerecent)

#I also want to find some p-values between certain comparisons ----  

country_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, country, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(country_distinct, "country_delete_list.csv")

abun_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, abundance, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(abun_distinct, "abundance_delete_list.csv")

alpha_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, S.obs, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(alpha_distinct, "alpha_delete_list.csv")

test <- merge_labels %>% 
  filter(!genus %in% c("Acutalibacter",
                       "Aeromonas",
                       "CAG-353",
                       "CAG-590",
                       "CAG-632",
                       "Cedecea",
                       "Clostridium",
                       "Clostridium_A",
                       "Clostridium_E",
                       "Eisenbergiella",
                       "Enterococcus",
                       "HGM12650",
                       "HGM13222",
                       "Intestinibacter",
                       "Intestinimonas",
                       "Lactobacillus",
                       "Lactococcus_A",
                       "Mitsuokella",
                       "Phyllobacterium",
                       "Plesiomonas",
                       "Proteus",
                       "SFJ001",
                       "UBA7597",
                       "UMGS1071",
                       "UMGS1241",
                       "UMGS1370",
                       "UMGS172",
                       "UMGS1889",
                       "UMGS973",
                       "WG-1",
                       "Bacteroides_F",
                       "CAG-217",
                       "CAG-533",
                       "CAG-877",
                       "Eikenella",
                       "Enterobacter_D",
                       "F082",
                       "Marvinbryantia",
                       "UBA1740",
                       "UBA2883",
                       "UBA9414",
                       "UMGS1810",
                       "CAG-81",
                       "Catenibacterium",
                       "Citrobacter_B",
                       "Evtepia",
                       "Haemophilus_A",
                       "MGYG000002683",
                       "Pseudomonas_E",
                       "UBA11524",
                       "UBA737",
                       "UMGS874",
                       "MGYG000000674",
                       "OF09-33XD",
                       "Pantoea",
                       "PeH17",
                       "MGYG000004445",
                       "Pseudescherichia",
                       "Tidjanibacter",
                       "UMGS1251",
                       "Aggregatibacter",
                       "Dialister",
                       "Eggerthella",
                       "OM05-12",
                       "Amulumruptor",
                       "Leclercia",
                       "Scardovia",
                       "Mesosutterella",
                       "Serratia",
                       "VUNI01",
                       "Raoultella",
                       "Turicimonas",
                       "CAG-45",
                       "CAG-41",
                       "Muribaculum",
                       "MGYG000004418",
                       "Methanobacterium",
                       "Bacillus_O",
                       "HGM08974",
                       "Nosocomiicoccus",
                       "MGYG000003416",
                       "Rubeoparvulum",
                       "UMGS1260",
                       "Winkia",
                       "Bacillus_BD",
                       "Epilithonimonas",
                       "Spiro-01",
                       "Niallia",
                       "28L",
                       "Lactococcus_A",
                       "HGM05376",
                       "Enterovibrio",
                       "Parascardovia",
                       "MGYG000004418",
                       "UBA7488",
                       "Massilibacterium",
                       "HGM08974",
                       "UBA7862",
                       "Thalassobacillus",
                       "Aliarcobacter",
                       "MGYG000000407",
                       "Flaviflexus",
                       "Anaerosalibacter",
                       "UBA3375",
                       "MGYG000004418",
                       "Methanomassiliicoccus",
                       "Jonquetella",
                       "Epilithonimonas",
                       "Fructilactobacillus",
                       "Methanobacterium",
                       "UBA4877",
                       "Mycoplasmopsis",
                       "MGYG000000407",
                       "Niallia",
                       "Exiguobacterium_A",
                       "MGYG000003235",
                       "Nosocomiicoccus",
                       "Winkia",
                       "MGYG000003772",
                       "MGYG000004144",
                       "Metamycoplasma",
                       "Numidum",
                       "Aliarcobacter",
                       "UBA7862",
                       "Epilithonimonas",
                       "MGYG000003562",
                       "Ureaplasma",
                       "UBA7488",
                       "UMGS1260",
                       "Psychrobacillus",
                       "Methanobacterium",
                       "Methanocorpusculum",
                       "Acetomicrobium",
                       "UBA1103",
                       "Helicobacter_F",
                       "MGYG000004418",
                       "Rubeoparvulum",
                       "MGYG000003562",
                       "Alloscardovia",
                       "Arachnia",
                       "Fructilactobacillus",
                       "Gleimia",
                       "ISO4-G1",
                       "MGYG000003809",
                       "Pradoshia",
                       "UBA1103",
                       "UBA7741",
                       "UMGS743",
                       "Ureaplasma",
                       "Weizmannia",
                       "MGYG000000407",
                       "MGYG000004144",
                       "UBA4877",
                       "Aliarcobacter",
                       "Epilithonimonas",
                       "Lagierella",
                       "MGYG000003562",
                       "Parascardovia",
                       "Thalassobacillus",
                       "Bacillus_BD",
                       "Chryseobacterium",
                       "Olegusella",
                       "Paenisporosarcina",
                       "Anaerosalibacter",
                       "Moellerella",
                       "Planococcus",
                       "Vitreoscilla",
                       "KA00274",
                       "Listeria",
                       "MGYG000002053",
                       "Turicimonas",
                       "Garciella",
                       "HGM05376",
                       "Helicobacter_F",
                       "28L",
                       "Flavobacterium",
                       "Desulfitobacterium",
                       "Acetomicrobium",
                       "MGYG000004418",
                       "Mycoplasmopsis",
                       "Exiguobacterium_A",
                       "Metamycoplasma",
                       "Methanobacterium",
                       "UMGS1260",
                       "UBA3375",
                       "Acetomicrobium",
                       "Bacillus_O",
                       "Epilithonimonas",
                       "HGM10890",
                       "MGYG000000286",
                       "MGYG000003416",
                       "MGYG000004606",
                       "Massilibacterium",
                       "Methanocorpusculum",
                       "Pseudomonas_B",
                       "Staphylococcus_A",
                       "Timonella",
                       "UBA7862",
                       "Clostridium_AA",
                       "Flaviflexus",
                       "MGYG000003389",
                       "Nigerium",
                       "Numidum",
                       "Paenibacillus_F",
                       "Psychrobacillus",
                       "Rikenella",
                       "SFEB01",
                       "Winkia",
                       "Capnocytophaga",
                       "Helicobacter_F",
                       "Laribacter",
                       "Methanobacterium",
                       "Moellerella",
                       "Mt11",
                       "Alloprevotella",
                       "Cupriavidus",
                       "Lactococcus_A",
                       "Rubeoparvulum",
                       "UMGS268",
                       "W0P29-029",
                       "Gracilibacillus",
                       "Nitratidesulfovibrio",
                       "Amedibacillus",
                       "Arachnia",
                       "UBA2804",
                       "MGYG000004443",
                       "Kocuria",
                       "Tidjanibacter",
                       "Turicimonas",
                       "UBA5809",
                       "UBA3375",
                       "HGM08974",
                       "Metamycoplasma",
                       "Nosocomiicoccus",
                       "UBA1103",
                       "Jonquetella",
                       "Lapidilactobacillus",
                       "HGM10611",
                       "MGYG000004606",
                       "Helicobacter_F",
                       "Exiguobacterium_A",
                       "Ureaplasma",
                       "ZJ304",
                       "Chryseobacterium",
                       "Methanomassiliicoccus",
                       "Lapidilactobacillus",
                       "Paenibacillus_F",
                       "Mycoplasmopsis",
                       "Aliarcobacter",
                       "HGM10890",
                       "Metamycoplasma",
                       "Bacillus_BD",
                       "Ureaplasma",
                       "Paenisporosarcina",
                       "Peribacillus",
                       "Planococcus",
                       "Methanocorpusculum",
                       "MGYG000004418",
                       "UBA7862",
                       "Bulleidia",
                       "Butyricicoccus_A",
                       "CAG-274",
                       "CAG-307",
                       "CAG-533",
                       "CAG-831",
                       "CAG-977",
                       "Clostridioides",
                       "Clostridium_A",
                       "DNF00809",
                       "Dysgonomonas",
                       "Edwardsiella",
                       "Enorma",
                       "Eubacterium_J",
                       "Flavonifractor",
                       "Gemella",
                       "HGM12650",
                       "MGYG000004380",
                       "Marseille-P4683",
                       "Mobilibacterium",
                       "NSJ-50",
                       "Niameybacter",
                       "Pseudobutyrivibrio",
                       "Pseudoflavonifractor",
                       "RQCD01",
                       "RUG806",
                       "SFEB01",
                       "SFEL01",
                       "UBA1221",
                       "UBA1234",
                       "UBA1409",
                       "UBA7102",
                       "UBA7862",
                       "UMGS1585",
                       "UMGS1840",
                       "UMGS403",
                       "UMGS687",
                       "Acetivibrio_A",
                       "CAG-345",
                       "CAG-41",
                       "CAG-568",
                       "CAG-617",
                       "Clostridium_AA",
                       "Clostridium_J",
                       "Comamonas",
                       "F082",
                       "HGM11530",
                       "Lachnoanaerobaculum",
                       "Paenibacillus",
                       "Paramuribaculum",
                       "Pauljensenia",
                       "Peptoniphilus_B",
                       "Phil1",
                       "SFFH01",
                       "UBA1685",
                       "UBA1724",
                       "UBA2821",
                       "UBA7597",
                       "UMGS1071",
                       "UMGS1279",
                       "UMGS416",
                       "Bacillus_A",
                       "Bacteroides_F",
                       "CAG-273",
                       "CAG-56",
                       "CAG-877",
                       "Faecalibacillus",
                       "Intestinimonas",
                       "Lactobacillus",
                       "Marseille-P3106",
                       "Pluralibacter",
                       "Pseudoruminococcus",
                       "RUG410",
                       "UBA1206",
                       "UMGS1872",
                       "CAG-594",
                       "Campylobacter",
                       "Campylobacter_D",
                       "Eggerthella",
                       "Enterococcus",
                       "Lachnoclostridium",
                       "NK4A136",
                       "SFMI01",
                       "RUG131",
                       "UMGS172",
                       "Barnesiella",
                       "Cedecea",
                       "Frisingicoccus",
                       "Serratia",
                       "UBA3375",
                       "CAG-1435",
                       "Companilactobacillus",
                       "Holdemania",
                       "UBA3402",
                       "Clostridium_P",
                       "WG-1",
                       "CAG-603",
                       "Helicobacter_C",
                       "CAG-411",
                       "UBA1740",
                       "UMGS1338",
                       "Anaeroplasma",
                       "Kluyvera",
                       "Lactococcus",
                       "Mycoplasmopsis",
                       "MGYG000003772",
                       "Mobilicoccus",
                       "Jonquetella",
                       "UMGS1260",
                       "28L")) %>% 
  group_by(cog_category, genus) %>% 
  mutate(pvcountry =kruskal.test(tally ~ country)$p.value)


#okay let's plot ----

site_color <- (c(wes_palette("Darjeeling1", 4, type = c("discrete")), wes_palette("GrandBudapest2", 2, type = c("discrete")),  wes_palette("Royal2", 1, type = c("discrete")), wes_palette("Darjeeling2", 3, type = c("discrete")), wes_palette("Zissou1", 5, type = c("discrete")), wes_palette("Chevalier1", 4, type = c("discrete")), wes_palette("FantasticFox1", 5, type = c("discrete")), wes_palette("Moonrise1", 3, type = c("discrete"))))

show_col(site_color)

#plot everyone ----

merge_list_percent[[1]] %>% 
  filter(id %in% c(1:50)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_1.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(51:101)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +   
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_2.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(102:152)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_3.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(153:203)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_4.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(204:254)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_5.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(305:355)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_6.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(406:456)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +   
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_7.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(507:557)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_8.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(608:658)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_9.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(659:709)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15))  + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_10.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(710:760)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_11.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(761:811)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_12.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(812:862)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_13.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(863:913)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_14.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(914:964)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_15.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[1]] %>% 
  filter(id %in% c(965:1032)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_all_16.pdf", width = 40, height = 20, units = "cm")

#okay now just high relative abundance ----

merge_list_percent[[8]] %>% 
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_highrelabun_1.pdf", width = 40, height = 20, units = "cm")

#how about those same genera at low abundance, in comparison with those with very high abundance ----

low_abud <- merge_list_percent[[9]] %>% 
  filter(genus %in% c("Acetatifactor",
                      "CAG-83", 
                      "CAG-103", 
                      "WG-1", 
                      "CAG-177", 
                      "CAG-110",
                      "Agathobacter",
                      "CAG-303", 
                      "Alistipes",
                      "Bacteroides",
                      "Butyrivibrio_A ",
                      "Bifidobacterium",
                      "Dialister",
                      "Escherichia",
                      "F082",
                      "Faecalibacterium",
                      "Klebsiella",
                      "Phocaeicola",
                      "Prevotella",
                      "RC9",
                      "RUG115",
                      "Roseburia",
                      "SFDB01",
                      "Succinivibrio",
                      "Treponema_D")) %>% 
  select(genus, cog_category, percent_by_genus, label) %>% 
  rename(low = percent_by_genus)

high_abun <- merge_list_percent[[8]] %>% 
  select(genus, cog_category, percent_by_genus, label) %>% 
  rename(high = percent_by_genus)

merge(low_abud, high_abun, by = c("genus", "cog_category", "label"), keep = TRUE) %>% 
  pivot_longer(names_to = "abun", values_to = "percent_by_genus", low:high) %>%
  group_by(abun) %>%
  ggplot(aes(x=factor(abun, levels = unique(abun)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_abuncomp.pdf", width = 40, height = 20, units = "cm")

#okay now just super low relative abundance ----

merge_list_percent[[10]] %>% 
  filter(id %in% c(1:50)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_1.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(51:101)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +   
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_2.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(102:152)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_3.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(153:203)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_4.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(204:254)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_5.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(305:355)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_6.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(406:456)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +   
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_7.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(507:557)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_8.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(608:658)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_9.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(659:709)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15))  + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_10.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(710:760)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_11.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(761:811)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_12.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(812:862)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_13.pdf", width = 40, height = 20, units = "cm")

merge_list_percent[[10]] %>% 
  filter(id %in% c(863:913)) %>%
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) +
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/function_genus/cog_function_superlow_14.pdf", width = 40, height = 20, units = "cm")


#now let's do country comparisons ----

cam <- merge_list_percent[[2]] %>% 
  select(genus, cog_category, percent_by_genus, label, id) %>% 
  mutate(country = "Cameroon")

bot <- merge_list_percent[[3]] %>% 
  select(genus, cog_category, percent_by_genus, label, id) %>% 
  mutate(country = "Botswana")

tan <- merge_list_percent[[4]] %>% 
  select(genus, cog_category, percent_by_genus, label, id) %>% 
  mutate(country = "Tanzania")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(1:20)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_1.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(21:41)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_2.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(42:62)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_3.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(63:83)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_4.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(84:104)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_5.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(105:125)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_6.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(126:146)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_7.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(147:167)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_8.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(168:188)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_9.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(189:219)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_10.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(220:240)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_11.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(241:261)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_12.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(262:282)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_13.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(283:303)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_14.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(304:324)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_15.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(325:345)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_16.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(346:366)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_17.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(367:387)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_18.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(388:408)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_19.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(409:429)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_20.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(430:450)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_21.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(451:471)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_22.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(472:492)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_23.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(493:513)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_24.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(514:534)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_25.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(515:535)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_26.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(536:556)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_27.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(557:577)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_28.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(578:598)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_29.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(599:619)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_30.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(620:640)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_31.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(641:661)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_32.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(662:682)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_33.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(683:703)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_34.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(704:724)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_35.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(725:745)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_36.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(746:766)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_37.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(767:787)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_38.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(788:808)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_39.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(809:829)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_40.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(830:850)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_41.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(851:871)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_42.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(872:892)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_43.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(893:913)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_44.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(914:934)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_45.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(935:955)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_46.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(956:976)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_47.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(id %in% c(977:997)) %>%
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_country_48.pdf", width = 40, height = 20, units = "cm")

rbind(cam, bot, tan) %>% 
  filter(genus %in% c("Alistipes",
                      "Bacteroides",
                      "CAG-791",
                      "Klebsiella",
                      "Lachnospira",
                      "Odoribacter",
                      "RC9",
                      "Oribacterium",
                      "Sodaliphilus",
                      "Succinivibrio",
                      "Sutterella",
                      "UBA4372",
                      "V9D3004",
                      "Eubacterium_J",
                      "Lachnospira",
                      "Slackia_A",
                      "Treponema_D",
                      "UBA1740",
                      "UBA2882",
                      "UBA3282",
                      "UBA4644",
                      "UMGS1484",
                      "CAG-45",
                      "Butyricicoccus_A",
                      "CAG-110",
                      "CAG-238")) %>% 
  ggplot(aes(x=factor(country, levels = unique(country)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90)) 

ggsave("output/function_genus/cog_function_country_tof.pdf", width = 40, height = 20, units = "cm")



#and now alpha diversity comparisons ---- 
high <- merge_list_percent[[11]] %>% 
  select(genus, cog_category, percent_by_genus, label, id) %>% 
  mutate(alpha = "high")

low <- merge_list_percent[[12]] %>% 
  select(genus, cog_category, percent_by_genus, label, id) %>% 
  mutate(alpha = "low")

rbind(low, high) %>% 
  filter(id %in% c(1:20)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_1.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(21:41)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_2.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(42:62)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_3.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(63:83)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_4.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(84:104)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_5.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(105:125)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_6.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(126:146)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_7.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(147:167)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_8.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(168:188)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_9.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(189:219)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_10.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(220:240)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_11.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(241:261)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_12.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(262:282)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_13.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(283:303)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_14.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(304:324)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_15.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(325:345)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_16.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(346:366)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_17.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(367:387)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_18.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(388:408)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_19.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(409:429)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_20.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(430:450)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_21.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(451:471)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_22.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(472:492)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_23.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(493:513)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_24.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(514:534)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_25.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(515:535)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_26.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(536:556)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_27.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(557:577)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_28.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(578:598)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_29.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(599:619)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_30.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(620:640)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_31.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(641:661)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_32.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(662:682)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_33.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(683:703)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_34.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(704:724)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_35.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(725:745)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_36.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(746:766)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_37.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(767:787)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_38.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(788:808)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_39.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(809:829)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_40.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(830:850)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_41.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(851:871)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_42.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(872:892)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_43.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(893:913)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_44.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(914:934)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_45.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(935:955)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_46.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(956:976)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_47.pdf", width = 40, height = 20, units = "cm")

rbind(low, high) %>% 
  filter(id %in% c(977:997)) %>%
  ggplot(aes(x=factor(alpha, levels = unique(alpha)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 10)) +    
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key.height = unit(1, "null")) + 
  facet_grid(~ genus, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) + 
  theme(strip.text.x = element_text(size = 8, angle = 90))

ggsave("output/function_genus/cog_function_alpha_48.pdf", width = 40, height = 20, units = "cm")


#finally, the genera that are highly DA ---- 

merge_list_percent[[13]] %>% 
  ggplot(aes(x=factor(genus, levels = unique(genus)), y=percent_by_genus, fill=label)) + 
  geom_col(position="stack") + 
  theme_minimal() + theme(legend.key.height = unit(1, "null")) + 
  labs(x="Genus",y="Protein Function", fill = "COG Category") + 
  scale_fill_manual(values=as.vector(site_color), labels = function(x) str_wrap(x, width = 25)) +
  theme(axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(size = 15)) + 
  guides(fill = guide_legend(keywidth = 0.5, label.theme = element_text(size = 15))) + 
  guides(fill = guide_legend(label.theme = element_text(size = 10))) + 
  theme(legend.key.height = unit(1, "null")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  geom_text(aes(label=round(percent_by_genus, digits = 3)), size = 4, position = position_stack(vjust = 0.5), check_overlap = TRUE) 

ggsave("output/function_genus/cog_function_all.pdf", width = 40, height = 20, units = "cm")

