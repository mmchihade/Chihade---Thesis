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


#load data ----

data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE) %>%
  select(-record_id, -topmed_id, -collection_name, -sample_site, -longitude, -latitude, -date, -year, -is_adult, -age, -waist_circumference, -body_fat_percent, -bmi) %>% 
  rename(merge_id = lab_id) %>%
  rename(lab_id = field_id) %>%
  mutate(lab_id = toupper(lab_id))

#create a dataframe with all of the cog tally files, with the file name as a column so that we can match column to tally ----

merged_df <- list.files(path = "data/cog_tally_overall",
                        pattern = "\\.csv$", 
                        full.names = TRUE) %>% 
  set_names() %>%
  map_dfr(read_csv, .id = "file_name") 

#we need the cog categories too! ----

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#let's modify our merged file to make it more  -----

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

#I want to know the total number of reads per cog category ----

cog_count <- full_tally %>%
  group_by(cog_category) %>% 
  summarise(count=sum(tally))

write_csv(cog_count, "data/cog_count.csv")

#now let's find the fraction of each taxa's proteome that has each function ----

w <- full_tally %>% 
  group_by(genus) %>% 
  mutate(sum_b = sum(tally, na.rm =TRUE, group_by=TRUE)) %>% 
  ungroup() %>%
  group_by(genus, cog_category) %>% 
  mutate(sum_bc = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  ungroup() %>%
  mutate(w = sum_bc/sum_b) %>% 
  select(lab_id, genus, cog_category, w)

#let's also calculate the relative abundance of each taxa based on these numbers----

a1 <- full_tally %>% 
  group_by(lab_id) %>% 
  summarise(sum_a = sum(tally, na.rm =TRUE, group_by=TRUE))

  
a2 <- full_tally %>% 
  group_by(lab_id, genus) %>% 
  summarise(sum_t = sum(tally, na.rm =TRUE, group_by=TRUE))

a <- merge(a1, a2, by="lab_id") %>% 
  mutate(a = sum_t/sum_a)

#finally, let's find the fraction of each individual's GM that has each function due to each taxa ----

f <- merge(a, w, by = c("lab_id", "genus")) %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  mutate(f = w*a/(sum(w*a)))
  

#umm, I am going to try to make a combined dataset, let's see if R can handle it ----

merge <- merge(f, data_grp, by = "lab_id")

merge_labels <- merge(merge, cog, by = "cog_category") 


#now I want to divide this big guy up into cog categories -----

split_df <- split(merge_labels, with(merge_labels, as.factor(cog_category)))

A <- split_df$A 
B <- split_df$B
C <- split_df$C
D <- split_df$D
E <- split_df$E
F <- split_df$F 
G <- split_df$G
H <- split_df$H
I <- split_df$I 
J <- split_df$J 
K <- split_df$K 
L <- split_df$L
M <- split_df$M 
N <- split_df$N
O <- split_df$O
P <- split_df$P
Q <- split_df$Q
S <- split_df$S
T <- split_df$T
U <- split_df$U 
V <- split_df$V 
W <- split_df$W
X <- split_df$X
Y <- split_df$Y
Z <- split_df$Z

#I would love to find p values for each genus between all the factors I will then average by ---- 

#To find the p values here, I need to make sure that the genera have tallys in each of the categories present.

#To make sure I have the best data possible, I am going to filter out the very low prevalence genera and treat them separately 

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

sex_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, sex, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(sex_distinct, "sex_delete_list.csv")

region_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, region, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(region_distinct, "region_delete_list.csv")

ethnicity_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, ethnicity, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(ethnicity_distinct, "ethnicity_delete_list.csv")

studygroup_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, study_group, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(studygroup_distinct, "studygroup_delete_list.csv")

countrygroup_distinct <- merge_labels %>% 
  group_by(genus, cog_category) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(n < 65) %>% 
  select(cog_category, genus, country_group, n) %>%  
  distinct() %>% 
  arrange(n, genus, cog_category) %>%
  group_by(genus, cog_category) %>% 
  mutate(num = n()) %>%
  filter(num == 1) %>% 
  arrange(cog_category)

write_csv(countrygroup_distinct, "countrygroup_delete_list.csv")

#all of these lists show specific genera that are associated with just one category in each grouping, so i know to filter them out for stats

#finding each p_value and sd is a little tricky though, because I need to cut out all of the genera that aren't present in more than one individual: 

#but here we go! 

#to find the pvalues, I need to create functions to use PERMANOVA for each of the groupings

permanova_country <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ country, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fcountry" = "F") %>% 
    rename("pvcountry" = "Pr(>F)")
}

permanova_sex <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ sex, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fsex" = "F") %>% 
    rename("pvsex" = "Pr(>F)")
}

permanova_countrygroup <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ country_group, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fcountrygroup" = "F") %>% 
    rename("pvcountrygroup" = "Pr(>F)")
}

permanova_studygroup <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ study_group, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fstudygroup" = "F") %>% 
    rename("pvstudygroup" = "Pr(>F)")
}

permanova_ethnicity <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ ethnicity, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fethnicity" = "F") %>% 
    rename("pvethnicity" = "Pr(>F)")
}

permanova_region <- function(i) { 
  cut <- i %>% 
    ungroup() %>%
    select(f) 
  
  matrix <- dist(cut, method = "euc") 
  
  values <- adonis2(matrix ~ region, i, method = "euc") 
  
  values %>% 
    na.omit() %>%
    select("F", "Pr(>F)") %>% 
    rename("Fregion" = "F") %>% 
    rename("pvregion" = "Pr(>F)")
}

A_country <- A %>% 
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
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
                       "UMGS1124",
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
                       "Muribaculum"))

split_A_country <- split(A_country, with(A_country, as.factor(genus)))

pv_A_country <- lapply(split_A_country, permanova_country)

pv_A_country <- bind_rows(pv_A_country, .id = "genus")

A_sex <- A %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acutalibacter",
                       "Angelakisella",
                       "CAG-353",
                       "Clostridium_A",
                       "Coprococcus",
                       "Eubacterium_I",
                       "HGM12650",
                       "HGM13222",
                       "Lactobacillus",
                       "Lactococcus_A",
                       "Mitsuokella",
                       "Neobittarella",
                       "Olsenella_E",
                       "Phyllobacterium",
                       "Plesiomonas",
                       "Pseudomonas",
                       "SFFH01",
                       "UBA7597",
                       "UMGS1071",
                       "UMGS1124",
                       "UMGS1241",
                       "UMGS1370",
                       "UMGS1649",
                       "UMGS172",
                       "UMGS973",
                       "Veillonella_A",
                       "CAG-1193",
                       "CAG-217",
                       "Eikenella",
                       "Enterobacter_D",
                       "UBA1740",
                       "UBA2883",
                       "UMGS1217",
                       "Citrobacter_B",
                       "Pseudomonas_E",
                       "RUG11247",
                       "UBA1067",
                       "CAG-1031",
                       "PeH17",
                       "UMGS874",
                       "Kosakonia",
                       "Tidjanibacter")) 

split_A_sex <- split(A_sex, with(A_sex, as.factor(genus)))

pv_A_sex <- lapply(split_A_sex, permanova_sex)

pv_A_sex <- bind_rows(pv_A_sex, .id = "genus")

A_studygroup <- A %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Angelakisella",
                       "CAG-632",
                       "CAG-822",
                       "Clostridium",
                       "Clostridium_E",
                       "Coprococcus",
                       "Enterococcus",
                       "Eubacterium_I",
                       "Firm-07",
                       "Flavonifractor",
                       "HGM12650",
                       "Intestinibacter",
                       "Intestinimonas",
                       "Neobittarella",
                       "Olsenella_E",
                       "Plesiomonas",
                       "UBA7597",
                       "UMGS1124",
                       "UMGS1370",
                       "UMGS1649",
                       "UMGS973",
                       "Veillonella_A",
                       "WG-1",
                       "CAG-877",
                       "F082",
                       "Ruminococcus_A",
                       "UBA1740",
                       "Acinetobacter",
                       "UBA7160")) 

split_A_studygroup <- split(A_studygroup, with(A_studygroup, as.factor(genus)))

pv_A_studygroup <- lapply(split_A_studygroup, permanova_studygroup)

pv_A_studygroup <- bind_rows(pv_A_studygroup, .id = "genus")

A_region <- A %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Aeromonas",
                       "CAG-590",
                       "Clostridium",
                       "Clostridium_E",
                       "HGM12650",
                       "Intestinibacter",
                       "Lactococcus_A",
                       "Mitsuokella",
                       "Plesiomonas",
                       "UMGS1071",
                       "UMGS1124",
                       "UMGS1370",
                       "UMGS973",
                       "WG-1",
                       "Bacteroides_F",
                       "F082",
                       "CAG-81",
                       "Tidjanibacter",
                       "Dialister",
                       "Amulumruptor")) 

split_A_region <- split(A_region, with(A_region, as.factor(genus)))

pv_A_region <- lapply(split_A_region, permanova_region)

pv_A_region <- bind_rows(pv_A_region, .id = "genus")

A_ethnicity <- A %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Clostridium",
                       "Clostridium_E",
                       "Intestinibacter",
                       "Plesiomonas",
                       "UMGS1124",
                       "UMGS973",
                       "WG-1",
                       "F082",
                       "Acinetobacter")) 

split_A_ethnicity <- split(A_ethnicity, with(A_ethnicity, as.factor(genus)))

pv_A_ethnicity <- lapply(split_A_ethnicity, permanova_ethnicity)

pv_A_ethnicity <- bind_rows(pv_A_ethnicity, .id = "genus")

A_countrygroup <- A %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("CAG-632",
                       "Clostridium",
                       "Clostridium_E",
                       "Enterococcus",
                       "HGM12650",
                       "Intestinibacter",
                       "Intestinimonas",
                       "Plesiomonas",
                       "UBA7597",
                       "UMGS1124",
                       "UMGS1370",
                       "UMGS973",
                       "WG-1",
                       "CAG-877",
                       "F082",
                       "UBA1740")) 

split_A_countrygroup <- split(A_countrygroup, with(A_countrygroup, as.factor(genus)))

pv_A_countrygroup <- lapply(split_A_countrygroup, permanova_countrygroup)

pv_A_countrygroup <- bind_rows(pv_A_countrygroup, .id = "genus")


sd_A <- A %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "A", 
                   sd=sd(f))

pv_A <- merge(pv_A_country, pv_A_sex, by = c("genus"), all=TRUE)

pv_A <- merge(pv_A, pv_A_countrygroup, by = c("genus"), all=TRUE)

pv_A <- merge(pv_A, pv_A_ethnicity, by = c("genus"), all=TRUE) 

pv_A <- merge(pv_A, pv_A_region, by = c("genus"), all=TRUE)

pv_A <- merge(pv_A, pv_A_studygroup, by = c("genus"), all=TRUE)

stats_A <- merge(pv_A, sd_A, by = c("genus"), all=TRUE) %>% 
  drop_na(sd)

B_country <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("51-20",
                       "Alistipes_A",
                       "Atlantibacter",
                       "Blastococcus",
                       "CACZQA01",
                       "CAG-632",
                       "Citrobacter",
                       "Enterococcus_D",
                       "Finegoldia",
                       "Kluyvera",
                       "Lactobacillus",
                       "MGYG000004461",
                       "Negativicoccus",
                       "Phytobacter",
                       "Pseudomonas",
                       "RUG115",
                       "RUG420",
                       "Rubneribacter",
                       "SFTJ01",
                       "Salmonella",
                       "UBA11490",
                       "UBA1740",
                       "UBA5905",
                       "UBA6984",
                       "UCG-010",
                       "UMGS1202",
                       "WG-1",
                       "Afipia",
                       "CAAEEV01",
                       "Campylobacter",
                       "Clostridioides",
                       "GCA-900066755",
                       "MGYG000000574",
                       "OF09-33XD",
                       "Pantoea",
                       "Pseudomonas_E",
                       "RF16",
                       "SFHK01",
                       "UBA1259",
                       "UBA1691",
                       "UBA4334",
                       "Bacteroides_F",
                       "CAG-411",
                       "Clostridium_AP",
                       "Dorea_B",
                       "Enterococcus_A",
                       "Harryflintia",
                       "Hungatella_A",
                       "JAAYBB01",
                       "MGYG000004850",
                       "Olegusella",
                       "CAG-196",
                       "CAG-484",
                       "Clostridium_J",
                       "Neofamilia",
                       "Turicimonas",
                       "UMGS2037",
                       "Caproiciproducens",
                       "Raoultella",
                       "MGYG000004171",
                       "Ruminococcus_B",
                       "UBA1712",
                       "Cetobacterium_A",
                       "Muribaculum",
                       "Extibacter",
                       "HGM12814",
                       "Paraprevotella",
                       "CAG-45"))

split_B_country <- split(B_country, with(B_country, as.factor(genus)))

pv_B_country <- lapply(split_B_country, permanova_country)

pv_B_country <- bind_rows(pv_B_country, .id = "genus")

B_sex <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("51-20",
                       "Alistipes_A",
                       "Atlantibacter",
                       "CACZQA01",
                       "CAG-1031",
                       "CAG-349",
                       "Enterococcus_D",
                       "Finegoldia",
                       "HGM12998",
                       "Lactobacillus",
                       "MGYG000004461",
                       "Phytobacter",
                       "RUG115",
                       "RUG420",
                       "Rubneribacter",
                       "SFTJ01",
                       "Salmonella",
                       "UBA5905",
                       "UMGS1202",
                       "WG-1",
                       "MGYG000000574",
                       "Pantoea",
                       "Pseudomonas_E",
                       "RF16",
                       "Varibaculum",
                       "Bacteroides_F",
                       "MGYG000000456",
                       "QFNR01",
                       "Bittarella",
                       "Slackia",
                       "UMGS2037",
                       "CAG-460",
                       "MGYG000002758",
                       "UMGS1783",
                       "Fusobacterium_C",
                       "Porphyromonas"))

split_B_sex <- split(B_sex, with(B_sex, as.factor(genus)))

pv_B_sex <- lapply(split_B_sex, permanova_sex)

pv_B_sex <- bind_rows(pv_B_sex, .id = "genus")

B_studygroup <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alistipes_A",
                       "Atlantibacter",
                       "CACZQA01",
                       "CAG-1031",
                       "Christensenella",
                       "Citrobacter",
                       "Enterococcus_D",
                       "GCA-900066495",
                       "HGM11507",
                       "HGM12998",
                       "Phytobacter",
                       "Pseudomonas",
                       "QAKW01",
                       "RUG11130",
                       "Rubneribacter",
                       "SFTJ01",
                       "UBA1740",
                       "UBA5905",
                       "UBA6984",
                       "CAG-475",
                       "GCA-900066755",
                       "MGYG000000574",
                       "UMGS856",
                       "Raoultibacter",
                       "CAG-877"))

split_B_studygroup <- split(B_studygroup, with(B_studygroup, as.factor(genus)))

pv_B_studygroup <- lapply(split_B_studygroup, permanova_studygroup)

pv_B_studygroup <- bind_rows(pv_B_studygroup, .id = "genus")

B_region <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alistipes_A",
                       "Atlantibacter",
                       "CACZQA01",
                       "Pseudomonas",
                       "RUG420",
                       "SFTJ01",
                       "UBA5905",
                       "UCG-010",
                       "Pseudomonas_E",
                       "UBA1691",
                       "Bacteroides_F",
                       "Clostridium_AP",
                       "Enterococcus_A",
                       "Harryflintia",
                       "MGYG000004850",
                       "CAG-196",
                       "Raoultella",
                       "MGYG000004171")) 

split_B_region <- split(B_region, with(B_region, as.factor(genus)))

pv_B_region <- lapply(split_B_region, permanova_region)

pv_B_region <- bind_rows(pv_B_region, .id = "genus")

B_ethnicity <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alistipes_A",
                       "Atlantibacter",
                       "CACZQA01",
                       "Pseudomonas",
                       "SFTJ01",
                       "UBA5905"))


split_B_ethnicity <- split(B_ethnicity, with(B_ethnicity, as.factor(genus)))

pv_B_ethnicity <- lapply(split_B_ethnicity, permanova_ethnicity)

pv_B_ethnicity <- bind_rows(pv_B_ethnicity, .id = "genus")

B_countrygroup <- B %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alistipes_A",
                       "Atlantibacter",
                       "CACZQA01",
                       "Citrobacter",
                       "Enterococcus_D",
                       "Phytobacter",
                       "Pseudomonas",
                       "Rubneribacter",
                       "SFTJ01",
                       "UBA1740",
                       "UBA5905",
                       "UBA6984",
                       "GCA-900066755",
                       "MGYG000000574"))

split_B_countrygroup <- split(B_countrygroup, with(B_countrygroup, as.factor(genus)))

pv_B_countrygroup <- lapply(split_B_countrygroup, permanova_countrygroup)

pv_B_countrygroup <- bind_rows(pv_B_countrygroup, .id = "genus")

sd_B <- B %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "B", 
                   sd=sd(f))

pv_B <- merge(pv_B_country, pv_B_sex, by = c("genus"), all=TRUE)

pv_B <- merge(pv_B , pv_B_countrygroup, by = c("genus"), all=TRUE)

pv_B <- merge(pv_B , pv_B_ethnicity, by = c("genus"), all=TRUE) 

pv_B <- merge(pv_B , pv_B_region, by = c("genus"), all=TRUE)

pv_B <- merge(pv_B , pv_B_studygroup, by = c("genus"), all=TRUE)

stats_B <- merge(pv_B , sd_B , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


C_country <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418",
                       "Methanobacterium"))

split_C_country <- split(C_country, with(C_country, as.factor(genus)))

pv_C_country <- lapply(split_C_country, permanova_country)

pv_C_country <- bind_rows(pv_C_country, .id = "genus")

C_sex <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma"))

split_C_sex <- split(C_sex, with(C_sex, as.factor(genus)))

pv_C_sex <- lapply(split_C_sex, permanova_sex)

pv_C_sex <- bind_rows(pv_C_sex, .id = "genus")

C_studygroup <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_C_studygroup <- split(C_studygroup, with(C_studygroup, as.factor(genus)))

pv_C_studygroup <- lapply(split_C_studygroup, permanova_studygroup)

pv_C_studygroup <- bind_rows(pv_C_studygroup, .id = "genus")

C_region <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_C_region <- split(C_region, with(C_region, as.factor(genus)))

pv_C_region <- lapply(split_C_region, permanova_region)

pv_C_region <- bind_rows(pv_C_region, .id = "genus")

C_ethnicity <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_C_ethnicity <- split(C_ethnicity, with(C_ethnicity, as.factor(genus)))

pv_C_ethnicity <- lapply(split_C_ethnicity, permanova_ethnicity)

pv_C_ethnicity <- bind_rows(pv_C_ethnicity, .id = "genus")

C_countrygroup <- C %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_C_countrygroup <- split(C_countrygroup, with(C_countrygroup, as.factor(genus)))

pv_C_countrygroup <- lapply(split_C_countrygroup, permanova_countrygroup)

pv_C_countrygroup <- bind_rows(pv_C_countrygroup, .id = "genus")

sd_C <- C %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "C", 
                   sd=sd(f))

pv_C <- merge(pv_C_country, pv_C_sex, by = c("genus"), all=TRUE)

pv_C <- merge(pv_C , pv_C_countrygroup, by = c("genus"), all=TRUE)

pv_C <- merge(pv_C , pv_C_ethnicity, by = c("genus"), all=TRUE) 

pv_C <- merge(pv_C , pv_C_region, by = c("genus"), all=TRUE)

pv_C <- merge(pv_C , pv_C_studygroup, by = c("genus"), all=TRUE)

stats_C <- merge(pv_C , sd_C , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


D_country <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bacillus_O",
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
                       "Lactococcus_A",
                       "HGM05376",
                       "Enterovibrio",
                       "Parascardovia"))

split_D_country <- split(D_country, with(D_country, as.factor(genus)))

pv_D_country <- lapply(split_D_country, permanova_country)

pv_D_country <- bind_rows(pv_D_country, .id = "genus")


D_sex <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Dermabacter",
                       "HGM08974",
                       "Haloferax",
                       "Nesterenkonia",
                       "Risungbinella",
                       "Garciella"))

split_D_sex <- split(D_sex, with(D_sex, as.factor(genus)))

pv_D_sex <- lapply(split_D_sex, permanova_sex)

pv_D_sex <- bind_rows(pv_D_sex, .id = "genus")

D_studygroup <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
                       "Bacillus_O",
                       "Risungbinella"))

split_D_studygroup <- split(D_studygroup, with(D_studygroup, as.factor(genus)))

pv_D_studygroup <- lapply(split_D_studygroup, permanova_studygroup)

pv_D_studygroup <- bind_rows(pv_D_studygroup, .id = "genus")

D_region <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_D_region <- split(D_region, with(D_region, as.factor(genus)))

pv_D_region <- lapply(split_D_region, permanova_region)

pv_D_region <- bind_rows(pv_D_region, .id = "genus")

D_ethnicity <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_D_ethnicity <- split(D_ethnicity, with(D_ethnicity, as.factor(genus)))

pv_D_ethnicity <- lapply(split_D_ethnicity, permanova_ethnicity)

pv_D_ethnicity <- bind_rows(pv_D_ethnicity, .id = "genus")

D_countrygroup <- D %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bacillus_O")) 

split_D_countrygroup <- split(D_countrygroup, with(D_countrygroup, as.factor(genus)))

pv_D_countrygroup <- lapply(split_D_countrygroup, permanova_countrygroup)

pv_D_countrygroup <- bind_rows(pv_D_countrygroup, .id = "genus")


sd_D <- D %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "D", 
                   sd=sd(f))

pv_D <- merge(pv_D_country, pv_D_sex, by = c("genus"), all=TRUE)

pv_D <- merge(pv_D , pv_D_countrygroup, by = c("genus"), all=TRUE)

pv_D <- merge(pv_D , pv_D_ethnicity, by = c("genus"), all=TRUE) 

pv_D <- merge(pv_D , pv_D_region, by = c("genus"), all=TRUE)

pv_D <- merge(pv_D , pv_D_studygroup, by = c("genus"), all=TRUE)

stats_D <- merge(pv_D , sd_D , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



E_country <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418",
                       "UBA7488")) 

split_E_country <- split(E_country, with(E_country, as.factor(genus)))

pv_E_country <- lapply(split_E_country, permanova_country)

pv_E_country <- bind_rows(pv_E_country, .id = "genus")

E_sex <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma"))

split_E_sex <- split(E_sex, with(E_sex, as.factor(genus)))

pv_E_sex <- lapply(split_E_sex, permanova_sex)

pv_E_sex <- bind_rows(pv_E_sex, .id = "genus")

E_studygroup <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418")) 

split_E_studygroup <- split(E_studygroup, with(E_studygroup, as.factor(genus)))

pv_E_studygroup <- lapply(split_E_studygroup, permanova_studygroup)

pv_E_studygroup <- bind_rows(pv_E_studygroup, .id = "genus")

E_region <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_E_region <- split(E_region, with(E_region, as.factor(genus)))

pv_E_region <- lapply(split_E_region, permanova_region)

pv_E_region <- bind_rows(pv_E_region, .id = "genus")

E_ethnicity <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418")) 

split_E_ethnicity <- split(E_ethnicity, with(E_ethnicity, as.factor(genus)))

pv_E_ethnicity <- lapply(split_E_ethnicity, permanova_ethnicity)

pv_E_ethnicity <- bind_rows(pv_E_ethnicity, .id = "genus")

E_countrygroup <- E %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418")) 

split_E_countrygroup <- split(E_countrygroup, with(E_countrygroup, as.factor(genus)))

pv_E_countrygroup <- lapply(split_E_countrygroup, permanova_countrygroup)

pv_E_countrygroup <- bind_rows(pv_E_countrygroup, .id = "genus")


sd_E <- E %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "E", 
                   sd=sd(f))

pv_E <- merge(pv_E_country, pv_E_sex, by = c("genus"), all=TRUE)

pv_E <- merge(pv_E , pv_E_countrygroup, by = c("genus"), all=TRUE)

pv_E <- merge(pv_E , pv_E_ethnicity, by = c("genus"), all=TRUE) 

pv_E <- merge(pv_E , pv_E_region, by = c("genus"), all=TRUE)

pv_E <- merge(pv_E , pv_E_studygroup, by = c("genus"), all=TRUE)

stats_E <- merge(pv_E , sd_E , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



F_country <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("HGM08974",
                       "UBA7862",
                       "Thalassobacillus",
                       "Aliarcobacter",
                       "MGYG000000407",
                       "Flaviflexus",
                       "Anaerosalibacter",
                       "UBA3375")) 

split_F_country <- split(F_country, with(F_country, as.factor(genus)))

pv_F_country <- lapply(split_F_country, permanova_country)

pv_F_country <- bind_rows(pv_F_country, .id = "genus")

F_sex <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanomassiliicoccus",
                       "Metamycoplasma")) 

split_F_sex <- split(F_sex, with(F_sex, as.factor(genus)))

pv_F_sex <- lapply(split_F_sex, permanova_sex)

pv_F_sex <- bind_rows(pv_F_sex, .id = "genus")

F_studygroup <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("UBA2834"))

split_F_studygroup <- split(F_studygroup, with(F_studygroup, as.factor(genus)))

pv_F_studygroup <- lapply(split_F_studygroup, permanova_studygroup)

pv_F_studygroup <- bind_rows(pv_F_studygroup, .id = "genus")

F_region <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_F_region <- split(F_region, with(F_region, as.factor(genus)))

pv_F_region <- lapply(split_F_region, permanova_region)

pv_F_region <- bind_rows(pv_F_region, .id = "genus")

F_ethnicity <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_F_ethnicity <- split(F_ethnicity, with(F_ethnicity, as.factor(genus)))

pv_F_ethnicity <- lapply(split_F_ethnicity, permanova_ethnicity)

pv_F_ethnicity <- bind_rows(pv_F_ethnicity, .id = "genus")

F_countrygroup <- F %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_F_countrygroup <- split(F_countrygroup, with(F_countrygroup, as.factor(genus)))

pv_F_countrygroup <- lapply(split_F_countrygroup, permanova_countrygroup)

pv_F_countrygroup <- bind_rows(pv_F_countrygroup, .id = "genus")


sd_F <- F %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "F", 
                   sd=sd(f))

pv_F <- merge(pv_F_country, pv_F_sex, by = c("genus"), all=TRUE)

pv_F <- merge(pv_F , pv_F_countrygroup, by = c("genus"), all=TRUE)

pv_F <- merge(pv_F , pv_F_ethnicity, by = c("genus"), all=TRUE) 

pv_F <- merge(pv_F , pv_F_region, by = c("genus"), all=TRUE)

pv_F <- merge(pv_F , pv_F_studygroup, by = c("genus"), all=TRUE)

stats_F <- merge(pv_F , sd_F , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



G_country <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418",
                       "Methanomassiliicoccus",
                       "Jonquetella")) 

split_G_country <- split(G_country, with(G_country, as.factor(genus)))

pv_G_country <- lapply(split_G_country, permanova_country)

pv_G_country <- bind_rows(pv_G_country, .id = "genus")

G_sex <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_G_sex <- split(G_sex, with(G_sex, as.factor(genus)))

pv_G_sex <- lapply(split_G_sex, permanova_sex)

pv_G_sex <- bind_rows(pv_G_sex, .id = "genus")

G_studygroup <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_G_studygroup <- split(G_studygroup, with(G_studygroup, as.factor(genus)))

pv_G_studygroup <- lapply(split_G_studygroup, permanova_studygroup)

pv_G_studygroup <- bind_rows(pv_G_studygroup, .id = "genus")

G_region <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_G_region <- split(G_region, with(G_region, as.factor(genus)))

pv_G_region <- lapply(split_G_region, permanova_region)

pv_G_region <- bind_rows(pv_G_region, .id = "genus")

G_ethnicity <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_G_ethnicity <- split(G_ethnicity, with(G_ethnicity, as.factor(genus)))

pv_G_ethnicity <- lapply(split_G_ethnicity, permanova_ethnicity)

pv_G_ethnicity <- bind_rows(pv_G_ethnicity, .id = "genus")

G_countrygroup <- G %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_G_countrygroup <- split(G_countrygroup, with(G_countrygroup, as.factor(genus)))

pv_G_countrygroup <- lapply(split_G_countrygroup, permanova_countrygroup)

pv_G_countrygroup <- bind_rows(pv_G_countrygroup, .id = "genus") 


sd_G <- G %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "G", 
                   sd=sd(f))

pv_G <- merge(pv_G_country, pv_G_sex, by = c("genus"), all=TRUE)

pv_G <- merge(pv_G , pv_G_countrygroup, by = c("genus"), all=TRUE)

pv_G <- merge(pv_G , pv_G_ethnicity, by = c("genus"), all=TRUE) 

pv_G <- merge(pv_G , pv_G_region, by = c("genus"), all=TRUE)

pv_G <- merge(pv_G , pv_G_studygroup, by = c("genus"), all=TRUE)

stats_G <- merge(pv_G , sd_G , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



H_country <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Epilithonimonas",
                       "Fructilactobacillus",
                       "UBA4877",
                       "Mycoplasmopsis",
                       "Niallia")) 

split_H_country <- split(H_country, with(H_country, as.factor(genus)))

pv_H_country <- lapply(split_H_country, permanova_country)

pv_H_country <- bind_rows(pv_H_country, .id = "genus")

H_sex <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Epilithonimonas")) 

split_H_sex <- split(H_sex, with(H_sex, as.factor(genus)))

pv_H_sex <- lapply(split_H_sex, permanova_sex)

pv_H_sex <- bind_rows(pv_H_sex, .id = "genus")

H_studygroup <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_H_studygroup <- split(H_studygroup, with(H_studygroup, as.factor(genus)))

pv_H_studygroup <- lapply(split_H_studygroup, permanova_studygroup)

pv_H_studygroup <- bind_rows(pv_H_studygroup, .id = "genus")

H_region <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_H_region <- split(H_region, with(H_region, as.factor(genus)))

pv_H_region <- lapply(split_H_region, permanova_region)

pv_H_region <- bind_rows(pv_H_region, .id = "genus")

H_ethnicity <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_H_ethnicity <- split(H_ethnicity, with(H_ethnicity, as.factor(genus)))

pv_H_ethnicity <- lapply(split_H_ethnicity, permanova_ethnicity)

pv_H_ethnicity <- bind_rows(pv_H_ethnicity, .id = "genus")

H_countrygroup <- H %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_H_countrygroup <- split(H_countrygroup, with(H_countrygroup, as.factor(genus)))

pv_H_countrygroup <- lapply(split_H_countrygroup, permanova_countrygroup)

pv_H_countrygroup <- bind_rows(pv_H_countrygroup, .id = "genus") 



sd_H <- H %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "H", 
                   sd=sd(f))

pv_H <- merge(pv_H_country, pv_H_sex, by = c("genus"), all=TRUE)

pv_H <- merge(pv_H , pv_H_countrygroup, by = c("genus"), all=TRUE)

pv_H <- merge(pv_H , pv_H_ethnicity, by = c("genus"), all=TRUE) 

pv_H <- merge(pv_H , pv_H_region, by = c("genus"), all=TRUE)

pv_H <- merge(pv_H , pv_H_studygroup, by = c("genus"), all=TRUE)

stats_H <- merge(pv_H , sd_H , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



I_country <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Exiguobacterium_A",
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
                       "Psychrobacillus")) 

split_I_country <- split(I_country, with(I_country, as.factor(genus)))

pv_I_country <- lapply(split_I_country, permanova_country)

pv_I_country <- bind_rows(pv_I_country, .id = "genus")


I_sex <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("HGM10611",
                       "Nosocomiicoccus",
                       "Gleimia",
                       "MGYG000004144",
                       "Metamycoplasma",
                       "HGM10890",
                       "Actinobaculum",
                       "Dermabacter",
                       "Epilithonimonas",
                       "Haloferax"))

split_I_sex <- split(I_sex, with(I_sex, as.factor(genus)))

pv_I_sex <- lapply(split_I_sex, permanova_sex)

pv_I_sex <- bind_rows(pv_I_sex, .id = "genus")

I_studygroup <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Exiguobacterium_A",
                       "Nosocomiicoccus",
                       "Bacillus_O",
                       "Gleimia",
                       "Numidum",
                       "HGM10890",
                       "Pelistega",
                       "UBA2834",
                       "Halorubrum")) 

split_I_studygroup <- split(I_studygroup, with(I_studygroup, as.factor(genus)))

pv_I_studygroup <- lapply(split_I_studygroup, permanova_studygroup)

pv_I_studygroup <- bind_rows(pv_I_studygroup, .id = "genus")

I_region <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Nosocomiicoccus")) 

split_I_region <- split(I_region, with(I_region, as.factor(genus)))

pv_I_region <- lapply(split_I_region, permanova_region)

pv_I_region <- bind_rows(pv_I_region, .id = "genus")

I_ethnicity <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Gleimia")) 

split_I_ethnicity <- split(I_ethnicity, with(I_ethnicity, as.factor(genus)))

pv_I_ethnicity <- lapply(split_I_ethnicity, permanova_ethnicity)

pv_I_ethnicity <- bind_rows(pv_I_ethnicity, .id = "genus")

I_countrygroup <- I %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Exiguobacterium_A",
                       "Nosocomiicoccus",
                       "Numidum")) 

split_I_countrygroup <- split(I_countrygroup, with(I_countrygroup, as.factor(genus)))

pv_I_countrygroup <- lapply(split_I_countrygroup, permanova_countrygroup)

pv_I_countrygroup <- bind_rows(pv_I_countrygroup, .id = "genus")

sd_I <- I %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "I", 
                   sd=sd(f))

pv_I <- merge(pv_I_country, pv_I_sex, by = c("genus"), all=TRUE)

pv_I <- merge(pv_I , pv_I_countrygroup, by = c("genus"), all=TRUE)

pv_I <- merge(pv_I , pv_I_ethnicity, by = c("genus"), all=TRUE) 

pv_I <- merge(pv_I , pv_I_region, by = c("genus"), all=TRUE)

pv_I <- merge(pv_I , pv_I_studygroup, by = c("genus"), all=TRUE)

stats_I <- merge(pv_I , sd_I , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


J_country <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_J_country <- split(J_country, with(J_country, as.factor(genus)))

pv_J_country <- lapply(split_J_country, permanova_country)

pv_J_country <- bind_rows(pv_J_country, .id = "genus")

J_sex <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_J_sex <- split(J_sex, with(J_sex, as.factor(genus)))

pv_J_sex <- lapply(split_J_sex, permanova_sex)

pv_J_sex <- bind_rows(pv_J_sex, .id = "genus")

J_studygroup <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_J_studygroup <- split(J_studygroup, with(J_studygroup, as.factor(genus)))

pv_J_studygroup <- lapply(split_J_studygroup, permanova_studygroup)

pv_J_studygroup <- bind_rows(pv_J_studygroup, .id = "genus")

J_region <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_J_region <- split(J_region, with(J_region, as.factor(genus)))

pv_J_region <- lapply(split_J_region, permanova_region)

pv_J_region <- bind_rows(pv_J_region, .id = "genus")

J_ethnicity <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_J_ethnicity <- split(J_ethnicity, with(J_ethnicity, as.factor(genus)))

pv_J_ethnicity <- lapply(split_J_ethnicity, permanova_ethnicity)

pv_J_ethnicity <- bind_rows(pv_J_ethnicity, .id = "genus")

J_countrygroup <- J %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_J_countrygroup <- split(J_countrygroup, with(J_countrygroup, as.factor(genus)))

pv_J_countrygroup <- lapply(split_J_countrygroup, permanova_countrygroup)

pv_J_countrygroup <- bind_rows(pv_J_countrygroup, .id = "genus") 




sd_J <- J %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "J", 
                   sd=sd(f))

pv_J <- merge(pv_J_country, pv_J_sex, by = c("genus"), all=TRUE)

pv_J <- merge(pv_J , pv_J_countrygroup, by = c("genus"), all=TRUE)

pv_J <- merge(pv_J , pv_J_ethnicity, by = c("genus"), all=TRUE) 

pv_J <- merge(pv_J , pv_J_region, by = c("genus"), all=TRUE)

pv_J <- merge(pv_J , pv_J_studygroup, by = c("genus"), all=TRUE)

stats_J <- merge(pv_J , sd_J , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



K_country <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanobacterium",
                       "Methanocorpusculum",
                       "Acetomicrobium",
                       "UBA1103",
                       "Helicobacter_F")) 

split_K_country <- split(K_country, with(K_country, as.factor(genus)))

pv_K_country <- lapply(split_K_country, permanova_country)

pv_K_country <- bind_rows(pv_K_country, .id = "genus")

K_sex <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_K_sex <- split(K_sex, with(K_sex, as.factor(genus)))

pv_K_sex <- lapply(split_K_sex, permanova_sex)

pv_K_sex <- bind_rows(pv_K_sex, .id = "genus")

K_studygroup <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_K_studygroup <- split(K_studygroup, with(K_studygroup, as.factor(genus)))

pv_K_studygroup <- lapply(split_K_studygroup, permanova_studygroup)

pv_K_studygroup <- bind_rows(pv_K_studygroup, .id = "genus")

K_region <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_K_region <- split(K_region, with(K_region, as.factor(genus)))

pv_K_region <- lapply(split_K_region, permanova_region)

pv_K_region <- bind_rows(pv_K_region, .id = "genus")

K_ethnicity <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_K_ethnicity <- split(K_ethnicity, with(K_ethnicity, as.factor(genus)))

pv_K_ethnicity <- lapply(split_K_ethnicity, permanova_ethnicity)

pv_K_ethnicity <- bind_rows(pv_K_ethnicity, .id = "genus")

K_countrygroup <- K %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_K_countrygroup <- split(K_countrygroup, with(K_countrygroup, as.factor(genus)))

pv_K_countrygroup <- lapply(split_K_countrygroup, permanova_countrygroup)

pv_K_countrygroup <- bind_rows(pv_K_countrygroup, .id = "genus") 



sd_K <- K %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "K", 
                   sd=sd(f))

pv_K <- merge(pv_K_country, pv_K_sex, by = c("genus"), all=TRUE)

pv_K <- merge(pv_K , pv_K_countrygroup, by = c("genus"), all=TRUE)

pv_K <- merge(pv_K , pv_K_ethnicity, by = c("genus"), all=TRUE) 

pv_K <- merge(pv_K , pv_K_region, by = c("genus"), all=TRUE)

pv_K <- merge(pv_K , pv_K_studygroup, by = c("genus"), all=TRUE)

stats_K <- merge(pv_K , sd_K , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


L_country <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418")) 

split_L_country <- split(L_country, with(L_country, as.factor(genus)))

pv_L_country <- lapply(split_L_country, permanova_country)

pv_L_country <- bind_rows(pv_L_country, .id = "genus")

L_sex <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_L_sex <- split(L_sex, with(L_sex, as.factor(genus)))

pv_L_sex <- lapply(split_L_sex, permanova_sex)

pv_L_sex <- bind_rows(pv_L_sex, .id = "genus")

L_studygroup <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_L_studygroup <- split(L_studygroup, with(L_studygroup, as.factor(genus)))

pv_L_studygroup <- lapply(split_L_studygroup, permanova_studygroup)

pv_L_studygroup <- bind_rows(pv_L_studygroup, .id = "genus")

L_region <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_L_region <- split(L_region, with(L_region, as.factor(genus)))

pv_L_region <- lapply(split_L_region, permanova_region)

pv_L_region <- bind_rows(pv_L_region, .id = "genus")

L_ethnicity <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_L_ethnicity <- split(L_ethnicity, with(L_ethnicity, as.factor(genus)))

pv_L_ethnicity <- lapply(split_L_ethnicity, permanova_ethnicity)

pv_L_ethnicity <- bind_rows(pv_L_ethnicity, .id = "genus")

L_countrygroup <- L %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_L_countrygroup <- split(L_countrygroup, with(L_countrygroup, as.factor(genus)))

pv_L_countrygroup <- lapply(split_L_countrygroup, permanova_countrygroup)

pv_L_countrygroup <- bind_rows(pv_L_countrygroup, .id = "genus") 




sd_L <- L %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "L", 
                   sd=sd(f))

pv_L <- merge(pv_L_country, pv_L_sex, by = c("genus"), all=TRUE)

pv_L <- merge(pv_L , pv_L_countrygroup, by = c("genus"), all=TRUE)

pv_L <- merge(pv_L , pv_L_ethnicity, by = c("genus"), all=TRUE) 

pv_L <- merge(pv_L , pv_L_region, by = c("genus"), all=TRUE)

pv_L <- merge(pv_L , pv_L_studygroup, by = c("genus"), all=TRUE)

stats_L <- merge(pv_L , sd_L , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


M_country <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Rubeoparvulum",
                       "MGYG000003562")) 

split_M_country <- split(M_country, with(M_country, as.factor(genus)))

pv_M_country <- lapply(split_M_country, permanova_country)

pv_M_country <- bind_rows(pv_M_country, .id = "genus")

M_sex <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_M_sex <- split(M_sex, with(M_sex, as.factor(genus)))

pv_M_sex <- lapply(split_M_sex, permanova_sex)

pv_M_sex <- bind_rows(pv_M_sex, .id = "genus")

M_studygroup <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_M_studygroup <- split(M_studygroup, with(M_studygroup, as.factor(genus)))

pv_M_studygroup <- lapply(split_M_studygroup, permanova_studygroup)

pv_M_studygroup <- bind_rows(pv_M_studygroup, .id = "genus")

M_region <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_M_region <- split(M_region, with(M_region, as.factor(genus)))

pv_M_region <- lapply(split_M_region, permanova_region)

pv_M_region <- bind_rows(pv_M_region, .id = "genus")

M_ethnicity <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_M_ethnicity <- split(M_ethnicity, with(M_ethnicity, as.factor(genus)))

pv_M_ethnicity <- lapply(split_M_ethnicity, permanova_ethnicity)

pv_M_ethnicity <- bind_rows(pv_M_ethnicity, .id = "genus")

M_countrygroup <- M %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_M_countrygroup <- split(M_countrygroup, with(M_countrygroup, as.factor(genus)))

pv_M_countrygroup <- lapply(split_M_countrygroup, permanova_countrygroup)

pv_M_countrygroup <- bind_rows(pv_M_countrygroup, .id = "genus")


sd_M <- M %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "M", 
                   sd=sd(f))

pv_M <- merge(pv_M_country, pv_M_sex, by = c("genus"), all=TRUE)

pv_M <- merge(pv_M , pv_M_countrygroup, by = c("genus"), all=TRUE)

pv_M <- merge(pv_M , pv_M_ethnicity, by = c("genus"), all=TRUE) 

pv_M <- merge(pv_M , pv_M_region, by = c("genus"), all=TRUE)

pv_M <- merge(pv_M , pv_M_studygroup, by = c("genus"), all=TRUE)

stats_M <- merge(pv_M , sd_M , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



N_country <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alloscardovia",
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
                       "CAG-274",
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
                       "Flavobacterium")) 

split_N_country <- split(N_country, with(N_country, as.factor(genus)))

pv_N_country <- lapply(split_N_country, permanova_country)

pv_N_country <- bind_rows(pv_N_country, .id = "genus")

N_sex <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Alloscardovia",
                       "Arachnia",
                       "Dermabacter",
                       "Mycobacterium",
                       "Pelistega",
                       "Pradoshia",
                       "UBA1103",
                       "UMGS743",
                       "Lactococcus_A",
                       "MGYG000004606",
                       "Weizmannia",
                       "Peptococcus",
                       "Winkia",
                       "Fannyhessea")) 

split_N_sex <- split(N_sex, with(N_sex, as.factor(genus)))

pv_N_sex <- lapply(split_N_sex, permanova_sex)

pv_N_sex <- bind_rows(pv_N_sex, .id = "genus")

N_studygroup <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Arachnia",
                       "Fructilactobacillus",
                       "ISO4-G1",
                       "MGYG000003809",
                       "Mycobacterium",
                       "Pradoshia",
                       "RUG100",
                       "ZJ304",
                       "Amedibacillus",
                       "MGYG000003807",
                       "Microvirga",
                       "Aliarcobacter",
                       "MGYG000000286",
                       "Thalassobacillus")) 

split_N_studygroup <- split(N_studygroup, with(N_studygroup, as.factor(genus)))

pv_N_studygroup <- lapply(split_N_studygroup, permanova_studygroup)

pv_N_studygroup <- bind_rows(pv_N_studygroup, .id = "genus")

N_region <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Fructilactobacillus",
                       "Gleimia",
                       "ISO4-G1",
                       "MGYG000003809",
                       "Pradoshia",
                       "UBA7741",
                       "UMGS743",
                       "Thalassobacillus"))

split_N_region <- split(N_region, with(N_region, as.factor(genus)))

pv_N_region <- lapply(split_N_region, permanova_region)

pv_N_region <- bind_rows(pv_N_region, .id = "genus")

N_ethnicity <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Fructilactobacillus",
                       "ISO4-G1",
                       "MGYG000003809",
                       "Pradoshia",
                       "Thalassobacillus")) 


split_N_ethnicity <- split(N_ethnicity, with(N_ethnicity, as.factor(genus)))

pv_N_ethnicity <- lapply(split_N_ethnicity, permanova_ethnicity)

pv_N_ethnicity <- bind_rows(pv_N_ethnicity, .id = "genus")

N_countrygroup <- N %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Arachnia",
                       "Fructilactobacillus",
                       "ISO4-G1",
                       "MGYG000003809",
                       "Pradoshia",
                       "Aliarcobacter",
                       "Thalassobacillus")) 


split_N_countrygroup <- split(N_countrygroup, with(N_countrygroup, as.factor(genus)))

pv_N_countrygroup <- lapply(split_N_countrygroup, permanova_countrygroup)

pv_N_countrygroup <- bind_rows(pv_N_countrygroup, .id = "genus")

sd_N <- N %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "N", 
                   sd=sd(f))

pv_N <- merge(pv_N_country, pv_N_sex, by = c("genus"), all=TRUE)

pv_N <- merge(pv_N , pv_N_countrygroup, by = c("genus"), all=TRUE)

pv_N <- merge(pv_N , pv_N_ethnicity, by = c("genus"), all=TRUE) 

pv_N <- merge(pv_N , pv_N_region, by = c("genus"), all=TRUE)

pv_N <- merge(pv_N , pv_N_studygroup, by = c("genus"), all=TRUE)

stats_N <- merge(pv_N , sd_N , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



O_country <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
                       "MGYG000004418",
                       "Mycoplasmopsis")) 

split_O_country <- split(O_country, with(O_country, as.factor(genus)))

pv_O_country <- lapply(split_O_country, permanova_country)

pv_O_country <- bind_rows(pv_O_country, .id = "genus")

O_sex <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("MGYG000004418")) 

split_O_sex <- split(O_sex, with(O_sex, as.factor(genus)))

pv_O_sex <- lapply(split_O_sex, permanova_sex)

pv_O_sex <- bind_rows(pv_O_sex, .id = "genus")

O_studygroup <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_O_studygroup <- split(O_studygroup, with(O_studygroup, as.factor(genus)))

pv_O_studygroup <- lapply(split_O_studygroup, permanova_studygroup)

pv_O_studygroup <- bind_rows(pv_O_studygroup, .id = "genus")

O_region <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_O_region <- split(O_region, with(O_region, as.factor(genus)))

pv_O_region <- lapply(split_O_region, permanova_region)

pv_O_region <- bind_rows(pv_O_region, .id = "genus")

O_ethnicity <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_O_ethnicity <- split(O_ethnicity, with(O_ethnicity, as.factor(genus)))

pv_O_ethnicity <- lapply(split_O_ethnicity, permanova_ethnicity)

pv_O_ethnicity <- bind_rows(pv_O_ethnicity, .id = "genus")

O_countrygroup <- O %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_O_countrygroup <- split(O_countrygroup, with(O_countrygroup, as.factor(genus)))

pv_O_countrygroup <- lapply(split_O_countrygroup, permanova_countrygroup)

pv_O_countrygroup <- bind_rows(pv_O_countrygroup, .id = "genus")


sd_O <- O %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "O", 
                   sd=sd(f))

pv_O <- merge(pv_O_country, pv_O_sex, by = c("genus"), all=TRUE)

pv_O <- merge(pv_O , pv_O_countrygroup, by = c("genus"), all=TRUE)

pv_O <- merge(pv_O , pv_O_ethnicity, by = c("genus"), all=TRUE) 

pv_O <- merge(pv_O , pv_O_region, by = c("genus"), all=TRUE)

pv_O <- merge(pv_O , pv_O_studygroup, by = c("genus"), all=TRUE)

stats_O <- merge(pv_O , sd_O , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



P_country <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma",
                       "Methanobacterium",
                       "UMGS1260",
                       "UBA3375"))

split_P_country <- split(P_country, with(P_country, as.factor(genus)))

pv_P_country <- lapply(split_P_country, permanova_country)

pv_P_country <- bind_rows(pv_P_country, .id = "genus")

P_sex <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_P_sex <- split(P_sex, with(P_sex, as.factor(genus)))

pv_P_sex <- lapply(split_P_sex, permanova_sex)

pv_P_sex <- bind_rows(pv_P_sex, .id = "genus")

P_studygroup <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("UBA1103",
                       "UBA2834"))

split_P_studygroup <- split(P_studygroup, with(P_studygroup, as.factor(genus)))

pv_P_studygroup <- lapply(split_P_studygroup, permanova_studygroup)

pv_P_studygroup <- bind_rows(pv_P_studygroup, .id = "genus")

P_region <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_P_region <- split(P_region, with(P_region, as.factor(genus)))

pv_P_region <- lapply(split_P_region, permanova_region)

pv_P_region <- bind_rows(pv_P_region, .id = "genus")

P_ethnicity <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_P_ethnicity <- split(P_ethnicity, with(P_ethnicity, as.factor(genus)))

pv_P_ethnicity <- lapply(split_P_ethnicity, permanova_ethnicity)

pv_P_ethnicity <- bind_rows(pv_P_ethnicity, .id = "genus")

P_countrygroup <- P %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_P_countrygroup <- split(P_countrygroup, with(P_countrygroup, as.factor(genus)))

pv_P_countrygroup <- lapply(split_P_countrygroup, permanova_countrygroup)

pv_P_countrygroup <- bind_rows(pv_P_countrygroup, .id = "genus") 

sd_P <- P %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "P", 
                   sd=sd(f))

pv_P <- merge(pv_P_country, pv_P_sex, by = c("genus"), all=TRUE)

pv_P <- merge(pv_P , pv_P_countrygroup, by = c("genus"), all=TRUE)

pv_P <- merge(pv_P , pv_P_ethnicity, by = c("genus"), all=TRUE) 

pv_P <- merge(pv_P , pv_P_region, by = c("genus"), all=TRUE)

pv_P <- merge(pv_P , pv_P_studygroup, by = c("genus"), all=TRUE)

stats_P <- merge(pv_P , sd_P , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



Q_country <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
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
                       "UBA3375")) 

split_Q_country <- split(Q_country, with(Q_country, as.factor(genus)))

pv_Q_country <- lapply(split_Q_country, permanova_country)

pv_Q_country <- bind_rows(pv_Q_country, .id = "genus")

Q_sex <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
                       "Epilithonimonas",
                       "Gleimia",
                       "HGM10890",
                       "HGM11377",
                       "MGYG000000286",
                       "MGYG000000407",
                       "MGYG000003416",
                       "MGYG000004606",
                       "Methanocorpusculum",
                       "Pseudomonas_B",
                       "Robertmurraya",
                       "Staphylococcus_A",
                       "Timonella",
                       "UBA7862",
                       "Vitreoscilla",
                       "Companilactobacillus",
                       "Filifactor",
                       "SFEB01",
                       "SFWF01",
                       "Winkia",
                       "MGYG000004537",
                       "Nigerium",
                       "Arachnia",
                       "QFNR01"))

split_Q_sex <- split(Q_sex, with(Q_sex, as.factor(genus)))

pv_Q_sex <- lapply(split_Q_sex, permanova_sex)

pv_Q_sex <- bind_rows(pv_Q_sex, .id = "genus")

Q_studygroup <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
                       "Gleimia",
                       "HGM10890",
                       "MGYG000000286",
                       "MGYG000000407",
                       "Methanocorpusculum",
                       "Pseudomonas_B",
                       "Staphylococcus_A",
                       "UBA4877",
                       "Filifactor",
                       "Niallia",
                       "Numidum",
                       "UMGS1707",
                       "Winkia",
                       "Schleiferilactobacillus")) 

split_Q_studygroup <- split(Q_studygroup, with(Q_studygroup, as.factor(genus)))

pv_Q_studygroup <- lapply(split_Q_studygroup, permanova_studygroup)

pv_Q_studygroup <- bind_rows(pv_Q_studygroup, .id = "genus")

Q_region <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Epilithonimonas",
                       "MGYG000000286",
                       "Methanocorpusculum",
                       "Staphylococcus_A",
                       "Winkia",
                       "Helicobacter_F")) 

split_Q_region <- split(Q_region, with(Q_region, as.factor(genus)))

pv_Q_region <- lapply(split_Q_region, permanova_region)

pv_Q_region <- bind_rows(pv_Q_region, .id = "genus")

Q_ethnicity <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Gleimia",
                       "HGM10890",
                       "MGYG000000286",
                       "Methanocorpusculum",
                       "Staphylococcus_A",
                       "Winkia")) 

split_Q_ethnicity <- split(Q_ethnicity, with(Q_ethnicity, as.factor(genus)))

pv_Q_ethnicity <- lapply(split_Q_ethnicity, permanova_ethnicity)

pv_Q_ethnicity <- bind_rows(pv_Q_ethnicity, .id = "genus")

Q_countrygroup <- Q %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Acetomicrobium",
                       "HGM10890",
                       "MGYG000000286",
                       "Methanocorpusculum",
                       "Pseudomonas_B",
                       "Staphylococcus_A",
                       "Numidum",
                       "Winkia")) 

split_Q_countrygroup <- split(Q_countrygroup, with(Q_countrygroup, as.factor(genus)))

pv_Q_countrygroup <- lapply(split_Q_countrygroup, permanova_countrygroup)

pv_Q_countrygroup <- bind_rows(pv_Q_countrygroup, .id = "genus")


sd_Q <- Q %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "Q", 
                   sd=sd(f))

pv_Q <- merge(pv_Q_country, pv_Q_sex, by = c("genus"), all=TRUE)

pv_Q <- merge(pv_Q , pv_Q_countrygroup, by = c("genus"), all=TRUE)

pv_Q <- merge(pv_Q , pv_Q_ethnicity, by = c("genus"), all=TRUE) 

pv_Q <- merge(pv_Q , pv_Q_region, by = c("genus"), all=TRUE)

pv_Q <- merge(pv_Q , pv_Q_studygroup, by = c("genus"), all=TRUE)

stats_Q <- merge(pv_Q , sd_Q , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


S_country <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_country <- split(S_country, with(S_country, as.factor(genus)))

pv_S_country <- lapply(split_S_country, permanova_country)

pv_S_country <- bind_rows(pv_S_country, .id = "genus")

S_sex <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_sex <- split(S_sex, with(S_sex, as.factor(genus)))

pv_S_sex <- lapply(split_S_sex, permanova_sex)

pv_S_sex <- bind_rows(pv_S_sex, .id = "genus")

S_studygroup <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_studygroup <- split(S_studygroup, with(S_studygroup, as.factor(genus)))

pv_S_studygroup <- lapply(split_S_studygroup, permanova_studygroup)

pv_S_studygroup <- bind_rows(pv_S_studygroup, .id = "genus")

S_region <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_region <- split(S_region, with(S_region, as.factor(genus)))

pv_S_region <- lapply(split_S_region, permanova_region)

pv_S_region <- bind_rows(pv_S_region, .id = "genus")

S_ethnicity <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_ethnicity <- split(S_ethnicity, with(S_ethnicity, as.factor(genus)))

pv_S_ethnicity <- lapply(split_S_ethnicity, permanova_ethnicity)

pv_S_ethnicity <- bind_rows(pv_S_ethnicity, .id = "genus")

S_countrygroup <- S %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_S_countrygroup <- split(S_countrygroup, with(S_countrygroup, as.factor(genus)))

pv_S_countrygroup <- lapply(split_S_countrygroup, permanova_countrygroup)

pv_S_countrygroup <- bind_rows(pv_S_countrygroup, .id = "genus")



sd_S <- S %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "S", 
                   sd=sd(f))

pv_S <- merge(pv_S_country, pv_S_sex, by = c("genus"), all=TRUE)

pv_S <- merge(pv_S , pv_S_countrygroup, by = c("genus"), all=TRUE)

pv_S <- merge(pv_S , pv_S_ethnicity, by = c("genus"), all=TRUE) 

pv_S <- merge(pv_S , pv_S_region, by = c("genus"), all=TRUE)

pv_S <- merge(pv_S , pv_S_studygroup, by = c("genus"), all=TRUE)

stats_S <- merge(pv_S , sd_S , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



T_country <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma",
                       "Nosocomiicoccus",
                       "UBA1103",
                       "Jonquetella",
                       "Lapidilactobacillus",
                       "HGM10611",
                       "MGYG000004606",
                       "Helicobacter_F",
                       "Exiguobacterium_A",
                       "Ureaplasma",
                       "ZJ304")) 

split_T_country <- split(T_country, with(T_country, as.factor(genus)))

pv_T_country <- lapply(split_T_country, permanova_country)

pv_T_country <- bind_rows(pv_T_country, .id = "genus")

T_sex <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma",
                       "Mycoplasmopsis",
                       "Nosocomiicoccus",
                       "UBA1103",
                       "Jonquetella",
                       "Chryseobacterium")) 

split_T_sex <- split(T_sex, with(T_sex, as.factor(genus)))

pv_T_sex <- lapply(split_T_sex, permanova_sex)

pv_T_sex <- bind_rows(pv_T_sex, .id = "genus")

T_studygroup <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Mycoplasmopsis")) 

split_T_studygroup <- split(T_studygroup, with(T_studygroup, as.factor(genus)))

pv_T_studygroup <- lapply(split_T_studygroup, permanova_studygroup)

pv_T_studygroup <- bind_rows(pv_T_studygroup, .id = "genus")

T_region <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Metamycoplasma")) 

split_T_region <- split(T_region, with(T_region, as.factor(genus)))

pv_T_region <- lapply(split_T_region, permanova_region)

pv_T_region <- bind_rows(pv_T_region, .id = "genus")

T_ethnicity <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_T_ethnicity <- split(T_ethnicity, with(T_ethnicity, as.factor(genus)))

pv_T_ethnicity <- lapply(split_T_ethnicity, permanova_ethnicity)

pv_T_ethnicity <- bind_rows(pv_T_ethnicity, .id = "genus")

T_countrygroup <- T %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_T_countrygroup <- split(T_countrygroup, with(T_countrygroup, as.factor(genus)))

pv_T_countrygroup <- lapply(split_T_countrygroup, permanova_countrygroup)

pv_T_countrygroup <- bind_rows(pv_T_countrygroup, .id = "genus")



sd_T <- T %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "T", 
                   sd=sd(f))

pv_T <- merge(pv_T_country, pv_T_sex, by = c("genus"), all=TRUE)

pv_T <- merge(pv_T , pv_T_countrygroup, by = c("genus"), all=TRUE)

pv_T <- merge(pv_T , pv_T_ethnicity, by = c("genus"), all=TRUE) 

pv_T <- merge(pv_T , pv_T_region, by = c("genus"), all=TRUE)

pv_T <- merge(pv_T , pv_T_studygroup, by = c("genus"), all=TRUE)

stats_T <- merge(pv_T , sd_T , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



U_country <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Chryseobacterium",
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
                       "Planococcus")) 

split_U_country <- split(U_country, with(U_country, as.factor(genus)))

pv_U_country <- lapply(split_U_country, permanova_country)

pv_U_country <- bind_rows(pv_U_country, .id = "genus")

U_sex <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Risungbinella",
                       "Atopobium",
                       "MGYG000003772",
                       "Metamycoplasma"))

split_U_sex <- split(U_sex, with(U_sex, as.factor(genus)))

pv_U_sex <- lapply(split_U_sex, permanova_sex)

pv_U_sex <- bind_rows(pv_U_sex, .id = "genus")

U_studygroup <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Paenibacillus_F",
                       "Nesterenkonia")) 

split_U_studygroup <- split(U_studygroup, with(U_studygroup, as.factor(genus)))

pv_U_studygroup <- lapply(split_U_studygroup, permanova_studygroup)

pv_U_studygroup <- bind_rows(pv_U_studygroup, .id = "genus")

U_region <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_U_region <- split(U_region, with(U_region, as.factor(genus)))

pv_U_region <- lapply(split_U_region, permanova_region)

pv_U_region <- bind_rows(pv_U_region, .id = "genus")

U_ethnicity <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1)

split_U_ethnicity <- split(U_ethnicity, with(U_ethnicity, as.factor(genus)))

pv_U_ethnicity <- lapply(split_U_ethnicity, permanova_ethnicity)

pv_U_ethnicity <- bind_rows(pv_U_ethnicity, .id = "genus")

U_countrygroup <- U %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Paenibacillus_F")) 

split_U_countrygroup <- split(U_countrygroup, with(U_countrygroup, as.factor(genus)))

pv_U_countrygroup <- lapply(split_U_countrygroup, permanova_countrygroup)

pv_U_countrygroup <- bind_rows(pv_U_countrygroup, .id = "genus")


sd_U <- U %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "U", 
                   sd=sd(f))

pv_U <- merge(pv_U_country, pv_U_sex, by = c("genus"), all=TRUE)

pv_U <- merge(pv_U , pv_U_countrygroup, by = c("genus"), all=TRUE)

pv_U <- merge(pv_U , pv_U_ethnicity, by = c("genus"), all=TRUE) 

pv_U <- merge(pv_U , pv_U_region, by = c("genus"), all=TRUE)

pv_U <- merge(pv_U , pv_U_studygroup, by = c("genus"), all=TRUE)

stats_U <- merge(pv_U , sd_U , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



V_country <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum",
                       "MGYG000004418",
                       "UBA7862")) 

split_V_country <- split(V_country, with(V_country, as.factor(genus)))

pv_V_country <- lapply(split_V_country, permanova_country)

pv_V_country <- bind_rows(pv_V_country, .id = "genus")

V_sex <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_V_sex <- split(V_sex, with(V_sex, as.factor(genus)))

pv_V_sex <- lapply(split_V_sex, permanova_sex)

pv_V_sex <- bind_rows(pv_V_sex, .id = "genus")

V_studygroup <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_V_studygroup <- split(V_studygroup, with(V_studygroup, as.factor(genus)))

pv_V_studygroup <- lapply(split_V_studygroup, permanova_studygroup)

pv_V_studygroup <- bind_rows(pv_V_studygroup, .id = "genus")

V_region <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_V_region <- split(V_region, with(V_region, as.factor(genus)))

pv_V_region <- lapply(split_V_region, permanova_region)

pv_V_region <- bind_rows(pv_V_region, .id = "genus")

V_ethnicity <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_V_ethnicity <- split(V_ethnicity, with(V_ethnicity, as.factor(genus)))

pv_V_ethnicity <- lapply(split_V_ethnicity, permanova_ethnicity)

pv_V_ethnicity <- bind_rows(pv_V_ethnicity, .id = "genus")

V_countrygroup <- V %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Methanocorpusculum")) 

split_V_countrygroup <- split(V_countrygroup, with(V_countrygroup, as.factor(genus)))

pv_V_countrygroup <- lapply(split_V_countrygroup, permanova_countrygroup)

pv_V_countrygroup <- bind_rows(pv_V_countrygroup, .id = "genus")


sd_V <- V %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "V", 
                   sd=sd(f))

pv_V <- merge(pv_V_country, pv_V_sex, by = c("genus"), all=TRUE)

pv_V <- merge(pv_V , pv_V_countrygroup, by = c("genus"), all=TRUE)

pv_V <- merge(pv_V , pv_V_ethnicity, by = c("genus"), all=TRUE) 

pv_V <- merge(pv_V , pv_V_region, by = c("genus"), all=TRUE)

pv_V <- merge(pv_V , pv_V_studygroup, by = c("genus"), all=TRUE)

stats_V <- merge(pv_V , sd_V , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)

W_country <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bulleidia",
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
                       "Lactococcus")) 

split_W_country <- split(W_country, with(W_country, as.factor(genus)))

pv_W_country <- lapply(split_W_country, permanova_country)

pv_W_country <- bind_rows(pv_W_country, .id = "genus")

W_sex <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("AF33-28",
                       "Bulleidia",
                       "CAG-533",
                       "CAG-831",
                       "DNF00809",
                       "Dysgonomonas",
                       "Edwardsiella",
                       "Enorma",
                       "Eubacterium_J",
                       "Levyella",
                       "MGYG000002075",
                       "Niameybacter",
                       "Pseudobutyrivibrio",
                       "Pseudoflavonifractor",
                       "RQCD01",
                       "RUG705",
                       "SFEB01",
                       "UBA11774",
                       "UBA1234",
                       "UMGS1585",
                       "UMGS1840",
                       "CAG-1427",
                       "CAG-345",
                       "CAG-914",
                       "HGM11530",
                       "MGYG000003772",
                       "Pauljensenia",
                       "SFFH01",
                       "UBA1685",
                       "UBA7597",
                       "UMGS1600",
                       "UMGS1663",
                       "UMGS416",
                       "An23",
                       "Coprobacillus",
                       "Ligilactobacillus",
                       "UMGS1872",
                       "SFMI01")) 

split_W_sex <- split(W_sex, with(W_sex, as.factor(genus)))

pv_W_sex <- lapply(split_W_sex, permanova_sex)

pv_W_sex <- bind_rows(pv_W_sex, .id = "genus")

W_studygroup <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("AF33-28",
                       "Anaerocolumna",
                       "Bulleidia",
                       "CAG-274",
                       "CAG-710",
                       "CAG-977",
                       "Edwardsiella",
                       "Enorma",
                       "Flavonifractor",
                       "HGM11575",
                       "Mobilibacterium",
                       "Pseudoflavonifractor",
                       "RUG806",
                       "SFEL01",
                       "UBA11774",
                       "UBA1234",
                       "UBA7862",
                       "UMGS1441",
                       "UMGS1840",
                       "CAG-568",
                       "CAG-617",
                       "F082",
                       "HGM11530",
                       "NK3B98",
                       "UMGS1071",
                       "UMGS416",
                       "Fusobacterium_A",
                       "UBA1206",
                       "Borkfalkia",
                       "Leclercia",
                       "Fibrobacter_A")) 

split_W_studygroup <- split(W_studygroup, with(W_studygroup, as.factor(genus)))

pv_W_studygroup <- lapply(split_W_studygroup, permanova_studygroup)

pv_W_studygroup <- bind_rows(pv_W_studygroup, .id = "genus")

W_region <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bulleidia",
                       "CAG-274",
                       "CAG-307",
                       "CAG-977",
                       "Dysgonomonas",
                       "Enorma",
                       "RQCD01",
                       "RUG806",
                       "SFEL01",
                       "UBA7862",
                       "UMGS1585",
                       "UMGS403",
                       "CAG-568",
                       "Clostridium_AA",
                       "F082",
                       "HGM11530",
                       "UBA2821",
                       "UMGS1071",
                       "UMGS416",
                       "Faecalibacillus",
                       "UBA1206",
                       "Holdemania")) 

split_W_region <- split(W_region, with(W_region, as.factor(genus)))

pv_W_region <- lapply(split_W_region, permanova_region)

pv_W_region <- bind_rows(pv_W_region, .id = "genus")

W_ethnicity <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bulleidia",
                       "CAG-274",
                       "CAG-977",
                       "Enorma",
                       "HGM11575",
                       "RUG806",
                       "SFEL01",
                       "UBA7862",
                       "F082",
                       "HGM11530",
                       "UMGS1071",
                       "UMGS416",
                       "UBA1206")) 

split_W_ethnicity <- split(W_ethnicity, with(W_ethnicity, as.factor(genus)))

pv_W_ethnicity <- lapply(split_W_ethnicity, permanova_ethnicity)

pv_W_ethnicity <- bind_rows(pv_W_ethnicity, .id = "genus")

W_countrygroup <- W %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Bulleidia",
                       "CAG-274",
                       "CAG-977",
                       "Edwardsiella",
                       "Enorma",
                       "Flavonifractor",
                       "Mobilibacterium",
                       "Pseudoflavonifractor",
                       "RUG806",
                       "SFEL01",
                       "UBA1234",
                       "UBA7862",
                       "UMGS1840",
                       "CAG-568",
                       "CAG-617",
                       "F082",
                       "HGM11530",
                       "UMGS1071",
                       "UMGS416",
                       "UBA1206")) 

split_W_countrygroup <- split(W_countrygroup, with(W_countrygroup, as.factor(genus)))

pv_W_countrygroup <- lapply(split_W_countrygroup, permanova_countrygroup)

pv_W_countrygroup <- bind_rows(pv_W_countrygroup, .id = "genus")


sd_W <- W %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "W", 
                   sd=sd(f))

pv_W <- merge(pv_W_country, pv_W_sex, by = c("genus"), all=TRUE)

pv_W <- merge(pv_W , pv_W_countrygroup, by = c("genus"), all=TRUE)

pv_W <- merge(pv_W , pv_W_ethnicity, by = c("genus"), all=TRUE) 

pv_W <- merge(pv_W , pv_W_region, by = c("genus"), all=TRUE)

pv_W <- merge(pv_W , pv_W_studygroup, by = c("genus"), all=TRUE)

stats_W <- merge(pv_W , sd_W , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


X_country <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Mycoplasmopsis",
                       "MGYG000003772",
                       "Mobilicoccus",
                       "Jonquetella",
                       "UMGS1260",
                       "28L"))

split_X_country <- split(X_country, with(X_country, as.factor(genus)))

pv_X_country <- lapply(split_X_country, permanova_country)

pv_X_country <- bind_rows(pv_X_country, .id = "genus")

X_sex <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Timonella", 
                       "MGYG000003772")) 

split_X_sex <- split(X_sex, with(X_sex, as.factor(genus)))

pv_X_sex <- lapply(split_X_sex, permanova_sex)

pv_X_sex <- bind_rows(pv_X_sex, .id = "genus")

X_studygroup <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_X_studygroup <- split(X_studygroup, with(X_studygroup, as.factor(genus)))

pv_X_studygroup <- lapply(split_X_studygroup, permanova_studygroup)

pv_X_studygroup <- bind_rows(pv_X_studygroup, .id = "genus")

X_region <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_X_region <- split(X_region, with(X_region, as.factor(genus)))

pv_X_region <- lapply(split_X_region, permanova_region)

pv_X_region <- bind_rows(pv_X_region, .id = "genus")

X_ethnicity <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_X_ethnicity <- split(X_ethnicity, with(X_ethnicity, as.factor(genus)))

pv_X_ethnicity <- lapply(split_X_ethnicity, permanova_ethnicity)

pv_X_ethnicity <- bind_rows(pv_X_ethnicity, .id = "genus")

X_countrygroup <- X %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_X_countrygroup <- split(X_countrygroup, with(X_countrygroup, as.factor(genus)))

pv_X_countrygroup <- lapply(split_X_countrygroup, permanova_countrygroup)

pv_X_countrygroup <- bind_rows(pv_X_countrygroup, .id = "genus")



sd_X <- X %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "X", 
                   sd=sd(f))

pv_X <- merge(pv_X_country, pv_X_sex, by = c("genus"), all=TRUE)

pv_X <- merge(pv_X , pv_X_countrygroup, by = c("genus"), all=TRUE)

pv_X <- merge(pv_X , pv_X_ethnicity, by = c("genus"), all=TRUE) 

pv_X <- merge(pv_X , pv_X_region, by = c("genus"), all=TRUE)

pv_X <- merge(pv_X , pv_X_studygroup, by = c("genus"), all=TRUE)

stats_X <- merge(pv_X , sd_X , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)



Y_country <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_country <- split(Y_country, with(Y_country, as.factor(genus)))

pv_Y_country <- lapply(split_Y_country, permanova_country)

pv_Y_country <- bind_rows(pv_Y_country, .id = "genus")

Y_sex <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_sex <- split(Y_sex, with(Y_sex, as.factor(genus)))

pv_Y_sex <- lapply(split_Y_sex, permanova_sex)

pv_Y_sex <- bind_rows(pv_Y_sex, .id = "genus")

Y_studygroup <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_studygroup <- split(Y_studygroup, with(Y_studygroup, as.factor(genus)))

pv_Y_studygroup <- lapply(split_Y_studygroup, permanova_studygroup)

pv_Y_studygroup <- bind_rows(pv_Y_studygroup, .id = "genus")

Y_region <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_region <- split(Y_region, with(Y_region, as.factor(genus)))

pv_Y_region <- lapply(split_Y_region, permanova_region)

pv_Y_region <- bind_rows(pv_Y_region, .id = "genus")

Y_ethnicity <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_ethnicity <- split(Y_ethnicity, with(Y_ethnicity, as.factor(genus)))

pv_Y_ethnicity <- lapply(split_Y_ethnicity, permanova_ethnicity)

pv_Y_ethnicity <- bind_rows(pv_Y_ethnicity, .id = "genus")

Y_countrygroup <- Y %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) 

split_Y_countrygroup <- split(Y_countrygroup, with(Y_countrygroup, as.factor(genus)))

pv_Y_countrygroup <- lapply(split_Y_countrygroup, permanova_countrygroup)

pv_Y_countrygroup <- bind_rows(pv_Y_countrygroup, .id = "genus")




sd_Y <- Y %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "Y", 
                   sd=sd(f))

pv_Y <- merge(pv_Y_country, pv_Y_sex, by = c("genus"), all=TRUE)

pv_Y <- merge(pv_Y , pv_Y_countrygroup, by = c("genus"), all=TRUE)

pv_Y <- merge(pv_Y , pv_Y_ethnicity, by = c("genus"), all=TRUE) 

pv_Y <- merge(pv_Y , pv_Y_region, by = c("genus"), all=TRUE)

pv_Y <- merge(pv_Y , pv_Y_studygroup, by = c("genus"), all=TRUE)

stats_Y <- merge(pv_Y , sd_Y , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


Z_country <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("CAG-776",
                       "CAG-822",
                       "Cloacibacillus",
                       "Emergencia",
                       "Eubacterium",
                       "Faecalicoccus",
                       "HGM12669",
                       "Lachnoanaerobaculum",
                       "Lachnoclostridium",
                       "Lacrimispora",
                       "Leclercia",
                       "Longibaculum",
                       "MGYG000000728",
                       "MGYG000000747",
                       "MGYG000002010",
                       "MGYG000002695",
                       "MGYG000003809",
                       "MGYG000003970",
                       "MGYG000004033",
                       "MGYG000004429",
                       "Merdimonas",
                       "Paenibacillus_A",
                       "Rothia",
                       "Salmonella",
                       "UBA4675",
                       "UMGS1975",
                       "UMGS599",
                       "UMGS731",
                       "Bacillus",
                       "CAAFZY01",
                       "CAG-484",
                       "Caecibacter",
                       "Citrobacter_A",
                       "Clostridium_AP",
                       "Enterococcus_B",
                       "Frisingicoccus",
                       "HGM11386",
                       "HGM11808",
                       "Johnsonella",
                       "MGYG000000387",
                       "MGYG000004485",
                       "Methanomassiliicoccus_A",
                       "Raoultibacter",
                       "Ruminococcus_B",
                       "UBA11549",
                       "UBA2804",
                       "UBA4855",
                       "UMGS1004",
                       "Eggerthella",
                       "Enterococcus",
                       "HGM13010",
                       "MGYG000002719",
                       "MGYG000004172",
                       "MGYG000004380",
                       "Neobacillus",
                       "GCA-900066755",
                       "HGM10836",
                       "Lachnoclostridium_A",
                       "Peptacetobacter",
                       "Tidjanibacter",
                       "CAG-462",
                       "MGYG000002696",
                       "Selenomonas_A",
                       "UMGS1601",
                       "Ligilactobacillus",
                       "Parasutterella",
                       "Turicibacter",
                       "Pygmaiobacter",
                       "Schaedlerella",
                       "MGYG000004445",
                       "Fusobacterium",
                       "Acidaminococcus",
                       "UMGS1839",
                       "Kluyvera",
                       "Muribaculum",
                       "UBA3375",
                       "Morganella",
                       "CAG-274")) 

split_Z_country <- split(Z_country, with(Z_country, as.factor(genus)))

pv_Z_country <- lapply(split_Z_country, permanova_country)

pv_Z_country <- bind_rows(pv_Z_country, .id = "genus")

Z_sex <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("CAG-776",
                       "Cloacibacillus",
                       "Emergencia",
                       "Lachnoclostridium",
                       "Lacrimispora",
                       "Leclercia",
                       "MGYG000000747",
                       "MGYG000002010",
                       "MGYG000003809",
                       "MGYG000003970",
                       "MGYG000004033",
                       "Merdimonas",
                       "Paenibacillus_A",
                       "Rothia",
                       "Salmonella",
                       "UBA1375",
                       "UBA1390",
                       "UBA4675",
                       "UBA5416",
                       "UBA9414",
                       "UMGS1851",
                       "UMGS1975",
                       "CAAFZY01",
                       "CAG-1031",
                       "CAG-484",
                       "Clostridium_AP",
                       "Enterococcus_B",
                       "Frisingicoccus",
                       "Johnsonella",
                       "QALW01",
                       "UBA2804",
                       "UBA6382",
                       "UMGS1908",
                       "Slackia_A",
                       "TWA4",
                       "UMGS1241",
                       "HGM11501",
                       "Lachnoclostridium_A",
                       "Tidjanibacter",
                       "UBA11471",
                       "UBA1232",
                       "RUG14215")) 

split_Z_sex <- split(Z_sex, with(Z_sex, as.factor(genus)))

pv_Z_sex <- lapply(split_Z_sex, permanova_sex)

pv_Z_sex <- bind_rows(pv_Z_sex, .id = "genus")

Z_studygroup <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Cloacibacillus",
                       "HGM12669",
                       "Leclercia",
                       "Longibaculum",
                       "MGYG000000350",
                       "MGYG000000747",
                       "MGYG000002010",
                       "MGYG000002695",
                       "MGYG000003945",
                       "MGYG000004243",
                       "MGYG000004429",
                       "Paenibacillus_A",
                       "Salmonella",
                       "UBA1375",
                       "UBA1390",
                       "UMGS599",
                       "UMGS731",
                       "UMGS946",
                       "CAG-484",
                       "Caecibacter",
                       "Providencia",
                       "Ruminococcus_B",
                       "UBA11549",
                       "UBA6382",
                       "AF33-28",
                       "Plesiomonas", 
                       "Emergencia", 
                       "Eubacterium")) 

split_Z_studygroup <- split(Z_studygroup, with(Z_studygroup, as.factor(genus)))

pv_Z_studygroup <- lapply(split_Z_studygroup, permanova_studygroup)

pv_Z_studygroup <- bind_rows(pv_Z_studygroup, .id = "genus")

Z_region <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Emergencia",
                       "HGM12669",
                       "Longibaculum",
                       "MGYG000000728",
                       "MGYG000002010",
                       "MGYG000002695",
                       "MGYG000003809",
                       "MGYG000004429",
                       "Salmonella",
                       "UBA4675",
                       "UMGS599",
                       "Caecibacter",
                       "Citrobacter_A",
                       "Methanomassiliicoccus_A",
                       "UBA11549",
                       "MGYG000004172",
                       "Tidjanibacter",
                       "Parasutterella")) 

split_Z_region <- split(Z_region, with(Z_region, as.factor(genus)))

pv_Z_region <- lapply(split_Z_region, permanova_region)

pv_Z_region <- bind_rows(pv_Z_region, .id = "genus")

Z_ethnicity <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Emergencia",
                       "HGM12669",
                       "Longibaculum",
                       "MGYG000002010",
                       "MGYG000002695",
                       "MGYG000004429",
                       "Salmonella",
                       "UBA1390",
                       "UMGS731",
                       "Caecibacter",
                       "UBA11549"))

split_Z_ethnicity <- split(Z_ethnicity, with(Z_ethnicity, as.factor(genus)))

pv_Z_ethnicity <- lapply(split_Z_ethnicity, permanova_ethnicity)

pv_Z_ethnicity <- bind_rows(pv_Z_ethnicity, .id = "genus")

Z_countrygroup <- Z %>%
  group_by(genus) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>%
  filter(!genus %in% c("Cloacibacillus",
                       "Emergencia",
                       "Eubacterium",
                       "HGM12669",
                       "Leclercia",
                       "Longibaculum",
                       "MGYG000000747",
                       "MGYG000002010",
                       "MGYG000002695",
                       "MGYG000004429",
                       "Paenibacillus_A",
                       "Salmonella",
                       "UMGS599",
                       "UMGS731",
                       "CAG-484",
                       "Caecibacter",
                       "Ruminococcus_B",
                       "UBA11549")) 

split_Z_countrygroup <- split(Z_countrygroup, with(Z_countrygroup, as.factor(genus)))

pv_Z_countrygroup <- lapply(split_Z_countrygroup, permanova_countrygroup)

pv_Z_countrygroup <- bind_rows(pv_Z_countrygroup, .id = "genus")

sd_Z <- Z %>%
  group_by(genus) %>%
  dplyr::summarize(cog_category = "Z", 
                   sd=sd(f))

pv_Z <- merge(pv_Z_country, pv_Z_sex, by = c("genus"), all=TRUE)

pv_Z <- merge(pv_Z , pv_Z_countrygroup, by = c("genus"), all=TRUE)

pv_Z <- merge(pv_Z , pv_Z_ethnicity, by = c("genus"), all=TRUE) 

pv_Z <- merge(pv_Z , pv_Z_region, by = c("genus"), all=TRUE)

pv_Z <- merge(pv_Z , pv_Z_studygroup, by = c("genus"), all=TRUE)

stats_Z <- merge(pv_Z , sd_Z , by = c("genus"), all=TRUE) %>% 
  drop_na(sd)


df_list <- list(stats_A, stats_B, stats_C, stats_D, stats_E, stats_F, stats_G, stats_H, stats_I, stats_J, stats_K, stats_L, stats_M, stats_N, stats_O, stats_P, stats_Q, stats_S, stats_T, stats_U, stats_V, stats_W, stats_X, stats_Y, stats_Z)      

#merge all data frames together ----

stats_all <- df_list %>% bind_rows() 

#I want to adjust the p values to make sure they really represent significance

stats_all$pvcountry <- p.adjust(stats_all$pvcountry, method = "BH")

stats_all$pvcountrygroup <- p.adjust(stats_all$pvcountrygroup, method = "BH") 

stats_all$region <- p.adjust(stats_all$pvregion, method = "BH")

stats_all$pvsex <- p.adjust(stats_all$pvsex, method = "BH")

stats_all$pvethnicity <- p.adjust(stats_all$pvethnicity, method = "BH")

stats_all$pvstudygroup <- p.adjust(stats_all$pvstudygroup, method = "BH")

write_csv(stats_all, "data/cog_genus_variance_values.csv") 

#stats_all <- read_csv("data/cog_genus_variance_values.csv")

#wait I also want my ranks for DA for all genera, let's see if we can add that here: 

stats_all <- stats_all %>% 
  mutate(taxon = tolower(genus))

DArank <- read.csv("output/DAranks/everyone.csv", header=TRUE) %>% 
  filter(grepl('g_', taxon)) %>%
  mutate(across(everything(), gsub, pattern = "g_", replacement = "")) %>% 
  mutate(na_sum = as.numeric(na_sum)) %>%
  mutate(DArank = 55-na_sum) %>% 
  select(taxon, DArank)

stats_all <- merge(stats_all, DArank, by = "taxon") %>% 
  mutate(DArank = as.numeric(DArank)) 

#and also let's add the the number of zeros for each genera so we can not look at super low abundance taxa: 

otu <- read.csv("data/OTU_tables/zerocount_kraken2_READS_sg_2024_decontam_G.csv", header=TRUE) %>% 
  rename(taxon = genus) %>%
  mutate(across(everything(), gsub, pattern = "g_", replacement = "")) %>% 
  select(taxon, sum_na, low_value_count, low_value_count)

stats_all <- merge(stats_all, otu, by="taxon") %>% 
  select(-taxon) %>% 
  ungroup() %>% 
  mutate(low_value_count = as.numeric(low_value_count))

#I also want to make two types of plots to answer two main questions ---- 

#first, does f vary across individuals for DA taxa? ----

q1_A <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count) %>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>% 
  filter(sum_na < 1 & low_value_count < 15)

q1_A %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'A', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/A.pdf", width = 40, height = 20, units = "cm")

q1_B <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_B %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'B', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/B.pdf", width = 40, height = 20, units = "cm")

q1_C <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_C %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'C', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/C.pdf", width = 40, height = 20, units = "cm")

q1_D <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_D %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'D', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/D.pdf", width = 40, height = 20, units = "cm")

q1_E <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_E %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'E', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/E.pdf", width = 40, height = 20, units = "cm")

q1_F <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_F %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'F', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/F.pdf", width = 40, height = 20, units = "cm")

q1_G <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_G %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'G', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/G.pdf", width = 40, height = 20, units = "cm")

q1_H <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)


q1_H %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'H', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/H.pdf", width = 40, height = 20, units = "cm")

q1_I <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_I %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'I', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/I.pdf", width = 40, height = 20, units = "cm")


q1_J <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_J %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'J', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/J.pdf", width = 40, height = 20, units = "cm")

q1_K <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_K %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'K', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/K.pdf", width = 40, height = 20, units = "cm")

q1_L <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_L %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'L', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/L.pdf", width = 40, height = 20, units = "cm")

q1_M <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_M %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'M', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/M.pdf", width = 40, height = 20, units = "cm")


q1_N <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_N %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'N', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/N.pdf", width = 40, height = 20, units = "cm")

q1_O <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_O %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'O', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/O.pdf", width = 40, height = 20, units = "cm")

q1_P <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_P %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'P', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/O.pdf", width = 40, height = 20, units = "cm")


q1_Q <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_Q %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'Q', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/Q.pdf", width = 40, height = 20, units = "cm")

q1_S <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_S %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'S', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/S.pdf", width = 40, height = 20, units = "cm")

q1_T <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_T %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'T', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/T.pdf", width = 40, height = 20, units = "cm")

q1_U <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_U %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'U', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/U.pdf", width = 40, height = 20, units = "cm")

q1_V <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_V %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'V', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/V.pdf", width = 40, height = 20, units = "cm")

q1_W <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_W %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'W', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/W.pdf", width = 40, height = 20, units = "cm")

q1_X <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_X %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'X', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/X.pdf", width = 40, height = 20, units = "cm")

q1_Y <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_Y %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'Y', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/Y.pdf", width = 40, height = 20, units = "cm")

q1_Z <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, sd, DArank, sum_na, low_value_count)%>% 
  arrange(-sd) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q1_Z %>% 
  ggplot(aes(x=fct_reorder(genus, -sd), y=sd, fill = DArank)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title= 'Z', x= 'Genus', y= 'Standard derivation', fill = "number of groups with differential abudance ") + 
  guides(fill = guide_colourbar(position = "top")) +
  theme(plot.title = element_text(hjust = 0, vjust = -5)) + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q1/Z.pdf", width = 40, height = 20, units = "cm")

#second, does f vary between groups of vary between groups of individuals, for each cog category and each taxa? 

q2_A_country <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y=Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Acountry.pdf", width = 40, height = 20, units = "cm")

q2_A_studygroup <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Astudygroup.pdf", width = 40, height = 20, units = "cm")

q2_A_ethnicity <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Aethnicity.pdf", width = 40, height = 20, units = "cm")

q2_A_region <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Aregion.pdf", width = 40, height = 20, units = "cm")

q2_A_sex <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Asex.pdf", width = 40, height = 20, units = "cm")

q2_A_countrygroup <- stats_all %>% 
  filter(cog_category == "A") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_A_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'A - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Acountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_B_country <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bcountry.pdf", width = 40, height = 20, units = "cm")

q2_B_studygroup <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_B_ethnicity <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bethnicity.pdf", width = 40, height = 20, units = "cm")

q2_B_region <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bregion.pdf", width = 40, height = 20, units = "cm")

q2_B_sex <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bsex.pdf", width = 40, height = 20, units = "cm")

q2_B_countrygroup <- stats_all %>% 
  filter(cog_category == "B") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_B_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'B - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Bcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_C_country <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ccountry.pdf", width = 40, height = 20, units = "cm")

q2_C_studygroup <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Cstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_C_ethnicity <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Cethnicity.pdf", width = 40, height = 20, units = "cm")

q2_C_region <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Cregion.pdf", width = 40, height = 20, units = "cm")

q2_C_sex <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Csex.pdf", width = 40, height = 20, units = "cm")

q2_C_countrygroup <- stats_all %>% 
  filter(cog_category == "C") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_C_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'C - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ccountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_D_country <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dcountry.pdf", width = 40, height = 20, units = "cm")

q2_D_studygroup <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_D_ethnicity <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dethnicity.pdf", width = 40, height = 20, units = "cm")

q2_D_region <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dregion.pdf", width = 40, height = 20, units = "cm")

q2_D_sex <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dsex.pdf", width = 40, height = 20, units = "cm")

q2_D_countrygroup <- stats_all %>% 
  filter(cog_category == "D") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_D_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'D - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Dcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_E_country <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ecountry.pdf", width = 40, height = 20, units = "cm")

q2_E_studygroup <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Estudygroup.pdf", width = 40, height = 20, units = "cm")

q2_E_ethnicity <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Eethnicity.pdf", width = 40, height = 20, units = "cm")

q2_E_region <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Eregion.pdf", width = 40, height = 20, units = "cm")

q2_E_sex <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Esex.pdf", width = 40, height = 20, units = "cm")

q2_E_countrygroup <- stats_all %>% 
  filter(cog_category == "E") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_E_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'E - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ecountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_F_country <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fcountry.pdf", width = 40, height = 20, units = "cm")

q2_F_studygroup <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_F_ethnicity <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fethnicity.pdf", width = 40, height = 20, units = "cm")

q2_F_region <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fregion.pdf", width = 40, height = 20, units = "cm")

q2_F_sex <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fsex.pdf", width = 40, height = 20, units = "cm")

q2_F_countrygroup <- stats_all %>% 
  filter(cog_category == "F") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_F_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'F - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Fcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_G_country <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gcountry.pdf", width = 40, height = 20, units = "cm")

q2_G_studygroup <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_G_ethnicity <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gethnicity.pdf", width = 40, height = 20, units = "cm")

q2_G_region <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gregion.pdf", width = 40, height = 20, units = "cm")

q2_G_sex <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gsex.pdf", width = 40, height = 20, units = "cm")

q2_G_countrygroup <- stats_all %>% 
  filter(cog_category == "G") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_G_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'G - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Gcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_H_country <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hcountry.pdf", width = 40, height = 20, units = "cm")

q2_H_studygroup <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_H_ethnicity <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hethnicity.pdf", width = 40, height = 20, units = "cm")

q2_H_region <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hregion.pdf", width = 40, height = 20, units = "cm")

q2_H_sex <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hsex.pdf", width = 40, height = 20, units = "cm")

q2_H_countrygroup <- stats_all %>% 
  filter(cog_category == "H") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_H_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'H - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Hcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_I_country <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Icountry.pdf", width = 40, height = 20, units = "cm")

q2_I_studygroup <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Istudygroup.pdf", width = 40, height = 20, units = "cm")

q2_I_ethnicity <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Iethnicity.pdf", width = 40, height = 20, units = "cm")

q2_I_region <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Iregion.pdf", width = 40, height = 20, units = "cm")

q2_I_sex <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Isex.pdf", width = 40, height = 20, units = "cm")

q2_I_countrygroup <- stats_all %>% 
  filter(cog_category == "I") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_I_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'I - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Icountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_J_country <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jcountry.pdf", width = 40, height = 20, units = "cm")

q2_J_studygroup <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_J_ethnicity <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jethnicity.pdf", width = 40, height = 20, units = "cm")

q2_J_region <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jregion.pdf", width = 40, height = 20, units = "cm")

q2_J_sex <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jsex.pdf", width = 40, height = 20, units = "cm")

q2_J_countrygroup <- stats_all %>% 
  filter(cog_category == "J") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_J_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'J - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Jcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_K_country <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Kcountry.pdf", width = 40, height = 20, units = "cm")

q2_K_studygroup <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Kstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_K_ethnicity <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Kethnicity.pdf", width = 40, height = 20, units = "cm")

q2_K_region <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Kregion.pdf", width = 40, height = 20, units = "cm")

q2_K_sex <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ksex.pdf", width = 40, height = 20, units = "cm")

q2_K_countrygroup <- stats_all %>% 
  filter(cog_category == "K") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_K_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'K - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Kcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_L_country <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lcountry.pdf", width = 40, height = 20, units = "cm")

q2_L_studygroup <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_L_ethnicity <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lethnicity.pdf", width = 40, height = 20, units = "cm")

q2_L_region <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lregion.pdf", width = 40, height = 20, units = "cm")

q2_L_sex <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lsex.pdf", width = 40, height = 20, units = "cm")

q2_L_countrygroup <- stats_all %>% 
  filter(cog_category == "L") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_L_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'L - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Lcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_M_country <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Mcountry.pdf", width = 40, height = 20, units = "cm")

q2_M_studygroup <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Mstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_M_ethnicity <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Methnicity.pdf", width = 40, height = 20, units = "cm")

q2_M_region <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Mregion.pdf", width = 40, height = 20, units = "cm")

q2_M_sex <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Msex.pdf", width = 40, height = 20, units = "cm")

q2_M_countrygroup <- stats_all %>% 
  filter(cog_category == "M") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_M_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'M - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Mcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_N_country <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ncountry.pdf", width = 40, height = 20, units = "cm")

q2_N_studygroup <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Nstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_N_ethnicity <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Nethnicity.pdf", width = 40, height = 20, units = "cm")

q2_N_region <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Nregion.pdf", width = 40, height = 20, units = "cm")

q2_N_sex <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Nsex.pdf", width = 40, height = 20, units = "cm")

q2_N_countrygroup <- stats_all %>% 
  filter(cog_category == "N") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_N_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'N - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ncountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_O_country <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ocountry.pdf", width = 40, height = 20, units = "cm")

q2_O_studygroup <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ostudygroup.pdf", width = 40, height = 20, units = "cm")

q2_O_ethnicity <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Oethnicity.pdf", width = 40, height = 20, units = "cm")

q2_O_region <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Oregion.pdf", width = 40, height = 20, units = "cm")

q2_O_sex <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Osex.pdf", width = 40, height = 20, units = "cm")

q2_O_countrygroup <- stats_all %>% 
  filter(cog_category == "O") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_O_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'O - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ocountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_P_country <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Pcountry.pdf", width = 40, height = 20, units = "cm")

q2_P_studygroup <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Pstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_P_ethnicity <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Pethnicity.pdf", width = 40, height = 20, units = "cm")

q2_P_region <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Pregion.pdf", width = 40, height = 20, units = "cm")

q2_P_sex <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Psex.pdf", width = 40, height = 20, units = "cm")

q2_P_countrygroup <- stats_all %>% 
  filter(cog_category == "P") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_P_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'P - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Pcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_Q_country <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qcountry.pdf", width = 40, height = 20, units = "cm")

q2_Q_studygroup <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_Q_ethnicity <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qethnicity.pdf", width = 40, height = 20, units = "cm")

q2_Q_region <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qregion.pdf", width = 40, height = 20, units = "cm")

q2_Q_sex <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qsex.pdf", width = 40, height = 20, units = "cm")

q2_Q_countrygroup <- stats_all %>% 
  filter(cog_category == "Q") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Q_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Q - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Qcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_S_country <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'S - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Scountry.pdf", width = 40, height = 20, units = "cm")

q2_S_studygroup <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'S - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Sstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_S_ethnicity <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c("red","green")) + 
  labs(title= 'S - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Sethnicity.pdf", width = 40, height = 20, units = "cm")

q2_S_region <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c("red","green")) + 
  labs(title= 'S - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Sregion.pdf", width = 40, height = 20, units = "cm")

q2_S_sex <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'S - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ssex.pdf", width = 40, height = 20, units = "cm")

q2_S_countrygroup <- stats_all %>% 
  filter(cog_category == "S") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_S_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c("red","green")) + 
  labs(title= 'S - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Scountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_T_country <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tcountry.pdf", width = 40, height = 20, units = "cm")

q2_T_studygroup <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_T_ethnicity <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tethnicity.pdf", width = 40, height = 20, units = "cm")

q2_T_region <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tregion.pdf", width = 40, height = 20, units = "cm")

q2_T_sex <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tsex.pdf", width = 40, height = 20, units = "cm")

q2_T_countrygroup <- stats_all %>% 
  filter(cog_category == "T") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_T_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'T - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Tcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_U_country <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ucountry.pdf", width = 40, height = 20, units = "cm")

q2_U_studygroup <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ustudygroup.pdf", width = 40, height = 20, units = "cm")

q2_U_ethnicity <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Uethnicity.pdf", width = 40, height = 20, units = "cm")

q2_U_region <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Uregion.pdf", width = 40, height = 20, units = "cm")

q2_U_sex <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Usex.pdf", width = 40, height = 20, units = "cm")

q2_U_countrygroup <- stats_all %>% 
  filter(cog_category == "U") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_U_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'U - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ucountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_V_country <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vcountry.pdf", width = 40, height = 20, units = "cm")

q2_V_studygroup <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_V_ethnicity <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vethnicity.pdf", width = 40, height = 20, units = "cm")

q2_V_region <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vregion.pdf", width = 40, height = 20, units = "cm")

q2_V_sex <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vsex.pdf", width = 40, height = 20, units = "cm")

q2_V_countrygroup <- stats_all %>% 
  filter(cog_category == "V") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_V_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'V - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Vcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_W_country <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wcountry.pdf", width = 40, height = 20, units = "cm")

q2_W_studygroup <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_W_ethnicity <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wethnicity.pdf", width = 40, height = 20, units = "cm")

q2_W_region <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wregion.pdf", width = 40, height = 20, units = "cm")

q2_W_sex <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wsex.pdf", width = 40, height = 20, units = "cm")

q2_W_countrygroup <- stats_all %>% 
  filter(cog_category == "W") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_W_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'W - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Wcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_X_country <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xcountry.pdf", width = 40, height = 20, units = "cm")

q2_X_studygroup <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_X_ethnicity <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xethnicity.pdf", width = 40, height = 20, units = "cm")

q2_X_region <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xregion.pdf", width = 40, height = 20, units = "cm")

q2_X_sex <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xsex.pdf", width = 40, height = 20, units = "cm")

q2_X_countrygroup <- stats_all %>% 
  filter(cog_category == "X") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_X_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'X - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Xcountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_Y_country <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ycountry.pdf", width = 40, height = 20, units = "cm")

q2_Y_studygroup <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ystudygroup.pdf", width = 40, height = 20, units = "cm")

q2_Y_ethnicity <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Yethnicity.pdf", width = 40, height = 20, units = "cm")

q2_Y_region <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Yregion.pdf", width = 40, height = 20, units = "cm")

q2_Y_sex <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ysex.pdf", width = 40, height = 20, units = "cm")

q2_Y_countrygroup <- stats_all %>% 
  filter(cog_category == "Y") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Y_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Y - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Ycountrygroup.pdf", width = 40, height = 20, units = "cm")

q2_Z_country <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvcountry, Fcountry, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountry) %>% 
  arrange(pvcountry) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_country %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountry), y= Fcountry, alpha = DArank, fill = pvcountry < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Country', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zcountry.pdf", width = 40, height = 20, units = "cm")

q2_Z_studygroup <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvstudygroup, Fstudygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvstudygroup) %>% 
  arrange(pvstudygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_studygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fstudygroup), y= Fstudygroup, alpha = DArank, fill = pvstudygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Study Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zstudygroup.pdf", width = 40, height = 20, units = "cm")

q2_Z_ethnicity <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvethnicity, Fethnicity,DArank, sum_na, low_value_count)%>% 
  drop_na(pvethnicity) %>% 
  arrange(pvethnicity) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_ethnicity %>% 
  ggplot(aes(x=fct_reorder(genus, -Fethnicity), y= Fethnicity, alpha = DArank, fill = pvethnicity < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Ethnicity', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zethnicity.pdf", width = 40, height = 20, units = "cm")

q2_Z_region <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvregion, Fregion, DArank, sum_na, low_value_count)%>% 
  drop_na(pvregion) %>% 
  arrange(pvregion) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_region %>% 
  ggplot(aes(x=fct_reorder(genus, -Fregion), y= Fregion, alpha = DArank, fill = pvregion < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Region', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zregion.pdf", width = 40, height = 20, units = "cm")

q2_Z_sex <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvsex, Fsex, DArank, sum_na, low_value_count)%>% 
  drop_na(pvsex) %>% 
  arrange(pvsex) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_sex %>% 
  ggplot(aes(x=fct_reorder(genus, -Fsex), y= Fsex, alpha = DArank, fill = pvsex < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Sex', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zsex.pdf", width = 40, height = 20, units = "cm")

q2_Z_countrygroup <- stats_all %>% 
  filter(cog_category == "Z") %>%
  select(genus, cog_category, pvcountrygroup, Fcountrygroup, DArank, sum_na, low_value_count)%>% 
  drop_na(pvcountrygroup) %>% 
  arrange(pvcountrygroup) %>% 
  filter(DArank %in% c(8:18) | genus %in% c("Prevotella", "Bacteroides")) %>%    filter(sum_na < 1 & low_value_count < 15)

q2_Z_countrygroup %>% 
  ggplot(aes(x=fct_reorder(genus, -Fcountrygroup), y= Fcountrygroup, alpha = DArank, fill = pvcountrygroup < 0.05)) +
  geom_col() + 
  scale_fill_manual(values = c( "red", "green")) + 
  labs(title= 'Z - Country Group', x= 'Genus', y= 'F statistic', alpha = "number of groups with differential abudance", fill = "significant p value") +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(plot.title = element_text(hjust = 0.8, vjust = -20)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=20,face="bold"), title = element_text(size=22,face="bold"), plot.tag = element_text(size = 20))

ggsave("output/q2/Zcountrygroup.pdf", width = 40, height = 20, units = "cm")

#okay, I also want to add a averaged value for each of these cog category: 

stat_average <- stats_all %>% 
  select(genus, cog_category, Fcountry, pvcountry, Fsex, pvsex, Fstudygroup, pvstudygroup, Fethnicity, pvethnicity, Fregion, pvregion, Fcountrygroup, pvcountrygroup, sd) %>% 
  pivot_longer(Fcountry:sd, names_to = "grouping", values_to = "stat_value") %>% 
  group_by(cog_category, grouping) %>% 
  mutate(mean = mean(stat_value, group_by=T, na.rm =TRUE)) %>% 
  select(-genus, -stat_value) %>% 
  distinct()

write_csv(stat_average, "data/cog_genus_variance_means.csv") 

