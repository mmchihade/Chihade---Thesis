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


#load data ----

data_grp <- read.csv("data/clean_simple_metadata_tishkofflab_2024.csv", header=TRUE) %>%
  select(-record_id, -lab_id, -topmed_id, -collection_name, -sample_site, -longitude, -latitude, -date, -year, -is_adult, -age, -waist_circumference, -body_fat_percent, -bmi) %>% 
  rename(id_lab = field_id) %>%
  mutate(id_lab = toupper(id_lab))

#create a dataframe with all of the cog tally files, with the file name as a column so that we can match column to tally

merged_df <- list.files(path = "data/cog_tally_overall",
                        pattern = "\\.csv$", 
                        full.names = TRUE) %>% 
  set_names() %>%
  map_dfr(read_csv, .id = "file_name") 

#we need the cog categories too! 

cog <- read.csv("data/cog_categories.csv", header=TRUE)

#now let's clean the data up  

cog_tally <- clean_names(merged_df) %>% 
  rename(lab_id = file_name) %>% 
  select(-genus) %>% 
  mutate(across(everything(), gsub, pattern = "data/cog_tally_overall/", replacement = "")) %>%
  mutate(across(everything(), gsub, pattern = ".csv", replacement = "")) %>%
  mutate(n = as.numeric(n)) %>%
  group_by(lab_id, cog_category) %>% 
  mutate(tally = sum(n, na.rm =TRUE, group_by=TRUE)) %>% 
  select(-n) %>%
  distinct()

#I also have to replace all of the TZ codes with B codes so that I can divide by country and ethnicity 

tally <- cog_tally %>% 
  mutate(lab_id = recode(lab_id,
                         "TZBG067"="B607",
                         "TZBG069"="B609",
                         "TZBG070"="B610",
                         "TZBG072"="B612",
                         "TZBG073"="B613",
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
  rename(id_lab = lab_id) %>% 
  mutate(cog_category = recode(cog_category,
                               "?" = "S",
                               "!" = "S")) %>% 
  group_by(id_lab, cog_category) %>% 
  mutate(tally = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#I also want to separate all of the multiple letter columns. They will probably be useful later on, but it is hella overwhelming now
full_tally <- tally %>% 
  separate(cog_category, into = c("cog_1", "cog_2", "cog_3", "cog_4", "cog_5", "cog_6", "cog_7"), sep = "") %>%
  select(-cog_1) %>% 
  pivot_longer(names_to="cog_name", values_to = "cog_category", cog_2:cog_7)%>%
  na.omit() %>%
  select(-cog_name) %>%
  group_by(id_lab, cog_category) %>% 
  mutate(tally = sum(tally, na.rm =TRUE, group_by=TRUE)) %>%
  distinct()

#umm, I am going to try to make a combined dataset, let's see if R can handle it ----

merge <- merge(full_tally, data_grp, by = "id_lab")

merge_labels <- merge(merge, cog, by = "cog_category")

#yay!!!!!

#let's find p values for all our factors of interest ----

#first I need to create separate dataframes for each cog category 

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
P <- split_df$P
R <- split_df$R
S <- split_df$S
T <- split_df$T
U <- split_df$U 
V <- split_df$V 
W <- split_df$W
X <- split_df$X
Y <- split_df$Y
Z <- split_df$Z

pv_a <- A %>% 
  dplyr::summarize(cog_category = "A",
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = A))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_c <- C %>% 
  dplyr::summarize(cog_category = "C",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = C))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_d <- D %>% 
  dplyr::summarize(cog_category = "D",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = D))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_e <- E %>% 
  dplyr::summarize(cog_category = "E",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = E))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_f <- F %>% 
  dplyr::summarize(cog_category = "F",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = F))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_g <- G %>% 
  dplyr::summarize(cog_category = "G",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = G))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_h <- H %>% 
  dplyr::summarize(cog_category = "H",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = H))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_i <- I %>% 
  dplyr::summarize(cog_category = "I",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = I))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_j <- J %>% 
  dplyr::summarize(cog_category = "J",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = J))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_k <- K %>% 
  dplyr::summarize(cog_category = "K",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = K))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_l <- L %>% 
  dplyr::summarize(cog_category = "L",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = L))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_m <- M %>% 
  dplyr::summarize(cog_category = "M",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = M))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_n <- N %>% 
  dplyr::summarize(cog_category = "N",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = N))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_o <- O %>% 
  dplyr::summarize(cog_category = "O",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = O))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_p <- P %>% 
  dplyr::summarize(cog_category = "P",
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = P))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_q <- Q %>% 
  dplyr::summarize(cog_category = "Q",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = Q))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_s <- S %>% 
  dplyr::summarize(cog_category = "S",
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = S))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_t <- T %>% 
  dplyr::summarize(cog_category = "T",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = T))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_u <- U %>% 
  dplyr::summarize(cog_category = "U",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = U))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_v <- V %>% 
  dplyr::summarize(cog_category = "V",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = V))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_w <- W %>% 
  dplyr::summarize(cog_category = "W",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = W))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))

pv_x <- X %>% 
  dplyr::summarize(cog_category = "X",
                   
                   pvcountry =kruskal.test(tally ~ country)$p.value,
                   pvsex =summary(aov(tally ~ sex, data = X))[[1]][["Pr(>F)"]][[1]],
                   pvethnicity =kruskal.test(tally ~ ethnicity)$p.value,
                   pvregion =kruskal.test(tally ~ region)$p.value, 
                   pvstudy_group =kruskal.test(tally ~ study_group)$p.value,
                   pvcountry_group =kruskal.test(tally ~ country_group)$p.value, 
                   stdev = sd(tally))



df_list <- list(pv_a, pv_c, pv_d, pv_e, pv_f, pv_g, pv_h, pv_i, pv_j, pv_k, pv_l, pv_m, pv_n, pv_o, pv_q, pv_p, pv_s, pv_t, pv_u, pv_v, pv_w, pv_x)      

#merge all data frames together

pv_all <- df_list %>% bind_rows() 

pv_all$pvcountry <- p.adjust(pv_all$pvcountry, method = "BH")

pv_all$pvcountry_group <- p.adjust(pv_all$pvcountry_group, method = "BH") 

pv_all$pvregion <- p.adjust(pv_all$pvregion, method = "BH")

pv_all$pvsex <- p.adjust(pv_all$pvsex, method = "BH")

pv_all$pvethnicity <- p.adjust(pv_all$pvethnicity, method = "BH")

pv_all$pvstudy_group <- p.adjust(pv_all$pvstudy_group, method = "BH")

write_csv(pv_all, "data/cog_variance_values.csv")








