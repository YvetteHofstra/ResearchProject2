# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#               Research Project 2, Data Analysis
#       The impact of salt tolerance on pollinator behaviour:
#      effects of salinity on floral traits of Medicago sativa
#                         Year: 2026
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ---- Before staring ----

# Clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# Restore the correct versions of the used packages
renv::restore() 

# Load the packages that are needed for this project
library(car)
library(chemodiv)
library(corrr)
library(emmeans)
library(factoextra)
library(glmmTMB)
library(gtsummary)
library(MASS)
library(multcomp)
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately

readr::read_csv # This makes a tibble instead of a table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Loading the data ----

Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

Flowers_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=166596564&single=true&output=csv")

Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

Observations_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=2034963164&single=true&output=csv")

Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=224800747&single=true&output=csv") 

## Replace all NAs in the entire data frame with 0 ##
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

## Show the type of data R will treat it as (numeric, character, factor, etc.) ##
str(Phenotype)
str(Nectar)
str(Flowers)
str(Flowers_2)
str(Flowering_date)
str(Repotting)
str(Observations)
str(Soil)

## Make numeric where needed ##
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)
Flowers_2$Number_Inflorescences <- as.numeric(Flowers_2$Number_Inflorescences)
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)

## Make the yes and no into factor ##
Repotting$Nodules_present <- factor(
  Repotting$Nodules_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

Repotting$Seeds_present <- factor(
  Repotting$Seeds_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

## Correct the observations (June/July) ##

# Filter the data to only include flowering plants and replace NA with 0 for arthropod counts after filtering flowering plants. This is important because we only want to consider arthropod visits to flowering plants, and any NA values in the arthropod counts should be treated as zero visits.
Observations_clean <- Observations_2 %>%
  filter(Flowering == 1)

arth_cols <- c(
  "Aglais_io", 	"Andrena", "Anthomyiidae",	"Apis_mellifera", "Autographa_gamma",	"Bibio_marci",	"Bombus_lapidarius", "Bombus_pascuorum", "Bombus_terrestris", "Celastrina_argiolus", "Closterotomus_norwegicus", "Cynomya_mortuorum", "Empis_tessellata", "Larinioides_cornutus", "Oedemera",	"Pieris_rapae", "Pollenia",	"Polyommatus_icarus", "Pseudovadonia_livida",	"Rhagonycha_fulva",	"Sarcophagidae",	"Stenichneumon_culpator",	"Syrphidae",	"Tetanocera", "Thereva_nobilitata",	"Thymelicus_lineola",	"Vanessa_atalanta", "Zygaena_filipendulae"
)

Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ as.numeric(as.character(.))))

# Replace NA with 0 
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ replace_na(., 0)))

# Pivot longer
Observations_long <- Observations_clean %>%
  pivot_longer(
    cols = all_of(arth_cols),
    names_to = "Visitor",
    values_to = "Count"
  )

print(table(Observations_long$Count, useNA = "ifany"))

arth_summary <- Observations_long %>%
  group_by(Cultivar, Visitor, Treatment_worded) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

print(summary(arth_summary$Count))


# ---- Use for model significance testing ----

# For models without interactions (Type II LRT)
# Anova(best_model, type = "II")

# For models with interactions (Type III LRT)
# Anova(best_model, type = "III")

# When more than 2 levels, use emmeans to get pairwise comparisons
# emmeans(best_model, pairwise ~ Cultivar, type = "response")

# When there is an interaction
# emmeans(best_model, pairwise ~ Cultivar | Treatment_worded)


# ---- Phenotype models ----

## Number of flowers ##

# Candidate GLMs
m1 <- glm.nb(Number_flowers ~ Cultivar,
             data = Phenotype)
m2 <- glm.nb(Number_flowers ~ Treatment_worded,
             data = Phenotype)
m3 <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded,
             data = Phenotype)
m4 <- glm.nb(Number_flowers ~ Cultivar * Treatment_worded,
             data = Phenotype)

# Candidate GLMMs
m5 <- glmmTMB(Number_flowers ~ Cultivar + (1|Block),
              family = nbinom2, data = Phenotype)
m6 <- glmmTMB(Number_flowers ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)
m7 <- glmmTMB(Number_flowers ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)
m8 <- glmmTMB(Number_flowers ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")

# Post hoc comparisons (only if significant)
emmeans(best_model, pairwise ~ Cultivar, type = "response")


## Number of inflorescences (March) ##

# Candidate GLMs
m1 <- glm.nb(Number_inflorescences ~ Cultivar,
             data = Phenotype)
m2 <- glm.nb(Number_inflorescences ~ Treatment_worded,
             data = Phenotype)
m3 <- glm.nb(Number_inflorescences ~ Cultivar + Treatment_worded,
             data = Phenotype)
m4 <- glm.nb(Number_inflorescences ~ Cultivar * Treatment_worded,
             data = Phenotype)

# Candidate GLMMs
m5 <- glmmTMB(Number_inflorescences ~ Cultivar + (1|Block),
              family = nbinom2, data = Phenotype)
m6 <- glmmTMB(Number_inflorescences ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)
m7 <- glmmTMB(Number_inflorescences ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)
m8 <- glmmTMB(Number_inflorescences ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Phenotype)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")

# Post hoc comparisons (only if significant)
emmeans(best_model, pairwise ~ Cultivar, type = "response")


## Number of inflorescences (June/July) ##

# Candidate GLMs
m1 <- glm.nb(Number_Inflorescences ~ Cultivar,
             data = Flowers_2)
m2 <- glm.nb(Number_Inflorescences ~ Treatment_worded,
             data = Flowers_2)
m3 <- glm.nb(Number_Inflorescences ~ Cultivar + Treatment_worded,
             data = Flowers_2)
m4 <- glm.nb(Number_Inflorescences ~ Cultivar * Treatment_worded,
             data = Flowers_2)

# Candidate GLMMs
m5 <- glmmTMB(Number_Inflorescences ~ Cultivar + (1|Block),
              family = nbinom2, data = Flowers_2)
m6 <- glmmTMB(Number_Inflorescences ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)
m7 <- glmmTMB(Number_Inflorescences ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)
m8 <- glmmTMB(Number_Inflorescences ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Nectar models ----

## Nectar collection all plants ##

# Candidate GLMs
m1 <- glm.nb(Microliter ~ Cultivar, 
             data = Nectar)
m2 <- glm.nb(Microliter ~ Treatment_worded, 
             data = Nectar)
m3 <- glm.nb(Microliter ~ Cultivar + Treatment_worded,
             data = Nectar)
m4 <- glm.nb(Microliter ~ Cultivar * Treatment_worded,
             data = Nectar)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")


## Nectar collection fully harvested plants ##

# First remove them, using their Plant_IDs
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))

# Candidate GLMs
m1 <- glm.nb(Microliter ~ Cultivar, 
             data = Nectar_cleaned)
m2 <- glm.nb(Microliter ~ Treatment_worded, 
             data = Nectar_cleaned)
m3 <- glm.nb(Microliter ~ Cultivar + Treatment_worded,
             data = Nectar_cleaned)
m4 <- glm.nb(Microliter ~ Cultivar * Treatment_worded,
             data = Nectar_cleaned)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Flower (round) models ----

## Average length (June/July only) ##

# Candidate GLMs
m1 <- glm.nb(Average_Inflorescence_Length ~ Cultivar,
             data = Flowers_2)
m2 <- glm.nb(Average_Inflorescence_Length ~ Treatment_worded,
             data = Flowers_2)
m3 <- glm.nb(Average_Inflorescence_Length ~ Cultivar + Treatment_worded,
             data = Flowers_2)
m4 <- glm.nb(Average_Inflorescence_Length ~ Cultivar * Treatment_worded,
             data = Flowers_2)

# Candidate GLMMs
m5 <- glmmTMB(Average_Inflorescence_Length ~ Cultivar + (1|Block),
              family = nbinom2, data = Flowers_2)
m6 <- glmmTMB(Average_Inflorescence_Length ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)
m7 <- glmmTMB(Average_Inflorescence_Length ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)
m8 <- glmmTMB(Average_Inflorescence_Length ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Flowers_2)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")


## Number of inflorescences closest to nectar collection (April only) ##

# This may be informative for significance between nectar and number of inflorescences

# Candidate GLMs
m1 <- glm.nb(Number_Inflorescences ~ Cultivar,
             data = Flowers)
m2 <- glm.nb(Number_Inflorescences ~ Treatment_worded,
             data = Flowers)
m3 <- glm.nb(Number_Inflorescences ~ Cultivar + Treatment_worded,
             data = Flowers)
m4 <- glm.nb(Number_Inflorescences ~ Cultivar * Treatment_worded,
             data = Flowers)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")

# Post hoc comparisons (only if significant)
emmeans(best_model, pairwise ~ Cultivar, type = "response")

## Average length (April only) ##

# Candidate GLMs
m1 <- glm.nb(Average_Inflorescence_Length ~ Cultivar,
             data = Flowers)
m2 <- glm.nb(Average_Inflorescence_Length ~ Treatment_worded,
             data = Flowers)
m3 <- glm.nb(Average_Inflorescence_Length ~ Cultivar + Treatment_worded,
             data = Flowers)
m4 <- glm.nb(Average_Inflorescence_Length ~ Cultivar * Treatment_worded,
             data = Flowers)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Flower (date) models ----

## Flowering date field experiment ##

# Candidate GLMs
m1 <- glm.nb(Date_numbered ~ Cultivar,
             data = Flowering_date)
m2 <- glm.nb(Date_numbered ~ Treatment_worded,
             data = Flowering_date)
m3 <- glm.nb(Date_numbered ~ Cultivar + Treatment_worded,
             data = Flowering_date)
m4 <- glm.nb(Date_numbered ~ Cultivar * Treatment_worded,
             data = Flowering_date)

# Candidate GLMMs
m5 <- glmmTMB(Date_numbered ~ Cultivar + (1|Block),
              family = nbinom2, data = Flowering_date)
m6 <- glmmTMB(Date_numbered ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Flowering_date)
m7 <- glmmTMB(Date_numbered ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Flowering_date)
m8 <- glmmTMB(Date_numbered ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Flowering_date)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Repotting/seed (pod) models ----

## Seed presence ##

# Candidate GLMs
m1 <- glm.nb(Seeds_present ~ Cultivar,
             data = Nectar)
m2 <- glm.nb(Seeds_present ~ Treatment_worded,
             data = Nectar)
m3 <- glm.nb(Seeds_present ~ Cultivar + Treatment_worded,
             data = Nectar)
m4 <- glm.nb(Seeds_present ~ Cultivar * Treatment_worded,
             data = Nectar)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## Seed pod abundance ##

# Candidate GLMs
m1 <- glm.nb(Seed_pod_abundance ~ Cultivar,
             data = Nectar)
m2 <- glm.nb(Seed_pod_abundance ~ Treatment_worded,
             data = Nectar)
m3 <- glm.nb(Seed_pod_abundance ~ Cultivar + Treatment_worded,
             data = Nectar)
m4 <- glm.nb(Seed_pod_abundance ~ Cultivar * Treatment_worded,
             data = Nectar)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")


## Seed abundance ##

# Candidate GLMs
m1 <- glm.nb(Seed_abundance ~ Cultivar,
             data = Nectar)
m2 <- glm.nb(Seed_abundance ~ Treatment_worded,
             data = Nectar)
m3 <- glm.nb(Seed_abundance ~ Cultivar + Treatment_worded,
             data = Nectar)
m4 <- glm.nb(Seed_abundance ~ Cultivar * Treatment_worded,
             data = Nectar)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")

# Post hoc comparisons (only if significant)
emmeans(best_model, pairwise ~ Cultivar, type = "response")


## Seed weight ##

# Candidate GLMs
m1 <- glm.nb(Seed_weight ~ Cultivar,
             data = Nectar)
m2 <- glm.nb(Seed_weight ~ Treatment_worded,
             data = Nectar)
m3 <- glm.nb(Seed_weight ~ Cultivar + Treatment_worded,
             data = Nectar)
m4 <- glm.nb(Seed_weight ~ Cultivar * Treatment_worded,
             data = Nectar)

# Test the models using AIC
AIC(m1, m2, m3, m4)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## Nodule abundance ##

# Candidate GLMs
m1 <- glm.nb(Nodule_abundance ~ Cultivar,
             data = Repotting)
m2 <- glm.nb(Nodule_abundance ~ Treatment_worded,
             data = Repotting)
m3 <- glm.nb(Nodule_abundance ~ Cultivar + Treatment_worded,
             data = Repotting)
m4 <- glm.nb(Nodule_abundance ~ Cultivar * Treatment_worded,
             data = Repotting)

# Candidate GLMMs
m5 <- glmmTMB(Nodule_abundance ~ Cultivar + (1|Block),
              family = nbinom2, data = Repotting)
m6 <- glmmTMB(Nodule_abundance ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m7 <- glmmTMB(Nodule_abundance ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m8 <- glmmTMB(Nodule_abundance ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m1

# Test fixed effects
car::Anova(best_model, type = "II")


## Nodule shape ##

# Candidate GLMs
m1 <- glm.nb(Nodule_shape ~ Cultivar,
             data = Repotting)
m2 <- glm.nb(Nodule_shape ~ Treatment_worded,
             data = Repotting)
m3 <- glm.nb(Nodule_shape ~ Cultivar + Treatment_worded,
             data = Repotting)
m4 <- glm.nb(Nodule_shape ~ Cultivar * Treatment_worded,
             data = Repotting)

# Candidate GLMMs
m5 <- glmmTMB(Nodule_shape ~ Cultivar + (1|Block),
              family = nbinom2, data = Repotting)
m6 <- glmmTMB(Nodule_shape ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m7 <- glmmTMB(Nodule_shape ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m8 <- glmmTMB(Nodule_shape ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## Nodule size ##

# Candidate GLMs
m1 <- glm.nb(Nodule_size ~ Cultivar,
             data = Repotting)
m2 <- glm.nb(Nodule_size ~ Treatment_worded,
             data = Repotting)
m3 <- glm.nb(Nodule_size ~ Cultivar + Treatment_worded,
             data = Repotting)
m4 <- glm.nb(Nodule_size ~ Cultivar * Treatment_worded,
             data = Repotting)

# Candidate GLMMs
m5 <- glmmTMB(Nodule_size ~ Cultivar + (1|Block),
              family = nbinom2, data = Repotting)
m6 <- glmmTMB(Nodule_size ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m7 <- glmmTMB(Nodule_size ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m8 <- glmmTMB(Nodule_size ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## Root abundance ##

# Candidate GLMs
m1 <- glm.nb(Root_abundance ~ Cultivar,
             data = Repotting)
m2 <- glm.nb(Root_abundance ~ Treatment_worded,
             data = Repotting)
m3 <- glm.nb(Root_abundance ~ Cultivar + Treatment_worded,
             data = Repotting)
m4 <- glm.nb(Root_abundance ~ Cultivar * Treatment_worded,
             data = Repotting)

# Candidate GLMMs
m5 <- glmmTMB(Root_abundance ~ Cultivar + (1|Block),
              family = nbinom2, data = Repotting)
m6 <- glmmTMB(Root_abundance ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m7 <- glmmTMB(Root_abundance ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)
m8 <- glmmTMB(Root_abundance ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Repotting)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Soil models (First measurements, test if treatment works) ----

## Wetness (%) of soil ##

# Candidate GLMs
m1 <- glm.nb(Wet ~ Cultivar,
             data = Soil)
m2 <- glm.nb(Wet ~ Treatment_worded,
             data = Soil)
m3 <- glm.nb(Wet ~ Cultivar + Treatment_worded,
             data = Soil)
m4 <- glm.nb(Wet ~ Cultivar * Treatment_worded,
             data = Soil)

# Candidate GLMMs
m5 <- glmmTMB(Wet ~ Cultivar + (1|Block),
              family = nbinom2, data = Soil)
m6 <- glmmTMB(Wet ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m7 <- glmmTMB(Wet ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m8 <- glmmTMB(Wet ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## ECp ##

# Candidate GLMs
m1 <- glm.nb(ECp ~ Cultivar,
             data = Soil)
m2 <- glm.nb(ECp ~ Treatment_worded,
             data = Soil)
m3 <- glm.nb(ECp ~ Cultivar + Treatment_worded,
             data = Soil)
m4 <- glm.nb(ECp ~ Cultivar * Treatment_worded,
             data = Soil)

# Candidate GLMMs
m5 <- glmmTMB(ECp ~ Cultivar + (1|Block),
              family = nbinom2, data = Soil)
m6 <- glmmTMB(ECp ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m7 <- glmmTMB(ECp ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m8 <- glmmTMB(ECp ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m6

# Test fixed effects
car::Anova(best_model, type = "II")


## Eb ##

# Candidate GLMs
m1 <- glm.nb(Eb ~ Cultivar,
             data = Soil)
m2 <- glm.nb(Eb ~ Treatment_worded,
             data = Soil)
m3 <- glm.nb(Eb ~ Cultivar + Treatment_worded,
             data = Soil)
m4 <- glm.nb(Eb ~ Cultivar * Treatment_worded,
             data = Soil)

# Candidate GLMMs
m5 <- glmmTMB(Eb ~ Cultivar + (1|Block),
              family = nbinom2, data = Soil)
m6 <- glmmTMB(Eb ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m7 <- glmmTMB(Eb ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m8 <- glmmTMB(Eb ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


## ECb ##

# Candidate GLMs
m1 <- glm.nb(ECb ~ Cultivar,
             data = Soil)
m2 <- glm.nb(ECb ~ Treatment_worded,
             data = Soil)
m3 <- glm.nb(ECb ~ Cultivar + Treatment_worded,
             data = Soil)
m4 <- glm.nb(ECb ~ Cultivar * Treatment_worded,
             data = Soil)

# Candidate GLMMs
m5 <- glmmTMB(ECb ~ Cultivar + (1|Block),
              family = nbinom2, data = Soil)
m6 <- glmmTMB(ECb ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m7 <- glmmTMB(ECb ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)
m8 <- glmmTMB(ECb ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Soil)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8)

# Best model
best_model <- m2

# Test fixed effects
car::Anova(best_model, type = "II")


# ---- Observation models ----

## Total arthropods (June/July) ##

# Candidate GLMs
m1 <- glm.nb(Total_arthropods ~ Cultivar,
             data = Observations_2)
m2 <- glm.nb(Total_arthropods ~ Treatment_worded,
             data = Observations_2)
m3 <- glm.nb(Total_arthropods ~ Cultivar + Treatment_worded,
             data = Observations_2)
m4 <- glm.nb(Total_arthropods ~ Cultivar * Treatment_worded,
             data = Observations_2)

# Candidate GLMMs
m5 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Block),
              family = nbinom2, data = Observations_2)
m6 <- glmmTMB(Total_arthropods ~ Treatment_worded + (1|Block),
              family = nbinom2, data = Observations_2)
m7 <- glmmTMB(Total_arthropods ~ Cultivar + Treatment_worded + (1|Block),
              family = nbinom2, data = Observations_2)
m8 <- glmmTMB(Total_arthropods ~ Cultivar * Treatment_worded + (1|Block),
              family = nbinom2, data = Observations_2)
m9 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Plant_ID),
              family = nbinom2, data = Observations_2)
m10 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Block) + (1|Plant_ID),
               family = nbinom2, data = Observations_2)
m11 <- glmmTMB(Total_arthropods ~ Treatment_worded + (1|Plant_ID),
               family = nbinom2, data = Observations_2)
m12 <- glmmTMB(Total_arthropods ~ Treatment_worded + (1|Block) + (1|Plant_ID),
               family = nbinom2, data = Observations_2)
m13 <- glmmTMB(Total_arthropods ~ Cultivar + Treatment_worded + (1|Plant_ID),
               family = nbinom2, data = Observations_2)
m14 <- glmmTMB(Total_arthropods ~ Cultivar + Treatment_worded + (1|Block) + (1|Plant_ID),
               family = nbinom2, data = Observations_2)
m15 <- glmmTMB(Total_arthropods ~ Cultivar * Treatment_worded + (1|Block) + (1|Plant_ID),
               family = nbinom2, data = Observations_2)

# Test the models using AIC
AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15)

# Best model
best_model <- m15

# Test fixed effects
car::Anova(best_model, type = "II")

# Post hoc comparisons (only if significant)
emmeans(best_model, pairwise ~ Cultivar, type = "response")
emmeans(best_model, pairwise ~ Cultivar | Treatment_worded, type = "response")


# ---- Combined / other models ----

# Combine the data frames into one data frame, using the common column "Plant_ID" to be able to work with all data in one data frame.
Combined_data <- Phenotype %>%
  #  left_join(Nectar, by = "Plant_ID") %>%
  left_join(Nectar_cleaned, by = "Plant_ID") %>%
  left_join(Flowers, by = "Plant_ID") %>%
  left_join(Repotting, by = "Plant_ID") %>%
  left_join(Observations, by = "Plant_ID") %>%
  left_join(Soil, by = "Plant_ID") %>%
  left_join(Flowering_date, by = "Plant_ID")





m1 <- glm.nb(Microliter ~ Number_Inflorescences, data = Combined_data)
m2 <- glm.nb(Microliter ~ Number_Inflorescences + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Microliter ~ Number_Inflorescences + Cultivar, data = Combined_data)
m4 <- glm.nb(Microliter ~ Number_Inflorescences + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Microliter ~ Average_Inflorescence_Length, data = Combined_data)
m6 <- glm.nb(Microliter ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
m7 <- glm.nb(Microliter ~ Average_Inflorescence_Length + Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m6 to be preferred, lowest AIC.

Model <- glm.nb(Microliter ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
car::Anova(Model , type = "III")





m1 <- glm.nb(Microliter ~ Total_arthropods, data = Combined_data)
m2 <- glm.nb(Microliter ~ Total_arthropods + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Microliter ~ Total_arthropods + Cultivar, data = Combined_data)
m4 <- glm.nb(Microliter ~ Total_arthropods + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Microliter ~ Total_arthropods + Cultivar * Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)

AIC(m1, m2, m3, m4, m5)
# Overall this shows m3 to be preferred, lowest AIC. But not significant.

Model <- glm.nb(Microliter ~ Total_arthropods + Cultivar, data = Combined_data)
car::Anova(Model , type = "III")




m1 <- glm.nb(Number_flowers ~ Cultivar, data = Combined_data)
m2 <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Number_flowers ~ Nodules_present, data = Combined_data)
m4 <- glm.nb(Number_flowers ~ Treatment_worded, data = Combined_data)
m5 <- glm.nb(Number_flowers ~ Root_abundance, data = Combined_data)
m6 <- glm.nb(Number_flowers ~ Seeds_present, data = Combined_data)
m7 <- glm.nb(Number_flowers ~ Seed_pod_abundance, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m2 to be preferred, lowest AIC.

Model <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Combined_data)
car::Anova(Model , type = "III")
































