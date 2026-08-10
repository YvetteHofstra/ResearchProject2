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
library(nnet)
library(ordinal)
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

Repotting$Nodule_abundance  <- as.factor(Repotting$Nodule_abundance)
Repotting$Nodule_shape <- as.factor(Repotting$Nodule_shape)
Repotting$Nodule_size  <- as.factor(Repotting$Nodule_size)
Repotting$Root_abundance <- as.factor(Repotting$Root_abundance)


## Correct the observations (June/July) ##

# Filter the data to only include flowering plants and replace NA with 0 for arthropod counts after filtering flowering plants. This is important because we only want to consider arthropod visits to flowering plants, and any NA values in the arthropod counts should be treated as zero visits.
Observations_clean <- Observations_2 %>%
  filter(Flowering == 1)

# All arthropod taxa observed and given a name to make it easier to work with
arth_cols <- c(
  "Aglais_io", 	"Andrena", "Anthomyiidae",	"Apis_mellifera", "Autographa_gamma",	"Bibio_marci",	"Bombus_lapidarius", "Bombus_pascuorum", "Bombus_terrestris", "Celastrina_argiolus", "Closterotomus_norwegicus", "Cynomya_mortuorum", "Empis_tessellata", "Larinioides_cornutus", "Oedemera",	"Pieris_rapae", "Pollenia",	"Polyommatus_icarus", "Pseudovadonia_livida",	"Rhagonycha_fulva",	"Sarcophagidae",	"Stenichneumon_culpator",	"Syrphidae",	"Tetanocera", "Thereva_nobilitata",	"Thymelicus_lineola",	"Vanessa_atalanta", "Zygaena_filipendulae"
)

# Use the new name for all taxa to make all numeric
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ as.numeric(as.character(.))))

# Replace NA with 0 
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ replace_na(., 0)))

Plant_arthropods <- Observations_clean %>%
  group_by(Plant_ID, Cultivar, Treatment_worded) %>%
  summarise(
    Total_arthropods = sum(Total_arthropods, na.rm = TRUE),
    .groups = "drop"
  )

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

# Add a taxa summary to identify if diversity is influenced by cultivar or treatment
Orders_per_plant <- Observations_long %>%
  group_by(Plant_ID,Cultivar, Treatment_worded,Visitor,Block) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Plant_ID, Cultivar, Treatment_worded,Block) %>%
  summarise(NumTaxa = sum(Count > 0), .groups = "drop")


# ----  ----




