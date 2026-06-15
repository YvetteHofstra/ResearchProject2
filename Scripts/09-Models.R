# Restore the correct versions of the used packages
renv::restore() 

# Clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# When a package is used for the first time, also add to lockfile
# renv::snapshot() # This shows a new library downloaded in manuscript
# Load the packages that are needed for this project
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately 
# load the required packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(chemodiv)
library(corrr)
library(factoextra)
library(ggpubr)
library(purrr)
library(tibble)
library(MASS)
library(multcomp)
library(emmeans)

readr::read_csv # This makes a tibble instead of a table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Loading the data ----

Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTR4kOxATM525UmU7895FgLFgjHlL2RJ_Cgtb5fepWR-vRVZpwzLF3OIc4ZtvtTDQge1iUkyZY5W8Se/pub?gid=1854960858&single=true&output=csv") 

# Replace all NAs in the entire data frame with 0
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Phenotype)
str(Nectar)
str(Flowers)
str(Flowering_date)
str(Repotting)
str(Observations)

# Make numeric instead of integer
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)

# Make the yes and no into factor
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


# ---- Phenotype models ----

# Is the number of flowers significantly different between cultivars
m1 <- glm.nb(Number_flowers ~ Cultivar * Treatment_worded, data = Phenotype)
m2 <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Phenotype)
m3 <- glm.nb(Number_flowers ~ Cultivar, data = Phenotype)
m4 <- glm.nb(Number_flowers ~ Treatment_worded, data = Phenotype)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Is the number of inflorescences significantly different between cultivars
m1 <- glm.nb(Number_inflorescences ~ Cultivar * Treatment_worded, data = Phenotype)
m2 <- glm.nb(Number_inflorescences ~ Cultivar + Treatment_worded, data = Phenotype)
m3 <- glm.nb(Number_inflorescences ~ Cultivar, data = Phenotype)
m4 <- glm.nb(Number_inflorescences ~ Treatment_worded, data = Phenotype)

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Nectar models ----

# Is the amount of nectar significantly different between cultivars
m1 <- glm.nb(Filled_until_mm ~ Cultivar * Treatment_worded, data = Nectar)
m2 <- glm.nb(Filled_until_mm ~ Cultivar + Treatment_worded, data = Nectar)
m3 <- glm.nb(Filled_until_mm ~ Cultivar, data = Nectar)
m4 <- glm.nb(Filled_until_mm ~ Treatment_worded, data = Nectar)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Now remove the plants that were not completely tried 
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))

m1 <- glm.nb(Filled_until_mm ~ Cultivar * Treatment_worded, data = Nectar_cleaned )
m2 <- glm.nb(Filled_until_mm ~ Cultivar + Treatment_worded, data = Nectar_cleaned )
m3 <- glm.nb(Filled_until_mm ~ Cultivar, data = Nectar_cleaned )
m4 <- glm.nb(Filled_until_mm ~ Treatment_worded, data = Nectar_cleaned )

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Flower (round) models ----

m1 <- glm.nb(Average_Inflorescence_Length ~ Cultivar * Treatment_worded, data = Flowers)
m2 <- glm.nb(Average_Inflorescence_Length ~ Cultivar + Treatment_worded, data = Flowers)
m3 <- glm.nb(Average_Inflorescence_Length ~ Cultivar, data = Flowers)
m4 <- glm.nb(Average_Inflorescence_Length ~ Treatment_worded, data = Flowers)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Now the number of inflorescences with the counting closest to nectar (could later be combined perhaps to figure out significance between nectar and number of inflorescences)
m1 <- glm.nb(Number_Inflorescences ~ Cultivar * Treatment_worded, data = Flowers)
m2 <- glm.nb(Number_Inflorescences ~ Cultivar + Treatment_worded, data = Flowers)
m3 <- glm.nb(Number_Inflorescences ~ Cultivar, data = Flowers)
m4 <- glm.nb(Number_Inflorescences ~ Treatment_worded, data = Flowers)

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Flower (date) models ----

m1 <- glm.nb(Date ~ Cultivar * Treatment_worded, data = Flowering_date)
m2 <- glm.nb(Date ~ Cultivar + Treatment_worded, data = Flowering_date)
m3 <- glm.nb(Date ~ Cultivar, data = Flowering_date)
m4 <- glm.nb(Date ~ Treatment_worded, data = Flowering_date)

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Repotting models ----

# Seed pod number
m1 <- glm.nb(Seed_pod_number ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_pod_number ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_pod_number ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_pod_number ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Seed pod weight
m1 <- glm.nb(Seed_pod_weight ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_pod_weight ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_pod_weight ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_pod_weight ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Seed number
m1 <- glm.nb(Seed_number ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_number ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_number ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_number ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Seed weight
m1 <- glm.nb(Seed_weight ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_weight ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_weight ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_weight ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Seed germination?
m1 <- glm.nb(Seed_germination ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_germination ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_germination ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_germination ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Soil models ----

# Wetness (%) of soil
m1 <- glm.nb(Wet ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(Wet ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(Wet ~ Cultivar, data = Soil)
m4 <- glm.nb(Wet ~ Treatment_worded, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# ECp
m1 <- glm.nb(ECp ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(ECp ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(ECp ~ Cultivar, data = Soil)
m4 <- glm.nb(ECp ~ Treatment_worded, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# Eb
m1 <- glm.nb(Eb ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(Eb ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(Eb ~ Cultivar, data = Soil)
m4 <- glm.nb(Eb ~ Treatment_worded, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

# ECb
m1 <- glm.nb(ECb ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(ECb ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(ECb ~ Cultivar, data = Soil)
m4 <- glm.nb(ECb ~ Treatment_worded, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)


# ---- Observation models ----

# Can the total arthropods (pollinators) be explained by cultivar or treatment? What about the block?
m1 <- glm.nb(Total_arthropods ~ Cultivar * Treatment_worded, data = Observations)
m2 <- glm.nb(Total_arthropods ~ Cultivar + Treatment_worded, data = Observations)
m3 <- glm.nb(Total_arthropods ~ Cultivar, data = Observations)
m4 <- glm.nb(Total_arthropods ~ Treatment_worded, data = Observations)
m5 <- glm.nb(Total_arthropods ~ Block, data = Observations)
m6 <- glm.nb(Total_arthropods ~ Block + Cultivar, data = Observations)
m7 <- glm.nb(Total_arthropods ~ Block + Treatment_worded, data = Observations)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)


# ---- Combined / other models ----
















