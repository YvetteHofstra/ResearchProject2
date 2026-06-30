# Restore the correct versions of the used packages
renv::restore() 

# Clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# When a package is used for the first time, also add to lockfile
# renv::snapshot() # This shows a new library downloaded in manuscript
# Load the packages that are needed for this project
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately 
# load the required packages
library(broom)
library(chemodiv)
library(corrr)
library(cowplot)
library(DHARMa)
library(dplyr)
library(emmeans)
library(factoextra)
library(FSA)
library(ggeffects)
library(ggh4x)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggthemes)
library(glmmTMB)
library(grid)
library(lmtest)
library(lme4)
library(MASS)
library(Matrix)
library(multcomp)
library(multcompView)
library(patchwork)
library(performance)
library(pscl)
library(purrr)
library(RColorBrewer)
library(scales)
library(sjPlot)
library(readxl)
library(tibble)
library(tidyr)
library(tidyverse)
library(vegan)

readr::read_csv # This makes a tibble instead of table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Exploratory Analysis ----

# Load the data you want to use
# Load the data file from google drive 
Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTR4kOxATM525UmU7895FgLFgjHlL2RJ_Cgtb5fepWR-vRVZpwzLF3OIc4ZtvtTDQge1iUkyZY5W8Se/pub?gid=1854960858&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Soil)

# Make numeric where needed
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)

# Make an exploratory graph
ggplot(Soil, aes(x = Cultivar, y = ECp, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "ECp (mS⋅m⁻¹)",
       fill = "Treatment") +
  theme_minimal() 
nobs(Soil)
# Save
# ggsave("Graphs/Soil_salinity.png", width = 8, height = 6, dpi = 300)
# ggsave("Graphs/Soil_ECp_presentation.png",width = 11,height = 8,dpi = 600,units = "in")

ggplot(Soil, aes(x = Cultivar, y = ECb, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "ECb (mS⋅m⁻¹)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Soil_ECb.png", width = 8, height = 6, dpi = 300)
# ggsave("Graphs/Soil_ECb_presentation.png",width = 11,height = 8,dpi = 600,units = "in")

ggplot(Soil, aes(x = Cultivar, y = Eb, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Soil salinity Medicago sativa",
       x = "Cultivar",
       y = "Eb",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Soil_Eb.png", width = 8, height = 6, dpi = 300)

# Now one with temperature
ggplot(Soil, aes(x = Cultivar, y = Tmp, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Soil temperature Medicago sativa",
       x = "Cultivar",
       y = "Temperature",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Soil_temp.png", width = 8, height = 6, dpi = 300)

# Now one with wetness
ggplot(Soil, aes(x = Cultivar, y = Wet, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Wetness (%)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Soil_wet.png", width = 8, height = 6, dpi = 300)
# ggsave("Graphs/Soil_wet_presentation.png",width = 11,height = 8,dpi = 600,units = "in")



















