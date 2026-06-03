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

readr::read_csv # This makes a tibble instead of table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Exploratory Analysis ----

# Load the data you want to use
# Load the data file from google drive 
Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Nectar)

# Make an exploratory graph
ggplot(Nectar, aes(x = Cultivar, y = Filled_until_mm)) +
  geom_point() +
  labs(title = "Nectar per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length.png", width = 8, height = 6, dpi = 300)

# Make it into a boxplot
ggplot(Nectar, aes(x = Cultivar, y = Filled_until_mm, fill = Cultivar)) +
  geom_boxplot() +
  labs(title = "Number of inflorescences per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_boxplot.png", width = 8, height = 6, dpi = 300)

# Now add the treatment
ggplot(Nectar, aes(x = Treatment_worded, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Treatment",
       y = "Nectar (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_treatment_boxplot.png", width = 8, height = 6, dpi = 300)

# Now divide the treatment to also show cultivars
ggplot(Nectar, aes(x = Treatment_worded, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  facet_wrap(~ Cultivar) +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Treatment",
       y = "Nectar (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_treatment_boxplot_cultivar.png", width = 8, height = 6, dpi = 300)

ggplot(Nectar, aes(x = Cultivar, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  facet_wrap(~ Treatment_worded) +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_cultivar_boxplot_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Nectar, aes(x = Cultivar, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_cultivar_boxplot_treatment_no_facet.png", width = 8, height = 6, dpi = 300)














