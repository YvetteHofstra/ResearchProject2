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
Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Phenotype)

# Make Number_inflorescences and Number_flowers numeric instead of integer
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)

# Make an exploratory graph
ggplot(Phenotype, aes(x = Cultivar, y = Number_inflorescences)) +
  geom_point() +
  labs(title = "Number of inflorescences per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Inflorescences.png", width = 8, height = 6, dpi = 300)

# Make the same but now with the amount of flowers
ggplot(Phenotype, aes(x = Cultivar, y = Number_flowers)) +
  geom_point() +
  labs(title = "Number of flowers per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Flowers.png", width = 8, height = 6, dpi = 300)

# Treatment with the amount of flowers
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_flowers)) +
  geom_point() +
  labs(title = "Number of flowers per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Treatment_with_flowers.png", width = 8, height = 6, dpi = 300)

# Treatment with the amount of inflorescences
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_inflorescences)) +
  geom_point() +
  labs(title = "Number of inflorescences per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Treatment_with_inflorescence.png", width = 8, height = 6, dpi = 300)










