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
Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

# Replace all NAs in the entire data frame with 0
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Observations)

# Make an exploratory graph
ggplot(Observations, aes(x = Time, y = Total_arthropods)) +
  geom_point() +
  labs(x = "Time",
       y = "Arthropods (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Arthropods_per_time_block.png", width = 8, height = 6, dpi = 300)

ggplot(Observations, aes(x = Block, y = Total_arthropods)) +
  geom_point() +
  labs(x = "Block",
       y = "Arthropods (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Arthropods_per__block.png", width = 8, height = 6, dpi = 300)

ggplot(Observations, aes(x = Cultivar, y = Total_arthropods)) +
  geom_point() +
  labs(x = "Cultivar",
       y = "Arthropods (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Arthropods_per_cultivar.png", width = 8, height = 6, dpi = 300)

ggplot(Observations, aes(x = Treatment_worded, y = Total_arthropods)) +
  geom_point() +
  labs(x = "Treatment",
       y = "Arthropods (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Arthropods_per_treatment.png", width = 8, height = 6, dpi = 300)



















