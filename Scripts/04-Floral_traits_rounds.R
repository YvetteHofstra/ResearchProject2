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
Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Flowers)

# Make Number_inflorescences numeric instead of integer
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)

# Make an exploratory graph
ggplot(Flowers, aes(x = Cultivar, y = Number_Inflorescences)) +
  geom_point() +
  labs(title = "Number of inflorescences per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescences_rounds.png", width = 8, height = 6, dpi = 300)

# Treatment with the amount of inflorescences
ggplot(Flowers, aes(x = Treatment_worded, y = Number_Inflorescences)) +
  geom_point() +
  labs(title = "Number of inflorescences per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_inflorescence_rounds.png", width = 8, height = 6, dpi = 300)

# Inflorescence per treatment and colored by cultivar
ggplot(Flowers, aes(x = Treatment_worded, y = Number_Inflorescences, color = Cultivar)) +
  geom_point() +
  labs(title = "Number of inflorescences per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_inflorescences_and_cultivar_rounds.png", width = 8, height = 6, dpi = 300)

# Now add the average inflorescence length per plant (so 36 still)
ggplot(Flowers, aes(x = Treatment_worded, y = Average_Inflorescence_Length, color = Cultivar)) +
  geom_point() +
  labs(title = "Length of inflorescences per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence length (mm)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_inflorescence_length_and_cultivar_rounds.png", width = 8, height = 6, dpi = 300)

# Now make a boxplot with this information, first try without treatment separation
ggplot(Flowers, aes(x = Cultivar, y = Average_Inflorescence_Length, fill = Cultivar)) +
  geom_boxplot() +
  labs(title = "Inflorescence length per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Inflorescence length (mm)") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescence_length_boxplot_rounds.png", width = 8, height = 6, dpi = 300)

# Now add the treatment
ggplot(Flowers, aes(x = Treatment_worded, y = Average_Inflorescence_Length, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Inflorescence length per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence length (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescence_length_treatment_boxplot_rounds.png", width = 8, height = 6, dpi = 300)

# Now divide the treatment to also show cultivars
ggplot(Flowers, aes(x = Treatment_worded, y = Average_Inflorescence_Length, fill = Treatment_worded)) +
  geom_boxplot() +
  facet_wrap(~ Cultivar) +
  labs(title = "Inflorescence length of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence length (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescence_length_treatment_boxplot_cultivar_rounds.png", width = 8, height = 6, dpi = 300)

# After the first round, the flowering dat was noted for each plant, add that data sheet and make some exploratory graphs
Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Flowering_date)

# Make numeric 
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Flowering_date$Flowered <- as.factor(Flowering_date$Flowered)

ggplot(Flowering_date, aes(x = Treatment_worded, y = Flowered, fill = Treatment_worded)) +
  geom_boxplot() +
  facet_wrap(~ Cultivar) +
  labs(title = "Flowering success of Medicago sativa",
       x = "Treatment",
       y = "Flowered",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Flowering_success.png", width = 8, height = 6, dpi = 300)

ggplot(Flowering_date, aes(x = Treatment_worded, y = Date_numbered, fill = Treatment_worded)) +
  geom_boxplot() +
  facet_wrap(~ Cultivar) +
  labs(title = "Flowering timing of Medicago sativa",
       x = "Treatment",
       y = "Flowering since",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Flowering_date_numbered_not_date_itself.png", width = 8, height = 6, dpi = 300)


























