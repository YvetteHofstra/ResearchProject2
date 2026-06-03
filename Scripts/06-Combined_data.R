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

readr::read_csv # This makes a tibble instead of table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Exploratory Analysis ----

# Load the data you want to use
# Load the data file from google drive 
Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

# Replace all NAs in the entire data frame with 0
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Phenotype)
str(Nectar)
str(Flowers)
str(Repotting)
str(Observations)

# Make numeric instead of integer
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)

# Combine the data frames into one data frame, using the common column "Plant_ID" to be able to work with all data in one data frame.
Combined_data <- Phenotype %>%
  left_join(Nectar, by = "Plant_ID") %>%
  left_join(Flowers, by = "Plant_ID") %>%
  left_join(Repotting, by = "Plant_ID") %>%
  left_join(Observations, by = "Plant_ID")

# Make an exploratory graph with nectar volume and number of inflorescences
ggplot(Combined_data, aes(x = Number_Inflorescences, y = Filled_until_mm)) +
  geom_point() +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_inflorescences_pointplot.png", width = 8, height = 6, dpi = 300)

# Now make it more concrete by adding treatment and cultivar
ggplot(Combined_data, aes(x = Number_Inflorescences, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Inflorescence (#)",
       color = "Treatment") +
  theme_minimal() +
  facet_wrap(~ Cultivar)
# Save
# ggsave("Graphs/Nectar_and_inflorescences_pointplot_treatment_cultivar_but_messy.png", width = 8, height = 6, dpi = 300)

# Try the point plot but now with a line
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_Inflorescences, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa",
       y = "Nectar volume (mm)",
       x = "Inflorescence (#)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_scatterplot_no_line_treatment.png", width = 8, height = 6, dpi = 300)

# Separate line for each treatment with se
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_Inflorescences, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Inflorescence (#)",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_inflorescence_scatterplot_line_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Separate line for each treatment without se
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_Inflorescences, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Number of inflorescences",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_inflorescences_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Do the same for # flowers
ggplot(Combined_data, aes(x = Number_flowers, y = Filled_until_mm)) +
  geom_point() +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_flowers_pointplot.png", width = 8, height = 6, dpi = 300)

# Now make it more concrete by adding treatment and cultivar
ggplot(Combined_data, aes(x = Number_flowers, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Flowers (#)",
       color = "Treatment") +
  theme_minimal() +
  facet_wrap(~ Cultivar)
# Save
# ggsave("Graphs/Nectar_and_flowers_pointplot_treatment_cultivar_but_messy.png", width = 8, height = 6, dpi = 300)

# Boxplot try, but does not make sense since both are counted and axes do not look nice when trying to fit it into a boxplot
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_flowers, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Medicago sativa",
       y = "Nectar volume (mm)",
       x = "Flowers (#)",
       fill = "Treatment") +
  theme_minimal() +
  facet_wrap(~ Cultivar)
# Save
# ggsave("Graphs/Nectar_and_flowers_boxplot_treatment_cultivar.png", width = 8, height = 6, dpi = 300)

# Try the point plot but now with a line
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_flowers, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa",
       y = "Nectar volume (mm)",
       x = "Flowers (#)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_flowers_scatterplot_no_line_treatment.png", width = 8, height = 6, dpi = 300)

# Separate line for each treatment with se
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_flowers, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Number of flowers",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_flowers_scatterplot_line_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Separate line for each treatment without se
ggplot(Combined_data, aes(x = Filled_until_mm, y = Number_flowers, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Number of flowers",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_flowers_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)




















