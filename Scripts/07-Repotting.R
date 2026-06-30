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
Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv") 

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Repotting)

# Make the yes and no into factor
Repotting$Nodules_present <- factor(
  Repotting$Nodules_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

Repotting$Nodule_abundance <- factor(
  Repotting$Nodule_abundance,
  levels = c("None", "Few", "Abundant")
)

Repotting$Nodule_shape <- factor(
  Repotting$Nodule_shape,
  levels = c("None", "Simple", "Both", "Complex")
)

Repotting$Seeds_present <- factor(
  Repotting$Seeds_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

# Make an exploratory graph
ggplot(Repotting, aes(x = Treatment_worded, fill = Nodules_present)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Nodule presence of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Nodule presence") +
  theme_minimal()
# Save
# ggsave("Graphs/Nodule_presence_per_cultivar_treatment.png", width = 8, height = 6, dpi = 300)

# Make the same but now with the seed 
ggplot(Repotting, aes(x = Treatment_worded, fill = Seeds_present)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Seed presence of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Seed presence") +
  theme_minimal()
# Save
# ggsave("Graphs/Seed_presence_per_cultivar_and_treatment.png", width = 8, height = 6, dpi = 300)

# Now make a graph with the root abundance 
ggplot(Repotting, aes(x = Treatment_worded, fill = Root_abundance)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Root abundance of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Root abundance") +
  theme_minimal()
# Save
# ggsave("Graphs/Root_abundance_per_cultivar_and_treatment.png", width = 8, height = 6, dpi = 300)


# Now make a plot for nodules to show in the presentation

# just the abundance of nodules between cultivars with the treatment on the x-axis and the abundance of nodules on the y-axis, with a fill color for the abundance of nodules.
ggplot(Repotting, aes(x = Treatment_worded, fill = Nodule_abundance)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(y = "Number of plants",
       x = "Treatment",
       fill = "Nodule abundance") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12)
  ) 
# ggsave("Graphs/Plants_w_Nodule_abundance.png", width = 12, height = 8, dpi = 300)

ggplot(Repotting, aes(x = Treatment_worded, fill = Nodule_shape)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(y = "Number of plants",
       x = "Treatment",
       fill = "Nodule complexity") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12)
  ) 
# ggsave("Graphs/Plants_w_Nodule_complexity.png", width = 12, height = 8, dpi = 300)







