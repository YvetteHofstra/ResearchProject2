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
# ggsave("Graphs/Inflorescences.png", width = 8, height = 6, dpi = 300)

# Make the same but now with the amount of flowers
ggplot(Phenotype, aes(x = Cultivar, y = Number_flowers)) +
  geom_point() +
  labs(title = "Number of flowers per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Flowers.png", width = 8, height = 6, dpi = 300)

# Treatment with the amount of flowers
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_flowers)) +
  geom_point() +
  labs(title = "Number of flowers per treatment of Medicago sativa",
       x = "Treatment",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_flowers.png", width = 8, height = 6, dpi = 300)

# Treatment with the amount of inflorescences
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_inflorescences)) +
  geom_point() +
  labs(title = "Number of inflorescences per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_inflorescence.png", width = 8, height = 6, dpi = 300)

# Inflorescence per treatment and colored by cultivar
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_inflorescences, color = Cultivar)) +
  geom_point() +
  labs(title = "Number of inflorescences per treatment of Medicago sativa",
       x = "Treatment",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_inflorescences_and_cultivar.png", width = 8, height = 6, dpi = 300)

# Flowers per treatment and colored by cultivar
ggplot(Phenotype, aes(x = Treatment_worded, y = Number_flowers, color = Cultivar)) +
  geom_point() +
  labs(title = "Number of flowers per treatment of Medicago sativa",
       x = "Treatment",
       y = "Flower (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Treatment_with_flowers_and_cultivar.png", width = 8, height = 6, dpi = 300)

# Try a more logical approach of cultivar with # flowers and colour by treatment
ggplot(Phenotype, aes(x = Cultivar, y = Number_flowers, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Number of flowers per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Flower (#)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Cultivar_with_flowers_divided_by_treatment.png", width = 8, height = 6, dpi = 300)

# Again but with inflorescences
ggplot(Phenotype, aes(x = Cultivar, y = Number_inflorescences, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Number of inflorescences",
       x = "Cultivar",
       y = "Inflorescence (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Cultivar_with_inflorescences_divided_by_treatment.png", width = 8, height = 6, dpi = 300)

# Try a boxplot
ggplot(Phenotype, aes(x = Cultivar, y = Number_flowers, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Number of flowers",
       x = "Cultivar",
       y = "Flower (#)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Cultivar_with_flowers_divided_by_treatment_boxplot.png", width = 8, height = 6, dpi = 300)

# Check if it makes sense to add significance to the plot 
model_flowers_nb <- glm.nb(Number_flowers ~ Cultivar * Treatment_worded,
                           data = Phenotype)
anova(model_flowers_nb)
# Only cultivar is significant, so do not add significance to previous plot
emm <- emmeans(model_flowers_nb, ~ Cultivar)
pairs(emm)
cld(emm)
# Flower number differed significantly among cultivars (negative binomial GLM, p = 0.022), but not between treatment

# Again but with inflorescences
ggplot(Phenotype, aes(x = Cultivar, y = Number_inflorescences, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Number of inflorescences",
       x = "Cultivar",
       y = "Inflorescence (#)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Cultivar_with_inflorescences_divided_by_treatment_boxplot.png", width = 8, height = 6, dpi = 300)

model_inflorescence_nb <- glm.nb(Number_inflorescences ~ Cultivar * Treatment_worded,
                                 data = Phenotype)
anova(model_inflorescence_nb)
# Only cultivar is significant, so do not add significance to previous plot
emm <- emmeans(model_inflorescence_nb, ~ Cultivar)
pairs(emm)
cld(emm)
# Inflorescence number differed significantly among cultivars (negative binomial GLM, p = < 0.05), but not between treatment

# Try to add the flower color in a plot. E.g. first the flower color per variety
ggplot(Phenotype, aes(x = Cultivar, fill = Flower_color_simple)) +
  geom_bar() +
  labs(title = "Flower color per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Count") +
  theme_minimal()
# Save
# ggsave("Graphs/Flower_color_per_cultivar.png", width = 8, height = 6, dpi = 300)

# Now make the colors correspond to the bars
ggplot(Phenotype, aes(x = Cultivar, fill = Flower_color_simple)) +
  geom_bar() +
  scale_fill_manual(values = c(
    "Yellow" = "yellow",
    "Purple" = "#9d07c5ff"
  )) +
  labs(title = "Flower color per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Count",
       fill = "Flower color") +
  theme_minimal()
# Save
# ggsave("Graphs/Flower_color_per_cultivar_with_actual_colors.png", width = 8, height = 6, dpi = 300)

# Now the treatment split to show the amount of plants per cultivar and per treatment which color the plant has
ggplot(Phenotype, aes(x = Cultivar, fill = Flower_color_simple)) +
  geom_bar() +
  facet_wrap(~ Treatment_worded) +
  scale_fill_manual(values = c(
    "Yellow" = "yellow",
    "Purple" = "#9d07c5ff"
  )) +
  labs(title = "Flower color per cultivar and treatment of Medicago sativa",
       x = "Cultivar",
       y = "Plants (#)",
       fill = "Flower color") +
  theme_minimal()

ggplot(Phenotype, aes(x = Treatment_worded, fill = Flower_color_simple)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  scale_fill_manual(values = c(
    "Yellow" = "yellow",
    "Purple" = "#9d07c5ff"
  )) +
  labs(title = "Flower color per cultivar of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Flower color") +
  theme_minimal()
# Save
# ggsave("Graphs/Flower_color_per_treatment_and_cultivar.png", width = 8, height = 6, dpi = 300)

# Now add growth direction of the plants to a plot
ggplot(Phenotype, aes(x = Treatment_worded, fill = Growth_direction)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Growth direction of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Growth direction") +
  theme_minimal()
# Save
# ggsave("Graphs/Growth_direction_per_cultivar_and_treatment.png", width = 8, height = 6, dpi = 300)

# Stem colour, change the colors to match the actual colors of the plants
ggplot(Phenotype, aes(x = Treatment_worded, fill = Stem_color)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  scale_fill_manual(values = c(
    "Green" = "#008300",
    "Intermediate" = "#DCB695",
    "Red" = "red"
  )) +
  labs(title = "Stem colour of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Stem colour") +
  theme_minimal()
# Save
# ggsave("Graphs/Stem_colour_per_cultivar_and_treatment.png", width = 8, height = 6, dpi = 300)

# Nodules, but change to have the yes (1) and no (0) as labels instead of continuous
Phenotype$Nodules <- factor(
  Phenotype$Nodules,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

ggplot(Phenotype, aes(x = Treatment_worded, fill = Nodules)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Nodule presence of Medicago sativa",
       x = "Treatment",
       y = "Plants (#)",
       fill = "Nodule presence") +
  theme_minimal()
# Save
# ggsave("Graphs/Nodule_presence_per_cultivar_and_treatment.png", width = 8, height = 6, dpi = 300)
























