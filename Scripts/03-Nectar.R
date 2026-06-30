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
library(lme4)
library(lmerTest)

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
  labs(title = "Nectar per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_boxplot.png", width = 8, height = 6, dpi = 300)

# Make one with the microliter
ggplot(Nectar, aes(x = Cultivar, y = Microliter, fill = Cultivar)) +
  geom_boxplot() +
  labs(title = "Nectar per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (microliter)") +
  theme_minimal()

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

# The same but with microliter
ggplot(Nectar, aes(x = Treatment_worded, y = Microliter, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Treatment",
       y = "Nectar (microliter)",
       fill = "Treatment") +
  theme_minimal()

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

# Now combine the data points with a boxplot
ggplot(Nectar,
       aes(x = Cultivar,
           y = Filled_until_mm,
           fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75
              ),
              size = 2) +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_length_cultivar_boxplot_treatment_no_facet_points.png", width = 8, height = 6, dpi = 300)

# Make the same one but with microliter
ggplot(Nectar,
       aes(x = Cultivar,
           y = Microliter,
           fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75
              ),
              size = 2) +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (µL)",
       fill = "Treatment") +
  theme_minimal()

# Now remove the plants we did not fully harvest (as we had enough) 
# This means exclude A76, V64, V67, V68, V71, V73, C89, C85, C64, C68, C73, C76, C62
# Then we can make a plot to show: per plant x nectar, keep control and treatment
ggplot(Nectar, aes(x = Plant_ID, y = Filled_until_mm, fill = Plant_ID)) +
  geom_boxplot() +
  labs(title = "Nectar per plant of Medicago sativa",
       x = "Plants (#)",
       y = "Nectar (mm)") +
  theme_minimal()
# That is all plants

# Now delete the plants we did not fully harvest
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))
ggplot(Nectar_cleaned, aes(x = Plant_ID, y = Filled_until_mm, fill = Plant_ID)) +
  geom_boxplot() +
  labs(title = "Nectar per plant of Medicago sativa",
       x = "Plants (#)",
       y = "Nectar (mm)") +
  theme_minimal()

ggplot(Nectar_cleaned, aes(x = Cultivar, y = Microliter, fill = Cultivar)) +
  geom_boxplot() +
  labs(title = "Nectar per cultivar of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (microliter)") +
  theme_minimal()

# Check if there are significant differences
model <- glm.nb(Microliter ~ Cultivar,
                data = Nectar_cleaned)
anova(model)

model2 <- glm.nb(Microliter ~ Treatment_worded,
                 data = Nectar_cleaned)
anova(model2)

# Use the 'cleaned' to make the boxplot with treatment separated as well
ggplot(Nectar_cleaned,
       aes(x = Cultivar,
           y = Microliter,
           fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75
              ),
              size = 2) +
  labs(title = "Nectar per treatment of Medicago sativa",
       x = "Cultivar",
       y = "Nectar (microliter)",
       fill = "Treatment") +
  theme_minimal()

ggplot(Nectar_cleaned,aes(x = Cultivar, y = Microliter, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Nectar (µL)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Nectar_microliter_presentation.png", width = 12, height = 8, dpi = 300)

# Model
model <- lm(Microliter ~ Treatment_worded * Cultivar,
            data = Nectar_cleaned)

anova(model)
summary(model)


# Now we can work with seed abundance and weight rather than absence/presence of seed pods

# Not nice...
ggplot(Nectar, aes(fill = Treatment_worded, x = Seed_abundance)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(title = "Seed abundance of Medicago sativa",
       x = "Seed abundance",
       y = "Plants (#)",
       fill = "Treatment") +
  theme_minimal()

ggplot(Nectar,
       aes(x = Cultivar,
           y = Seed_abundance,
           fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Treatment_worded),
              position = position_jitterdodge(
                jitter.width = 0.1,
                dodge.width = 0.75
              ),
              size = 2) +
  labs(title = "Seed production of Medicago sativa",
       x = "Cultivar",
       y = "Seeds per plant",
       fill = "Treatment",
       color = "Treatment") +
  theme_minimal()

# Similar plot as the nectar volume one
ggplot(Nectar, aes(x = Cultivar,
                           y = Seed_abundance,
                           fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Treatment_worded),
              color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Number of seeds",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Seed_abundance_presentation.png", width = 12, height = 8, dpi = 300)

# Now with seed weight
ggplot(Nectar,
       aes(x = Seed_abundance,
           y = Seed_weight,
           color = Treatment_worded)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Cultivar) +
  labs(title = "Relationship between seed number and seed weight",
       x = "Seeds per plant",
       y = "Seed weight (g)") +
  theme_minimal()

ggplot(Nectar, aes(x = Cultivar,
                   y = Seed_weight,
                   fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Treatment_worded),
              color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Seed weight (g)",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Seed_weight_presentation.png", width = 12, height = 8, dpi = 300)

# With seed pods rather than seeds
ggplot(Nectar, aes(x = Cultivar,
                   y = Seed_pod_abundance,
                   fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Treatment_worded),
              color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Number of seed pods",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Seed_pods_plot.png", width = 12, height = 8, dpi = 300)














