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
Observations_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=2034963164&single=true&output=csv")

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
# ggsave("Graphs/Arthropods_per_block.png", width = 8, height = 6, dpi = 300)

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

ggplot(Observations, aes(x = Cultivar, y = Total_arthropods)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Arthropods (#)") +
  theme_minimal()

ggplot(Observations, aes(x = Treatment_worded, y = Total_arthropods)) +
  geom_boxplot() +
  labs(x = "Treatment",
       y = "Arthropods (#)") +
  theme_minimal()

ggplot(Observations, aes(x = Time, y = Total_arthropods, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa pollinators per time",
       x = "Time",
       y = "Arthropods (#)",
       color = "Treatment") +
  theme_minimal() +
  facet_wrap(~ Cultivar)

ggplot(Observations, aes(x = Time, y = Total_arthropods, color = Cultivar)) +
  geom_point() +
  labs(title = "Medicago sativa pollinators per time",
       x = "Time",
       y = "Arthropods (#)",
       color = "Cultivar") +
  theme_minimal() 

ggplot(Observations, aes(x = Time, y = Total_arthropods, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa pollinators per time",
       x = "Time",
       y = "Arthropods (#)",
       color = "Treatment") +
  theme_minimal() 

# Try and make some nice plots with arthropods
Observations_arthropods <- Observations |>
  pivot_longer(
    cols = c(Bibio_marci,
             Bombus_pascuorum,
             Empis_tessellata,
             Larinioides_cornutus),
    names_to = "Visitor",
    values_to = "Count"
  )

table(Observations_arthropods$Visitor)

arth_summary <- Observations_arthropods |>
  group_by(Cultivar, Visitor) |>
  summarise(Count = sum(Count), .groups = "drop")

# Then make the plot to show total arthropods (and composition) per cultivar
ggplot(arth_summary,
       aes(x = Cultivar,
           y = Count,
           fill = Visitor)) +
  geom_col() +
  labs(y = "Total arthropod visits"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_fill_discrete(
    labels = c(
      expression(italic("Bibio marci")),
      expression(italic("Bombus pascuorum")),
      expression(italic("Empis tessellata")),
      expression(italic("Larinioides cornutus"))
    )
  )
# Save
# ggsave("Graphs/Total_arthropod_cultivar.png", width = 8, height = 6, dpi = 300)

# Not yet included the june observations correctly


# First just show the round 1 observations as that is complete for the 36 plants.
arth_plant <- Observations %>%
  group_by(Plant_ID, Cultivar, Treatment_worded) %>%
  summarise(
    Total_visits = sum(Total_arthropods, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(arth_plant, aes(x = Cultivar,
                       y = Total_visits,
                       fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Total arthropod visits",
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
# ggsave("Graphs/Arthropods_Round_one_presentation.png", width = 12, height = 8, dpi = 300)

Observations_clean <- Observations_2 %>%
  filter(Flowering == 1)

arth_cols <- c(
  "Apis_mellifera", "Bibio_marci", "Bombus_pascuorum",
  "Bombus_terrestris", "Closterotomus_norwegicus",
  "Empis_tessellata", "Larinioides_cornutus",
  "Oedemera", "Stenichneumon_culpator", "Syrphidae"
)

Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ as.numeric(as.character(.))))

# Replace NA with 0 (ONLY after filtering flowering plants)
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ replace_na(., 0)))

# Pivot longer
Observations_long <- Observations_clean %>%
  pivot_longer(
    cols = all_of(arth_cols),
    names_to = "Visitor",
    values_to = "Count"
  )

print(table(Observations_long$Count, useNA = "ifany"))

arth_summary <- Observations_long %>%
  group_by(Cultivar, Visitor, Treatment_worded) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

print(summary(arth_summary$Count))

ggplot(arth_summary,
       aes(x = Cultivar,
           y = Count,
           fill = Visitor)) +
  geom_col() +
  labs(x = "Cultivar",
    y = "Total arthropod visits"
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_fill_discrete(
    labels = c(
      expression(italic("Apis mellifera")),
      expression(italic("Bibio marci")),
      expression(italic("Bombus pascuorum")),
      expression(italic("Bombus terrestris")),
      expression(italic("Closterotomus norwegicus")),
      expression(italic("Empis tessellata")),
      expression(italic("Larinioides cornutus")),
      expression(italic("Oedemera")),
      expression(italic("Stenichneumon culpator")),
      "Syrphidae"
    )
  )
# Save
# ggsave("Graphs/Arthropods_Round_two_presentation.png", width = 12, height = 8, dpi = 300)

ggplot(arth_summary, aes(x = Cultivar,
                       y = Count,
                       fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.75),
              size = 2) +
  labs(x = "Cultivar",
       y = "Total arthropod visits",
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
# ggsave("Graphs/Arthropods_Round_two.png", width = 12, height = 8, dpi = 300)







