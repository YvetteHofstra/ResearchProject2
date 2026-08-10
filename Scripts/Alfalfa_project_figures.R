# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#               Research Project 2, Data Analysis
#       The impact of salt tolerance on pollinator behaviour:
#      effects of salinity on floral traits of Medicago sativa
#                         Year: 2026
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- Before staring ----

# Clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# Restore the correct versions of the used packages
renv::restore() 

# Load the packages that are needed for this project
library(car)
library(chemodiv)
library(corrr)
library(emmeans)
library(factoextra)
library(glmmTMB)
library(gtsummary)
library(MASS)
library(multcomp)
library(nnet)
library(ordinal)
library(patchwork)
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately

readr::read_csv # This makes a tibble instead of a table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Loading the data ----

Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

Flowers_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=166596564&single=true&output=csv")

Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

Observations_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=2034963164&single=true&output=csv")

Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=224800747&single=true&output=csv") 

# We also need inflorescence abundance for April (before nectar collection) and June (before observations)
April <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

June <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=166596564&single=true&output=csv") 

## Replace all NAs in the entire data frame with 0 ##
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

## Make numeric where needed ##
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)
Flowers_2$Number_Inflorescences <- as.numeric(Flowers_2$Number_Inflorescences)
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)
April$Number_Inflorescences <- as.numeric(April$Number_Inflorescences)
June$Number_Inflorescences <- as.numeric(June$Number_Inflorescences)

## Make the yes and no into factor ##
Repotting$Nodules_present <- factor(
  Repotting$Nodules_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

Repotting$Seeds_present <- factor(
  Repotting$Seeds_present,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

Repotting$Nodule_abundance  <- as.factor(Repotting$Nodule_abundance)
Repotting$Nodule_shape <- as.factor(Repotting$Nodule_shape)
Repotting$Nodule_size  <- as.factor(Repotting$Nodule_size)
Repotting$Root_abundance <- as.factor(Repotting$Root_abundance)

## Give the factors the correct order and names for the figure ##
Repotting$Nodule_abundance <- factor(
  Repotting$Nodule_abundance,
  levels = c("None", "Few", "Abundant")
)

Repotting$Nodule_shape <- factor(
  Repotting$Nodule_shape,
  levels = c("None", "Simple", "Both", "Complex")
)

Repotting$Nodule_size <- factor(
  Repotting$Nodule_size,
  levels = c("None", "Small", "Both", "Big")
)

## Make a 'clean' nectar version, filter out not fully harvested plants
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", 
                          "V71", "V73", "C89", "C85",
                          "C64", "C68", "C73", "C76", "C62")) %>%
  dplyr::select(Plant_ID, Microliter, Treatment_worded, Cultivar)

## Flower measurements ##
Flowers_inside <- Flowers %>%
  dplyr::select(
    Plant_ID,
    Number_Inflorescences_inside = Number_Inflorescences
  )

Flowers_outside <- Flowers_2 %>%
  dplyr::select(
    Plant_ID,
    Number_Inflorescences_outside = Number_Inflorescences
  )

## Make a filter for soil measurement time point 1 and 2 ##
Soil_1 <- Soil %>% filter(Time_point == 1)
Soil_2 <- Soil %>% filter(Time_point == 2)

## Correct the observations (June/July) ##

# Filter the data to only include flowering plants and replace NA with 0 for arthropod counts after filtering flowering plants. This is important because we only want to consider arthropod visits to flowering plants, and any NA values in the arthropod counts should be treated as zero visits.
Observations_clean <- Observations_2 %>%
  filter(Flowering == 1)

# All arthropod taxa observed and given a name to make it easier to work with
arth_cols <- c(
  "Aglais_io", 	"Andrena", "Anthomyiidae",	"Apis_mellifera", "Autographa_gamma",	"Bibio_marci",	"Bombus_lapidarius", "Bombus_pascuorum", "Bombus_terrestris", "Celastrina_argiolus", "Closterotomus_norwegicus", "Cynomya_mortuorum", "Empis_tessellata", "Larinioides_cornutus", "Oedemera",	"Pieris_rapae", "Pollenia",	"Polyommatus_icarus", "Pseudovadonia_livida",	"Rhagonycha_fulva",	"Sarcophagidae",	"Stenichneumon_culpator",	"Syrphidae",	"Tetanocera", "Thereva_nobilitata",	"Thymelicus_lineola",	"Vanessa_atalanta", "Zygaena_filipendulae"
)

# Use the new name for all taxa to make all numeric
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ as.numeric(as.character(.))))

# Replace NA with 0 
Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ replace_na(., 0)))

Plant_arthropods <- Observations_clean %>%
  group_by(Plant_ID, Cultivar, Treatment_worded) %>%
  summarise(
    Total_arthropods = sum(Total_arthropods, na.rm = TRUE),
    .groups = "drop"
  )

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

# Combine to some groups, e.g. the bombus and butterfly ones to a "bee" and "butterfly" group, respectively, to make it easier to read the plot
Observations_long <- Observations_long %>%
  mutate(
    Group = case_when(
      Visitor %in% c("Andrena",
                     "Apis_mellifera",
                     "Bombus_lapidarius",
                     "Bombus_pascuorum",
                     "Bombus_terrestris") ~ "Bees",
      
      Visitor == "Syrphidae" ~ "Hoverflies",
      
      Visitor %in% c("Anthomyiidae",
                     "Bibio_marci",
                     "Cynomya_mortuorum",
                     "Empis_tessellata",
                     "Pollenia",
                     "Sarcophagidae",
                     "Tetanocera",
                     "Thereva_nobilitata") ~ "Non-syrphid flies",
      
      Visitor %in% c("Aglais_io",
                     "Autographa_gamma",
                     "Celastrina_argiolus",
                     "Pieris_rapae",
                     "Polyommatus_icarus",
                     "Thymelicus_lineola",
                     "Vanessa_atalanta",
                     "Zygaena_filipendulae") ~ "Butterflies & moths",
      
      Visitor %in% c("Oedemera",
                     "Pseudovadonia_livida",
                     "Rhagonycha_fulva") ~ "Beetles",
      
      Visitor == "Closterotomus_norwegicus" ~ "True bugs",
      
      Visitor == "Stenichneumon_culpator" ~ "Wasps",
      
      Visitor == "Larinioides_cornutus" ~ "Spiders",
      
      TRUE ~ "Other"
    )
  )

arth_groups <- Observations_long %>%
  group_by(Cultivar, Treatment_worded, Group) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

# Add a taxa summary to identify if diversity is influenced by cultivar or treatment
Orders_per_plant <- Observations_long %>%
  group_by(Plant_ID,Cultivar, Treatment_worded,Visitor,Block) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Plant_ID, Cultivar, Treatment_worded,Block) %>%
  summarise(NumTaxa = sum(Count > 0), .groups = "drop")

# Final combined dataset for plots
Combined_data2 <- Phenotype %>%
  dplyr::select(Plant_ID, Cultivar, Treatment_worded) %>%
  left_join(Nectar_cleaned, by = c("Plant_ID")) %>%
  left_join(Plant_arthropods, by = "Plant_ID") %>%
  left_join(Flowers_inside, by = "Plant_ID") %>%
  left_join(Flowers_outside, by = "Plant_ID") %>%
  mutate(
    Total_arthropods = replace_na(Total_arthropods, 0)
  )

Combined_data2 %>%
  summarise(
    NA_nectar = sum(is.na(Microliter)),
    NA_arthropods = sum(is.na(Total_arthropods)),
    NA_inside = sum(is.na(Number_Inflorescences_inside)),
    NA_outside = sum(is.na(Number_Inflorescences_outside))
  )

# For the combined plots, only plants with an actual nectar measurement AND a valid measurement
Combined_nectar <- Combined_data2 %>%
  filter(!is.na(Microliter),
         !is.na(Treatment_worded),
         Treatment_worded != "")


# ---- Soil measurements ----

## ECp ##

# First measuring moment
ggplot(Soil_1, aes(x = Cultivar, y = ECp, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "ECp (mS⋅m⁻¹)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# Second measuting moment
ggplot(Soil_2, aes(x = Cultivar, y = ECp, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "ECp (mS⋅m⁻¹)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Eb ##

ggplot(Soil, aes(x = Cultivar, y = Eb, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Soil salinity Medicago sativa",
       x = "Cultivar",
       y = "Eb",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## ECb ##

ggplot(Soil_1, aes(x = Cultivar, y = ECb, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "ECb (mS⋅m⁻¹)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Wetness ##

ggplot(Soil, aes(x = Cultivar, y = Wet, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Wetness (%)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Temperature ##

ggplot(Soil, aes(x = Cultivar, y = Tmp, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Soil temperature Medicago sativa",
       x = "Cultivar",
       y = "Temperature",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )


# ---- Vegetative traits ---- 

## Nodule abundance ##

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

## Nodule complexity ##

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

## Nodule size ## 

ggplot(Repotting, aes(x = Treatment_worded, fill = Nodule_size)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(y = "Number of plants",
       x = "Treatment",
       fill = "Nodule size") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12)
  ) 

## Root abundance ##

ggplot(Repotting, aes(x = Treatment_worded, fill = Root_abundance)) +
  geom_bar() +
  facet_wrap(~ Cultivar) +
  labs(x = "Treatment",
       y = "Number of plants",
       fill = "Root abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12)
  ) 

## Stem colour ##

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
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )


# ---- Floral traits ----

## Number of flowers ##

ggplot(Phenotype, aes(x = Cultivar, y = Number_flowers, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Number of flowers",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Number of inflorescences ##

# April
ggplot(April, aes(x = Cultivar, y = Number_Inflorescences, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Number of inflorescences",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# June
ggplot(June, aes(x = Cultivar, y = Number_Inflorescences, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Number of inflorescences",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# Combined 
Inflorescence_April <- ggplot(April, aes(x = Cultivar, y = Number_Inflorescences, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Number of inflorescences",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

Inflorescence_June <- ggplot(June, aes(x = Cultivar, y = Number_Inflorescences, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Number of inflorescences",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

combined_plot <- Inflorescence_April + Inflorescence_June + 
  plot_layout(ncol = 2, widths = c(1, 1))

## Inflorescence size ##

# April
ggplot(Flowers, aes(x = Cultivar, y = Average_Inflorescence_Length, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Inflorescence length (mm)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# June
ggplot(Flowers_2, aes(x = Cultivar, y = Average_Inflorescence_Length, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(x = "Cultivar",
       y = "Inflorescence length (mm)",
       fill = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Nectar production ##

# All plants
ggplot(Nectar, aes(x = Cultivar, y = Microliter, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.75), size = 2) +
  labs(x = "Cultivar",
       y = "Nectar (µL)",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# Only fully harvested
ggplot(Nectar_cleaned, aes(x = Cultivar, y = Microliter, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black",
              position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.75), size = 2) +
  labs(x = "Cultivar",
       y = "Nectar (µL)",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )


# ---- Reproductive traits ---- 

## Seed abundance ##

ggplot(Nectar, aes(x = Cultivar, y = Seed_abundance, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Treatment_worded),
              color = "black",
              position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.75), size = 2) +
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

## Seed weight ##

ggplot(Nectar, aes(x = Cultivar, y = Seed_weight, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Treatment_worded),
              color = "black",
              position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.75), size = 2) +
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


# ---- Observations ----

## Arthropod abundance ##

ggplot(Plant_arthropods, aes(x = Cultivar, y = Total_arthropods, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.75), size = 2) +
  labs(x = "Cultivar",
       y = "Total arthropod visits per plant",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Number of taxa ##

# Absolute
ggplot(Orders_per_plant, aes(x = Cultivar, y = NumTaxa, fill = Treatment_worded)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.75), size = 2) +
  labs(x = "Cultivar",
       y = "Number of arthropod taxa",
       fill = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

# Relative (grouped)
ggplot(arth_groups,
       aes(x = Treatment_worded, y = Count, fill = Group)) +
  geom_col(position = "fill") +
  facet_wrap(~Cultivar) +
  labs(
    x = "Treatment",
    y = "Relative visitor composition",
    fill = "Visitor group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(
    values = c(
      "Bees" = "#E69F00",
      "Hoverflies" = "#D55E00",
      "Non-syrphid flies" = "#56B4E9",
      "Butterflies & moths" = "#CC79A7",
      "Beetles" = "#009E73",
      "True bugs" = "#7A9A01",
      "Wasps" = "#A6761D",
      "Spiders" = "#7F7F7F"
    )
  )

# Relative (ungrouped)
ggplot(arth_groups,
       aes(x = Treatment_worded, y = Count, fill = Group)) +
  geom_col(position = "fill") +
  facet_wrap(~Cultivar) +
  labs(
    x = "Treatment",
    y = "Relative visitor composition",
    fill = "Visitor group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) 


# ---- Combined analyses ----

## Nectar volume with inflorescence abundance (April) ##

ggplot(Combined_nectar, aes(x = Number_Inflorescences_inside, y = Microliter, colour = Treatment_worded)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.75), size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Number of inflorescences",
       y = "Nectar (µL)",
       colour = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Inflorescence abundance (June) with arthropod abundance ##

ggplot(Combined_data2, aes(x = Number_Inflorescences_outside, y = Total_arthropods, colour = Treatment_worded)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.75), size = 2) +
  geom_smooth(method = "glm",
              method.args = list(family = "poisson"),
              se = FALSE) +
  labs(x = "Number of inflorescences",
       y = "Total arthropod visits",
       colour = "Treatment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

## Nectar volume with arthropod abundance ##

ggplot(Combined_nectar, aes(x = Microliter, y = Total_arthropods, colour = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "glm",
              method.args = list(family = "poisson"),
              se = FALSE) +
  labs(x = "Nectar (µL)",
       y = "Total arthropod visits",
       colour = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

