# Restore the correct versions of the used packages
renv::restore() 

# Clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# When a package is used for the first time, also add to lockfile
# renv::snapshot() # This shows a new library downloaded in manuscript
# Load the packages that are needed for this project
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately 
# load the required packages
library(ggplot2)
library(dplyr)
library(chemodiv)
library(corrr)
library(factoextra)
library(ggpubr)
library(purrr)
library(tibble)
library(tidyr)
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

Flowers_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=166596564&single=true&output=csv")

Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1067776784&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

Observations_2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=2034963164&single=true&output=csv")

Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=224800747&single=true&output=csv") 

## Replace all NAs in the entire data frame with 0 ##
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

## Show the type of data R will treat it as (numeric, character, factor, etc.) ##
str(Phenotype)
str(Nectar)
str(Flowers)
str(Flowers_2)
str(Flowering_date)
str(Repotting)
str(Observations)
str(Soil)

## Make numeric where needed ##
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)
Flowers_2$Number_Inflorescences <- as.numeric(Flowers_2$Number_Inflorescences)
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)

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


## Correct the observations (June/July) ##

# Filter the data to only include flowering plants and replace NA with 0 for arthropod counts after filtering flowering plants. This is important because we only want to consider arthropod visits to flowering plants, and any NA values in the arthropod counts should be treated as zero visits.
Observations_clean <- Observations_2 %>%
  filter(Flowering == 1)

arth_cols <- c(
  "Aglais_io", "Andrena", "Anthomyiidae", "Apis_mellifera",
  "Autographa_gamma", "Bibio_marci", "Bombus_lapidarius",
  "Bombus_pascuorum", "Bombus_terrestris", "Celastrina_argiolus",
  "Closterotomus_norwegicus", "Cynomya_mortuorum",
  "Empis_tessellata", "Larinioides_cornutus", "Oedemera",
  "Pieris_rapae", "Pollenia", "Polyommatus_icarus",
  "Pseudovadonia_livida", "Rhagonycha_fulva", "Sarcophagidae",
  "Stenichneumon_culpator", "Syrphidae", "Tetanocera",
  "Thereva_nobilitata", "Thymelicus_lineola",
  "Vanessa_atalanta", "Zygaena_filipendulae"
)

Observations_clean <- Observations_clean %>%
  mutate(across(all_of(arth_cols), ~ as.numeric(as.character(.)))) %>%
  mutate(across(all_of(arth_cols), ~ replace_na(., 0)))

# Total visits per plant
Plant_arthropods <- Observations_clean %>%
  group_by(Plant_ID) %>%
  summarise(
    Total_arthropods = sum(Total_arthropods, na.rm = TRUE),
    .groups = "drop"
  )

# Nectar dataset
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", 
                          "V71", "V73", "C89", "C85",
                          "C64", "C68", "C73", "C76", "C62")) %>%
  dplyr::select(Plant_ID, Microliter)

# Flower measurements
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

# Final combined dataset for plots/models
Combined_data2 <- Phenotype %>%
  dplyr::select(Plant_ID, Cultivar, Treatment_worded) %>%
  left_join(Nectar_cleaned, by = "Plant_ID") %>%
  left_join(Plant_arthropods, by = "Plant_ID") %>%
  left_join(Flowers_inside, by = "Plant_ID") %>%
  left_join(Flowers_outside, by = "Plant_ID")



### First combined graphs ### 

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

# Colour the treatment
ggplot(Combined_data, aes(x = Number_Inflorescences, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa",
       x = "Nectar volume (mm)",
       y = "Inflorescence (#)",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_inflorescences_pointplot_treatment.png", width = 8, height = 6, dpi = 600)

# Now make it more concrete by adding treatment and cultivar
ggplot(Combined_data, aes(x = Number_Inflorescences, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa inflorescence x nectar",
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
  labs(title = "Medicago sativa inflorescence x nectar",
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
  labs(title = "Medicago sativa inflorescence x nectar",
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
  labs(title = "Medicago sativa inflorescence x nectar",
       x = "Nectar volume (mm)",
       y = "Number of inflorescences",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_inflorescences_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Other way, now with microliter
ggplot(Combined_data, aes(x = Number_Inflorescences, y = Microliter, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Nectar (µL)",
       x = "Number of inflorescences",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Nectar_and_inflorescences_treatment_lines.png", width = 12, height = 8, dpi = 300)

# Other way around for the axes + now with microliter
ggplot(Combined_data, aes(x = Microliter, y = Number_Inflorescences, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Nectar (µL)",
       y = "Number of inflorescences",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Nectar_with_inflorescences_treatment_lines.png", width = 12, height = 8, dpi = 300)

# Remove plants that still could have more nectar extracted
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))

Combined_data_nectar <- Nectar_cleaned %>%
  left_join(Combined_data, by = "Plant_ID")

names(Combined_data_nectar)

ggplot(Combined_data_nectar, aes(x = Number_Inflorescences, y = Microliter.y, color = Treatment_worded.y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Nectar (µL)",
       x = "Number of inflorescences",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Filtered_nectar_with_inflorescences_treatment_lines.png", width = 12, height = 8, dpi = 300)

model <- lm(Microliter.y ~ Number_Inflorescences * Treatment_worded.y,
            data = Combined_data_nectar)
anova(model)

# Do it for # flowers
ggplot(Combined_data, aes(x = Number_flowers, y = Filled_until_mm)) +
  geom_point() +
  labs(title = "Medicago sativa flowers x nectar",
       x = "Nectar volume (mm)",
       y = "Flowers (#)") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_flowers_pointplot.png", width = 8, height = 6, dpi = 300)

# Now make it more concrete by adding treatment and cultivar
ggplot(Combined_data, aes(x = Number_flowers, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa flowers x nectar",
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
  labs(title = "Medicago sativa flowers x nectar",
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
  labs(title = "Medicago sativa flowers x nectar",
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
  labs(title = "Medicago sativa flowers x nectar",
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
  labs(title = "Medicago sativa flowers x nectar",
       x = "Nectar volume (mm)",
       y = "Number of flowers",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Nectar_and_flowers_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Try to combine inflorescence length and nectar (does length correlate to nectar?)
ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Treatment_worded)) +
  geom_point() +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_no_line_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Cultivar)) +
  geom_point() +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Cultivar") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_no_line_cultivar.png", width = 8, height = 6, dpi = 300)

# Try the point plot but now with a line
ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_line_SE_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_line__no_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Make the same ones, but than with cultivar
ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Treatment") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_line_SE_Cultivar.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Filled_until_mm, y = Average_Inflorescence_Length, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       x = "Nectar volume (mm)",
       y = "Inflorescence length (mm)",
       color = "Cultivar") +
  theme_minimal() 
# Save
# ggsave("Graphs/Nectar_and_inflorescence_length_scatterplot_line__no_SE_Cultivar.png", width = 8, height = 6, dpi = 300)

# Also see how it looks when nectar is on the y-axis
ggplot(Combined_data, aes(y = Filled_until_mm, x = Average_Inflorescence_Length, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       y = "Nectar volume (mm)",
       x = "Inflorescence length (mm)",
       color = "Treatment_worded") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescence_length_with_nectar_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(y = Filled_until_mm, x = Average_Inflorescence_Length, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence length x nectar",
       y = "Nectar volume (mm)",
       x = "Inflorescence length (mm)",
       color = "Cultivar") +
  theme_minimal() 
# Save
# ggsave("Graphs/Inflorescence_length_with_nectar_scatterplot_line_no_SE_Cultivar.png", width = 8, height = 6, dpi = 300)

# Try some models combining nectar and inflorescence # and inflorescence length
m1 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Cultivar * Treatment_worded, data = Combined_data)
m2 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Cultivar, data = Combined_data)
m4 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Filled_until_mm ~ Number_Inflorescences, data = Combined_data)
m6 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Cultivar, data = Combined_data)
m7 <- glm.nb(Filled_until_mm ~ Number_Inflorescences + Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m2 to be preferred, lowest AIC.Not significantly different from m1.

m1 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Cultivar * Treatment_worded, data = Combined_data)
m2 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
m4 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length, data = Combined_data)
m6 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
m7 <- glm.nb(Filled_until_mm ~ Average_Inflorescence_Length + Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m1 to be preferred, lowest AIC.

# Less logical as flowers were counted quite long before nectar was collected
m1 <- glm.nb(Filled_until_mm ~ Number_flowers + Cultivar * Treatment_worded, data = Combined_data)
m2 <- glm.nb(Filled_until_mm ~ Number_flowers + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Filled_until_mm ~ Number_flowers + Cultivar, data = Combined_data)
m4 <- glm.nb(Filled_until_mm ~ Number_flowers + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Filled_until_mm ~ Number_flowers, data = Combined_data)
m6 <- glm.nb(Filled_until_mm ~ Number_flowers + Cultivar, data = Combined_data)
m7 <- glm.nb(Filled_until_mm ~ Number_flowers + Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m1 to be preferred, lowest AIC. But m1 and m2 are very similar in AIC.


# Nectar with the total arthropod visitors
ggplot(Combined_data, aes(x = Total_arthropods, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa arthropods x nectar",
       x = "Total arthropods",
       y = "Nectar volume (mm)",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_nectar_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

# Work with april round and arthropod total
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
names(Observations_arthropods)

arth_summary <- Observations_arthropods |>
  group_by(Cultivar, Visitor, Plant_ID) |>
  summarise(Count = sum(Count), .groups = "drop")

names(arth_summary)

Combined_data_visits <- arth_summary %>%
  left_join(Combined_data, by = "Plant_ID")
names(Combined_data_visits)

ggplot(Combined_data_visits, aes(x = Microliter, y = Total_arthropods, color = Treatment_worded.y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Nectar (µL)",
       y = "Total arthropod visits",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# Save
# ggsave("Graphs/Nectar_with_visitors_treatment_lines.png", width = 12, height = 8, dpi = 300)

ggplot(Combined_data_visits, aes(y = Microliter, x = Total_arthropods, color = Treatment_worded.y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Nectar (µL)",
       x = "Total arthropod visits",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# ggsave("Graphs/Nectar_and_visitors_treatment_lines.png", width = 12, height = 8, dpi = 300)

model_visits <- lm(Total_arthropods ~ Microliter * Treatment_worded,
                   data = Combined_data_visits)

anova(model_visits)

# Now visitation with inflorescences
ggplot(Combined_data_visits, aes(y = Number_Inflorescences, x = Total_arthropods, color = Treatment_worded.y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Number of inflorescences",
       x = "Total arthropod visits",
       color = "Treatment") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
# ggsave("Graphs/Inflorescences_and_visitors_treatment_lines.png", width = 12, height = 8, dpi = 300)

model_visits <- lm(Total_arthropods ~ Number_Inflorescences * Treatment_worded,
                   data = Combined_data_visits)

anova(model_visits)

# Also try to see cultivar differences
ggplot(Combined_data, aes(x = Total_arthropods, y = Filled_until_mm, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa arthropods x nectar",
       x = "Total arthropods",
       y = "Nectar volume (mm)",
       color = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_nectar_scatterplot_line_no_SE_cultivar.png", width = 8, height = 6, dpi = 300)

# Make the arthropods in long table format
# First inspect how it looks like
Combined_data[, c("Bibio_marci",
                  "Bombus_pascuorum",
                  "Empis_tessellata",
                  "Larinioides_cornutus")]

Combined_data_arthropods <- Combined_data |>
  pivot_longer(
    cols = c(Bibio_marci,
             Bombus_pascuorum,
             Empis_tessellata,
             Larinioides_cornutus),
    names_to = "Genus",
    values_to = "Count"
  )

table(Combined_data_arthropods$Genus)

arth_summary <- Combined_data_arthropods |>
  group_by(Cultivar, Genus) |>
  summarise(Count = sum(Count), .groups = "drop")

# Then make the plot to show total arthropods (and composition) per cultivar
ggplot(arth_summary,
       aes(x = Cultivar,
           y = Count,
           fill = Genus)) +
  geom_col() +
  labs(
    title = "Arthropod community composition by cultivar",
    y = "Total observations"
  ) +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_cultivar_pointplot_cultivar.png", width = 8, height = 6, dpi = 300)

# Try visitation with time block
ggplot(Combined_data, aes(x = Time, y = Total_arthropods, color = Treatment_worded)) +
  geom_point(size = 2) +
  labs(title = "Medicago sativa arthropod abundance",
       x = "Time block",
       y = "Total arthropods",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_time_block_scatterplot_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Time, y = Total_arthropods, color = Cultivar)) +
  geom_point(size = 2) +
  labs(title = "Medicago sativa arthropod abundance",
       x = "Time block",
       y = "Total arthropods",
       color = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_time_block_scatterplot_cultivar.png", width = 8, height = 6, dpi = 300)






# Not yet working

arth_time <- Combined_data_arthropods |>
  group_by(Time) |>
  summarise(Total = sum(Count), .groups = "drop")

ggplot(arth_time, aes(x = Time, y = Total_arthropods, color = Genus)) +
  geom_point(size = 2) +
  labs(title = "Medicago sativa",
       x = "Time block",
       y = "Total arthropods",
       color = "Genus") +
  theme_minimal()
# Save
# ggsave("Graphs/Total_arthropods_with_time_block_scatterplot_treatment.png", width = 8, height = 6, dpi = 300)











# Nectar with roots
ggplot(Combined_data, aes(x = Root_abundance, y = Filled_until_mm, color = Treatment_worded)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa roots x nectar",
       x = "Root abundance",
       y = "Nectar volume (mm)",
       color = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Root_abundance_with_nectar_scatterplot_line_no_SE_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Root_abundance, y = Filled_until_mm, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa roots x nectar",
       x = "Root abundance",
       y = "Nectar volume (mm)",
       color = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/Root_abundance_with_nectar_scatterplot_line_no_SE_cultivar.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Root_abundance, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Nectar volume by root abundance",
       x = "Root abundance",
       y = "Nectar volume (mm)",
       fill = "Treatment") +
  theme_minimal()
# Save
# ggsave("Graphs/Root_abundance_with_nectar_boxplot_treatment.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_data, aes(x = Root_abundance, y = Filled_until_mm, fill = Cultivar)) +
  geom_boxplot() +
  labs(title = "Nectar volume by root abundance",
       x = "Root abundance",
       y = "Nectar volume (mm)",
       fill = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/Root_abundance_with_nectar_boxplot_cultivar.png", width = 8, height = 6, dpi = 300)

# Now the inflorescences (the ones measured closest to plants outside) and visitors
obs_summary <- Observations |>
  group_by(Plant_ID) |>
  summarise(
    Total_arthropods = sum(Total_arthropods),
    .groups = "drop"
  )

Combined_flower_obs <- Flowers |>
  left_join(obs_summary, by = "Plant_ID")

ggplot(Combined_flower_obs, aes(x = Total_arthropods, y = Number_Inflorescences, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence x arthropod",
       y = "Inflorescence (#)",
       x = "Total arthropods",
       color = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/Inflorescence_number_arthropod_total_cultivars.png", width = 8, height = 6, dpi = 300)

ggplot(Combined_flower_obs, aes(y = Total_arthropods, x = Number_Inflorescences, color = Cultivar)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Medicago sativa inflorescence x arthropod",
       x = "Inflorescence (#)",
       y = "Total arthropods",
       color = "Cultivar") +
  theme_minimal()
# Save
# ggsave("Graphs/arthropod_total_Inflorescence_number_cultivars.png", width = 8, height = 6, dpi = 300)


# Combine the data frames into one data frame, using the common column "Plant_ID" to be able to work with all data in one data frame.
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))

Combined <- Phenotype %>%
  #  left_join(Nectar, by = "Plant_ID") %>%
  left_join(Nectar_cleaned, by = "Plant_ID") %>%
  left_join(Flowers, by = "Plant_ID") %>%
  left_join(Repotting, by = "Plant_ID") %>%
  left_join(Observations, by = "Plant_ID") %>%
  left_join(Soil, by = "Plant_ID") 

ggplot(Combined_data, aes(x = Root_abundance, y = Filled_until_mm, fill = Treatment_worded)) +
  geom_boxplot() +
  labs(title = "Nectar volume by root abundance",
       x = "Root abundance",
       y = "Nectar volume (microliter)",
       fill = "Treatment") +
  theme_minimal()



### Second try ###

########
# These combine the 3 variables needed: inflorescence, nectar, visits

ggplot(Combined_data2,
       aes(x = Number_Inflorescences_inside,
           y = Microliter,
           colour = Treatment_worded)) +
  geom_point(position = position_jitterdodge(
    jitter.width = 0.15,
    dodge.width = 0.75),
    size = 2) +
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
# Save
# ggsave("Graphs/Combined_plots_Round_2_Infl_Nectar.png", width = 6, height = 4, dpi = 300)


ggplot(Combined_data2,
       aes(x = Number_Inflorescences_outside,
           y = Total_arthropods,
           colour = Treatment_worded)) +
  geom_point(position = position_jitterdodge(
    jitter.width = 0.15,
    dodge.width = 0.75),
    size = 2) +
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
# Save
# ggsave("Graphs/Combined_plots_Round_2_Infl_visits.png", width = 6, height = 4, dpi = 300)

ggplot(Combined_data2,
       aes(x = Microliter,
           y = Total_arthropods,
           colour = Treatment_worded)) +
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
# Save
# ggsave("Graphs/Combined_plots_Round_2_Nectar_visits.png", width = 6, height = 4, dpi = 300)






































