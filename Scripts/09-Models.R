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
library(glmmTMB)
library(gtsummary)
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

readr::read_csv # This makes a tibble instead of a table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Loading the data ----

Phenotype <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=206224982&single=true&output=csv") 

Nectar <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv") 

Flowers <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSQgacYoLmN4V7eOLdZ4JMNue1B_q67kQHkpzGDWt3DCY8FyHVBW5Ml_TR2rViu7jViE_WXihuuZiRc/pub?gid=0&single=true&output=csv") 

Flowering_date <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=1460555223&single=true&output=csv") 

Repotting <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=217767211&single=true&output=csv")

Observations <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRNix7qqZS7cB-KkXmk4Yu7XNvI8uNFhS_ZCfTGwVIziLeXCzH-VlHzEzrndxrzLGgWUj-ssOHRmORV/pub?gid=1102638602&single=true&output=csv")

Soil <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrKk4lVr_GFFwaudVT_jG4tLL9LhCNixrmjzVfOHbsHk3y-3YA8C9dtlWfm4QyFoy9Xmhn2AQmr7SY/pub?gid=224800747&single=true&output=csv") 

# Replace all NAs in the entire data frame with 0
Observations <- Observations |>
  mutate_all(~ ifelse(is.na(.), 0, .))

# Show the type of data R will treat it as (numeric, character, factor, etc.)
str(Phenotype)
str(Nectar)
str(Flowers)
str(Flowering_date)
str(Repotting)
str(Observations)
str(Soil)

# Make numeric instead of integer
Phenotype$Number_inflorescences <- as.numeric(Phenotype$Number_inflorescences)
Phenotype$Number_flowers <- as.numeric(Phenotype$Number_flowers)
Flowers$Number_Inflorescences <- as.numeric(Flowers$Number_Inflorescences)
Flowering_date$Date_numbered <- as.numeric(Flowering_date$Date_numbered)
Soil$ECp <- as.numeric(Soil$ECp)
Soil$ECb <- as.numeric(Soil$ECb)

# Make the yes and no into factor
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


# ---- Phenotype models ----

# Is the number of flowers significantly different between cultivars
m1 <- glm.nb(Number_flowers ~ Cultivar * Treatment_worded, data = Phenotype)
m2 <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Phenotype)
m3 <- glm.nb(Number_flowers ~ Cultivar, data = Phenotype)
m4 <- glm.nb(Number_flowers ~ Treatment_worded, data = Phenotype)
m5 <- glmmTMB(Number_flowers ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)
m6 <- glmmTMB(Number_flowers ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)
m7 <- glmmTMB(Number_flowers ~ Cultivar + (1|Block), family = nbinom2, data = Phenotype)
m8 <- glmmTMB(Number_flowers ~ Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8) # m7 is best of these

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m3 to be preferred, lowest AIC. But not significant. Can just use the glm.nb as the random effect does not significantly improve the model.

Model <- glm.nb(Number_flowers ~ Cultivar, data = Phenotype)
car::Anova(Model , type = "III")
ModelglmmTMB <- glmmTMB(Number_flowers ~ Cultivar, family = nbinom2, data = Phenotype)
car::Anova(ModelglmmTMB , type = "III")

# Display the results with significant effects highlighted
tbl_regression(Model) %>%
  bold_p()
tbl_regression(ModelglmmTMB) %>%
  bold_p()

# Is the number of inflorescences significantly different between cultivars
m1 <- glm.nb(Number_inflorescences ~ Cultivar * Treatment_worded, data = Phenotype)
m2 <- glm.nb(Number_inflorescences ~ Cultivar + Treatment_worded, data = Phenotype)
m3 <- glm.nb(Number_inflorescences ~ Cultivar, data = Phenotype)
m4 <- glm.nb(Number_inflorescences ~ Treatment_worded, data = Phenotype)
m5 <- glmmTMB(Number_inflorescences ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)
m6 <- glmmTMB(Number_inflorescences ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)
m7 <- glmmTMB(Number_inflorescences ~ Cultivar + (1|Block), family = nbinom2, data = Phenotype)
m8 <- glmmTMB(Number_inflorescences ~ Treatment_worded + (1|Block), family = nbinom2, data = Phenotype)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8) # m7 is best of these

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m3 to be preferred, lowest AIC. No significant difference (AIC <2) but it is also the most parsimonious. 

Model <- glm.nb(Number_inflorescences ~ Cultivar, data = Phenotype)
car::Anova(Model , type = "III")

# ---- Nectar models ----

# Is the amount of nectar significantly different between cultivars
m1 <- glm.nb(Microliter ~ Cultivar * Treatment_worded, data = Nectar)
m2 <- glm.nb(Microliter ~ Cultivar + Treatment_worded, data = Nectar)
m3 <- glm.nb(Microliter ~ Cultivar, data = Nectar)
m4 <- glm.nb(Microliter ~ Treatment_worded, data = Nectar)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m3 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Microliter ~ Cultivar, data = Nectar)
car::Anova(Model , type = "III")

# Now remove the plants that were not completely tried 
Nectar_cleaned <- Nectar %>%
  filter(!Plant_ID %in% c("A76", "V64", "V67", "V68", "V71", "V73", "C89", "C85", "C64", "C68", "C73", "C76", "C62"))

m1 <- glm.nb(Microliter ~ Cultivar * Treatment_worded, data = Nectar_cleaned )
m2 <- glm.nb(Microliter ~ Cultivar + Treatment_worded, data = Nectar_cleaned )
m3 <- glm.nb(Microliter ~ Cultivar, data = Nectar_cleaned )
m4 <- glm.nb(Microliter ~ Treatment_worded, data = Nectar_cleaned )

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m4 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Microliter ~ Treatment_worded, data = Nectar_cleaned)
car::Anova(Model , type = "III")


# ---- Flower (round) models ----

m1 <- glm.nb(Average_Inflorescence_Length ~ Cultivar * Treatment_worded, data = Flowers)
m2 <- glm.nb(Average_Inflorescence_Length ~ Cultivar + Treatment_worded, data = Flowers)
m3 <- glm.nb(Average_Inflorescence_Length ~ Cultivar, data = Flowers)
m4 <- glm.nb(Average_Inflorescence_Length ~ Treatment_worded, data = Flowers)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m4 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Average_Inflorescence_Length ~ Treatment_worded, data = Flowers)
car::Anova(Model , type = "III")

# Now the number of inflorescences with the counting closest to nectar (could later be combined perhaps to figure out significance between nectar and number of inflorescences)
m1 <- glm.nb(Number_Inflorescences ~ Cultivar * Treatment_worded, data = Flowers)
m2 <- glm.nb(Number_Inflorescences ~ Cultivar + Treatment_worded, data = Flowers)
m3 <- glm.nb(Number_Inflorescences ~ Cultivar, data = Flowers)
m4 <- glm.nb(Number_Inflorescences ~ Treatment_worded, data = Flowers)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m3 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Number_Inflorescences ~ Cultivar, data = Flowers)
car::Anova(Model , type = "III")


# ---- Flower (date) models ----

m1 <- glm.nb(Date ~ Cultivar * Treatment_worded, data = Flowering_date)
m2 <- glm.nb(Date ~ Cultivar + Treatment_worded, data = Flowering_date)
m3 <- glm.nb(Date ~ Cultivar, data = Flowering_date)
m4 <- glm.nb(Date ~ Treatment_worded, data = Flowering_date)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m3 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Date ~ Cultivar, data = Flowering_date)
car::Anova(Model , type = "III")


# ---- Repotting models ----

# Seed presence
m1 <- glm.nb(Seeds_present ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seeds_present ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seeds_present ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seeds_present ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m3 to be preferred, lowest AIC. But not significantly.

Model <- glm.nb(Seeds_present ~ Cultivar, data = Repotting)
car::Anova(Model , type = "III")

# Seed pod number
m1 <- glm.nb(Seed_pod_abundance ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_pod_abundance ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_pod_abundance ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_pod_abundance ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m to be preferred, lowest AIC. 

Model <- glm.nb(Seeds_pod_abundance ~ , data = Repotting)
car::Anova(Model , type = "III")

# Seed number
m1 <- glm.nb(Seed_abundance ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_abundance ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_abundance ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_abundance ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m3 to be preferred, lowest AIC. But, not significant.

Model <- glm.nb(Seed_abundance ~ Cultivar, data = Repotting)
car::Anova(Model , type = "III")

# Seed weight
m1 <- glm.nb(Seed_weight ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_weight ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_weight ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_weight ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m4 to be preferred, lowest AIC.

Model <- glm.nb(Seed_weight ~ Cultivar, data = Repotting)
car::Anova(Model , type = "III")

# Seed germination?
m1 <- glm.nb(Seed_germination ~ Cultivar * Treatment_worded, data = Repotting)
m2 <- glm.nb(Seed_germination ~ Cultivar + Treatment_worded, data = Repotting)
m3 <- glm.nb(Seed_germination ~ Cultivar, data = Repotting)
m4 <- glm.nb(Seed_germination ~ Treatment_worded, data = Repotting)

anova(m1)
anova(m2)
anova(m3)
anova(m4)

AIC(m1, m2, m3, m4)
# Overall this shows m to be preferred, lowest AIC.

Model <- glm.nb(Seed_germination ~ , data = Repotting)
car::Anova(Model , type = "III")


# ---- Soil models ----

# Wetness (%) of soil
m1 <- glm.nb(Wet ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(Wet ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(Wet ~ Cultivar, data = Soil)
m4 <- glm.nb(Wet ~ Treatment_worded, data = Soil)
m5 <- glmmTMB(Wet ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m6 <- glmmTMB(Wet ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m7 <- glmmTMB(Wet ~ Cultivar + (1|Block), family = nbinom2, data = Soil)
m8 <- glmmTMB(Wet ~ Treatment_worded + (1|Block), family = nbinom2, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8) # m5 is preferred, lowest AIC. But not significant from m6

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m8 to be preferred, lowest AIC. But not significant, so go to most simple one: m4

Model <- glm.nb(Wet ~ Treatment_worded, data = Soil)
car::Anova(Model , type = "III")

# ECp
m1 <- glm.nb(ECp ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(ECp ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(ECp ~ Cultivar, data = Soil)
m4 <- glm.nb(ECp ~ Treatment_worded, data = Soil)
m5 <- glmmTMB(ECp ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m6 <- glmmTMB(ECp ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m7 <- glmmTMB(ECp ~ Cultivar + (1|Block), family = nbinom2, data = Soil)
m8 <- glmmTMB(ECp ~ Treatment_worded + (1|Block), family = nbinom2, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8)

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m8 to be preferred, lowest AIC. 

Model <- glmmTMB(ECp ~ Treatment_worded + (1|Block), family = nbinom2, data = Soil)
car::Anova(Model , type = "III")

# Eb
m1 <- glm.nb(Eb ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(Eb ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(Eb ~ Cultivar, data = Soil)
m4 <- glm.nb(Eb ~ Treatment_worded, data = Soil)
m5 <- glmmTMB(Eb ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m6 <- glmmTMB(Eb ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m7 <- glmmTMB(Eb ~ Cultivar + (1|Block), family = nbinom2, data = Soil)
m8 <- glmmTMB(Eb ~ Treatment_worded + (1|Block), family = nbinom2, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8)

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m4 to be preferred, lowest AIC. But not significant, however it is the more simple model.

Model <- glm.nb(Eb ~ Treatment_worded, data = Soil)
car::Anova(Model , type = "III")

# ECb
m1 <- glm.nb(ECb ~ Cultivar * Treatment_worded, data = Soil)
m2 <- glm.nb(ECb ~ Cultivar + Treatment_worded, data = Soil)
m3 <- glm.nb(ECb ~ Cultivar, data = Soil)
m4 <- glm.nb(ECb ~ Treatment_worded, data = Soil)
m5 <- glmmTMB(ECb ~ Cultivar * Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m6 <- glmmTMB(ECb ~ Cultivar + Treatment_worded + (1|Block), family = nbinom2, data = Soil)
m7 <- glmmTMB(ECb ~ Cultivar + (1|Block), family = nbinom2, data = Soil)
m8 <- glmmTMB(ECb ~ Treatment_worded + (1|Block), family = nbinom2, data = Soil)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5, m6, m7, m8)

AIC(m1, m2, m3, m4, m5, m6, m7, m8)
# Overall this shows m8 to be preferred, lowest AIC. But not significantly, so go for the simpler m4.

Model <- glm.nb(ECb ~ Treatment_worded, data = Soil)
car::Anova(Model , type = "III")


# ---- Observation models ----

# Can the total arthropods (pollinators) be explained by cultivar or treatment? What about the block?
m1 <- glm.nb(Total_arthropods ~ Cultivar * Treatment_worded, data = Observations)
m2 <- glm.nb(Total_arthropods ~ Cultivar + Treatment_worded, data = Observations)
m3 <- glm.nb(Total_arthropods ~ Cultivar, data = Observations)
m4 <- glm.nb(Total_arthropods ~ Treatment_worded, data = Observations)
m5 <- glm.nb(Total_arthropods ~ Block, data = Observations)
m6 <- glm.nb(Total_arthropods ~ Block + Cultivar, data = Observations)
m7 <- glm.nb(Total_arthropods ~ Block + Treatment_worded, data = Observations)
m8 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Block), family = nbinom2, data = Observations)
m9 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Time), family = nbinom2, data = Observations)
m10 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Time) + (1|Block), family = nbinom2, data = Observations)
m11 <- glmmTMB(Total_arthropods ~ Cultivar + (1|Time) + (1|Block) + (1|Treatment_worded), family = nbinom2, data = Observations)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)
anova(m8, m9, m10, m11)

AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)
# Overall this shows m1 to be preferred, lowest AIC.

Model <- glm.nb(Total_arthropods ~ Cultivar * Treatment_worded, data = Observations)
car::Anova(Model , type = "III")


# ---- Combined / other models ----

# Combine the data frames into one data frame, using the common column "Plant_ID" to be able to work with all data in one data frame.
Combined_data <- Phenotype %>%
  #  left_join(Nectar, by = "Plant_ID") %>%
  left_join(Nectar_cleaned, by = "Plant_ID") %>%
  left_join(Flowers, by = "Plant_ID") %>%
  left_join(Repotting, by = "Plant_ID") %>%
  left_join(Observations, by = "Plant_ID") %>%
  left_join(Soil, by = "Plant_ID") %>%
  left_join(Flowering_date, by = "Plant_ID")

m1 <- glm.nb(Microliter ~ Number_Inflorescences, data = Combined_data)
m2 <- glm.nb(Microliter ~ Number_Inflorescences + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Microliter ~ Number_Inflorescences + Cultivar, data = Combined_data)
m4 <- glm.nb(Microliter ~ Number_Inflorescences + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Microliter ~ Average_Inflorescence_Length, data = Combined_data)
m6 <- glm.nb(Microliter ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
m7 <- glm.nb(Microliter ~ Average_Inflorescence_Length + Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m6 to be preferred, lowest AIC.

Model <- glm.nb(Microliter ~ Average_Inflorescence_Length + Cultivar, data = Combined_data)
car::Anova(Model , type = "III")

m1 <- glm.nb(Microliter ~ Total_arthropods, data = Combined_data)
m2 <- glm.nb(Microliter ~ Total_arthropods + Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Microliter ~ Total_arthropods + Cultivar, data = Combined_data)
m4 <- glm.nb(Microliter ~ Total_arthropods + Treatment_worded, data = Combined_data)
m5 <- glm.nb(Microliter ~ Total_arthropods + Cultivar * Treatment_worded, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)

AIC(m1, m2, m3, m4, m5)
# Overall this shows m3 to be preferred, lowest AIC. But not significant.

Model <- glm.nb(Microliter ~ Total_arthropods + Cultivar, data = Combined_data)
car::Anova(Model , type = "III")

m1 <- glm.nb(Number_flowers ~ Cultivar, data = Combined_data)
m2 <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Combined_data)
m3 <- glm.nb(Number_flowers ~ Nodules_present, data = Combined_data)
m4 <- glm.nb(Number_flowers ~ Treatment_worded, data = Combined_data)
m5 <- glm.nb(Number_flowers ~ Root_abundance, data = Combined_data)
m6 <- glm.nb(Number_flowers ~ Seeds_present, data = Combined_data)
m7 <- glm.nb(Number_flowers ~ Seed_pod_abundance, data = Combined_data)

anova(m1)
anova(m2)
anova(m3)
anova(m4)
anova(m5)
anova(m6)
anova(m7)

AIC(m1, m2, m3, m4, m5, m6, m7)
# Overall this shows m2 to be preferred, lowest AIC.

Model <- glm.nb(Number_flowers ~ Cultivar + Treatment_worded, data = Combined_data)
car::Anova(Model , type = "III")






























