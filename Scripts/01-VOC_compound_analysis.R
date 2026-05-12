# ---- Before starting ----
# As GitHub is used, before running the code, make sure to ' PULL ' before getting started to get the latest version of the code and data
# Also go through the steps to save the code while working to avoid losing work
# Save (file) - Stage (box) - Commit (add what you changed) - Push (safe changes online)

# Make sure that the library is synchronized, this ensures that the code will always run as the same version is used of the packages
# First a lockfile needs to be made to store the versions of the packages used
# renv::init() -> done, now a lockfile is created

renv::restore() # can now be used to restore the versions of the specific packages previously used and loaded to work through the code

# clean the environment to avoid conflicts with other projects or names
rm(list = ls())

# When a package is used for the first time, also add to lockfile
# renv::snapshot() # This shows a new library downloaded in manuscript
# Load the packages that are needed for this project
library(tidyverse) # Includes ggplot2, dplyr, etc. can also add them separately 
# load the required packages
library(tidyverse)
library(ggplot2)
library(dplyr)

readr::read_csv # This makes a tibble instead of table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Exploratory Analysis ----

# Load the data you want to use
# Load the data file from google drive 
volatiles <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ3RbmLQrb8FdeBQgmGatpg5KtOOr4TqRqxVhJZwY8k312ufvbwdcagTyuHrxvOprSR95EtlhW4Oh5B/pub?gid=0&single=true&output=csv") 

# Replace all NAs in the entire data frame with 0
Volatiles <- volatiles |>
  mutate_all(~ ifelse(is.na(.), 0, .))

library(chemodiv)
library(corrr)
library(factoextra)
library(ggpubr)
library(purrr)
library(tibble)

# Rename the data sheet to something easier to work with, and make it a data frame
Vs <- Volatiles
Vs <- as.data.frame(Vs)
rownames(Vs) <- Vs$"...1"
new_names <- Vs$Sample
Vs <- Vs[,-c(1)] #cols
X <- as.matrix(Vs)

# Trying z-scores to scale the data before PCA, as the compounds have different ranges. This is important for PCA as it is sensitive to the scale of the data. Z-scores will standardize the data to make comparison between the compounds in the PCA possible
X_z <- scale(X)
colMeans(X_z)    # should be ~0
apply(X_z, 2, sd)  # should be 1

pca.vs <- prcomp(X_z, scale. = F, center = T)

# Create treatment/genotype factor so I can colour/shape the PCA by factors
treatment <- factor(c("s","u","u","s","s","s","u","u","u","u",
                      "s","s","s","s","s","u","u","s","u","u",
                      "s","s","u","u","s","u","u","u","s","s",
                      "u","u","s","s","s","u"))

genotype <- factor(c("A","A","A","A","A","A",
                     "A","A","A","A","A","A",
                     "C","C","C","C","C","C",
                     "C","C","C","C","C","C",
                     "V","V","V","V","V","V",
                     "V","V","V","V","V","V"))

treatment_levels <- levels(treatment)
pch_map <- setNames(c(18, 20), treatment_levels)  # triangle, square, circle
pch_items <- pch_map[as.character(treatment_levels)]
# To add the shapes for the PCA using the genus after creating the genus factor 
# Create a small shapes item

plot(pca.vs$x[,1], pca.vs$x[,2],
     col = genotype,           # color by genotype
     pch = pch_items,      
     cex = 0.9,
     xlab = paste0("PC1 (", round(summary(pca.vs)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca.vs)$importance[2,2]*100,1), "%)"),
     main = "Volatiles in salt treated alfalfa")


# Scree plot
fviz_eig(pca.vs)

plot(pca.vs$x[,1], pca.vs$x[,3],
     col = genotype,           # color by genotype
     pch = pch_items,      
     cex = 0.9,
     xlab = paste0("PC1 (", round(summary(pca.vs)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC3 (", round(summary(pca.vs)$importance[2,3]*100,1), "%)"),
     main = "Volatiles in salt treated alfalfa")

plot(pca.vs$x[,1], pca.vs$x[,4],
     col = genotype,           # color by genotype
     pch = pch_items,      
     cex = 0.9,
     xlab = paste0("PC1 (", round(summary(pca.vs)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC4 (", round(summary(pca.vs)$importance[2,4]*100,1), "%)"),
     main = "Volatiles in salt treated alfalfa")











