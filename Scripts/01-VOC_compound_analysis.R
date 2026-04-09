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
library(tidyverse) # This loads the tidyverse, it includes ggplot2, dplyr, etc. but can also add them separately 
# load the required packages
library(tidyverse)
library(ggplot2)
library(dplyr)

readr::read_csv # This makes a tibble instead of table, for every variable it stores what the type of variable is. It doesn't just stop at the length it can print, which is what table does.

# ---- Exploratory Analysis ----

# Load the data you want to use
# load the data file from google drive 
Volatiles <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ3RbmLQrb8FdeBQgmGatpg5KtOOr4TqRqxVhJZwY8k312ufvbwdcagTyuHrxvOprSR95EtlhW4Oh5B/pub?gid=0&single=true&output=csv") 






