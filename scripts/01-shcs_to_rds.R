library(tidyverse)
library(haven)
library(here)
library(bernr)


if (!dir.exists("processed")) {
  dir.create("processed")
} else {
  message("Processed-folder already exists.")
}

dir <- here("data", "2205stata/")
# Get all files in SHCS folder
files <- list.files(path = dir)
files_full <- str_c(dir, files)
names <- str_remove(files, ".dta")

# Read them into a list
all <- map(files_full, read_dta, encoding = "latin1")

# Set names to dta names
names(all) <-  names

# Zipp off all attributes using own function
all <- map(all, clean_dta)

# Write into processed folder
walk2(all, names, ~write_rds(.x, str_c("processed/", .y, ".rds")))

