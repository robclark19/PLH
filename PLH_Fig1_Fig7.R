

# Figures made by Rob Clark in revisions of Potato Leaf Hopper paper
# Rob cleaned data into a csv for better formatting in R
# April 3 2023

library("tidyverse")
library("janitor")
library("lubridate")
library("ggplot2")


# Figure 1 is a line plot of the 7 of the seasonal patterns for each site

season_dat <- read.csv("Fig 1 Data2023.csv") %>% clean_names()
str(season_dat)

# Replace "n.s." with NA
season_dat$orr_il[season_dat$orr_il == "n.s."] <- NA

# Convert variable to numeric
season_dat$orr_il <- as.numeric(season_dat$orr_il)

# Get Julian days


# Convert three-letter month codes to numeric months
season_dat$month <- match(toupper(season_dat$month), toupper(month.abb))

# Convert year, month, and day columns to a Date format
season_dat$date <- as.Date(paste(season_dat$year, season_dat$month, season_dat$day, sep = "-"), format = "%Y-%m-%d")

# Calculate Julian day using the yday() function from lubridate
season_dat$julian_day <- yday(season_dat$date)


# Make a kernal density plot for freeport_il

# Define a function that returns the kernel density estimate of the orr_il values
density_func <- approxfun(density(season_dat$orr_il))

# Plot the kernel density curve with orr_il values on the y-axis
ggplot(data_df, aes(x = orr_il, y = ..x..)) + 
  stat_function(fun = density_func, aes(colour = "Density")) +
  xlim(min(data_df$orr_il), max(data_df$orr_il)) + 
  ylim(min(data_df$julian_day_year), max(data_df$julian_day_year)) + 
  xlab("Orr, IL") + 
  ylab("Julian day year") +
  theme_classic()
