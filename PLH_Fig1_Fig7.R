

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


# Reshape the data frame into long format
season_long <- gather(season_dat, site, count, freeport_il:hancock_wi, factor_key=TRUE)

# Create sequential julian days
season_long$sequential_jd <- season_long$julian_day + ((season_long$year - min(season_long$year)) * 365) + as.numeric(format(as.Date(paste0(season_long$year, "-01-01")), "%j")) - 1
# not necessary if you fragment the axis by year



# Create a ggplot with 7 facets
# still borked. needs partitions along the x axis by year
ggplot(season_long, aes(x = sequential_jd , y = count, group = site)) +
  geom_line() +
  facet_wrap(~site, scales = "free_y") +
  scale_y_continuous(
    breaks = season_long$year,
    labels = season_long$year) +
  xlab("Column") +
  ylab("Value") +
  theme_bw()
