

# Figures made by Rob Clark in revisions of Potato Leaf Hopper paper
# Rob cleaned data into a csv for better formatting in R
# April 3 2023

library("tidyverse")
library("janitor")
library("lubridate")
library("ggplot2")
library("showtext")

font_add_google(name="Raleway", family="Raleway")
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

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


panel_labs = c("Freeport, IL",
               "Orr, IL",
               "Urbana-Champaign, IL",
               "Columbia City, IN",
               "Lafayette, IN",
               "Kanawha, IA",
               "Hancock, WI")
names(panel_labs) = c("freeport_il",
                      "orr_il",
                      "urbana_champaign_ii_il",
                      "columbia_city_in",
                      "lafayette_in",
                      "kanawha_ia",
                      "hancock_wi")

# Create a ggplot with 7 facets
# still borked. needs partitions along the x axis by year
base_pt = 12
ggplot(season_long, aes(x = date, y = count, group = site)) +
  geom_line(size=0.5) +
  facet_wrap(~site, scales = "free_y", labeller=labeller(site=panel_labs),ncol=1) +
  scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month", date_labels = "%Y") +
  xlab("Year") +
  ylab("Count") +
  theme_bw() + theme(
    text = element_text(family="Raleway"),
    panel.grid.major.x = element_line(linewidth=1.3,colour="lightgrey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(linewidth=1.3, colour="lightgrey"),
    axis.text.x = element_text(size=0.8*base_pt, color="black"),
    axis.title = element_text(size=base_pt),
    strip.text = element_text(size=0.8*base_pt, hjust=0.02, margin=margin(b=0,t=0)),
    strip.background = element_blank(),
    panel.border = element_rect(color="black",linewidth=0.8)
  )
ggsave(filename="timeseries.png",path="./figures/",
       width=12,height=30,units="cm",dpi=300,device=ragg::agg_png())
