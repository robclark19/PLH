# Notes ####
# Script for using Stepwise AIC for discovering a parsimonous mixed model with potato leaf hopper data
# Data from midwestern suction trap network provided by Mike Crossley and collaborators

# Tips from Rob on 3/2/2023
# All data is imported like prior scripts, but check to make sure I made the correct assumptions
# Model averaging in this context would be a serious effort, perhaps not yet warranted.
# model.avg() example code is retained for model 1.
# stepAIC() is flexible and the models are easily changed to include interactions
# You can force stepAIC() to keep certain independent variables with the "keep =" statement
# More data & time would be needed for lasso/ridge regression, so that's ruled out for now

# In a larger project we could write functions to automatically take the trimmed model and provide outputs
# For now the trimmed models have to be entered manually when inspecting output from stepAIC(model)$anova

# Libraries ####
library("tidyverse")
library("janitor")
library("lme4")
library("car")
library("multcomp")
library("ggplot2")
library("emmeans")
library("nlme")
library("MASS") #stepAIC()
library("glmmTMB")
library("MuMIn") # model.avg(). Rob has no experience with this approach so proceed with caution.
library("ggeffects") # produces model predictions for ggplot
library("ggpubr") # tools for making pub figures

# Import raw data ####

plh_dat = read.table('PLH_pheno_weather_land.txt',sep='\t',as.is=T,check.names=F,header=T) %>%
  clean_names() %>%
  mutate_at(vars(year,hoppers_week), factor) %>%
  mutate_at(vars(gdd_total), function(x) x/1000) # scale down growing degree days

# Trimming methods following original script ####
# Trim site*years with low sample coverage, Chase, & sites lacking first or last captures
threshold = 15 # cutoff

plh_dat <- plh_dat %>%
  subset(., n_week >= threshold) %>% 
  subset(., site != 'Chase, LA') %>% # Trim Chase, LA
  subset(., first < 213) %>%
  subset(., last > 213)


# Model 1 ####
# Effects on number of potato leafhoppers

# Fully specified model using all independent variables seen in original script

mod1_full <- glmmTMB(n_hoppers ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + tmean_summer,
                              family=nbinom2, # Tweedie would be better but that model doesn't converge
                              data=plh_dat) 

mod1_sa <- stepAIC(mod1_full, direction = "both")

mod1_sa$anova

# Step AIC runs all iterations of the model and then provides a 'conditional model'
# The conditional model is evaluated only on delta AIC (difference between original and final)

mod1_trim <- glmmTMB(n_hoppers ~ lat + p_beans, data = plh_dat)

Anova(mod1_trim)
summary(mod1_trim) # more or less looks like only beans matter

# Model averaging test case ####
# https://stats.stackexchange.com/questions/208724/interpreting-model-averaging-results-in-r
# https://search.r-project.org/CRAN/refmans/MuMIn/html/model.avg.html

# pick two models with similar AIC reported from stepAIC output
mod1_a<- glmmTMB(n_hoppers ~ lat + p_beans + n_week, data = plh_dat)
mod1_b <- glmmTMB(n_hoppers ~ lat + p_beans + p_potato, data = plh_dat)

# make a model list and then report summary of conditional models
mod_list <- list(mod1_a,mod1_b)
model.avg(mod_list) %>% summary()




# Model 2 ####
# Effects on date of first capture

mod2_full <- glmmTMB(first ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + tmean_summer,
                     family=gaussian(), #I assume these are Julian date values which are normally distributed
                     data=plh_dat) 

mod2_sa <- stepAIC(mod2_full, direction = "both")

mod2_sa$anova

mod2_trim <- glmmTMB(first ~ lat, data = plh_dat)

Anova(mod2_trim)
summary(mod2_trim)





# Model 3 ####
# Effects on date of last capture

mod3_full <- glmmTMB(last ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + tmean_summer,
                     family=gaussian(),
                     data=plh_dat) 

mod3_sa <- stepAIC(mod3_full, direction = "both")

mod3_sa$anova

mod3_trim <- glmmTMB(last ~ year + p_potato, data = plh_dat)

Anova(mod3_trim)
summary(mod3_trim)


# Model 4 ####
# Effects on the range between first and last

mod4_full <- glmmTMB(range ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + tmean_summer,
                     family=gaussian(),
                     data=plh_dat) 

mod4_sa <- stepAIC(mod4_full, direction = "both")

mod4_sa$anova

mod4_trim <- glmmTMB(range ~ year + lat + p_potato, data = plh_dat)

Anova(mod4_trim)
summary(mod4_trim)


# Example Data Viz ####

mod4_pred = ggpredict(mod4_trim, ci.lvl=NA)

mod4_plot_year = mod4_pred$year %>%
  ggplot(aes(x=x,y=predicted)) +
  geom_bar(stat="identity") + # pull points from raw data
  labs(x = "Year", y = "Date Range of Planthopper Observations")

mod4_plot_year

mod4_plot_lat = mod4_pred$lat %>%
  ggplot(aes(x=x,y=predicted)) +
  geom_point(data=plh_dat, aes(x=lat, y=range), size=0.5, alpha=0.5) + # pull points from raw data
  geom_line() +
  # scale_y_continuous(limits=c(40,160)) + # trim lower ones from viz if needed
  labs(x = "Latitude", y = "Date Range of Planthopper Observations")

mod4_plot_lat

mod4_plot_tater = mod4_pred$p_potato %>%
  ggplot(aes(x=x,y=predicted)) +
  geom_point(data=plh_dat, aes(x=p_potato, y=range), size=0.5, alpha=0.5) + # pull points from raw data
  geom_line() +
  labs(x = "Proportional coverage of nearby potato", y = "Date Range of Planthopper Observations")

mod4_plot_tater

# arrange and output as png for pubs
fig_1ab <- ggarrange(mod4_plot_year, mod4_plot_lat, mod4_plot_tater, labels = c("a", "b", "c"), nrow = 1,
                     common.legend = FALSE, widths = c(1, 1))

ggsave(filename = "./figures/fig_1abc.png", plot = fig_1ab , device = "png",
       width = 10, height = 4, units = "in", scale = 0.9)


