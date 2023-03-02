
# Script for using Stepwise AIC for discovering a parsimonous mixed model with potato leaf hopper data
# Data from midwestern suction trap network provided by Mike Crossley and collaborators

# Libraries
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

# Import raw data

plh_dat = read.table('PLH_pheno_weather_land.txt',sep='\t',as.is=T,check.names=F,header=T) %>%
  clean_names() %>%
  mutate_at(vars(year,hoppers_week), factor) %>%
  mutate_at(vars(gdd_total), function(x) x/1000) # scale down growing degree days

# Trimming methods following original script
# Trim site*years with low sample coverage, Chase, & sites lacking first or last captures
threshold = 15 # cutoff

plh_dat <- plh_dat %>%
  subset(., n_week >= threshold) %>% 
  subset(., site != 'Chase, LA') %>% # Trim Chase, LA
  subset(., first < 213) %>%
  subset(., last > 213)


# Model 1
# Effects on number of potato leafhoppers

# Fully specified model using all independent variables seen in original script

mod1_full <- glmmTMB(n_hoppers ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + gdd_total,
                              family=nbinom2(), # Tweedie would be better but that model doesn't converge
                              data=plh_dat) 

mod1_sa <- stepAIC(mod1_full, direction = "both")

mod1_sa$anova

# Step AIC runs all iterations of the model and then provides a 'conditional model'
# The conditional model is evaluated only on delta AIC (difference between original and final)

mod1_trim <- glmmTMB(n_hoppers ~ lat + p_beans, data = plh_dat)

Anova(mod1_trim)
summary(mod1_trim) # more or less looks like only beans matter




# Model 2
# Effects on date of first capture

mod2_full <- glmmTMB(first ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + gdd_total,
                     family=gaussian(), #I assume these are Julian date values which are normally distributed
                     data=plh_dat) 

mod2_sa <- stepAIC(mod2_full, direction = "both")

mod2_sa$anova

mod2_trim <- glmmTMB(first ~ lat, data = plh_dat)

Anova(mod2_trim)
summary(mod2_trim)





# Model 3
# Effects on date of last capture

mod3_full <- glmmTMB(last ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + gdd_total,
                     family=gaussian(),
                     data=plh_dat) 

mod3_sa <- stepAIC(mod3_full, direction = "both")

mod3_sa$anova

mod3_trim <- glmmTMB(last ~ year + p_potato, data = plh_dat)

Anova(mod3_trim)
summary(mod3_trim)


# Model 4
# Effects on the range between first and last

mod4_full <- glmmTMB(range ~ year + n_week + lat + p_alfalfa + p_potato + p_beans + gdd_total,
                     family=gaussian(),
                     data=plh_dat) 

mod4_sa <- stepAIC(mod4_full, direction = "both")

mod4_sa$anova

mod4_trim <- glmmTMB(range ~ year + lat + p_potato, data = plh_dat)

Anova(mod4_trim)
summary(mod4_trim)
