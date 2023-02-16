

######################################################################################
# Explain site*year-level variation in PLH aggregate activity-density and phenology

library(DHARMa)
library(glmmTMB)
library(performance)
library(ggplot2)
library(gridExtra)
library(sjPlot)

# For plotting significant effects
theme1 = theme_classic() +
	theme(axis.line = element_line(size=2)) +
	theme(axis.ticks = element_line(size=2)) + 
	theme(axis.ticks.length=unit(.25, "cm")) + 
	theme(axis.text = element_text(colour='black',size=20)) +
	theme(axis.title = element_text(colour='black',size=24)) + 
	theme(plot.margin = margin(t=15,r=5,b=5,l=5, unit="pt")) + 
	theme(plot.title = element_blank())


# Curated data
pheno = read.table('PLH_pheno_weather_land.txt',sep='\t',as.is=T,check.names=F,header=T)

# Trim site*years with low sample coverage
thresh = 15
keep = which(pheno$N.week>=thresh)
length(keep) / nrow(pheno) #retain 93% of data
pheno2 = pheno[keep,]
pheno2 = pheno2[which(pheno2$Site!='Chase, LA'),] # Trim Chase, LA
dim(pheno2)
plot(pheno2$Lon,pheno2$Lat) # Check sites retained after trimming

# Create trimmed datasets that remove site*years that lack expected first or last captures
pheno2.first = pheno2[which(pheno2$First<213),]; dim(pheno2.first)
pheno2.last = pheno2[which(pheno2$Last>213),]; dim(pheno2.last)
pheno2.range = pheno2[which(pheno2$First<213 & pheno2$Last>213),]; dim(pheno2.range)


###
# Model site-level differences in total PLH counts among years
# Format data
pheno2$Site. = as.factor(pheno2$Site)
pheno2$Year. = as.numeric(as.factor(pheno2$Year))
pheno2$Lat. = as.factor(pheno2$Lat)
pheno2$p.alfalfaZ = (pheno2$p.alfalfa - mean(pheno2$p.alfalfa,na.rm=T)) / sd(pheno2$p.alfalfa,na.rm=T)
pheno2$p.potatoZ = (pheno2$p.potato - mean(pheno2$p.potato,na.rm=T)) / sd(pheno2$p.potato,na.rm=T)
pheno2$p.beansZ = (pheno2$p.beans - mean(pheno2$p.beans,na.rm=T)) / sd(pheno2$p.beans,na.rm=T)
pheno2$LatZ = (pheno2$Lat - mean(pheno2$Lat)) / sd(pheno2$Lat)
pheno2$N.weekZ = (pheno2$N.week - mean(pheno2$N.week)) / sd(pheno2$N.week)
pheno2$GDD2 = pheno2$GDD.total/1000

mod1 = glmmTMB(N.hoppers ~ Year. + N.week + Lat + p.alfalfa + p.potato + p.beans, data=pheno2, family='nbinom2')
summary(mod1)
confint(mod1)
mod1b = glmmTMB(N.hoppers ~ Year. + N.week + GDD2 + Lat + p.alfalfa + p.potato + p.beans, data=pheno2, family='nbinom2')
summary(mod1b)
confint(mod1b)
#mod2 = glmmTMB(N.hoppers ~ Year. + N.week + p.alfalfa + p.potato + p.beans + (1|Lat.), data=pheno2, family='nbinom2')
#summary(mod2)
#mod3 = glmmTMB(N.hoppers ~ Year. + N.weekZ + LatZ + p.alfalfaZ + p.potatoZ + p.beansZ, data=pheno2, family='nbinom2')
#summary(mod3)
#mod4 = glmmTMB(N.hoppers ~ Year. + N.weekZ + p.alfalfaZ + p.potatoZ + p.beansZ + (1|Lat.), data=pheno2, family='nbinom2')
#summary(mod4)
#AIC(mod1,mod1b,mod2,mod3,mod4) #mod2 and mod4 lowest = 559 (mod1 and mod3 = 568)
#scaling covariates had no effect
#Year (positive - i.e., PLH increasing over time), Lat (negative - i.e., fewer PLH further north), and beans (positive - i.e., more PLH with more beans) effects are sig
#only Year (positive) is sig when Lat is included as random intercept

# Check colinearity
check_collinearity(mod1)
# Check for spatial autocorrelation
res = simulateResiduals(mod1)
res2 = recalculateResiduals(res, group=pheno2$Site)
testSpatialAutocorrelation(res2, x=aggregate(pheno2$Lon,list(pheno2$Site),mean)$x, y=aggregate(pheno2$Lat,list(pheno2$Site),mean)$x) #(none)

# Plot significant effects
v1 = plot_model(mod1b, type="pred", terms=c('Year.'), dot.size=5, line.size=2) + theme1 + xlab('Year') + ylab('Total PLH')
	ggsave(plot=v1,filename='Nhoppers vs year.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v2 = plot_model(mod1b, type="pred", terms=c('Lat'), dot.size=5, line.size=2) + theme1 + xlab('Latitude') + ylab('Total PLH')
	ggsave(plot=v2,filename='Nhoppers vs Lat.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v3 = plot_model(mod1b, type="pred", terms=c('p.beans'), dot.size=5, line.size=2) + theme1 + xlab('Prop. beans') + ylab('Total PLH')
	ggsave(plot=v3,filename='Nhoppers vs beans.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v4 = plot_model(mod1b, type="pred", terms=c('GDD2'), dot.size=5, line.size=2) + theme1 + xlab('GDD (/1000)') + ylab('Total PLH')
	ggsave(plot=v3,filename='Nhoppers vs GDD.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
# 2x2 multipanel figure
v.grid = grid.arrange(v1,v2,v3,v4,ncol=2)
ggsave(v.grid,filename='Nhoppers all.png',path='./plots/glmm/',width=10*4,height=10*2,units='cm',dpi=600)


###
# Model date of first capture
# Format data
pheno2.first$Site. = as.factor(pheno2.first$Site)
pheno2.first$Year. = as.numeric(as.factor(pheno2.first$Year))
pheno2.first$Lat. = as.factor(pheno2.first$Lat)
pheno2.first$p.alfalfaZ = (pheno2.first$p.alfalfa - mean(pheno2.first$p.alfalfa,na.rm=T)) / sd(pheno2.first$p.alfalfa,na.rm=T)
pheno2.first$p.potatoZ = (pheno2.first$p.potato - mean(pheno2.first$p.potato,na.rm=T)) / sd(pheno2.first$p.potato,na.rm=T)
pheno2.first$p.beansZ = (pheno2.first$p.beans - mean(pheno2.first$p.beans,na.rm=T)) / sd(pheno2.first$p.beans,na.rm=T)
pheno2.first$LatZ = (pheno2.first$Lat - mean(pheno2.first$Lat)) / sd(pheno2.first$Lat)
pheno2.first$N.weekZ = (pheno2.first$N.week - mean(pheno2.first$N.week)) / sd(pheno2.first$N.week)
pheno2.first$GDD.first.date2 = pheno2.first$GDD.first.date/1000

first.mod1 = glmmTMB(First ~ Year. + N.week + Lat + p.alfalfa + p.potato + p.beans, data=pheno2.first, family='nbinom2')
summary(first.mod1)
confint(first.mod1)
first.mod1b = glmmTMB(First ~ Year. + N.week + GDD.first.date2 + Lat + p.alfalfa + p.potato + p.beans, data=pheno2.first, family='nbinom2')
summary(first.mod1b)
#first.mod2 = glmmTMB(First ~ Year. + N.weekZ + p.alfalfaZ + p.potatoZ + p.beansZ + (1|Lat.), data=pheno2.first, family='nbinom2')
#summary(first.mod2)
#first.mod3 = glmmTMB(First ~ Year. + N.weekZ + LatZ + p.alfalfaZ + p.potatoZ + p.beansZ, data=pheno2.first, family='nbinom2')
#summary(first.mod3)
#AIC(first.mod1,first.mod2,first.mod3) #mod1 is best = 600 (mod2 = 606)
#model with Lat as a fixed effect found sig positive effect of Lat
#i.e., first captures were later further north

# Check colinearity
check_collinearity(first.mod1)
# Check for spatial autocorrelation
res = simulateResiduals(first.mod1)
res2 = recalculateResiduals(res, group=pheno2.first$Site)
testSpatialAutocorrelation(res2, x=aggregate(pheno2.first$Lon,list(pheno2.first$Site),mean)$x, y=aggregate(pheno2.first$Lat,list(pheno2.first$Site),mean)$x) #(none)

# Plot significant effects
v1 = plot_model(first.mod1b, type="pred", terms=c('Lat'), dot.size=5, line.size=2) + theme1 + xlab('Latitude') + ylab('Day of 1st detection')
	ggsave(plot=v1,filename='First vs Lat.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v2 = plot_model(first.mod1b, type="pred", terms=c('GDD.first.date2'), dot.size=5, line.size=2) + theme1 + xlab('GDD (/1000)') + ylab('Day of 1st detection')
	ggsave(plot=v1,filename='First vs GDD.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v.grid = grid.arrange(v1,v2,ncol=2)
ggsave(v.grid,filename='First all.png',path='./plots/glmm/',width=10*2,height=10,units='cm',dpi=600)


###
# Model date of last capture
# Format data
pheno2.last$Site. = as.factor(pheno2.last$Site)
pheno2.last$Year. = as.numeric(as.factor(pheno2.last$Year))
pheno2.last$Lat. = as.factor(pheno2.last$Lat)
pheno2.last$p.alfalfaZ = (pheno2.last$p.alfalfa - mean(pheno2.last$p.alfalfa,na.rm=T)) / sd(pheno2.last$p.alfalfa,na.rm=T)
pheno2.last$p.potatoZ = (pheno2.last$p.potato - mean(pheno2.last$p.potato,na.rm=T)) / sd(pheno2.last$p.potato,na.rm=T)
pheno2.last$p.beansZ = (pheno2.last$p.beans - mean(pheno2.last$p.beans,na.rm=T)) / sd(pheno2.last$p.beans,na.rm=T)
pheno2.last$LatZ = (pheno2.last$Lat - mean(pheno2.last$Lat)) / sd(pheno2.last$Lat)
pheno2.last$N.weekZ = (pheno2.last$N.week - mean(pheno2.last$N.week)) / sd(pheno2.last$N.week)
pheno2.last$GDD.last.date2 = pheno2.last$GDD.last.date/1000

last.mod1 = glmmTMB(Last ~ Year. + N.week + Lat + p.alfalfa + p.potato + p.beans, data=pheno2.last, family='nbinom2')
summary(last.mod1)
confint(last.mod1)
#last.mod1b = glmmTMB(Last ~ Year. + N.week + GDD.last.date2 + Lat + p.alfalfa + p.potato + p.beans, data=pheno2.last, family='nbinom2')
#summary(last.mod1b)
#last.mod2 = glmmTMB(Last ~ Year. + N.weekZ + p.alfalfaZ + p.potatoZ + p.beansZ + (1|Lat.), data=pheno2.last, family='nbinom2')
#summary(last.mod2)
#last.mod3 = glmmTMB(Last ~ Year. + N.weekZ + LatZ + p.alfalfaZ + p.potatoZ + p.beansZ, data=pheno2.last, family='nbinom2')
#summary(last.mod3)
#AIC(last.mod1,last.mod2,last.mod3) #mod1 is best = 545 (mod2 = 550)
#model with Lat as a fixed effect found sig positive effect of Year and sig negative effect of Lat
#i.e., last captures are getting later each year
#i.e., last captures are earlier where there is more potato

# Check colinearity
check_collinearity(last.mod1)
# Check for spatial autocorrelation
res = simulateResiduals(last.mod1)
res2 = recalculateResiduals(res, group=pheno2.last$Site)
testSpatialAutocorrelation(res2, x=aggregate(pheno2.last$Lon,list(pheno2.last$Site),mean)$x, y=aggregate(pheno2.last$Lat,list(pheno2.last$Site),mean)$x) #(none)

# Plot significatn effects
v1 = plot_model(last.mod1, type="pred", terms=c('Year.'), dot.size=5, line.size=2) + theme1 + xlab('Year') + ylab('Day of last detection')
	ggsave(plot=v1,filename='Last vs year.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v2 = plot_model(last.mod1, type="pred", terms=c('p.potato'), dot.size=5, line.size=2) + theme1 + xlab('Prop. potato') + ylab('Day of last detection')
	ggsave(plot=v2,filename='Last vs potato.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
# 2x2 multipanel figure
v.grid = grid.arrange(v1,v2,ncol=2)
ggsave(v.grid,filename='Last all.png',path='./plots/glmm/',width=10*3,height=10*2,units='cm',dpi=600)


# Model date range of captures
# Format data
pheno2.range$Site. = as.factor(pheno2.last$Site)
pheno2.range$Year. = as.numeric(as.factor(pheno2.range$Year))
pheno2.range$Lat. = as.factor(pheno2.range$Lat)
pheno2.range$p.alfalfaZ = (pheno2.range$p.alfalfa - mean(pheno2.range$p.alfalfa,na.rm=T)) / sd(pheno2.range$p.alfalfa,na.rm=T)
pheno2.range$p.potatoZ = (pheno2.range$p.potato - mean(pheno2.range$p.potato,na.rm=T)) / sd(pheno2.range$p.potato,na.rm=T)
pheno2.range$p.beansZ = (pheno2.range$p.beans - mean(pheno2.range$p.beans,na.rm=T)) / sd(pheno2.range$p.beans,na.rm=T)
pheno2.range$LatZ = (pheno2.range$Lat - mean(pheno2.range$Lat)) / sd(pheno2.range$Lat)
pheno2.range$N.weekZ = (pheno2.range$N.week - mean(pheno2.range$N.week)) / sd(pheno2.range$N.week)

range.mod1 = glmmTMB(Range ~ Year. + N.week + Lat + p.alfalfa + p.potato + p.beans, data=pheno2.range, family='nbinom2')
summary(range.mod1)
confint(range.mod1)
#range.mod2 = glmmTMB(Range ~ Year. + N.weekZ + p.alfalfaZ + p.potatoZ + p.beansZ + (1|Lat.), data=pheno2.range, family='nbinom2')
#summary(range.mod2)
#range.mod3 = glmmTMB(Range ~ Year. + N.weekZ + LatZ + p.alfalfaZ + p.potatoZ + p.beansZ, data=pheno2.range, family='nbinom2')
#summary(range.mod3)
#AIC(range.mod1,range.mod2,range.mod3) #mod1 is best = 715 (mod2 = 718)
#Both models found sig positive effect of Year and sig negative effect of potato
#i.e., duration of captures is getting longer each year
#i.e., duration of captures is shorter where there is more potato

# Check colinearity
check_collinearity(range.mod1)
# Check for spatial autocorrelation
res = simulateResiduals(range.mod1)
res2 = recalculateResiduals(res, group=pheno2.range$Site)
testSpatialAutocorrelation(res2, x=aggregate(pheno2.range$Lon,list(pheno2.range$Site),mean)$x, y=aggregate(pheno2.range$Lat,list(pheno2.range$Site),mean)$x) #(none)

# Plot significatn effects
v1 = plot_model(range.mod1, type="pred", terms=c('Year.'), dot.size=5, line.size=2) + theme1 + xlab('Year') + ylab('Duration of detections')
	ggsave(plot=v1,filename='Last vs year.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
v2 = plot_model(range.mod1, type="pred", terms=c('p.potato'), dot.size=5, line.size=2) + theme1 + xlab('Prop. potato') + ylab('Duration of detections')
	ggsave(plot=v2,filename='Last vs potato.png',path='./plots/glmm/',width=10,height=10,units='cm',dpi=600)
# 2x2 multipanel figure
v.grid = grid.arrange(v1,v2,ncol=2)
ggsave(v.grid,filename='Duration all.png',path='./plots/glmm/',width=10*3,height=10*2,units='cm',dpi=600)


