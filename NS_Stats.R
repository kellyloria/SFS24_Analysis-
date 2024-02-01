#' Script for statistical analysis SFS 2024 conference 
#' 
#' @description investigation of nearshore metabolism drivers 
#' @param lake 
#' 
#' @return 
#' @export 
#' 

##==============================================================================
## Created  02/01/2024 by KAL
#===============================================================================
library(tidyverse)
library(lubridate)

# plotting packages:
library(ggplot2)
library(reshape2)
library(scales)

# stats packages 
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)

se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

##===========================================
## read data aggregated data for the project:
#============================================
# * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/

SFS_datQ <- read.csv("./SFS24_analysis_dat/SFS24_analysis_dat.csv") %>%
  mutate(date = as.Date(date))

str(SFS_datQ)
##===========================================
## Begin analysis 
#============================================
# adjust to whatever covariate rows we may want to include in a mixed model:
corplot <- SFS_datQ[,c(13,21:26)]
chart.Correlation(corplot, histogram=TRUE, pch=19)
# light and air temp variables shouldn't be in the same model 

##===========================================
## GPP glms
#============================================
hist(SFS_datQ$middle_GPP) # strong left skew
hist(log(SFS_datQ$middle_GPP+1)) # better
names(SFS_datQ)

# create auto-regressive term to account for prior GPP
SFS_datQ$lag_middle_GPP <- lag(SFS_datQ$middle_GPP)

# global model with everything 
gpp_mod <- lmer(log(middle_GPP+1) ~ 
                  scale(lag_middle_GPP) +
                  scale(ppt_mm)+ 
                  scale(tmax_C) + 
                  scale(light_mean) +
                  scale(windsp_cv) +
                  scale(flow_mean) +
                  scale(spc_mean) +
                  (1|site.x), data=SFS_datQ)

summary(gpp_mod)
vif(gpp_mod)
hist(residuals(gpp_mod))
r.squaredGLMM(gpp_mod)

# global model with everything for air temp 
gpp_mod1 <- lmer(log(middle_GPP+1) ~ 
                  scale(lag_middle_GPP) +
                  scale(ppt_mm)+ 
                  scale(tmax_C) + 
                  #scale(light_mean) +
                  scale(windsp_cv) +
                  scale(flow_mean) +
                  scale(spc_mean) +
                  (1|site.x), data=SFS_datQ)

summary(gpp_mod1)
vif(gpp_mod1)
hist(residuals(gpp_mod1))
r.squaredGLMM(gpp_mod1)

# global model with everything for light 
gpp_mod2 <- lmer(log(middle_GPP+1) ~ 
                  scale(lag_middle_GPP) +
                  scale(ppt_mm)+ 
                  #scale(tmax_C) + 
                  scale(light_mean) +
                  scale(windsp_cv) +
                  scale(flow_mean) +
                  scale(spc_mean) +
                  (1|site.x), data=SFS_datQ)

summary(gpp_mod2)
vif(gpp_mod2)
hist(residuals(gpp_mod2))
r.squaredGLMM(gpp_mod2)



gpp_mod3 <- lmer(log(middle_GPP+1) ~ 
                   scale(lag_middle_GPP) +
                   scale(ppt_mm)+ 
                   scale(light_cv) +
                   scale(windsp_cv) +
                   scale(spc_mean) +
                   (1|site.x), data=SFS_datQ)

summary(gpp_mod3)
vif(gpp_mod3)
hist(residuals(gpp_mod3))
r.squaredGLMM(gpp_mod3)


gpp_mod4 <- lmer(log(middle_GPP+1) ~ 
                   scale(lag_middle_GPP) +
                   scale(light_cv) +
                   scale(windsp_cv) +
                   scale(spc_mean) +
                   (1|site.x), data=SFS_datQ)

summary(gpp_mod4)
vif(gpp_mod4)
hist(residuals(gpp_mod4))
r.squaredGLMM(gpp_mod4)

# BEST Model #
gpp_mod5 <- lmer(log(middle_GPP+1) ~ 
                   scale(lag_middle_GPP) +
                   scale(ppt_mm)+ 
                   scale(light_cv) +
                   scale(windsp_cv) +
                   scale(flow_cv) +
                   scale(spc_cv) +
                   (1|site.x), data=SFS_datQ)

summary(gpp_mod5)
vif(gpp_mod5)
hist(residuals(gpp_mod5))
r.squaredGLMM(gpp_mod5)

# BEST Model #
gpp_mod6 <- lmer(log(middle_GPP+1) ~ 
                   scale(lag_middle_GPP) +
                   scale(light_cv) +
                   scale(windsp_cv) +
                   #scale(flow_cv) +
                   scale(spc_cv) +
                   (1|site.x), data=SFS_datQ)

summary(gpp_mod6)
vif(gpp_mod6)
hist(residuals(gpp_mod6))
r.squaredGLMM(gpp_mod6)

# quick model comparison: 
AIC(gpp_mod, gpp_mod1, gpp_mod2, gpp_mod3, gpp_mod4, gpp_mod5, gpp_mod6)



