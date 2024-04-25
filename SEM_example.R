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

# setwd("/Users/kellyloria/Documents/UNR/MSMmetab")
library(piecewiseSEM)

getwd()
# SFS_datQ <- read.csv("./SFS24_analysis_dat/SFS24_analysis_dat.csv") %>%
#   mutate(date = as.Date(date))

SFS_datQ <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/SFS24_data_T.rds")
summary(SFS_datQ)
str(SFS_datQ)

## create new df for complete GPP obs 
GPP_df <- SFS_datQ%>%
  dplyr::select(-middle_NEP, -upper_NEP, -lower_NEP, 
                -middle_ER, -upper_ER, -lower_ER)

ER_df <- SFS_datQ%>%
  dplyr::select(-middle_NEP, -upper_NEP, -lower_NEP, 
                -middle_GPP, -upper_GPP, -lower_GPP)

# 
GPP_df1<- drop_na(GPP_df$middle_GPP)
drop_na()
summary(GPP_df1)

## start with all your glms:
GPP_df1 <- GPP_df %>%
  mutate(lag_middle_GPP = lag(middle_GPP, default = NA))

# global model with everything 
library(glmmTMB)

hist(GPP_df1$middle_GPP)
hist(log(GPP_df1$middle_GPP+1))

gpp_fit1 <- glmmTMB(middle_GPP ~ lag(middle_GPP, 1) + (1|site), data = GPP_df1, family = gaussian)
summary(gpp_fit1)
residuals <- residuals(gpp_fit1)
hist(residuals)

gpp_fit1 <- glmmTMB(log(middle_GPP+1) ~ lag(log(middle_GPP+1), 1) + (1|site), data = GPP_df1, family = gaussian)
summary(gpp_fit1)
residuals <- residuals(gpp_fit1)
hist(residuals)

# Plot ACF and PACF of residuals
par(mfrow=c(2,1))
acf(residuals, main="ACF of Residuals")
pacf(residuals, main="PACF of Residuals")



fit2 <- glmmTMB(middle_GPP ~ lag(middle_GPP, 1) + 
                  lag(middle_GPP, 2)+
                  # lag(middle_GPP, 3) + 
                  # lag(middle_GPP, 4) + 
                  # lag(middle_GPP, 5) + 
                  (1|site), data = GPP_df1, family = gaussian)
summary(fit2)
residuals <- residuals(fit2)

# Plot ACF and PACF of residuals
#par(mfrow=c(2,1))
acf(residuals, main="ACF of Residuals")
pacf(residuals, main="PACF of Residuals")



# gpp_auto_mod <- lmer(log(middle_GPP+1) ~ 
#                     scale(lag_middle_GPP) + 
#                     (1|site), data=GPP_df1)
# summary(gpp_auto_mod)
# hist(residuals(gpp_auto_mod))

names(GPP_df1)
library(PerformanceAnalytics)
Streamdf_cor <- GPP_df1[, c(8, 16, 25, 28,30)]
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

gpp_auto_mod <- lmer(log(middle_GPP+1) ~ 
                       scale(lag_middle_GPP) + 
                       (1|site), data=GPP_df1)

model <- psem(lm((middle_GPP) ~ (lag_middle_GPP) + (flow_mean) + Stream_temp, data=GPP_df1),
              lm(flow_mean ~  (ppt_mm), data=GPP_df1),
              lm(Stream_temp ~  (stream_SPC), data=GPP_df1)) 


summary(model)
        



# Create model with correlation structure
GPP.sem <- psem(
  spring.kelp <- gls(spring.kelp ~ wave * prev.kelp + reef.cover,
                     correlation = corCAR1(form = ~ YEAR | SITE / TRANSECT),
                     data = byrnes),
  kelp <- gls(summer.kelp ~ wave * prev.kelp + reef.cover + spring.kelp,
              correlation = corCAR1(form = ~ YEAR | SITE / TRANSECT),
              data = byrnes),
  richness <- gls(richness ~ summer.kelp + prev.kelp + reef.cover +
                    spring.kelp,
                  correlation = corCAR1(form = ~ YEAR | SITE / TRANSECT),
                  data = byrnes),
  linkdensity <- gls(linkdensity ~ richness + summer.kelp + prev.kelp +
                       reef.cover + spring.kelp,
                     correlation = corCAR1(form = ~ YEAR | SITE / TRANSECT),
                     data = byrnes)
)





#####-----Advances in piecewise estimation of path models-----------------------

# EFI/ESA Seminar
# Date: 05 December 2022
# Author: Jon Lefcheck
# Contact: LefcheckJ@si.edu

# Load required libraries
# devtools::install_github("jslefche/piecewiseSEM@devel") # version 2.3.0
library(piecewiseSEM)
library(lavaan)

# Load Keeley data set
data(keeley)

# Examine Keeley data
head(keeley)

#####-----Fit structural equation model to Keeley data--------------------------

# Fit simple/multiple regressions using `lm`

# Break down component regressions
abiotic_model <- lm(abiotic ~ distance, data = keeley)
summary(abiotic_model)
hetero_model <- lm(hetero ~ distance, data = keeley)
summary(hetero_model)
richness_model <- lm(rich ~ abiotic + hetero, data = keeley)
summary(richness_model)
# Use the `psem` function to create the SEM
model <- psem(abiotic_model, hetero_model, richness_model)

# Look at object
model # Sometimes it works better to specifiy models in line, with called df

# Step 1: conduct tests of directed separation
# Establish the basis set & evaluate independence claims
# Missing path #1:
dsep1 <- lm(abiotic ~ hetero + distance, data = keeley)
summary(dsep1)
# Missing path #2:
dsep2 <- lm(rich ~ distance + abiotic + hetero, data = keeley)
summary(dsep2)

# Get P-values
P1 <- summary(dsep1)$coefficients[2, "Pr(>|t|)"]
P2 <- summary(dsep2)$coefficients[2, "Pr(>|t|)"]

# Construct C-statistic
C <- -2 * (log(P1) + log(P2))

C

fisherC(model) # says model doesn't match associations of data  -
# dsep2 model path between distance and richness might drive model fit 

# Compare to chi-squared distribution with 2*2 degrees of freedom
1 - pchisq(C, 4) # P < 0.05 == poor fit!

# Can use `dsep` function to perform the tests automagically
dSep(model)

# Can use `fisherC` function to evaluate claims
fisherC(model)

# The relationship between rich and distance is significant
# Re-introduce to the model
model2 <- update(model, rich ~ abiotic + hetero + distance)

model2

dSep(model2) # only 1 claim now
fisherC(model2) # P > 0.05 == model fits well!

# Get coefficients from good-fitting SEM
coefs(model2)

# Standardized estimates are in units of standard deviations of the mean
#     unit less - so can compare standardized estimates of all paths 
#     So you can use to multiply indirect effects 
# Can be directly compared even though initial units are very different

# Plot SEM with standardized coefficients
plot(model2)

# Use `summary` function to get all information at once
summary(model2)

##### Chi-squared goodness of fit

# Evaluate fit using Chi-squared

# Sum Log-likelihoods from original model
LL_sem <- sum(sapply(model2[-4], logLik))

# Fit saturated model
saturated_model <- psem(
  lm(abiotic ~ distance + hetero, data = keeley),
  lm(hetero ~ distance, data = keeley),
  lm(rich ~ abiotic + hetero + distance, data = keeley)
)

# Sum Log-likelihoods from saturated model
LL_sat <- sum(sapply(saturated_model[-4], logLik))

Chi_sq <- -2 * (LL_sem - LL_sat)

Chi_sq

# Compare to chi-squared distribution with 1 df (one additional estimated 
# parameter in saturated model)
1 - pchisq(Chi_sq, 1) # P > 0.05 == good fit!

fisherC(model2) # slightly different than P-value from Fisher's C

# Re-fit in lavaan
library(lavaan)

form <- '
abiotic ~ distance 
hetero ~ distance
rich ~ abiotic + hetero + distance
'

sem(form, keeley) # same P-value!

# Can we test whether the model with the distance -> rich path 
# is statistically better?
AIC(model, model2)

anova(model, model2) # Chi-square difference test

#####-----Extensions to non-linear models---------------------------------------

# Re-fit Keeley example with GLM
model3 <- psem(
  lm(abiotic ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  glm(rich ~ abiotic + hetero + distance, family = poisson(link = "log"), 
      data = keeley)
)

# Get summary
summary(model3)

# Compare with SEM of just LM 
anova(model2, model3) # GLM actually is less likely!

# Imagine that distance -> hetero relationship is truly nonlinear
# Re-fit SEM using generalized additive model (GAM)
library(mgcv)

model4 <- psem(
  lm(abiotic ~ distance, data = keeley),
  gam(hetero ~ s(distance), data = keeley),
  glm(rich ~ abiotic + hetero + distance, family = poisson(link = "log"), 
      data = keeley)
)

# Get summary 
summary(model4)

# Formally compare
anova(model2, model4) # Again, linear model is best





##########
# global model with everything 
library(glmmTMB)

############################
##make composite for stream 
stream_model <- lm(middle_GPP ~ Stream_temp + stream_flow + stream_SPC, data = GPP_df2)
# Get loadings (regression coefficients)
beta_Stream_temp <- summary(stream_model)$coefficients["Stream_temp", "Estimate"]
beta_stream_flow <- summary(stream_model)$coefficients["stream_flow", "Estimate"]
beta_stream_SPC <- summary(stream_model)$coefficients["stream_SPC", "Estimate"]


GPP_df2$composite <- with(GPP_df2, beta_Stream_temp * Stream_temp +
                    beta_stream_flow * stream_flow +
                    beta_stream_SPC * stream_SPC)

summary(lm(middle_GPP ~ composite, data = GPP_df2))




hist(GPP_df2$middle_GPP)
hist(log(GPP_df2$middle_GPP+1))
hist(GPP_df2$lake_par_int)

# 
# gpp_fit1 <- glmmTMB((middle_GPP) ~ lag((middle_GPP), 1) + (1|site), data = GPP_df2, family = gaussian)
# summary(gpp_fit1)
# residuals <- residuals(gpp_fit1)
# hist(residuals)

gpp_fit1 <- glmmTMB(middle_GPP ~ lag((middle_GPP), 1), data = GPP_df1, family = gaussian)
summary(gpp_fit1)
residuals <- residuals(gpp_fit1)
hist(residuals)
# Plot ACF and PACF of residuals
acf(residuals, main="ACF of Residuals")
pacf(residuals, main="PACF of Residuals")
# 
# hist(GPP_df1$lake_temp)

##### WORKING HERE ######

GPP_df1$yearf<- as.factor(year(GPP_df1$date))
GPP_df1$GPP <- log(GPP_df1$middle_GPP+1)
GPP_df1$GPP_lag <- log(GPP_df1$lag_middle_GPP+1)

hist(GPP_df1$GPP)
names(GPP_df1)

unique(GPP_df1$shore)
unique(GPP_df1$site)

hist(GPP_df1$flow_mean)
hist(GPP_df1$flow_mean)
hist(log(GPP_df1$stream_flow)+1)
hist(GPP_df1$ppt_mm)
summary((GPP_df1$ppt_mm))

GPP_df1$precip_bi <- ifelse(GPP_df1$ppt_mm > 0, 1, 0)

hist(log(GPP_df1$ppt_mm+1))

GPP_df1$log_ppmt <- log(GPP_df1$ppt_mm+1)
GPP_df1$log_streamflow <- (log(GPP_df1$flow_mean)+1)
names(GPP_df1)

Streamdf_cor <- GPP_df1[, c("Stream_temp", "stream_flow", "stream_SPC", "log_streamflow", "log_ppmt")]

# 
# GPP_df_scale <- GPP_df1 %>%
#   scale(7:26, 28:43, 46)
# 
# GPP_df1[, c(7:26, 28:43, 46)] <- scale(GPP_df1[, c(7:26, 28:43, 46)])
# summary(GPP_df1)

hist((GPP_df1$lake_par))

GPP_df1<- GPP_df1%>%
  drop_na(site)

GPP_df1<- GPP_df1%>%
  drop_na(GPP)






library(PerformanceAnalytics)
names(GPP_df1)
Streamdf_cor <- GPP_df1[, c(17:26)]
Streamdf_cor <- GPP_df1[, c("Stream_temp", "ppt_mm", "stream_flow", "stream_SPC")]
Streamdf_cor <- GPP_df1[, c("Stream_temp", "stream_flow", "stream_SPC")]
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)


lakedf_cor <- GPP_df1[, c("lake_wtemp", "light_mean", "tmean_C", "windsp_mean", "flow_mean", "stream_temp", "spc_mean")]
chart.Correlation(lakedf_cor, histogram=TRUE, pch=19)




hist(GPP_df1$lake_wtemp)

# Remove the first row from the dataframe
lake_temp_mod <- lmer(lake_wtemp ~ 
                        windsp_mean +
                        precip_bi +
                        flow_mean +
                        light_mean +
                        spc_mean + 
                        (1|site), data = GPP_df1)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))
vif(lake_temp_mod)


lake_temp_mod2 <- lmer(lake_wtemp ~ 
                        windsp_mean +
                        precip_bi +
                        flow_mean +
                        light_mean +
                        spc_mean + 
                        (1|site), data = GPP_df1)
summary(lake_temp_mod2)
hist(residuals(lake_temp_mod2))
vif(lake_temp_mod2)



lake_temp_mod3 <- lmer(lake_wtemp ~ 
                         windsp_mean +
                         precip_bi +
                         flow_mean +
                         tmean_C +
                         spc_mean + 
                         (1|site), data = GPP_df1)
summary(lake_temp_mod3)
hist(residuals(lake_temp_mod3))
vif(lake_temp_mod3)

AIC(lake_temp_mod, lake_temp_mod2, lake_temp_mod3)


## Look at best fit lake temp model:
lake_temp_mod
hist(GPP_df1$lake_temp)
G_plot <- ggplot(data = GPP_df1, aes(x= (windsp_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot

G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot


G_plot <- ggplot(data = GPP_df1, aes(x= (tmean_C), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot


G_plot <- ggplot(data = GPP_df1, aes(x= (spc_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot

G_plot <- ggplot(data = GPP_df1, aes(x= (precip_bi), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot

### flow_mean# ##


## 
lakelight_cor <- GPP_df1[c("lake_par_int", "light_cv", "flow_mean", "precip_bi", "lake_wspeed", "lake_depth") ]
chart.Correlation(lakelight_cor, histogram=TRUE, pch=19)
hist(GPP_df1$lake_par_int)

lake_light_mod <- lmer(lake_par_int ~ 
                         flow_mean +
                         precip_bi +
                         windsp_mean +
                         lake_depth + (1|site), data = GPP_df1)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)


# DAG chapter in rethinking -
GPP_df1$yearf<- as.factor(year(GPP_df1$date))

gpp_auto1 <- lmer(middle_GPP ~  lag_middle_GPP + 
                    (lake_wtemp) + 
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto1)
residuals <- residuals(gpp_auto1)

gpp_auto2 <- lmer(middle_GPP ~  lag_middle_GPP + 
                    (lake_wtemp) + 
                    precip_bi +
                    windsp_mean +
                    flow_mean +
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto2)
residuals <- residuals(gpp_auto2)
vif(gpp_auto2)
hist(residuals(gpp_auto2))
r.squaredGLMM(gpp_auto2)


gpp_auto3 <- lmer(middle_GPP ~  lag_middle_GPP + 
                    (lake_wtemp) + 
                    precip_bi +
                    lake_wspeed +
                    #flow_mean +
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto3)
residuals <- residuals(gpp_auto3)
vif(gpp_auto3)
hist(residuals(gpp_auto3))

gpp_auto4 <- lmer(middle_GPP ~  lag_middle_GPP + 
                    (lake_wtemp) + 
                    precip_bi +
                    #lake_wspeed +
                    #flow_mean +
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto4)
residuals <- residuals(gpp_auto4)
vif(gpp_auto4)
hist(residuals(gpp_auto4))

AIC(gpp_auto1, gpp_auto2, gpp_auto3, gpp_auto4)

# Use the `psem` function to create the SEM
model <- psem(
  gpp_auto1,
  lake_temp_mod3,  
  lake_light_mod)

# Look at object
model # Sometimes it works better to specify models in line, with called df
summary(model)
plot(model)
fisherC(model) # says model

model2 <- psem(
  gpp_auto2,
  lake_temp_mod3,  
  lake_light_mod)

summary(model2)
plot(model2)

AIC(model, model2)

# how to interpret conditional independence statements ?
# could be worth collapsing the weekly level then creating auto regressive structures - re-watch re-thinking lecture 

### model 2 ###
lakeT_cor <- GPP_df1[, c("lake_temp","Stream_temp", "windsp_cv", "ppt_mm", "light_mean", "tmax_C", "tmean_C", "vpdmax_hPa")]
chart.Correlation(lakeT_cor, histogram=TRUE, pch=19)




# Remove the first row from the dataframe
lake_temp_mod1 <- lmer(lake_temp ~ 
                        windsp_cv +
                        ppt_mm +
                        Stream_temp +
                        #light_mean + 
                        #stream_SPC + 
                        (1|site) + (1|shore), data = GPP_df1)
summary(lake_temp_mod1)
hist(residuals(lake_temp_mod1))
vif(lake_temp_mod1)



# Remove the first row from the dataframe
lake_temp_mod2 <- lmer(lake_temp ~ 
                         windsp_cv +
                         ppt_mm +
                         tmean_C +
                         #light_mean + 
                         #stream_SPC + 
                         (1|site) + (1|shore), data = GPP_df1)
summary(lake_temp_mod2)
hist(residuals(lake_temp_mod2))
vif(lake_temp_mod2)


stream_temp_mod <- lmer(Stream_temp ~ 
                          stream_flow + 
                          stream_SPC + 
                          (1|site) + (1|shore), data = GPP_df1)
summary(stream_temp_mod)
hist(residuals(stream_temp_mod))
vif(stream_temp_mod)


## 
lakelight_cor <- GPP_df1[, c(8, 7)]
chart.Correlation(lakelight_cor, histogram=TRUE, pch=19)

lake_light_mod <- lmer(lake_par_int ~ 
                         stream_flow +
                         ppt_mm +
                         # stream_SPC + 
                         lake_z+
                         windsp_cv + (1|site) +(1|shore), data = GPP_df1)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)


# DAG chapter in rethinking -
GPP_df1$yearf<- as.factor(year(GPP_df1$date))
gpp_auto1 <- lmer(middle_GPP ~  lag_middle_GPP + (lake_temp) + (lake_par_int) +(1|site), data = GPP_df1)
summary(gpp_auto1)
residuals <- residuals(gpp_auto1)

gpp_auto2 <- lmer(GPP ~  (GPP_lag) + 
                    (lake_temp) + 
                    (lake_par_int) +
                    #(windsp_cv) +
                    #ppt_mm +
                    #(Stream_temp) +
                    #(stream_flow) + 
                    (1|site)+(1|shore), data = GPP_df1)

summary(gpp_auto2)
residuals <- residuals(gpp_auto2)
vif(gpp_auto2)
hist(residuals(gpp_auto2))

AIC(gpp_auto2_all, gpp_auto2_1, gpp_auto2_2, gpp_auto2_3, gpp_auto2_4)

# Use the `psem` function to create the SEM
model2 <- psem(
  #gpp_mod,
  gpp_auto2,
  stream_temp_mod,
  #lake_temp_mod1,
  lake_temp_mod2,
  lake_light_mod)

# Look at object
model2 # Sometimes it works better to specify models in line, with called df
summary(model2)
plot(model2)
fisherC(model2) # says model


model3 <- psem(
  #gpp_mod,
  gpp_auto2,
  stream_temp_mod,
  lake_temp_mod1,
  #lake_temp_mod2,
  lake_light_mod)

# Look at object
model3 # Sometimes it works better to specify models in line, with called df
summary(model3)
plot(model3)
fisherC(model3) # says model


AIC(model3, model2, model)










##########
###########
## ER ###
##########
ER_df$middle_ER_lag <- lag(ER_df$middle_ER)

ER_df_df2 <- ER_df%>%
  dplyr::select(site, middle_ER, middle_ER_lag, lake_temp, lake_par_int, 
                ppt_mm, windsp_cv, lake_z,
                Stream_temp, stream_flow, stream_SPC)


ER_df_df2<- na.omit(ER_df_df2)
hist(ER_df_df2$middle_ER)


##### WORKING HERE ######
# Remove the first row from the dataframe
lake_temp_mod <- lmer(lake_temp ~ 
                        #windsp_cv +
                        #ppt_mm +
                        Stream_temp +
                        #stream_flow + 
                        stream_SPC + (1|site), data = ER_df_df2)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))
vif(lake_temp_mod)

## 
lake_light_mod <- lmer(lake_par_int ~ 
                         lake_z +
                         stream_flow + 
                         #stream_SPC + 
                         windsp_cv + (1|site), data = ER_df_df2)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)

# this needs work 
ER_mod <- lmer(middle_ER ~ (lake_temp) + (lake_par_int) + (1|site), data = ER_df_df2)
summary(ER_mod)
hist(residuals(ER_mod))


ER_auto1 <- lmer(middle_ER ~  middle_ER_lag + (lake_temp) + (lake_par_int) +(1|site), data = ER_df_df2)
summary(ER_auto1)
residuals <- residuals(ER_auto1)

# Use the `psem` function to create the SEM
ER_model <- psem(
  ER_auto1,
  lake_temp_mod,  
  lake_light_mod)

# Look at object
ER_model # Sometimes it works better to specify models in line, with called df
summary(ER_model)
plot(ER_model)
fisherC(ER_model) # says model
