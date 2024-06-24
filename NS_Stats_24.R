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

getwd()
SFS_datQ <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_analysis_dat/SFS24_analysis_dat.rds")
summary(SFS_datQ)
str(SFS_datQ)

SFS_datQ<- SFS_datQ%>%
  drop_na(site)

unique(SFS_datQ$shore)
unique(SFS_datQ$site)


##===========================================
## useful variable creation
SFS_datQ1 <- SFS_datQ %>%
  mutate(
    position = case_when(
      site %in% c("BWNS1","GBNS1", "SSNS1", "SHNS1") ~ "north",
      site %in% c("BWNS2","GBNS2", "SSNS2", "SHNS2" ) ~ "center",
      site %in% c("BWNS3","GBNS3", "SHNS3") ~ "south",
      TRUE ~ as.character(site)
    )
  )

SFS_datQ1$precip_bi <- ifelse(SFS_datQ1$ppt_mm > 0, 1, 0)
SFS_datQ1$log_ppmt <- log(SFS_datQ1$ppt_mm+1)
SFS_datQ1$log_streamflow <- (log(SFS_datQ1$flow_mean+1))

## create column for year and doy 
SFS_datQ1$year <- year(SFS_datQ1$date)
SFS_datQ1$yday <- yday(SFS_datQ1$date)
SFS_datQ1$week <- week(SFS_datQ1$date)


SFS_datQ2 <- SFS_datQ1 %>%
  filter(middle_GPP<30 &  middle_ER>-30)
summary(SFS_datQ2)
##===========================================
## weekly aggreation  
SFS_week_dat <- SFS_datQ2 %>%
  group_by(site, shore, week, year, position) %>%
  summarise(
    middle_GPP=mean(middle_GPP, na.rm = TRUE),
    lower_GPP=mean(lower_GPP, na.rm = TRUE),
    upper_GPP=mean(upper_GPP, na.rm = TRUE),
    middle_ER=mean(middle_ER, na.rm = TRUE),
    lower_ER=mean(lower_ER, na.rm = TRUE),
    upper_ER=mean(upper_ER, na.rm = TRUE),
    ## Weather 
    tmean_C=mean(tmean_C, na.rm = TRUE),
    light_mean=mean(light_mean, na.rm = TRUE),
    windsp_mean=mean(windsp_mean, na.rm = TRUE),
    ppt_mm=mean(ppt_mm, na.rm = TRUE),
    log_ppmt=mean(log_ppmt, na.rm = TRUE),
    ppmt_sum=sum(precip_bi, na.rm = TRUE),
    ## lake quality 
    lake_tempC=mean(lake_tempC, na.rm = TRUE),
    lake_DO=mean(lake_DO, na.rm = TRUE),
    lake_SPC=mean(lake_SPC, na.rm = TRUE),
    Kd_fill=mean(Kd_fill, na.rm = TRUE),
    par_int_3m=mean(par_int_3m, na.rm = TRUE),
    real_NS_depth=mean(real_NS_depth, na.rm = TRUE),
    ## Stream quality 
    flow_mean_m=mean(flow_mean, na.rm = TRUE),
    log_streamflow=mean(log_streamflow, na.rm = TRUE),
    flow_sum=sum(flow_mean * 86400, na.rm = TRUE),
    stream_SPC=mean(stream_SPC, na.rm = TRUE),
    stream_DO=mean(stream_DO, na.rm = TRUE),
    stream_temp=mean(stream_temp, na.rm = TRUE))
    
summary(SFS_week_dat)
##===========================================
SFS_week_dat_BWt <- SFS_week_dat%>%
  filter(shore=="BW")

summary(SFS_week_dat_BWt)

SFS_week_dat_GBt <- SFS_week_dat%>%
  filter(shore=="GB")

summary(SFS_week_dat_GBt)

(292.0-44.065)/292.0
##===========================================



## 
## create a new df for complete GPP obs 
GPP_df <- SFS_datQ2%>%
  dplyr::select( -middle_ER, -upper_ER, -lower_ER)

GPP_df<- GPP_df%>%
  drop_na(middle_GPP)

ER_df <- SFS_datQ2%>%
  dplyr::select(-middle_GPP, -upper_GPP, -lower_GPP)

ER_df<- ER_df%>%
  drop_na(middle_ER)


ER_df$ER <- c(ER_df$middle_ER *-1)
ER_df$ER_low <- c(ER_df$lower_ER *-1)
ER_df$ER_up <- c(ER_df$upper_ER *-1)
#GPP_df1<- na.omit(GPP_df)
summary(GPP_df)
summary(ER_df)
##===========================================
## Begin analysis 
#============================================
# lake condition across sites
# When is GPP highest? lowest? 
# Test for different in Metab across shores 

BW_dat1 <- SFS_datQ %>%
  filter(shore=="BW" & date>as.Date("2023-06-01") & date<as.Date("2023-09-15"))


BW_dat1 <- BW_dat %>%
  filter(date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))
unique(BW_dat$shore)
round(mean(na.omit(BW_dat1$middle_GPP)),2)
round(se(na.omit(BW_dat1$middle_GPP)),2)
range(na.omit(BW_dat$middle_GPP))
round(mean(na.omit(BW_dat1$middle_ER)),2)
round(se(na.omit(BW_dat1$middle_ER)),2)
round(se(na.omit(BW_dat$middle_ER)),2)
range(na.omit(BW_dat$middle_ER))

round(mean(na.omit(BW_dat$middle_NEP)),2)
round(se(na.omit(BW_dat$middle_NEP)),2)

summary(BW_dat)
max_index <- max(BW_dat$middle_GPP)
print(BW_dat[max_index, ])
min_index <- which.min((BW_dat$middle_GPP))
print(BW_dat[min_index, ])


SS_dat <- SFS_datQ %>%
  filter(shore=="SS")
SS_dat1 <- SS_dat %>%
  filter(date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))
unique(SS_dat$shore)
round(mean(na.omit(SS_dat1$middle_GPP)),2)
round(se(na.omit(SS_dat1$middle_GPP)),2)
round(mean(na.omit(SS_dat1$middle_ER)),2)
round(se(na.omit(SS_dat1$middle_ER)),2)
summary(SS_dat)


GB_dat1 <- SFS_datQ %>%
  filter(shore=="GB"  & date>as.Date("2023-06-01") & date<as.Date("2023-09-15"))

GB_dat <- SFS_datQ %>%
  filter(shore=="GB")
GB_dat1 <- GB_dat %>%
  filter(date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))
unique(GB_dat$shore)
round(mean(na.omit(GB_dat1$middle_GPP)),2)
round(se(na.omit(GB_dat1$middle_GPP)),2)
range(na.omit(GB_dat1$middle_GPP))
round(mean(na.omit(GB_dat1$middle_ER)),2)
round(se(na.omit(GB_dat1$middle_ER)),2)
summary(GB_dat)
range(na.omit(GB_dat$middle_ER))



round(mean(na.omit(GB_dat1$middle_GPP)),2)
round(se(na.omit(GB_dat1$middle_GPP)),2)
range(na.omit(GB_dat1$middle_GPP))
round(mean(na.omit(GB_dat1$middle_ER)),2)
round(se(na.omit(GB_dat1$middle_ER)),2)


SH_dat <- SFS_datQ %>%
  filter(shore=="SH")

SH_dat1 <- SH_dat %>%
filter(date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))

unique(SH_dat$shore)
round(mean(na.omit(SH_dat1$middle_GPP)),2)
round(se(na.omit(SH_dat1$middle_GPP)),2)
round(mean(na.omit(SH_dat1$middle_ER)),2)
round(se(na.omit(SH_dat1$middle_ER)),2)
summary(SH_dat)

#============================================
# Figure out auto regressive structure:
#============================================
library(glmmTMB)

hist(GPP_df$middle_GPP)
hist(log(GPP_df$middle_GPP+1))

gpp_fit1 <- glmmTMB((middle_GPP) ~ 
                      #lag(middle_GPP, 1) + 
                      lag(middle_GPP, 2) + 
                      lag(middle_GPP, 4) + 
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit1)
residuals <- residuals(gpp_fit1)
hist(residuals)
# Plot ACF and PACF of residuals
acf(residuals, main="ACF of Residuals")
pacf(residuals, main="PACF of Residuals")


gpp_fit2 <- glmmTMB((middle_GPP) ~  lag(middle_GPP, 2) +
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit2)
residuals2 <- residuals(gpp_fit2)
hist(residuals2)
# Plot ACF and PACF of residuals
acf(residuals2, main="ACF of Residuals")
pacf(residuals2, main="PACF of Residuals")


gpp_fit3 <- glmmTMB((middle_GPP) ~  lag(middle_GPP, 1) +
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit3)
residuals2 <- residuals(gpp_fit3)
hist(residuals2)
# Plot ACF and PACF of residuals
acf(residuals2, main="ACF of Residuals")
pacf(residuals2, main="PACF of Residuals")
AIC(gpp_fit1, gpp_fit2, gpp_fit3)

##===========================================
## GPP glms for stream influence
#============================================
# subset for shore BW and GB
GPP_df_stream <- GPP_df %>%
  filter(shore=="BW" |shore=="GB")

GPP_df_stream$middle_GPP_lag <- lag(GPP_df_stream$middle_GPP)


ER_df_stream <- ER_df %>%
  filter(shore=="BW" |shore=="GB")

ER_df_stream$middle_GPP_lag <- lag(ER_df_stream$middle_ER)



##===========================================
## check co-variance 
names(GPP_df_stream)

# scale variables
# GPP_df_stream[, c(7, 19:26)] <- scale(GPP_df_stream[, c(7, 19:26)])
Streamdf_cor <- scale(GPP_df_stream[, c(9:26)])
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

Streamdf_cor <- scale(GPP_df_stream[, c(6, 12:14, 22:26)])
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

## ER ##
ER_df_stream[, c(7, 19:26)] <- scale(ER_df_stream[, c(7, 19:26)])
Streamdf_cor <- ER_df_stream[, c(7, 19:26)]
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

##==========================
### MODEL ###
## NULL:
gpp_stream_null <- glmmTMB((middle_GPP) ~ (1|shore / site) +(1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_null)
residuals3 <- residuals(gpp_stream_null)
hist(residuals3)

##===========================================
## Predictors for stream 
gpp_stream_mod <- glmmTMB((middle_GPP) ~ 
                             scale(stream_temp)  + 
                             scale(flow_sum)+
                             scale(ppmt_sum)+
                             (1|shore / site) +(1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)


gpp_stream_mod1 <- glmmTMB((middle_GPP) ~ 
                            scale(stream_temp)  + 
                            scale(flow_sum)+
                            scale(log_ppmt)+
                            (1|shore / site) +(1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod1)
residuals3 <- residuals(gpp_stream_mod1)
hist(residuals)
r.squaredGLMM(gpp_stream_mod1)


gpp_stream_mod2 <- glmmTMB((middle_GPP) ~ 
                             scale(stream_temp)  + 
                             scale(flow_mean)+
                             scale(log_ppmt)+
                             (1|shore / site) +(1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod2)
residuals3 <- residuals(gpp_stream_mod2)
hist(residuals)
r.squaredGLMM(gpp_stream_mod2)

AIC(gpp_stream_null, gpp_stream_mod, gpp_stream_mod1, gpp_stream_mod2)


##===========================================
## Predictors for lake
lakedf_cor <- scale(GPP_df_stream[, c(6, 9:20)])
chart.Correlation(lakedf_cor, histogram=TRUE, pch=19)
## ALL light correlations are negative -- maybe photo inhibition 

# incoming light is 0.73 correlated to air temp, and 0.94 correlated to light at depth , and 0.52 to lake temp

gpp_stream_mod <- glmmTMB((middle_GPP) ~ 
                            scale(light_mean)  + 
                            scale(lake_tempC)+ # maybe too correlated 
                            scale(real_NS_depth)+
                            scale(windsp_mean)+
                            scale(ppmt_sum) +
                             (1|shore / site) +(1|year), data = GPP_df, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)



gpp_stream_mod <- glmmTMB((middle_GPP) ~ 
                            scale(light_mean)  + 
                           # scale(lake_tempC)+ # maybe too correlated 
                            scale(real_NS_depth)+
                            scale(windsp_mean)+
                            scale(ppmt_sum) +
                            (1|shore / site) +(1|year), data = GPP_df, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)



gpp_stream_mod <- glmmTMB((middle_GPP) ~ 
                            scale(light_mean)  + 
                            # scale(lake_tempC)+ # maybe too correlated 
                            scale(real_NS_depth)+
                            scale(windsp_mean)+
                            #scale(ppmt_sum) +
                            (1|shore / site) +(1|year), data = GPP_df, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)


gpp_stream_mod <- glmmTMB((middle_GPP) ~ 
                            scale(light_mean)  + 
                            # scale(lake_tempC)+ # maybe too correlated 
                            scale(real_NS_depth)+
                            #scale(windsp_mean)+
                            #scale(ppmt_sum) +
                            (1|shore / site) +(1|year), data = GPP_df, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals)
r.squaredGLMM(gpp_stream_mod)

##===========================================
## Predictors for stream -- lmer
gpp_stream_lmer_mod1 <- lmer(log(middle_GPP+1) ~ 
                             scale(stream_temp)  + 
                             scale(flow_sum)+
                             scale(log_ppmt)+
                             (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_lmer_mod1)
residuals3 <- residuals(gpp_stream_lmer_mod1)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod1)
vif(gpp_stream_lmer_mod1)

gpp_stream_lmer_mod2 <- lmer(log(middle_GPP+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_mean_m)+
                               scale(log_ppmt)+
                               (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_lmer_mod2)
residuals3 <- residuals(gpp_stream_lmer_mod2)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod2)
vif(gpp_stream_lmer_mod2)


gpp_stream_lmer_mod3 <- lmer(log(middle_GPP+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_mean_m)+
                               scale(ppmt_sum)+
                               (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_lmer_mod3)
residuals3 <- residuals(gpp_stream_lmer_mod3)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod3)
vif(gpp_stream_lmer_mod3)



gpp_stream_lmer_mod4 <- lmer(log(middle_GPP+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_sum)+
                               #scale(ppmt_sum)+
                               (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_lmer_mod4)
residuals3 <- residuals(gpp_stream_lmer_mod4)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod4)
vif(gpp_stream_lmer_mod4)

AIC(gpp_stream_mod_null, gpp_stream_lmer_mod1, gpp_stream_lmer_mod2, gpp_stream_lmer_mod3, gpp_stream_lmer_mod4)

# best fit model:
# gpp_stream_lmer_mod4


##===========================================
## Stats wrap up:

summary(gpp_stream_lmer_mod4)
gpp_stream_lmer_mod4

GPP_df_stream_BW <- GPP_df_stream%>%
  filter(shore=="BW") 

summary(GPP_df_stream_GB)

gpp_stream_lmer_mod4 <- lmer(log(middle_GPP+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_sum)+
                               #scale(ppmt_sum)+
                               (1| site) +(1|year), data = GPP_df_stream_BW)
summary(gpp_stream_lmer_mod4)
residuals3 <- residuals(gpp_stream_lmer_mod4)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod4)
vif(gpp_stream_lmer_mod4)

GPP_df_stream_GB <- GPP_df_stream%>%
  filter(shore=="GB") 

gpp_stream_lmer_mod4 <- lmer(log(middle_GPP+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_sum)+
                               #scale(ppmt_sum)+
                               (1| site) +(1|year), data = GPP_df_stream_GB)
summary(gpp_stream_lmer_mod4)
residuals3 <- residuals(gpp_stream_lmer_mod4)
hist(residuals3)
r.squaredGLMM(gpp_stream_lmer_mod4)
vif(gpp_stream_lmer_mod4)

GPP_df_stream_BW1 <- GPP_df_stream%>%
  filter(shore=="BW" & site=="BWNS1") 
summary(GPP_df_stream_BW1)

GPP_df_stream_BW2 <- GPP_df_stream%>%
  filter(shore=="BW" & site=="BWNS2") 
summary(GPP_df_stream_BW2)

GPP_df_stream_BW3 <- GPP_df_stream%>%
  filter(shore=="BW" & site=="BWNS3") 
summary(GPP_df_stream_BW3)

GPP_df_stream_GB1 <- GPP_df_stream%>%
  filter(shore=="GB" & site=="GBNS1") 
summary(GPP_df_stream_GB1)

GPP_df_stream_GB2 <- GPP_df_stream%>%
  filter(shore=="GB" & site=="GBNS2") 
summary(GPP_df_stream_GB2)

GPP_df_stream_GB3 <- GPP_df_stream%>%
  filter(shore=="GB" & site=="GBNS3") 
summary(GPP_df_stream_GB3)

(11624.6-462.9)/11624.6

(5.424- 4.40778)/5.424
(4.40778-5.424)/4.40778

##===========================================
## Predictors for lake -- lmer
gpp_stream_mod_null <- lmer(log(middle_GPP+1) ~ 
                              (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_null)
residuals3 <- residuals(gpp_stream_mod_null)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_null)

gpp_stream_mod_lmer <- lmer(log(middle_GPP+1) ~ 
                            scale(light_mean)  + 
                            scale(lake_tempC)+ # maybe too correlated 
                            scale(real_NS_depth)+
                            scale(windsp_mean)+
                            scale(ppmt_sum) +
                            (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer)
residuals3 <- residuals(gpp_stream_mod_lmer)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer)
vif(gpp_stream_mod_lmer)


gpp_stream_mod_lmer1 <- lmer(log(middle_GPP+1) ~ 
                             # scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer1)
residuals3 <- residuals(gpp_stream_mod_lmer1)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer1)
vif(gpp_stream_mod_lmer1)


gpp_stream_mod_lmer1 <- lmer(log(middle_GPP+1) ~ 
                               # scale(light_mean)  + 
                               scale(lake_tempC)+ # maybe too correlated 
                               #scale(real_NS_depth)+
                               scale(windsp_mean)+
                               scale(ppmt_sum) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer1)
residuals3 <- residuals(gpp_stream_mod_lmer1)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer1)
vif(gpp_stream_mod_lmer1)



gpp_stream_mod_lmer2 <- lmer(log(middle_GPP+1) ~ 
                               # scale(light_mean)  + 
                               scale(lake_tempC)+ # maybe too correlated 
                               #scale(real_NS_depth)+
                               scale(windsp_mean)+
                              # scale(ppmt_sum) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer2)
residuals3 <- residuals(gpp_stream_mod_lmer1)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer2)
vif(gpp_stream_mod_lmer2)

gpp_stream_mod_lmer3 <- lmer(log(middle_GPP+1) ~ 
                               scale(light_mean)  + 
                               #scale(lake_tempC)+ # maybe too correlated 
                               scale(real_NS_depth)+
                               scale(windsp_mean)+
                              scale(ppmt_sum) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer3)
residuals3 <- residuals(gpp_stream_mod_lmer3)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer3)
vif(gpp_stream_mod_lmer3)

gpp_stream_mod_lmer4 <- lmer(log(middle_GPP+1) ~ 
                               scale(light_mean)  + 
                               #scale(lake_tempC)+ # maybe too correlated 
                               scale(real_NS_depth)+
                               scale(windsp_mean)+
                               #scale(ppmt_sum) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer4)
residuals3 <- residuals(gpp_stream_mod_lmer4)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer4)
vif(gpp_stream_mod_lmer4)


gpp_stream_mod_lmer5 <- lmer(log(middle_GPP+1) ~ 
                               scale(light_mean)  + 
                               #scale(lake_tempC)+ # maybe too correlated 
                               scale(real_NS_depth)+
                               scale(windsp_mean)+
                               scale(ppt_mm) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer5)
residuals3 <- residuals(gpp_stream_mod_lmer5)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer5)
vif(gpp_stream_mod_lmer5)



gpp_stream_mod_lmer6 <- lmer(log(middle_GPP+1) ~ 
                               scale(light_mean)  + 
                               #scale(lake_tempC)+ # maybe too correlated 
                               scale(real_NS_depth)+
                               scale(windsp_mean)+
                               scale(log_ppmt) +
                               (1|shore / site) +(1|year), data = GPP_df)
summary(gpp_stream_mod_lmer6)
residuals3 <- residuals(gpp_stream_mod_lmer6)
hist(residuals)
r.squaredGLMM(gpp_stream_mod_lmer6)
vif(gpp_stream_mod_lmer6)

AIC(gpp_stream_mod_null, gpp_stream_mod_lmer, gpp_stream_mod_lmer1, gpp_stream_mod_lmer2, gpp_stream_mod_lmer3,
    gpp_stream_mod_lmer4, gpp_stream_mod_lmer5, gpp_stream_mod_lmer6)

# BEST Model:
# gpp_stream_mod_lmer4 




##===========================================
## Predictors for stream -- lmer
hist(ER_df_stream$ER)
hist(log(ER_df_stream$ER+1))


er_stream_mod_null <- lmer(log(ER+1) ~ 
                             (1|shore / site) +(1|year), data = ER_df_stream)
summary(er_stream_mod_null)
residuals3 <- residuals(er_stream_mod_null)
hist(residuals)
r.squaredGLMM(er_stream_mod_null)


er_stream_lmer_mod1 <- lmer(log(ER+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_sum)+
                               scale(log_ppmt)+
                               (1|shore / site) +(1|year), data = ER_df_stream)
summary(er_stream_lmer_mod1)
residuals3 <- residuals(er_stream_lmer_mod1)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_mod1)
vif(er_stream_lmer_mod1)

er_stream_lmer_mod2 <- lmer(log(ER+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_mean_m)+
                               scale(log_ppmt)+
                               (1|shore / site) +(1|year), data = ER_df_stream)
summary(er_stream_lmer_mod2)
residuals3 <- residuals(er_stream_lmer_mod2)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_mod2)
vif(er_stream_lmer_mod2)


er_stream_lmer_mod3 <- lmer(log(ER+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_mean_m)+
                               scale(ppmt_sum)+
                               (1|shore / site) +(1|year), data = ER_df_stream)
summary(er_stream_lmer_mod3)
residuals3 <- residuals(er_stream_lmer_mod3)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_mod3)
vif(er_stream_lmer_mod3)



er_stream_lmer_mod4 <- lmer(log(ER+1) ~ 
                               scale(stream_temp)  + 
                               scale(flow_sum)+
                               #scale(ppmt_sum)+
                               (1|shore / site) +(1|year), data = ER_df_stream)
summary(er_stream_lmer_mod4)
residuals3 <- residuals(er_stream_lmer_mod4)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_mod4)
vif(er_stream_lmer_mod4)

AIC(er_stream_mod_null, er_stream_lmer_mod1, er_stream_lmer_mod2, er_stream_lmer_mod3, er_stream_lmer_mod4)

# best fit model:
# gpp_stream_lmer_mod4
ER_df_stream_BW <- ER_df_stream%>%
  filter(shore=="BW") 

er_stream_lmer_modt1 <- lmer(log(ER+1) ~ 
                              scale(stream_temp)  + 
                              scale(flow_sum)+
                              #scale(ppmt_sum)+
                              (1|site) +(1|year), data = ER_df_stream_BW)
summary(er_stream_lmer_modt1)
residuals3 <- residuals(er_stream_lmer_modt1)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_modt1)
vif(er_stream_lmer_modt1)


ER_df_stream_GB <- ER_df_stream%>%
  filter(shore=="GB") 

er_stream_lmer_modt <- lmer(log(ER+1) ~ 
                              scale(stream_temp)  + 
                              scale(flow_sum)+
                             # scale(ppmt_sum)+
                              (1|site) +(1|year), data = ER_df_stream_GB)
summary(er_stream_lmer_modt)
residuals3 <- residuals(er_stream_lmer_modt)
hist(residuals3)
r.squaredGLMM(er_stream_lmer_modt)
vif(er_stream_lmer_modt)



##===========================================
## Predictors for lake -- lmer
er_stream_mod_null <- lmer(log(ER+1) ~ 
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_null)
residuals3 <- residuals(er_stream_mod_null)
hist(residuals)
r.squaredGLMM(er_stream_mod_null)

er_stream_mod_lmer <- lmer(log(ER+1) ~ 
                              scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer)
residuals3 <- residuals(er_stream_mod_lmer)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer)
vif(er_stream_mod_lmer)


er_stream_mod_lmer1 <- lmer(log(ER+1) ~ 
                             #scale(light_mean)  + 
                             scale(lake_tempC)+ # maybe too correlated 
                             scale(real_NS_depth)+
                             scale(windsp_mean)+
                             scale(ppmt_sum) +
                             (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer1)
residuals3 <- residuals(er_stream_mod_lmer1)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer1)
vif(er_stream_mod_lmer1)


er_stream_mod_lmer2 <- lmer(log(ER+1) ~ 
                              #scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              #scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer2)
residuals3 <- residuals(er_stream_mod_lmer2)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer2)
vif(er_stream_mod_lmer2)


er_stream_mod_lmer3 <- lmer(log(ER+1) ~ 
                              #scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              #scale(real_NS_depth)+
                              scale(windsp_mean)+
                              #scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer3)
residuals3 <- residuals(er_stream_mod_lmer3)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer3)
vif(er_stream_mod_lmer3)


er_stream_mod_lmer4 <- lmer(log(ER+1) ~ 
                              scale(light_mean)  + 
                              #scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer4)
residuals3 <- residuals(er_stream_mod_lmer4)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer4)
vif(er_stream_mod_lmer4)

er_stream_mod_lmer5 <- lmer(log(ER+1) ~ 
                              scale(light_mean)  + 
                              #scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              #scale(windsp_mean)+
                              scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer5)
residuals3 <- residuals(er_stream_mod_lmer5)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer5)
vif(er_stream_mod_lmer5)

er_stream_mod_lmer6 <- lmer(log(ER+1) ~ 
                              #scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(ppt_mm) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer6)
residuals3 <- residuals(er_stream_mod_lmer6)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer6)
vif(er_stream_mod_lmer6)



er_stream_mod_lmer7 <- lmer(log(ER+1) ~ 
                              #scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              scale(real_NS_depth)+
                              scale(windsp_mean)+
                              scale(log_ppmt) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer7)
residuals3 <- residuals(er_stream_mod_lmer7)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer7)
vif(er_stream_mod_lmer7)

er_stream_mod_lmer8 <- lmer(log(ER+1) ~ 
                              #scale(light_mean)  + 
                              scale(lake_tempC)+ # maybe too correlated 
                              #scale(real_NS_depth)+
                              #scale(windsp_mean)+
                              #scale(ppmt_sum) +
                              (1|shore / site) +(1|year), data = ER_df)
summary(er_stream_mod_lmer8)
residuals3 <- residuals(er_stream_mod_lmer8)
hist(residuals)
r.squaredGLMM(er_stream_mod_lmer8)
vif(er_stream_mod_lmer8)


AIC(er_stream_mod_null, er_stream_mod_lmer, er_stream_mod_lmer1, er_stream_mod_lmer2, er_stream_mod_lmer3,
    er_stream_mod_lmer4, er_stream_mod_lmer5, er_stream_mod_lmer6, er_stream_mod_lmer7,er_stream_mod_lmer8)

# best fit model so far
# er_stream_mod_lmer8



# Trendline
fitted_values <- predict(gpp_stream_lmer_mod4, newdata = GPP_df_stream)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(flow_sum = GPP_df_stream$flow_sum, middle_GPP = fitted_values)

# Plot
G_plot <- ggplot(data = GPP_df_stream, aes(x = flow_sum, y = log(middle_GPP+1))) +
  geom_point(size = 1.5, alpha = 0.6, aes(shape=position, color = shore)) +
  geom_errorbar(aes(ymin = log(lower_GPP + 1), ymax = log(upper_GPP + 1), color = shore), width = 0, alpha = 0.05) +
  geom_smooth(data = plot_data, aes(y = (middle_GPP)), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression(log(GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(cumulative~flow~m^3)) + 
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_flowplot_a.png",sep=""), width=4,height=3,dpi=300)


# Trendline
fitted_values <- predict(gpp_stream_lmer_mod4, newdata = GPP_df_stream)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(stream_temp = GPP_df_stream$stream_temp, middle_GPP = fitted_values)

# Plot
G_plot2 <- ggplot(data = GPP_df_stream, aes(x = stream_temp, y = log(middle_GPP+1))) +
  geom_point(size = 1.5, alpha = 0.6, aes(shape=position, color = shore)) +
  geom_errorbar(aes(ymin = log(lower_GPP + 1), ymax = log(upper_GPP + 1), color = shore), width = 0, alpha = 0.05) +
  geom_smooth(data = plot_data, aes(y = (middle_GPP)), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression(log(GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(stream~temp.~C)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 20.0, by = 4.0),
                     limits = c(0, 20.0)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.50, by = 0.9),
                     limits = c(0, 4.50)) +  
  
  theme_bw() 
G_plot2




# Trendline
fitted_values <- predict(er_stream_lmer_mod4, newdata = ER_df_stream)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(flow_sum = ER_df_stream$flow_sum, ER = fitted_values)

# Plot
E_plot <- ggplot(data = ER_df_stream, aes(x = flow_sum, y = log(ER+1))) +
  geom_point(size = 1.5, alpha = 0.6, aes(shape=position, color = shore)) +
  geom_errorbar(aes(ymin = log(ER_low + 1), ymax = log(ER_up + 1), color = shore), width = 0, alpha = 0.05) +  # Color error bars by shore and make them transparent
  geom_smooth(data = plot_data, aes(y = ER), method = "lm", se = F,  color = "grey50") + 
  #geom_smooth(method = "lm", se = F, aes(color = shore)) +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression(log(ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(cumulative~flow~m^3)) + 
  theme_bw() 
E_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_flowplot_a.png",sep=""), width=4,height=3,dpi=300)


# Trendline
fitted_values <- predict(er_stream_lmer_mod4, newdata = GPP_df_stream)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(stream_temp = ER_df_stream$stream_temp, ER = fitted_values)

# Plot
E_plot2 <- ggplot(data = ER_df_stream, aes(x = stream_temp, y = log(ER+1))) +
  geom_point(size = 1.5, alpha = 0.6, aes(shape=position, color = shore)) +
  geom_errorbar(aes(ymin = log(ER_low + 1), ymax = log(ER_up + 1), color = shore), width = 0, alpha = 0.05) +  # Color error bars by shore and make them transparent
  geom_smooth(data = plot_data, aes(y = ER), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression(log(ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(stream~temp.~C)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 20.0, by = 4.0),
                     limits = c(0, 20.0)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.50, by = 0.9),
                     limits = c(0, 4.50)) +  
  
  theme_bw() 
E_plot2




stream_metab_grid <- ggarrange(G_plot,
                               G_plot2,
                       E_plot,
                       E_plot2,
                       ncol = 2, nrow = 2,
                       common.legend = TRUE, 
                       labels = c("A", "B", "C", "D"),
                       legend = "bottom")
###
# ggsave(plot = stream_metab_grid, filename = paste("./SFS24_Analysis/figures/stream_metab_grid_plot.png",sep=""),width=6.2,height=6.35,dpi=300)













# Filter out missing values from the original data frame
ER_df_stream_complete <- na.omit(ER_df_stream_week)
# Extract fitted values for the complete observations
fitted_values <- predict(ER_stream_mod4a, newdata = ER_df_stream_complete)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(flow_sum = ER_df_stream_complete$flow_sum, middle_ER = fitted_values)

# Plot
E_plot <- ggplot(data = ER_df_stream_week, aes(x = flow_sum, y = log(middle_ER+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
E_plot

# ggsave(plot = E_plot, filename = paste("./SFS24_Analysis/figures/ER_flowplot_a.png",sep=""), width=4,height=3,dpi=300)


# Plot
E_plot <- ggplot(data = ER_df_stream_week, aes(x = stream_temp, y = log(middle_ER+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
E_plot

# ggsave(plot = E_plot, filename = paste("./SFS24_Analysis/figures/ER_streamtemp_a.png",sep=""), width=4,height=3,dpi=300)

# Plot
E_plot <- ggplot(data = ER_df_stream_week, aes(x = precip_sum, y = log(middle_ER+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
E_plot

# ggsave(plot = E_plot, filename = paste("./SFS24_Analysis/figures/ER_precipevent_a.png",sep=""), width=4,height=3,dpi=300)

# Plot
G_plot <- ggplot(data = ER_df_stream_week, aes(x = windsp_mean, y = log(middle_ER+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/ER_windsp_a.png",sep=""), width=4,height=3,dpi=300)

##################
####################

gpp_auto <- lmer(log(middle_ER+1) ~  log(middle_ER_lag+1) + scale(windsp_mean) + (1|shore / site) +(1|year), 
                 data = ER_df_week)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)


G_plot_windw <- ggplot(data = ER_df_week, aes(x= (windsp_mean), y=log(middle_ER+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + #geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~wind~sp~m~s^-1)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_windw



gpp_auto <- lmer(log(middle_ER+1) ~  log(middle_ER_lag+1) + scale(tmean_C) + (1|shore / site) +(1|year), 
                 data = ER_df_week)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot_temp <- ggplot(data = ER_df_week, aes(x= (tmean_C), y=log(middle_ER+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + #geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_temp


gpp_auto <- lmer(log(middle_ER+1) ~  log(middle_ER_lag+1) + scale(lake_wtemp) + (1|shore / site) +(1|year), 
                 data = ER_df_week)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot_ltempw <- ggplot(data = ER_df_week, aes(x= (lake_wtemp), y=log(middle_ER+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~lake~temp~C)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_ltempw


gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(precip_bi) + (1|shore / site) +(1|year), 
                 data = GPP_df1)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)


gpp_auto <- lmer(log(middle_ER+1) ~  log(middle_ER_lag+1) + scale(precip_sum) + (1|shore / site) +(1|year), 
                 data = ER_df_week)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)



G_plot_precipw <- ggplot(data = ER_df_week, aes(x=as.factor(precip_sum), y=log(middle_ER+1), fill=shore, color = shore))+
  # geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  # Add box plot for reference
  labs(y=expression(log(mean~weekly~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Sum~weekly~precip~events)) + 
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_precipw

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_temp_w.png",sep=""),width=6,height=8,dpi=300)



gpp_autow <- lmer(log(middle_ER+1) ~  log(middle_ER_lag+1) + scale(lake_par_int) + (1|shore / site) +(1|year), data = ER_df_week%>%filter(shore=="BW"|shore=="GB"))
summary(gpp_autow)
hist(residuals(gpp_autow))
r.squaredGLMM(gpp_autow)

G_plot_lightw <- ggplot(data = ER_df_week%>%filter(shore=="BW"|shore=="GB"|shore=="SH"), aes(x= (lake_par_int), y=log(middle_ER+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + #geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~PAR~umol~m^-2~s^-1)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_lightw




G_plot_precipw
G_plot_ltempw
G_plot_windw
G_plot_lightw



metab23 <- ggarrange(G_plot_ltempw,
                     G_plot_windw,
                     G_plot_lightw,
                     G_plot_precipw,
                     ncol = 4, nrow = 1,
                     common.legend = TRUE, 
                     legend = "bottom")
###
ggsave(plot = metab23, filename = paste("./SFS24_Analysis/figures/ER1_weekly_plot.png",sep=""),width=12,height=4,dpi=300)



G_plot <- ggplot(data = GPP_df1, aes(x= (lake_par_int), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_par_w.png",sep=""),width=6,height=8,dpi=300)

##########
############
##########


# Chemistry models? 

chem_dat <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/NS_chem_dat_nh4_24.rds") %>% # NS_chem_dat_nh4_24.rds NS_chem_dat_24.rds
  filter(site!="BW10m" & site!="BW15m" & site!="BW20m" & 
           site!="GB10m" & site!="GB15m" &  site!="GB20m")

unique(chem_dat$site)
summary(chem_dat) 

### Stats ###
chem_dat_stream_sw <- chem_dat %>%
  filter(location=="stream" & substrate=="sw") %>%
  group_by(shore, location, date, substrate) %>%
  summarise(
    NO3_mgL_dl = mean(NO3_mgL_dl, na.rm=T),
    NH3_mgL_dl= mean(NH3_mgL_dl, na.rm=T),
    NH4_mgL_dl= mean(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl = mean(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl = mean(DOC_mgL_dl, na.rm=T),
    pH_infill= mean(pH_infill, na.rm=T))

summary(chem_dat_stream_sw)
chem_dat_stream_sw_gb <- chem_dat_stream_sw%>%
  filter(shore=="GB")
summary(chem_dat_stream_sw_gb)

chem_dat_stream_sw_BW <- chem_dat_stream_sw%>%
  filter(shore=="BW")
summary(chem_dat_stream_sw_BW)

chem_dat_stream_tab <- chem_dat %>%
  group_by(shore, location, substrate) %>%
  summarise(
    NO3_mgL_dl_m = mean(NO3_mgL_dl, na.rm=T),
    NO3_mgL_dl_min = min(NO3_mgL_dl, na.rm=T),
    NO3_mgL_dl_max = max(NO3_mgL_dl, na.rm=T),
    NH4_mgL_dl_m= mean(NH4_mgL_dl, na.rm=T),
    NH4_mgL_dl_min= min(NH4_mgL_dl, na.rm=T),
    NH4_mgL_dl_max= max(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl_m = mean(PO4_ugL_dl, na.rm=T),
    PO4_ugL_dl_min = min(PO4_ugL_dl, na.rm=T),
    PO4_ugL_dl_max = max(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl_m = mean(DOC_mgL_dl, na.rm=T),
    DOC_mgL_dl_min = min(DOC_mgL_dl, na.rm=T),
    DOC_mgL_dl_max = max(DOC_mgL_dl, na.rm=T),
    pH_infill_m= mean(pH_infill, na.rm=T),
    pH_infill_min= min(pH_infill, na.rm=T),
    pH_infill_max= max(pH_infill, na.rm=T))

# write.csv(chem_dat_stream_tab, file = "./Chem_table1.csv", row.names = TRUE)
###

chem_datQ <- chem_dat %>%
  group_by(shore, location, date, substrate) %>%
  summarise(
    NO3_mgL_dl = mean(NO3_mgL_dl, na.rm=T),
    NH3_mgL_dl= mean(NH3_mgL_dl, na.rm=T),
    NH4_mgL_dl= mean(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl = mean(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl = mean(DOC_mgL_dl, na.rm=T))


chem_dat_stats <- chem_dat %>%
  filter(location=="lake")%>%
  filter(date> as.Date("2023-01-01"))%>%
  group_by(shore, location) %>%
  summarise(
    NO3_mgL = mean(NO3_mgL_dl, na.rm=T),
    NH3_mgL= mean(NH3_mgL_dl, na.rm=T),
    NH4_mgL= mean(NH4_mgL_dl, na.rm=T),
    PO4 = mean(PO4_ugL_dl, na.rm=T),
    DOC_mgL = mean(DOC_mgL_dl, na.rm=T))

summary(chem_datQ)


(0.0363-0.0184)/0.0363
(0.0441-0.0207)/0.0441

(10.3-6.87)/10.3

(20.7-8.34)/20.7

(2.74-0.599)/2.74

(4.05-0.901)/4.05
# NH3_plot <- ggplot(data = chem_datQ%>%filter(location=="lake"), aes(x =date, y = NH3_mgL_dl, shape = substrate, color=shore)) +
#   geom_point()+
#   scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
#   theme_bw() +  geom_hline(yintercept = 0.002, color = "grey") +
#   labs(y=expression(NH[3]~mg~L^-1), x=NULL) #+ facet_grid(substrate~.)
# 

NH4_plot <- ggplot(data = chem_datQ%>%filter(location=="lake"), aes(x =date, y = NH4_mgL_dl, shape = substrate, color=shore)) +
  geom_point()+
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  geom_hline(yintercept = 0.002, color = "grey") +
  labs(y=expression(NH[4]~mg~L^+1), x=NULL) #+ facet_grid(substrate~.)


NO3_plot <- ggplot(data = chem_datQ%>%filter(location=="lake"), aes(x =date, y = NO3_mgL_dl, shape = substrate, color=shore)) +
  geom_point()+
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  geom_hline(yintercept = 0.002, color = "grey") +
  labs(y=expression(NO[3]~mg~L^-1), x=NULL) + facet_grid(substrate~.)


PO4_plot <- ggplot(data = chem_datQ%>%filter(location=="lake"), aes(x =date, y = PO4_ugL_dl, shape = substrate, color=shore)) +
  geom_point()+
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  geom_hline(yintercept = 0.002, color = "grey") +
  labs(y=expression(o-phos~ug~L^-1), x=NULL) + facet_grid(substrate~.)

######################################################
chem_dat_wide <- chem_datQ %>%
  filter(shore=="BW"|shore=="GB")%>%
  pivot_wider(
    id_cols = c("shore", "date", "substrate"),
    names_from = c("location"),
    values_from = c("NO3_mgL_dl", "NH4_mgL_dl", "PO4_ugL_dl", "DOC_mgL_dl")
  )

summary(chem_dat_wide)
names(chem_dat_wide)
stream_NO3 <- lmer(NO3_mgL_dl_stream~NO3_mgL_dl_lake + (1|shore), data=chem_dat_wide)
summary(stream_NO3)
r.squaredGLMM(stream_NO3)
hist(residuals(stream_NO3))


chem_dat_wide_datecheck <- chem_dat_wide%>%
  filter(site=="BWNS1" |site=="BWNS2" |site=="BWNS3"|
           site=="GBNS1" |site=="GBNS2" |site=="GBNS3")

str(chem_dat_wide)

dates <- unique(chem_dat_wide_datecheck$date)

### GLM for stream NO3 (either pw or sw) and lake NO3

str(as.data.frame(chem_datQ))

chem_datQ$week <- week(chem_datQ$date)
chem_datQ$year <- year(chem_datQ$date)
chem_datQ$weekyear <- paste(chem_datQ$year, chem_datQ$week, sep = "-")

chem_dat <-as.data.frame(chem_dat)


chem_dat_wide <- chem_datQ %>%
  filter(shore=="BW"|shore=="GB")%>%
  pivot_wider(
    id_cols = c("shore", "date", "substrate"),
    names_from = c("location"),
    values_from = c("NO3_mgL_dl", "NH3_mgL_dl", "PO4_ugL_dl", "DOC_mgL_dl")
  )

chem_dat_filtered <- chem_dat_wide %>%
  mutate(week_period = ceiling(week(date) / 3))

chem_dat_filtered <- chem_dat_filtered %>%
  mutate(grouping_var = paste0(year(date), "-", week_period))

chem_dat_wide_summary <- chem_dat_filtered %>%
  group_by(shore, substrate, grouping_var) %>%
  summarise(
    NO3_mgL_dl_stream = mean(NO3_mgL_dl_stream, na.rm = TRUE),
    NO3_mgL_dl_lake = mean(NO3_mgL_dl_lake, na.rm = TRUE),
    NO3_mgL_dl_interface = mean(NO3_mgL_dl_interface, na.rm = TRUE),
    NH4_mgL_dl_stream = mean(NH4_mgL_dl_stream, na.rm = TRUE),
    NH4_mgL_dl_lake = mean(NH4_mgL_dl_lake, na.rm = TRUE),
    NH4_mgL_dl_interface = mean(NH4_mgL_dl_interface, na.rm = TRUE),
    PO4_ugL_dl_stream = mean(PO4_ugL_dl_stream, na.rm = TRUE),
    PO4_ugL_dl_lake = mean(PO4_ugL_dl_lake, na.rm = TRUE),
    PO4_ugL_dl_interface = mean(PO4_ugL_dl_interface, na.rm = TRUE),
    DOC_mgL_dl_stream = mean(DOC_mgL_dl_stream, na.rm = TRUE),
    DOC_mgL_dl_lake = mean(DOC_mgL_dl_lake, na.rm = TRUE),
    DOC_mgL_dl_interface = mean(DOC_mgL_dl_interface, na.rm = TRUE),
    # Add other variables and summary statistics as needed
    .groups = "drop"
  )

NO3_plot <- ggplot(data = chem_dat_wide_summary%>%filter(NO3_mgL_dl_stream<0.15), aes(x =NO3_mgL_dl_lake,  y = NO3_mgL_dl_stream, color=shore, shape=substrate)) +
  geom_point(size=2, alpha=c(0.75)) + geom_smooth(method="lm", color="grey25", se=F, lty=1) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) + theme_bw() #+ facet_grid(substrate~.)

hist(chem_dat_wide_summary$NO3_mgL_dl_stream)
hist(log10(chem_dat_wide_summary$NO3_mgL_dl_stream +1))
stream_NO3 <- lmer(NO3_mgL_dl_stream~NO3_mgL_dl_lake + (1|shore), data=chem_dat_wide_summary)
summary(stream_NO3)

r.squaredGLMM(stream_NO3)
hist(residuals(stream_NO3))


######
#####

# Fit the linear mixed-effects model
stream_NO3 <- lmer(NO3_mgL_dl_lake ~ NO3_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary)
summary(stream_NO3)
r.squaredGLMM(stream_NO3)
hist(residuals(stream_NO3))
# Extract the fixed effects coefficients
coefficients <- fixef(stream_NO3)
# Extract the slope coefficient
slope <- coefficients["NO3_mgL_dl_stream"]
# Calculate the reciprocal of the slope
inverse_slope <- 1 / slope
# Create a data frame for plotting the trend line
trend_data <- data.frame(NO3_mgL_dl_stream = seq(min(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)),
                                                 max(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)), length.out = 100),
                         NO3_mgL_dl_lake = inverse_slope * seq(min(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)),
                                                               max(na.omit(chem_dat_wide_summary$NO3_mgL_dl_stream)), length.out = 100))

# Plot
NO3_plot <- ggplot(data = chem_dat_wide_summary, 
                   aes(y = NO3_mgL_dl_lake, x = NO3_mgL_dl_stream)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate, color=shore)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  xlim(0,0.2) +ylim(0,0.2)+
  geom_line(data = trend_data, aes(y = NO3_mgL_dl_lake, x = NO3_mgL_dl_stream), 
            color = "grey25", linetype = 1) +  # Add trend line
  geom_hline(yintercept = 0.003, color = "grey50") + 
  geom_vline(xintercept = 0.003, color = "grey50") +
  labs(y=expression(Lake~NO[3]^-1~mg~L^-1), x=expression(Stream~NO[3]^-1~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()


######
######
######


# Fit the linear mixed-effects model
stream_NH3 <- lmer(NH4_mgL_dl_lake ~ NH4_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary)
summary(stream_NH3)
r.squaredGLMM(stream_NH3)
hist(residuals(stream_NH3))
# Extract the fixed effects coefficients
coefficients <- fixef(stream_NH3)
# Extract the slope coefficient
slope <- coefficients["NH4_mgL_dl_stream"]
# Calculate the reciprocal of the slope
inverse_slope <- 1 / slope
# Create a data frame for plotting the trend line
trend_data <- data.frame(NH4_mgL_dl_stream = seq(min(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)),
                                                 max(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)), length.out = 100),
                         NH4_mgL_dl_lake = inverse_slope * seq(min(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)),
                                                               max(na.omit(chem_dat_wide_summary$NH4_mgL_dl_stream)), length.out = 100))


# Plot
NH3_plot <- ggplot(data = chem_dat_wide_summary, 
                   aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) + 
  geom_line(data = trend_data, aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream), 
            color = "grey25", linetype = 1) +  # Add trend line
  geom_hline(yintercept = 0.002, color = "grey50") + 
  geom_vline(xintercept = 0.002, color = "grey50") +
  labs(y = expression(Lake~NH[4]^+1~mg~L^-1), x = expression(Stream~NH[4]^+1~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  xlim(0,0.1) + ylim(0,0.1) +
  theme_bw()

NH3_plot



# # Plot
# NH3_plot <- ggplot(data = chem_dat_wide_summary, 
#                    aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream, color = shore)) +
#   geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) + xlim(0,0.1) +ylim(0,0.1)+
#   geom_line(data = trend_data, aes(y = NH4_mgL_dl_lake, x = NH4_mgL_dl_stream), 
#             color = "grey25", linetype = 1) +  # Add trend line
#   geom_hline(yintercept = 0.002, color = "grey50") + 
#   geom_vline(xintercept = 0.002, color = "grey50") +
#   labs(y=expression(Lake~NH[4]^+1~mg~L^-1), x=expression(Stream~NH[4]^+1~mg~L^-1)) +
#   scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
#   theme_bw()
# 
# NH3_plot


#####
#####
####


# Fit the linear mixed-effects model
stream_PO4 <- lmer(PO4_ugL_dl_lake ~ PO4_ugL_dl_stream + (1 | shore), data = chem_dat_wide_summary)
summary(stream_PO4)
# Create a data frame for plotting the trend line
# Plot
PO4_plot <- ggplot(data = chem_dat_wide_summary,
                   aes(x = PO4_ugL_dl_stream, y = PO4_ugL_dl_lake, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) +
  # geom_line(data = trend_data, aes(x = NO3_mgL_dl_lake, y = NO3_mgL_dl_stream), 
  #           color = "grey25", linetype = 1) +  # Add trend line
  labs(y=expression(Lake~ophos~ug~L^-1), x=expression(Stream~ophos~ug~L^-1)) +
  geom_hline(yintercept = 0.402, color = "grey50") + 
  geom_vline(xintercept = 0.402, color = "grey50") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 60, by = 12),
                     limits = c(0, 60)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 60, by = 12),
                     limits = c(0, 60)) +  
  theme_bw()
  



library(tidyverse)

# Remove rows with missing values
chem_dat_wide_summary_clean <- chem_dat_wide_summary %>%
  drop_na(DOC_mgL_dl_lake, DOC_mgL_dl_stream)

# Plot
DOC_plot <- ggplot(data = chem_dat_wide_summary_clean,
                   aes(y = DOC_mgL_dl_lake, x = DOC_mgL_dl_stream, color = shore)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) +
  geom_hline(yintercept = 0.25, color = "grey50") + 
  geom_vline(xintercept = 0.25, color = "grey50") +
  labs(y = expression(Lake~DOC~mg~L^-1), x = expression(Stream~DOC~mg~L^-1)) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.5, by = 0.9),
                     limits = c(0, 4.5)) +  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), 
                     breaks = seq(0, 4.5, by = 0.9),
                     limits = c(0, 4.5)) +  
  theme_bw()

DOC_plot



# Fit the linear mixed-effects model
stream_DOC <- lmer(DOC_mgL_dl_lake ~ DOC_mgL_dl_stream + (1 | shore), data = chem_dat_wide_summary_clean)
summary(stream_DOC)
# # Create a data frame for plotting the trend line
# # Plot
# DOC_plot <- ggplot(data = chem_dat_wide_summary,
#                    aes(y = DOC_mgL_dl_lake, x = DOC_mgL_dl_stream, color = shore)) +
#   geom_point(size = 2, alpha = 0.75, aes(shape = substrate)) +
#   geom_hline(yintercept = 0.25, color = "grey50") + 
#   geom_vline(xintercept = 0.25, color = "grey50") +
#   labs(y=expression(Lake~DOC~mg~L^-1), x=expression(Stream~DOC~mg~L^-1)) +
#   scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
#   scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  
#   xlim(0, 4.5) +ylim(0, 4.5)+
#   theme_bw()


r.squaredGLMM(stream_DOC)
hist(residuals(stream_DOC))




chem_grid <- ggarrange(NO3_plot,
                     NH3_plot,
                     DOC_plot,
                     PO4_plot,
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, 
                     labels = c("A", "B", "C", "D"),
                     legend = "bottom")
###
# ggsave(plot = chem_grid, filename = paste("./SFS24_Analysis/figures/Lake_streamChem_plot.png",sep=""),width=5,height=5.05,dpi=300)

#############



##===========================================
## SEM for the SFS presentation ##
names(GPP_df_stream)
GPP_df_stream 



lake_temp_mod <- lmer(lake_tempC ~ 
                        (ppmt_sum) +
                        #(tmean_C) +
                        #stream_temp +
                        (flow_sum) + 
                        (1|shore/site), data = GPP_df_stream)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))
vif(lake_temp_mod)

## 
lake_light_mod <- lmer(par_int_3m ~ 
                         (light_mean) + 
                         (ppmt_sum)+
                         #(windsp_mean)+
                         #(flow_sum)+
                         #scale(real_NS_depth)+
                         (1|shore/site), data = GPP_df_stream)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)

# this needs work 
###

GPP_mod_null <- lmer(middle_GPP ~ (middle_GPP_lag) + (1|shore/site), data = GPP_df_stream)
summary(GPP_mod_null)
hist(residuals(GPP_mod_null))

GPP_mod <- lmer(middle_GPP ~ (middle_GPP_lag) +
                  (lake_tempC) +
                  (par_int_3m) +
                  (1|shore/site), data = GPP_df_stream)
summary(GPP_mod)
hist(residuals(GPP_mod))
vif(GPP_mod)


GPP_mod1 <- lmer(middle_GPP ~ (middle_GPP_lag) +
                  (lake_tempC) +
                  (par_int_3m) +
                  #scale(windsp_mean) +
                  #scale(ppmt_sum) +
                  (1|shore/site), data = GPP_df_stream)
summary(GPP_mod1)
hist(residuals(GPP_mod1))
vif(GPP_mod1)


GPP_mod2 <- lmer(middle_GPP ~ (middle_GPP_lag) +
                   (lake_tempC) +
                   (par_int_3m) +
                   #(windsp_mean) +
                   (ppmt_sum) +
                   (flow_sum)+
                   (1|shore/site), data = GPP_df_stream)
summary(GPP_mod2)
hist(residuals(GPP_mod2))


GPP_mod3 <- lmer(middle_GPP ~ (middle_GPP_lag) +
                   (lake_tempC) +
                   (par_int_3m) +
                   #scale(windsp_mean) +
                   #scale(ppmt_sum) +
                   (flow_sum)+
                   (1|shore/site), data = GPP_df_stream)
summary(GPP_mod3)
hist(residuals(GPP_mod3))

GPP_mod3 <- lmer(middle_GPP ~ (middle_GPP_lag) +
                   (lake_tempC) +
                   (light_mean) +
                   (1|shore), data = GPP_df_stream)
summary(GPP_mod3)



GPP_mod4 <- lmer(middle_GPP ~ (middle_GPP_lag) +
                   (lake_tempC) +
                   (par_int_3m) +
                    windsp_mean +
                   (ppmt_sum) +
                   (flow_sum)+
                   (1|shore/site), data = GPP_df_stream)
summary(GPP_mod4)
hist(residuals(GPP_mod4))


AIC(GPP_mod_null, GPP_mod1, GPP_mod2, GPP_mod3)
library(piecewiseSEM)

# Use the `psem` function to create the SEM
GPP_psem <- psem(
  GPP_mod4,
  lake_temp_mod)
#, lake_light_mod)

# Look at object
GPP_psem # Sometimes it works better to specify models in line, with called df
summary(GPP_psem)
plot(GPP_psem)
fisherC(GPP_psem) # says model

#######
## maybe lets pull out one site at time. 

##===========================================
## SEM for the SFS presentation ##
names(GPP_df_stream)
GPP_df_stream 

GPP_df_BW <- GPP_df_stream %>%
  filter(shore=="BW")

## Issues with normality... 
hist(GPP_df_BW$lake_tempC)
hist(GPP_df_BW$middle_GPP)
hist(GPP_df_BW$par_int_3m)

GPP_df_BW$logGPP <- log(GPP_df_BW$middle_GPP+1)
hist(GPP_df_BW$logGPP)

GPP_df_BW$logGPPlag <- log(GPP_df_BW$middle_GPP_lag+1)
GPP_df_BW$loglake_tempC <- log(GPP_df_BW$lake_tempC)
hist(GPP_df_BW$loglake_tempC)


hist(GPP_df_BW$par_int_3m)
summary(GPP_df_BW)

GPP_df_BWQ <- GPP_df_BW%>%
  select(site, date, year, middle_GPP, middle_GPP_lag, logGPP, logGPPlag, lake_tempC, loglake_tempC, Kd_fill, par_int_3m, tmean_C, light_mean, flow_mean, log_streamflow, log_ppmt, windsp_mean)

GPP_df_BWQ<- drop_na(GPP_df_BWQ)

# ## let's scale everything other than gpp:
# GPP_BW_sem <- GPP_df_BW %>%
#   mutate(tmean_Cs = scale(tmean_C),
#          light_means = scale(light_mean),
#          Kd_fills = scale(Kd_fill),
#          par_int_3ms = scale(par_int_3m),
#          real_NS_depths = scale(real_NS_depth),
#          ppmt_sums = scale(ppmt_sum),
#          lake_tempCs = scale(lake_tempC), 
#          flow_mean_ms = scale(flow_mean_m),
#          flow_sums =scale(flow_sum),
#          stream_temps= scale(stream_temp)) %>%
#   dplyr::select(site, week, year, logGPP, logGPPlag, tmean_Cs, light_means, Kd_fills, par_int_3ms,
#                 real_NS_depths, ppmt_sums, lake_tempCs, flow_mean_ms, flow_sums, stream_temps)
#          
#          
#          
#          
#          
#          
#          
#          
#     
# 
# 
# lake_temp_mod <- lmer(loglake_tempC ~ 
#                         (ppmt_sum) +
#                         #(tmean_C) +
#                         #stream_temp +
#                         (flow_sum) + 
#                         (1|site), data = GPP_df_BW)
# summary(lake_temp_mod)
# hist(residuals(lake_temp_mod))
# vif(lake_temp_mod)
# 
# ## 
# lake_light_mod <- lmer(par_int_3m ~ 
#                          (light_mean) + 
#                          (ppmt_sum)+
#                          #(windsp_mean)+
#                          #(flow_sum)+
#                          #scale(real_NS_depth)+
#                          (1|site), data = GPP_df_BW)
# summary(lake_light_mod)
# hist(residuals(lake_light_mod))
# vif(lake_light_mod)
# 
# # this needs work 
# ###
# 
# GPP_mod_null <- lmer(log(middle_GPP+1) ~ 
#                        #log(middle_GPP_lag+1) + 
#                        (1|site), data = GPP_df_BW)
# summary(GPP_mod_null)
# hist(residuals(GPP_mod_null))
# 
# GPP_mod <- lmer(log(middle_GPP+1) ~ #(middle_GPP_lag) +
#                   (lake_tempC) +
#                   (par_int_3m) +
#                   (1|site), data = GPP_df_BW)
# summary(GPP_mod)
# hist(residuals(GPP_mod))
# vif(GPP_mod)
# 
# 
# GPP_mod1 <- lmer(log(middle_GPP+1) ~ 
#                    #(middle_GPP_lag) +
#                    (lake_tempC) +
#                    (par_int_3m) +
#                    scale(windsp_mean) +
#                    scale(ppmt_sum) +
#                    (1|site), data = GPP_df_BW)
# summary(GPP_mod1)
# hist(residuals(GPP_mod1))
# vif(GPP_mod1)
# 
# 
# GPP_mod2 <- lmer(log(middle_GPP+1) ~ 
#                    #(middle_GPP_lag) +
#                    (lake_tempC) +
#                    (par_int_3m) +
#                    #(windsp_mean) +
#                    (ppmt_sum) +
#                    (flow_sum)+
#                    (1|site), data = GPP_df_BW)
# summary(GPP_mod2)
# hist(residuals(GPP_mod2))
# 
# 
# GPP_mod3 <- lmer(log(middle_GPP+1) ~ 
#                    #(middle_GPP_lag) +
#                    (lake_tempC) +
#                    (par_int_3m) +
#                    #scale(windsp_mean) +
#                    #scale(ppmt_sum) +
#                    (flow_sum)+
#                    (1|site), data = GPP_df_BW)
# summary(GPP_mod3)
# hist(residuals(GPP_mod3))
# 
# GPP_mod3 <- lmer(middle_GPP ~ (middle_GPP_lag) +
#                    (lake_tempC) +
#                    (light_mean) +
#                    (1|shore), data = GPP_df_stream)
# summary(GPP_mod3)
# 
# 
# 
# GPP_mod4 <- lmer(middle_GPP ~ (middle_GPP_lag) +
#                    (lake_tempC) +
#                    (par_int_3m) +
#                    windsp_mean +
#                    (ppmt_sum) +
#                    (flow_sum)+
#                    (1|site), data = GPP_df_BW)
# summary(GPP_mod4)
# hist(residuals(GPP_mod4))
# 
# 
# AIC(GPP_mod_null, GPP_mod1, GPP_mod2, GPP_mod3, GPP_mod4)
# library(piecewiseSEM)
# 
# # Use the `psem` function to create the SEM
# GPP_psem <- psem(
#   GPP_mod4,
#   lake_temp_mod)
# #, lake_light_mod)
# 
# # Look at object
# GPP_psem # Sometimes it works better to specify models in line, with called df
# summary(GPP_psem)
# plot(GPP_psem)
# fisherC(GPP_psem)
# 
# 
# ######################
# #######################
# SFS_datQ2

SFS_datQ2_BW <- GPP_df_stream%>%
  filter(shore=="BW")


lake_temp_mod <- glmmTMB(loglake_tempC ~ 
                        #log_ppmt +
                        log_streamflow + 
                        (1|site), data = GPP_df_BWQ)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))


## 
lake_light_mod <- glmmTMB(par_int_3m ~ 
                         light_mean + 
                         log_ppmt+
                        # windsp_mean +
                         (1|site), data = GPP_df_BWQ)
summary(lake_light_mod)
hist(residuals(lake_light_mod))


GPP_mod1 <- glmmTMB(logGPP ~ logGPPlag +
                   loglake_tempC +
                   par_int_3m +
                  # windsp_mean +
                   #log_ppmt +
                   #log_streamflow+
                   (1|site), data = GPP_df_BWQ)
summary(GPP_mod1)
hist(residuals(GPP_mod1))

library(semEff)
# Use the `psem` function to create the SEM
GPP_psem <- psem(
  GPP_mod1,
  #lake_light_mod,
  lake_temp_mod)
#, lake_light_mod)

# Look at object
GPP_psem # Sometimes it works better to specify models in line, with called df
summary(GPP_psem)
plot(GPP_psem)
fisherC(GPP_psem)


shipley.sem.eff <- semEff(GPP_psem)


## GB ###

GPP_df_GB <- GPP_df_stream%>%
  filter(shore=="GB")


## Issues with normality... 
hist(GPP_df_GB$lake_tempC)
hist(GPP_df_GB$middle_GPP)
hist(GPP_df_GB$par_int_3m)

GPP_df_GB$logGPP <- log(GPP_df_GB$middle_GPP+1)
hist(GPP_df_GB$logGPP)

GPP_df_GB$logGPPlag <- log(GPP_df_GB$middle_GPP_lag+1)
GPP_df_GB$loglake_tempC <- log(GPP_df_GB$lake_tempC)
hist(GPP_df_GB$loglake_tempC)


hist(GPP_df_GB$par_int_3m)
summary(GPP_df_GB)

GPP_df_GBQ <- GPP_df_GB%>%
  select(site, date, year, middle_GPP, middle_GPP_lag, logGPP, logGPPlag, lake_tempC, loglake_tempC, Kd_fill, par_int_3m, tmean_C, light_mean, flow_mean, log_streamflow, log_ppmt, windsp_mean)

GPP_df_GB<- drop_na(GPP_df_GBQ)

GPP_df_GB <- GPP_df_GB[-1,]

lake_temp_mod <- glmmTMB(loglake_tempC ~ 
                           log_ppmt +
                           log_streamflow + 
                           (1|site), data = GPP_df_GB)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))

## 
lake_light_mod <- glmmTMB(par_int_3m ~ 
                            light_mean + 
                            log_ppmt+
                            # windsp_mean +
                            (1|site), data = GPP_df_GB)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(GPP_mod1)


GPP_mod1 <- glmmTMB(logGPP ~ logGPPlag +
                     # loglake_tempC +
                      #par_int_3m +
                      # windsp_mean +
                      log_ppmt +
                      log_streamflow+
                      (1|site), data = GPP_df_GB)
summary(GPP_mod1)
hist(residuals(GPP_mod1))

# Use the `psem` function to create the SEM
GPP_psem <- psem(
  GPP_mod1,
  #lake_light_mod,
  lake_temp_mod)
#, lake_light_mod)

# Look at object
GPP_psem # Sometimes it works better to specify models in line, with called df
summary(GPP_psem)
plot(GPP_psem)
fisherC(GPP_psem)