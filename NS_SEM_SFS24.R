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
library(glmmTMB)


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

############################################
####### create some new variables ##########
###########################################

GPP_df1$yearf<- as.factor(year(GPP_df1$date))
GPP_df1$year<- (year(GPP_df1$date))

GPP_df1$GPP <- log(GPP_df1$middle_GPP+1)
GPP_df1$GPP_lag <- log(GPP_df1$lag_middle_GPP+1)

hist(GPP_df1$GPP)
names(GPP_df1)

unique(GPP_df1$shore)
unique(GPP_df1$site)

hist(GPP_df1$flow_mean)
hist(GPP_df1$flow_mean)
hist(log(GPP_df1$flow_mean)+1)
hist(GPP_df1$ppt_mm)
summary((GPP_df1$ppt_mm))

GPP_df1$precip_bi <- ifelse(GPP_df1$ppt_mm > 0, 1, 0)

hist(log(GPP_df1$ppt_mm+1))

GPP_df1$log_ppmt <- log(GPP_df1$ppt_mm+1)
GPP_df1$log_streamflow <- (log(GPP_df1$flow_mean)+1)
hist(GPP_df1$log_streamflow)
names(GPP_df1)

# 
# GPP_df_scale <- GPP_df1 %>%
#   scale(7:26, 28:43, 46)
# 
# GPP_df1[, c(7:26, 28:43, 46)] <- scale(GPP_df1[, c(7:26, 28:43, 46)])
# summary(GPP_df1)

hist((GPP_df1$lake_par))

# GPP_df1<- GPP_df1%>%
#   drop_na(site)
# 
# GPP_df1<- GPP_df1%>%
#   drop_na(GPP)

########################################################
####### DAILY ##################
######################################

library(PerformanceAnalytics)
names(GPP_df1)
Streamdf_cor <- GPP_df1[, c(17:26)]
Streamdf_cor <- GPP_df1[, c("Stream_temp", "ppt_mm", "flow_mean", "stream_SPC")]
Streamdf_cor <- GPP_df1[, c("Stream_temp", "stream_flow", "stream_SPC")]
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)


lakedf_cor <- GPP_df_stream_week[, c("lake_wtemp", "lake_par_int", "light_mean", "tmean_C", "windsp_mean", "flow_mean", "flow_sum","stream_temp", "spc_mean", "precip_sum")]
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

lake_temp_mod4 <- lmer(lake_wtemp ~ 
                         windsp_mean +
                         precip_bi +
                         log_streamflow +
                         tmean_C  + 
                         (1|site), data = GPP_df1)
summary(lake_temp_mod4)
hist(residuals(lake_temp_mod4))
vif(lake_temp_mod4)
r.squaredGLMM(lake_temp_mod4)

AIC(lake_temp_mod, lake_temp_mod2, lake_temp_mod3, lake_temp_mod4)


## Look at best fit lake temp model:
lake_temp_mod3
hist(GPP_df1$lake_temp)
G_plot <- ggplot(data = GPP_df1, aes(x= (windsp_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_windsp.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()  + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_windsp.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x= (tmean_C), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_airtemp.png",sep=""),width=6,height=8,dpi=300)

G_plot <- ggplot(data = GPP_df1, aes(x= (spc_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_SPC.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x = factor(precip_bi), y = lake_wtemp, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_precip.png",sep=""),width=4,height=8,dpi=300)



## 
lakelight_cor <- GPP_df1[c("lake_par_int", "light_cv", "flow_mean", "precip_bi", "lake_wspeed", "lake_depth") ]
chart.Correlation(lakelight_cor, histogram=TRUE, pch=19)
hist(GPP_df1$lake_par_int)

lake_light_mod <- lmer(lake_par_int ~ 
                         flow_mean +
                         precip_bi +
                         windsp_mean + (1|site), data = GPP_df1)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)
r.squaredGLMM(lake_light_mod)




## Look at best fit lake temp model:
lake_temp_mod3
hist(GPP_df1$lake_temp)
G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_precip.png",sep=""),width=4,height=8,dpi=300)




G_plot <- ggplot(data = GPP_df1, aes(x= (windsp_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()  + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_windsp.png",sep=""),width=6,height=8,dpi=300)



G_plot <- ggplot(data = GPP_df1, aes(x= (spc_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot  + facet_grid(shore~.)

G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_flow.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x = factor(precip_bi), y = lake_par_int, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_precip.png",sep=""),width=4,height=8,dpi=300)






# DAG chapter in rethinking -
GPP_df1$yearf<- as.factor(year(GPP_df1$date))

gpp_auto1 <- lmer(GPP ~  GPP_lag + 
                    (lake_wtemp) + 
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto1)
residuals <- residuals(gpp_auto1)
hist(residuals)

gpp_auto2 <- lmer(GPP ~  GPP_lag + 
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



G_plot <- ggplot(data = GPP_df1, aes(x= (GPP_lag), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_GPPlag.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_flow.png",sep=""),width=6,height=8,dpi=300)



G_plot <- ggplot(data = GPP_df1, aes(x= (windsp_mean), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_windsp.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x = factor(precip_bi), y = GPP, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_precip.png",sep=""),width=4,height=8,dpi=300)


gpp_auto3 <- lmer(GPP ~  GPP_lag + 
                    (lake_wtemp) + 
                    precip_bi +
                    lake_wspeed +
                    #flow_mean +
                    (lake_par_int) + (1|site), data = GPP_df1)
summary(gpp_auto3)
residuals <- residuals(gpp_auto3)
vif(gpp_auto3)
hist(residuals(gpp_auto3))

gpp_auto4 <- lmer(GPP ~  GPP_lag + 
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
  gpp_auto2,
  lake_temp_mod3,  
  lake_light_mod)

# Look at object
model # Sometimes it works better to specify models in line, with called df
summary(model)
plot(model)
fisherC(model) # says model

model2 <- psem(
  gpp_auto1,
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


#################################################################################
#################
#### WEEKLY ####
################################################################################
GPP_df1$year<- year(GPP_df1$date)
GPP_df_week <- GPP_df1%>%
  group_by(site, shore, week, year, position) %>%
  summarise(middle_GPP=mean(middle_GPP, na.rm = TRUE),
            lower_GPP=mean(lower_GPP, na.rm = TRUE),
            upper_GPP=mean(upper_GPP, na.rm = TRUE),
            tmean_C=mean(tmean_C, na.rm = TRUE),
            light_mean=mean(light_mean, na.rm = TRUE),
            stream_temp=mean(stream_temp, na.rm = TRUE),
            spc_mean=mean(spc_mean, na.rm = TRUE),
            lake_do=mean(lake_do, na.rm = TRUE),
            lake_wtemp=mean(lake_wtemp, na.rm = TRUE),
            lake_par=mean(lake_par, na.rm = TRUE),
            lake_par_int=mean(lake_par_int, na.rm = TRUE),
            windsp_mean=mean(windsp_mean, na.rm = TRUE),
            precip_sum=sum(precip_bi, na.rm = TRUE),
            flow_sum=sum(flow_mean, na.rm = TRUE),
            flow_mean=mean(flow_mean, na.rm = TRUE))


GPP_df_week$precip_bi <- ifelse(GPP_df_week$precip_sum > 0, 1, 0)
GPP_df_week$middle_GPP_lag <- lag(GPP_df_week$middle_GPP)

## OVER LAP for weeks 22-31 

###########################################
# Remove the first row from the dataframe
###########################################
lake_temp_mod <- lmer(lake_wtemp ~ 
                        windsp_mean +
                        precip_sum +
                        flow_sum +
                        (1|shore), data = GPP_df_week)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))
vif(lake_temp_mod)

# Remove the first row from the dataframe
lake_temp_mod2 <- lmer(lake_wtemp ~ 
                        tmean_C +
                        precip_sum +
                        flow_sum +
                        (1|shore), data = GPP_df_week)
summary(lake_temp_mod2)
hist(residuals(lake_temp_mod2))
vif(lake_temp_mod2)


# Remove the first row from the dataframe
lake_temp_mod3 <- lmer(lake_wtemp ~ 
                         tmean_C +
                         precip_sum +
                         flow_sum +
                         (1|shore) + (1|year), data = GPP_df_week)
summary(lake_temp_mod3)
hist(residuals(lake_temp_mod3))
vif(lake_temp_mod3)

AIC(lake_temp_mod, lake_temp_mod2, lake_temp_mod3)

## Look at best fit lake temp model:
lake_temp_mod3
hist(GPP_df1$lake_temp)
G_plot <- ggplot(data = GPP_df_w, aes(x= (windsp_mean), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_WSP_week.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df_w, aes(x= log(flow_mean+1), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()  + facet_grid(shore~.)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_flow_w.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df_week, aes(x= flow_sum, y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_flow_w.png",sep=""),width=6,height=8,dpi=300)


ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_SPC.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df_w, aes(x = factor(precip_bi), y = lake_wtemp, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_point(alpha = 0.2) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_precip_w.png",sep=""),width=4,height=8,dpi=300)



## 
lakelight_cor <- GPP_df_w[c("lake_par_int", "light_mean", "flow_mean", "precip_bi", "precip_sum","windsp_mean") ]
chart.Correlation(lakelight_cor, histogram=TRUE, pch=19)
hist(GPP_df1$lake_par_int)

lake_light_mod <- lmer(lake_par_int ~ 
                         flow_mean +
                         precip_sum +
                         windsp_mean + (1|site), data = GPP_df_w)
summary(lake_light_mod)
hist(residuals(lake_light_mod))
vif(lake_light_mod)
r.squaredGLMM(lake_light_mod)

lake_light_mod2 <- lmer(lake_par_int ~ 
                         #flow_mean +
                         precip_sum +
                         windsp_mean + (1|site), data = GPP_df_w)
summary(lake_light_mod2)
hist(residuals(lake_light_mod2))
vif(lake_light_mod2)
r.squaredGLMM(lake_light_mod2)



## Look at best fit lake temp model:
lake_temp_mod3
hist(GPP_df1$lake_temp)
G_plot <- ggplot(data = GPP_df_w, aes(x= (precip_sum), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_laketemp_precip_w.png",sep=""),width=4,height=8,dpi=300)




G_plot <- ggplot(data = GPP_df_w, aes(x= (windsp_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()  + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_windsp_w.png",sep=""),width=6,height=8,dpi=300)



G_plot <- ggplot(data = GPP_df_w, aes(x= (spc_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
G_plot  + facet_grid(shore~.)

G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_flow.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x = factor(precip_bi), y = lake_par_int, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_light_precip.png",sep=""),width=4,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x= (GPP_lag), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_GPPlag.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x= (flow_mean), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_flow.png",sep=""),width=6,height=8,dpi=300)



G_plot <- ggplot(data = GPP_df1, aes(x= (windsp_mean), y=GPP, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_windsp.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df1, aes(x = factor(precip_bi), y = GPP, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.5) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_precip.png",sep=""),width=4,height=8,dpi=300)


gpp_auto1 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    scale(lake_wtemp) + 
                    scale(lake_par_int) + (1|site) + (1|year), data = GPP_df_week)
summary(gpp_auto1)
residuals <- residuals(gpp_auto1)
vif(gpp_auto1)
hist(residuals(gpp_auto1))

gpp_auto2 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    scale(lake_wtemp) + 
                    scale(precip_sum) +
                    scale(windsp_mean) +
                    scale(flow_sum) +
                    scale(lake_par_int) + (1|site)+(1|year), data = GPP_df_week)
summary(gpp_auto2)
residuals <- residuals(gpp_auto3)
vif(gpp_auto2)
hist(residuals(gpp_auto2))

gpp_auto3 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    scale(lake_wtemp) + 
                    scale(precip_sum) +
                    scale(windsp_mean) +
                    scale(flow_sum) +
                    scale(lake_par_int) + (1|site), data = GPP_df_week)
summary(gpp_auto3)
residuals <- residuals(gpp_auto3)
vif(gpp_auto3)
hist(residuals(gpp_auto3))

gpp_auto4 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    scale(lake_wtemp) + 
                    scale(precip_sum) +
                    scale(windsp_mean) +
                    scale(flow_mean) +
                    scale(lake_par_int) + (1|site), data = GPP_df_week)
summary(gpp_auto4)
residuals <- residuals(gpp_auto4)
vif(gpp_auto4)
hist(residuals(gpp_auto4))


gpp_auto5 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    scale(lake_wtemp) + 
                    scale(precip_sum) +
                    scale(windsp_mean) +
                    scale(flow_mean) +
                    scale(lake_par_int) + (1|site) +  (1|year), data = GPP_df_week)
summary(gpp_auto5)
residuals <- residuals(gpp_auto5)
vif(gpp_auto5)
hist(residuals(gpp_auto5))

AIC(gpp_auto1, gpp_auto2, gpp_auto3, gpp_auto4, gpp_auto5)


##############################


##############################



GPP_df_w1 <- GPP_df_w%>%
  drop_na(middle_GPP)

lake_light_mod2 <- lmer(lake_par_int ~ 
                          precip_sum +
                          windsp_mean + (1|shore), data = GPP_df_w1)
summary(lake_light_mod2)
hist(residuals(lake_light_mod2))
vif(lake_light_mod2)
r.squaredGLMM(lake_light_mod2)

G_plot <- ggplot(data = GPP_df_w1, aes(x= (windsp_mean), y=lake_par_int, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_par_windsp_w.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df_w1, aes(x = factor(precip_bi), y = lake_par_int, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_point(alpha = 0.4)+
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.2) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_par_precip_w.png",sep=""),width=4,height=8,dpi=300)


lake_temp_mod <- lmer(lake_wtemp ~ 
                        tmean_C +
                        precip_sum +
                        flow_sum +
                        #light_mean +
                        #spc_mean + 
                        (1|shore), data = GPP_df_w1)
summary(lake_temp_mod)
hist(residuals(lake_temp_mod))
vif(lake_temp_mod)



G_plot <- ggplot(data = GPP_df_w1, aes(x= (tmean_C), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_temp_atemp_w.png",sep=""),width=6,height=8,dpi=300)


G_plot <- ggplot(data = GPP_df_w1, aes(x = factor(precip_bi), y = lake_wtemp, fill = shore)) +
  geom_violin(alpha = 0.8) + 
  geom_point(alpha = 0.4)+
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.2) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_temp_precip_w.png",sep=""),width=4,height=8,dpi=300)



G_plot <- ggplot(data = GPP_df_w1%>%filter(shore=="BW"| shore=="GB"), aes(x= log(flow_sum+1), y=lake_wtemp, color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_temp_flow_w.png",sep=""),width=6,height=4,dpi=300)




gpp_auto3 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                   # tmean_C + 
                    lake_wtemp +
                    precip_sum +
                    windsp_mean +
                    flow_sum +
                    lake_par_int + (1|shore), data = GPP_df_w1)
summary(gpp_auto3)
residuals <- residuals(gpp_auto3)
vif(gpp_auto3)
hist(residuals(gpp_auto3))

G_plot <- ggplot(data = GPP_df_w1%>%filter(shore=="BW"| shore=="GB"), aes(x= log(flow_sum+1), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_flow_w.png",sep=""),width=6,height=4,dpi=300)



G_plot <- ggplot(data = GPP_df_week, aes(x = factor(precip_sum), y = log(middle_GPP+1), fill = shore)) +
  #geom_violin(alpha = 0.8) + 
  geom_point(alpha = 0.4)+
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.2) +  # Add box plot for reference
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  #labs(x = "Precipitation (Binary)", y = "Lake Water Temperature", title = "Distribution of Lake Water Temperature by Precipitation") +
  theme_bw() +
  facet_grid(shore ~ .)
G_plot
ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_precip_w.png",sep=""),width=4,height=8,dpi=300)


## report plot
gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(lake_par_int) + (1|shore / site) +(1|year), 
                 data = GPP_df1%>%filter(shore=="BW"|shore=="GB"|shore=="SH"))
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(light_mean) + (1|shore / site) +(1|year), 
                 data = GPP_df1%>%filter(shore=="BW"|shore=="GB"|shore=="SH"))
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot <- ggplot(data = GPP_df1%>%filter(shore=="BW"|shore=="GB"), aes(x= (lake_par_int), y=log(middle_GPP+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.5) + #geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_windsp_w.png",sep=""),width=6,height=8,dpi=300)


gpp_autow <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + scale(lake_par_int) + (1|shore / site) +(1|year), data = GPP_df_week%>%filter(shore=="BW"|shore=="GB"))
summary(gpp_autow)
hist(residuals(gpp_autow))
r.squaredGLMM(gpp_autow)

G_plot_lightw <- ggplot(data = GPP_df_week%>%filter(shore=="BW"|shore=="GB"|shore=="SH"), aes(x= (lake_par_int), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + #geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~PAR~umol~m^-2~s^-1)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_lightw

#ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_windsp_w.png",sep=""),width=6,height=8,dpi=300)


gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(windsp_mean) + (1|shore / site) +(1|year), 
                 data = GPP_df1)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot_windw <- ggplot(data = GPP_df_week, aes(x= (windsp_mean), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~wind~sp~m~s^-1)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_windw



gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(tmean_C) + (1|shore / site) +(1|year), 
                 data = GPP_df1)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot_temp <- ggplot(data = GPP_df1, aes(x= (tmean_C), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_temp


gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(lake_wtemp) + (1|shore / site) +(1|year), 
                 data = GPP_df1)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)

G_plot_ltempw <- ggplot(data = GPP_df_week, aes(x= (lake_wtemp), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  labs(y=expression(log(mean~weekly~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~weekly~lake~temp~C)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_ltempw


gpp_auto <- lmer(log(middle_GPP+1) ~  log(lag_middle_GPP+1) + scale(precip_bi) + (1|shore / site) +(1|year), 
                 data = GPP_df1)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)


gpp_auto <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + scale(precip_sum) + (1|shore / site) +(1|year), 
                 data = GPP_df_week)
summary(gpp_auto)
hist(residuals(gpp_auto))
r.squaredGLMM(gpp_auto)



G_plot_precipw <- ggplot(data = GPP_df_week, aes(x=as.factor(precip_sum), y=log(middle_GPP+1), fill=shore, color = shore))+
 # geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  # Add box plot for reference
  labs(y=expression(log(mean~weekly~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Sum~weekly~precip~events)) + 
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() #+ facet_grid(shore~.)
G_plot_precipw

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_temp_w.png",sep=""),width=6,height=8,dpi=300)


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
ggsave(plot = metab23, filename = paste("./SFS24_Analysis/figures/GPP1_weekly_plot.png",sep=""),width=12,height=4,dpi=300)



G_plot <- ggplot(data = GPP_df1, aes(x= (lake_par_int), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_par_w.png",sep=""),width=6,height=8,dpi=300)




G_plot <- ggplot(data = GPP_df_w1, aes(x=log(middle_GPP_lag+1), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(shore~.)
G_plot

ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/SEM_GPP_GPP_lag_w.png",sep=""),width=6,height=8,dpi=300)








# Use the `psem` function to create the SEM
model_a <- psem(
  gpp_auto3,
  lake_light_mod2,  
  lake_temp_mod)

# Look at object
model_a # Sometimes it works better to specify models in line, with called df
summary(model_a)
plot(model_a)
fisherC(model) # says model



gpp_auto2 <- lmer(log(middle_GPP+1) ~  log(middle_GPP_lag+1) + 
                    # tmean_C + 
                    lake_wtemp +
                    #precip_sum +
                    #windsp_mean +
                    flow_sum +
                    lake_par_int + (1|shore), data = GPP_df_w1)
summary(gpp_auto2)
residuals <- residuals(gpp_auto2)
vif(gpp_auto2)
hist(residuals(gpp_auto2))





# Use the `psem` function to create the SEM
model_b <- psem(
  gpp_auto2,
  lake_light_mod2,  
  lake_temp_mod)

# Look at object
model_b # Sometimes it works better to specify models in line, with called df
summary(model_b)
plot(model_b)
fisherC(model_b) # says model



model2 <- psem(
  gpp_auto1,
  lake_temp_mod3,  
  lake_light_mod)

summary(model2)
plot(model2)

AIC(model, model2)


######################
#####################

################################################################################
GPP_df1$year<- year(GPP_df1$date)


GPP_df_week$precip_bi <- ifelse(GPP_df_week$precip_sum > 0, 1, 0)
GPP_df_week$middle_GPP_lag <- lag(GPP_df_week$middle_GPP)


