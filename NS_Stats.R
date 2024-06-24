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


SFS_datQ1 <- SFS_datQ %>%
  mutate(
    position = case_when(
      site %in% c("BWNS1","GBNS1", "SSNS1", "SHNS1") ~ "north",
      site %in% c("BWNS2","GBNS2", "SSNS2", "SHNS2" ) ~ "center",
      site %in% c("BWNS3","GBNS3", "SHNS3") ~ "south",
      TRUE ~ as.character(site)
    )
  )



## create column for year and doy 
SFS_datQ1$year <- year(SFS_datQ1$date)
SFS_datQ1$yday <- yday(SFS_datQ1$date)

## creat new df for complete GPP obs 
GPP_df <- SFS_datQ1%>%
  dplyr::select(-middle_NEP, -upper_NEP, -lower_NEP, 
                -middle_ER, -upper_ER, -lower_ER)

GPP_df<- GPP_df%>%
  drop_na(middle_GPP)

ER_df <- SFS_datQ1%>%
  dplyr::select(-middle_NEP, -upper_NEP, -lower_NEP, 
                -middle_GPP, -upper_GPP, -lower_GPP)

ER_df<- ER_df%>%
  drop_na(middle_ER)

#GPP_df1<- na.omit(GPP_df)
summary(GPP_df)

##===========================================
## Begin analysis 
#============================================
# lake condition across sites



# When is GPP highest? lowest? 
# Test for different in Metab across shores 

BW_dat <- SFS_datQ %>%
  filter(shore=="BW")

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

summary(BW_dat1)
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

GB_dat <- SFS_datQ %>%
  filter(shore=="GB")
GB_dat1 <- GB_dat %>%
  filter(date> as.Date("2023-02-07") & date < as.Date("2023-09-06"))
unique(GB_dat$shore)
round(mean(na.omit(GB_dat$middle_GPP)),2)
round(se(na.omit(GB_dat$middle_GPP)),2)
range(na.omit(GB_dat$middle_GPP))
round(mean(na.omit(GB_dat$middle_ER)),2)
round(se(na.omit(GB_dat$middle_ER)),2)
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

gpp_fit1 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                      (1|shore/site), data = GPP_df, family = lognormal)
summary(gpp_fit1)
residuals <- residuals(gpp_fit1)
hist(residuals)
# Plot ACF and PACF of residuals
acf(residuals, main="ACF of Residuals")
pacf(residuals, main="PACF of Residuals")


gpp_fit2 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                      lag(middle_GPP, 2) +
                      (1|shore/site), data = GPP_df, family = gaussian)
summary(gpp_fit2)
residuals2 <- residuals(gpp_fit2)
hist(residuals2)
# Plot ACF and PACF of residuals
acf(residuals2, main="ACF of Residuals")
pacf(residuals2, main="PACF of Residuals")

AIC(gpp_fit1, gpp_fit2)

## 

# adjust to whatever covariate rows we may want to include in a mixed model:
corplot <- GPP_df[,c(13,21:26)]
chart.Correlation(corplot, histogram=TRUE, pch=19)
# light and air temp variables shouldn't be in the same model 


##===========================================
## GPP glms for stream influence
#============================================
# subset for shore BW and GB
GPP_df_stream <- GPP_df %>%
  filter(shore=="BW" |shore=="GB")

ER_df_stream <- ER_df %>%
  filter(shore=="BW" |shore=="GB")

##===========================================
## check co-variance 
library(PerformanceAnalytics)
names(GPP_df_stream)

# scale variables
# GPP_df_stream[, c(7, 19:26)] <- scale(GPP_df_stream[, c(7, 19:26)])
Streamdf_cor <- GPP_df_stream[, c(7, 19:26)]

chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

hist(GPP_df_stream$middle_GPP) # strong left skew
hist(log(GPP_df_stream$middle_GPP+1)) # better

## ER ##
ER_df_stream[, c(7, 19:26)] <- scale(ER_df_stream[, c(7, 19:26)])
Streamdf_cor <- ER_df_stream[, c(7, 19:26)]
chart.Correlation(Streamdf_cor, histogram=TRUE, pch=19)

GPP_df_stream$precip_bi <- ifelse(GPP_df_stream$ppt_mm > 0, 1, 0)

hist(log(GPP_df_stream$ppt_mm+1))

GPP_df_stream$log_ppmt <- log(GPP_df_stream$ppt_mm+1)
GPP_df_stream$log_streamflow <- (log(GPP_df_stream$flow_mean)+1)
names(GPP_df_stream)


ER_df_stream$middle_ER<- (-1*ER_df_stream$middle_ER)
hist(ER_df_stream$middle_ER)
ER_df_stream$ER<- log(ER_df_stream$middle_ER+1)
hist(ER_df_stream$ER)

ER_df_stream$precip_bi <- ifelse(ER_df_stream$ppt_mm > 0, 1, 0)

ER_df_stream$log_ppmt <- log(ER_df_stream$ppt_mm+1)
ER_df_stream$log_streamflow <- (log(ER_df_stream$flow_mean)+1)

ER_df_stream <- ER_df_stream %>%
  mutate(ER_lag = lag(ER, default = NA))


## 
lakedf_cor <- GPP_df_stream[, c("stream_temp", "windsp_mean", "flow_mean", "spc_mean", "log_streamflow", "log_ppmt", "precip_bi")]
chart.Correlation(lakedf_cor, histogram=TRUE, pch=19)

GPP_df_stream <- GPP_df_stream %>%
  mutate(lag_middle_GPP = lag(middle_GPP, default = NA))


GPP_df_stream$GPP <- log(GPP_df_stream$middle_GPP+1)

GPP_df_stream <- GPP_df_stream %>%
  mutate(GPP_lag = lag(GPP, default = NA))



##==========================
## GPP and stream influence 
## NULL:
gpp_stream_mod <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                      lag(middle_GPP, 2)+
                      # lag(middle_GPP, 3) +
                      #   lag(middle_GPP, 4) +
                      #   lag(middle_GPP, 5) +
                      #   lag(middle_GPP, 6) +
                      #   lag(middle_GPP, 7) +
                      (1|site), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals3)

AIC(gpp_stream_mod, gpp_stream_mod7)

acf(residuals3, main="ACF of Residuals")
pacf(residuals3, main="PACF of Residuals")


names(GPP_df_stream)

## Predictors
gpp_stream_mod2 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                           # lag(middle_GPP, 3) + 
                             log_streamflow  + 
                             stream_temp +
                            (1|site) + (1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod2)
residuals3 <- residuals(gpp_stream_mod2)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod2)

gpp_stream_mod1 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                             lag(middle_GPP, 2)+
                             # lag(middle_GPP, 3) + 
                             flow_mean + 
                             spc_mean + 
                             stream_temp +
                             + (1|shore / site) +(1|year), data = GPP_df_stream, family = lognormal)
summary(gpp_stream_mod1)
residuals3 <- residuals(gpp_stream_mod1)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod1)


## could be over predicting .... 
gpp_stream_mod_b <- lmer(GPP ~ GPP_lag +
                        flow_mean + 
                        spc_mean + 
                        stream_temp +
                          + (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod_b)
hist(residuals(gpp_stream_mod_b))
vif(gpp_stream_mod_b)
r.squaredGLMM(gpp_stream_mod_b)

AIC(gpp_stream_mod2, gpp_stream_mod1, gpp_stream_mod_b)


# lets look at individual responses from best fit model to see how they track:
summary(GPP_df_stream)
GGP_SPC_plot <- ggplot(data = GPP_df_stream, aes(x= (spc_mean), y=log(middle_GPP+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.5) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(.~position)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(spc_mean = GPP_df_stream_complete$spc_mean, middle_GPP = fitted_values)
# Plot
G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = spc_mean, y = log(middle_GPP), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = (middle_GPP*-0.001)), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 


# Filter out missing values from the original data frame
GPP_df_stream_complete <- na.omit(GPP_df_stream)

# Extract fitted values for the complete observations
fitted_values <- predict(gpp_stream_mod1, newdata = GPP_df_stream_complete)

# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(flow_mean = GPP_df_stream_complete$flow_mean, middle_GPP = fitted_values)

# 


G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = flow_mean, y = log(middle_GPP), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = middle_GPP), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

G_plot


# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GGP_SPC_plot_a.png",sep=""), width=5,height=4,dpi=300)



##########
######### REPORT FIGS
##########

gpp_stream_mod <- lmer(GPP ~ GPP_lag +
                         log(flow_mean+1) + 
                         (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod)
hist(residuals(gpp_stream_mod))
vif(gpp_stream_mod)
r.squaredGLMM(gpp_stream_mod)
###
GGP_SF_plot <- ggplot(data = GPP_df_stream, aes(x= log(flow_mean+1), y=log(middle_GPP+1), color = shore))+
  geom_point(size= 3, alpha = 0.25) + geom_smooth(method="lm", se=F, aes(color=shore)) +
  labs(y=expression(log(mean~daily~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(log(mean~flow+1)~m^3~s^-1~km^2)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
ggsave(plot = GGP_SF_plot, filename = paste("./SFS24_Analysis/figures/GGP_SF_plot_a.png",sep=""), width=5,height=4,dpi=300)



gpp_stream_mod <- lmer(GPP ~ GPP_lag + scale(stream_temp) + 
                         (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod)
hist(residuals(gpp_stream_mod))
vif(gpp_stream_mod)
r.squaredGLMM(gpp_stream_mod)

# Plot
GGP_ST_plot <- ggplot(data = GPP_df_stream, aes(x = stream_temp, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Stream~temp.~C)) + 
  geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_ST_plot, filename = paste("./SFS24_Analysis/figures/GGP_ST_plot.png",sep=""), width=5,height=4,dpi=300)


gpp_stream_mod <- lmer(GPP ~ GPP_lag + 
                         scale(lake_par_int) + 
                         (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod)
hist(residuals(gpp_stream_mod))
vif(gpp_stream_mod)
r.squaredGLMM(gpp_stream_mod)

GGP_PAR_plot <- ggplot(data = GPP_df_stream, aes(x = lake_par_int, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~PAR~umol~m^-2~s^-1)) + 
  #geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_PAR_plot, filename = paste("./SFS24_Analysis/figures/GGP_PAR_plot.png",sep=""), width=5,height=4,dpi=300)

##

gpp_stream_mod <- lmer(GPP ~ GPP_lag + scale(windsp_mean) + 
                         (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod)
hist(residuals(gpp_stream_mod))
vif(gpp_stream_mod)
r.squaredGLMM(gpp_stream_mod)

GGP_windsp_plot <- ggplot(data = GPP_df_stream, aes(x = windsp_mean, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~wind~speed~m~s^-1)) + 
  geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_PAR_plot, filename = paste("./SFS24_Analysis/figures/GGP_PAR_plot.png",sep=""), width=5,height=4,dpi=300)

## 
gpp_stream_mod <- lmer(GPP ~ GPP_lag + log(log_ppmt+1) + 
                         (1|shore / site) +(1|year), data = GPP_df_stream)
summary(gpp_stream_mod)
hist(residuals(gpp_stream_mod))
vif(gpp_stream_mod)
r.squaredGLMM(gpp_stream_mod)

GGP_ppmt_plot <- ggplot(data = GPP_df_stream, aes(x = log(log_ppmt+1), y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~GPP+1)~mmol~O[2]~m^-3~d^-1), x=expression(log(mean~precip+1)~mm)) + 
  #geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_ppmt_plot, filename = paste("./SFS24_Analysis/figures/GGP_ppmt_plot.png",sep=""), width=5,height=4,dpi=300)


#### 
metab23 <- ggarrange(GGP_SF_plot, 
                     GGP_ST_plot,
                     GGP_windsp_plot,
                     GGP_PAR_plot,
                     GGP_ppmt_plot,
                     ncol = 5, nrow = 1,
                     common.legend = TRUE, 
                     legend = "bottom")
###
ggsave(plot = metab23, filename = paste("./SFS24_Analysis/figures/GPP_dailystream_plot.png",sep=""),width=16,height=4,dpi=300)

########## 
########## report figures end daily stream GPP
##########

##==========================
## ER and stream influence 

## NULL:
ER_stream_mod <- glmmTMB((middle_ER) ~ lag(middle_ER, 1) + 
                            lag(middle_ER, 2)+
                            # lag(middle_GPP, 3) +
                            #   lag(middle_GPP, 4) +
                            #   lag(middle_GPP, 5) +
                            #   lag(middle_GPP, 6) +
                            #   lag(middle_GPP, 7) +
                            (1|site), data = ER_df_stream, family = lognormal)
summary(ER_stream_mod)
residuals3 <- residuals(ER_stream_mod)
hist(residuals3)
acf(residuals3, main="ACF of Residuals")
pacf(residuals3, main="PACF of Residuals")


## Predictors
ER_stream_mod1 <- glmmTMB((middle_ER) ~ lag(middle_ER, 1) + 
                             lag(middle_ER, 2)+
                            flow_mean + 
                            spc_mean + 
                            stream_temp +
                             (1|site) + (1|year), data = ER_df_stream, family = lognormal)
summary(ER_stream_mod1)
residuals3 <- residuals(ER_stream_mod1)
hist(residuals3)
r.squaredGLMM(ER_stream_mod1)




AIC(ER_stream_mod, ER_stream_mod1, ER_stream_mod2, ER_stream_mod3)


# lets look at individual responses from best fit model to see how they track:
summary(ER_df_stream)
ER_SPC_plot <- ggplot(data = ER_df_stream, aes(x= (spc_mean), y=log(middle_ER), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.6) + geom_smooth(method="lm", se=F, color="grey") +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + facet_grid(.~position)


# Filter out missing values from the original data frame
ER_df_stream_complete <- na.omit(ER_df_stream)

# Extract fitted values for the complete observations
fitted_values <- predict(ER_stream_mod1, newdata = ER_df_stream_complete)

# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(spc_mean = ER_df_stream_complete$spc_mean, middle_ER = fitted_values)

# Plot
ER_SPC_plot <- ggplot(data = ER_df_stream_complete, aes(x = spc_mean, y = log(middle_ER), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = middle_ER), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + 
  facet_grid(. ~ position)

ER_SPC_plot

ggsave(plot = ER_SPC_plot, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/ER_SPC_plot_a.png",sep=""),
       width=5,height=4,dpi=300)



########## 
########## report figures daily ER
##########


stream_mod <- lmer(ER ~ ER_lag +
                         log(flow_mean+1) + 
                         (1|shore / site) +(1|year), data = ER_df_stream)
summary(stream_mod)
hist(residuals(stream_mod))
vif(stream_mod)
r.squaredGLMM(stream_mod)
###
ER_SF_plot <- ggplot(data = ER_df_stream, aes(x= log(flow_mean+1), y=log(middle_ER+1), color = shore))+
  geom_point(size= 3, alpha = 0.25) + #geom_smooth(method="lm", se=F, aes(color=shore)) +
  labs(y=expression(log(mean~daily~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(log(mean~flow+1)~m^3~s^-1~km^2)) + 
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw()
ggsave(plot = ER_SF_plot, filename = paste("./SFS24_Analysis/figures/ER_SF_plot_a.png",sep=""), width=5,height=4,dpi=300)



stream_mod <- lmer(ER ~ ER_lag + 
                     scale(stream_temp) + 
                         (1|shore / site) +(1|year), data = ER_df_stream)
summary(stream_mod)
hist(residuals(stream_mod))
vif(stream_mod)
r.squaredGLMM(stream_mod)

# Plot
ER_ST_plot <- ggplot(data = ER_df_stream, aes(x = stream_temp, y =ER, color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Stream~temp.~C)) + 
 # geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = ER_ST_plot, filename = paste("./SFS24_Analysis/figures/ER_ST_plot.png",sep=""), width=5,height=4,dpi=300)


stream_mod <- lmer(ER ~ #ER_lag + 
                     scale(lake_par_int) + 
                         (1|shore / site) +(1|year), data = ER_df_stream)
summary(stream_mod)
hist(residuals(stream_mod))
vif(stream_mod)
r.squaredGLMM(stream_mod)

ER_PAR_plot <- ggplot(data = ER_df_stream, aes(x = lake_par_int, y =ER, color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~PAR~umol~m^-2~s^-1)) + 
  geom_smooth(method="lm", se=F, lty=1) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = ER_PAR_plot, filename = paste("./SFS24_Analysis/figures/ER_PAR_plot.png",sep=""), width=5,height=4,dpi=300)

##

stream_mod <- lmer(ER ~ER_lag + 
                     scale(windsp_mean) + 
                         (1|shore / site) +(1|year), data = ER_df_stream)
summary(stream_mod)
hist(residuals(stream_mod))
vif(stream_mod)
r.squaredGLMM(stream_mod)

ER_windsp_plot <- ggplot(data = ER_df_stream, aes(x = windsp_mean, y =ER, color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(Mean~wind~speed~m~s^-1)) + 
  #geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_PAR_plot, filename = paste("./SFS24_Analysis/figures/GGP_PAR_plot.png",sep=""), width=5,height=4,dpi=300)


stream_mod <- lmer(ER ~ #ER_lag + 
                     log(log_ppmt+1) + 
                         (1|shore / site) +(1|year), data = ER_df_stream)
summary(stream_mod)
hist(residuals(stream_mod))
vif(stream_mod)
r.squaredGLMM(stream_mod)

ER_ppmt_plot <- ggplot(data = ER_df_stream, aes(x = log(log_ppmt+1), y = ER, color = shore)) +
  geom_point(size = 3, alpha = 0.25, aes(shape=position)) +
  labs(y=expression(log(mean~daily~ER+1)~mmol~O[2]~m^-3~d^-1), x=expression(log(mean~precip+1)~mm)) + 
  #geom_smooth(method="lm", se=F, lty=2) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 

# ggsave(plot = GGP_ppmt_plot, filename = paste("./SFS24_Analysis/figures/GGP_ppmt_plot.png",sep=""), width=5,height=4,dpi=300)


#### 
metab23 <- ggarrange(ER_SF_plot, 
                     ER_ST_plot,
                     ER_windsp_plot,
                     ER_PAR_plot,
                     ER_ppmt_plot,
                     ncol = 5, nrow = 1,
                     common.legend = TRUE, 
                     legend = "bottom")
###
ggsave(plot = metab23, filename = paste("./SFS24_Analysis/figures/ER_dailystream_plot.png",sep=""),width=16,height=4,dpi=300)







##====================================
## GPP and stream influence overtime

GPP_stream_mod1 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                            scale(year) +
                            (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod1)
residuals3 <- residuals(GPP_stream_mod1)
hist(residuals3)
r.squaredGLMM(GPP_stream_mod1)

GPP_stream_mod2 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                             lag(middle_GPP, 2)+
                             scale(log_streamflow) + 
                             scale(year) +
                             (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod2)
residuals3 <- residuals(GPP_stream_mod2)
hist(residuals2)
r.squaredGLMM(GPP_stream_mod2)



GPP_stream_mod3 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                            scale(log_streamflow) + 
                            scale(log_streamflow * year)+
                            scale(year) +
                            (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod3)
residuals3 <- residuals(GPP_stream_mod2)
hist(residuals2)
r.squaredGLMM(GPP_stream_mod3)

AIC(GPP_stream_mod3, GPP_stream_mod2, GPP_stream_mod1)

GPP_stream_mod3 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                            log(stream_flow+1)  + 
                            scale(Stream_temp) +
                            scale(stream_SPC) +
                            scale(stream_SPC*year)+
                            scale(Stream_temp*year)+
                            scale(log(stream_flow+1) * year)+
                            scale(year) +
                            (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod3)
residuals3 <- residuals(GPP_stream_mod3)
hist(residuals3)
r.squaredGLMM(GPP_stream_mod3)


GPP_stream_mod4 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                            log(stream_flow+1)  + 
                            scale(Stream_temp) +
                            scale(Stream_temp*year)+
                            scale(log(stream_flow+1) * year)+
                            scale(year) +
                            (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod4)
residuals3 <- residuals(GPP_stream_mod4)
hist(residuals3)
r.squaredGLMM(GPP_stream_mod4)



GPP_stream_mod5 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                             lag(middle_GPP, 2)+
                             scale(Stream_temp) +
                             scale(Stream_temp*year)+
                             scale(year) +
                             (1|site), data = GPP_df_stream, family = lognormal)
summary(GPP_stream_mod5)
residuals3 <- residuals(GPP_stream_mod5)
hist(residuals3)
r.squaredGLMM(GPP_stream_mod5)

AIC(GPP_stream_mod1, GPP_stream_mod2, GPP_stream_mod3, GPP_stream_mod4, GPP_stream_mod5)

GPP_ST_plot_yr <- ggplot(data = GPP_df_stream, aes(x= log_streamflow, y=log(middle_GPP+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.6) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  facet_grid(position ~ year)

ggsave(plot = GPP_ST_plot_yr, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/GPP_STemp_year_plot.png",sep=""),
       width=10,height=8,dpi=300)


##====================================
## ER and stream influence overtime



ER_stream_mod1 <- glmmTMB((middle_ER) ~ lag(middle_ER, 1) + 
                            lag(middle_ER, 2)+
                            scale(year) +
                            (1|site), data = ER_df_stream, family = lognormal)
summary(ER_stream_mod1)
residuals3 <- residuals(ER_stream_mod1)
hist(residuals3)
r.squaredGLMM(ER_stream_mod1)




ER_stream_mod2 <- glmmTMB((middle_ER) ~ lag(middle_ER, 1) + 
                            lag(middle_ER, 2)+
                            scale(log_streamflow)+
                            scale(year) +
                            (1|site), data = ER_df_stream, family = lognormal)
summary(ER_stream_mod2)
residuals3 <- residuals(ER_stream_mod2)
hist(residuals2)
r.squaredGLMM(ER_stream_mod2)


ER_stream_mod3 <- glmmTMB((middle_ER) ~ lag(middle_ER, 1) + 
                            lag(middle_ER, 2)+
                            scale(log_streamflow)+
                            scale(year) +
                            scale(log_streamflow*year) +
                            (1|site), data = ER_df_stream, family = lognormal)
summary(ER_stream_mod3)
residuals3 <- residuals(ER_stream_mod3)
hist(residuals2)
r.squaredGLMM(ER_stream_mod3)



AIC(ER_stream_mod1, ER_stream_mod2, ER_stream_mod3)


predicted_middle_ER <- predict(ER_stream_mod2, type = "response")
plot(ER_df_stream$middle_ER ~ log(ER_df_stream$stream_flow + 1) * ER_df_stream$year, 
     xlab = "log(stream_flow + 1) * year", ylab = "middle_ER",
     main = "middle_ER vs. log(stream_flow + 1) * year",
     pch = 16, col = "blue")

ER_SF_plot_yr <- ggplot(data = ER_df_stream, aes(x=log_streamflow, y=log(middle_ER+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.6) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  facet_grid(position ~ year)

ggsave(plot = ER_SF_plot_yr, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/ER_flow_year_plot.png",sep=""),
       width=10,height=8,dpi=300)

ER_ST_plot_yr <- ggplot(data = ER_df_stream, aes(x= Stream_temp, y=log(middle_ER+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.6) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  facet_grid(position ~ year)


ggsave(plot = ER_ST_plot_yr, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/ER_STemp_year_plot.png",sep=""),
       width=10,height=8,dpi=300)

ER_SPC_plot_yr <- ggplot(data = ER_df_stream, aes(x= stream_SPC, y=log(middle_ER+1), color = shore, shape=position))+
  geom_point(size= 3, alpha = 0.6) + geom_smooth(method="lm", se=F) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +  facet_grid(position ~ year)

ggsave(plot = ER_SPC_plot_yr, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/ER_SPC_year_plot.png",sep=""),
       width=10,height=8,dpi=300)



##====================================
## SO what best explains gpp?
##====================================

library(Boruta)
library(mlbench)
library(randomForest)

## select data: ER_df
ER_df$precip_bi <- ifelse(ER_df$ppt_mm > 0, 1, 0)
ER_df$log_ppmt <- log(ER_df$ppt_mm+1)
ER_df$log_streamflow <- (log(ER_df$flow_mean)+1)
names(ER_df)
ER_df1<- ER_df%>%
  drop_na(lake_par_int)

ER_df2<- ER_df1%>%
  drop_na(middle_ER)
names(ER_df1)


## Select data for all four sites - no stream variables:
summary(ER_df2[,c(7:16,31:39,43,44)])

str(ER_df2[,c(7:16,31:39,43,44)])

boruta1 <- Boruta(ER_df2$middle_ER ~ ., data = ER_df2[,c(7:16,31:39,43,44)], doTrace = 2, maxRuns = 500)
print(boruta1)

summary(ER_df_stream[,c(7:26,28:39,43:45)])

png('./SFS24_Analysis/figures/ER_RF_nonstream_plot.png', width = 1000, height = 800, res = 150)
plot(boruta1, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta1$ImpHistory),function(i)
  boruta1$ImpHistory[is.finite(boruta1$ImpHistory[,i]),i])
names(lz) <- colnames(boruta1$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta1$ImpHistory), cex.axis = 0.45)
dev.off()



## select data: ER_df
GPP_df$precip_bi <- ifelse(GPP_df$ppt_mm > 0, 1, 0)
GPP_df$log_ppmt <- log(GPP_df$ppt_mm+1)
GPP_df$log_streamflow <- (log(GPP_df$flow_mean)+1)
names(GPP_df)
GPP_df1<- GPP_df%>%
  drop_na(lake_par_int)

GPP_df<- GPP_df%>%
  drop_na(lake_par_int)
names(ER_df1)

summary(GPP_df[,c(7:16,31:39,43,44)])

str(ER_df2[,c(7:16,31:39,43,44)])

boruta1 <- Boruta(GPP_df$middle_GPP ~ ., data = GPP_df[,c(7:16,31:39,43,44)], doTrace = 2, maxRuns = 500)
print(boruta1)

png('./SFS24_Analysis/figures/GPP_RF_nonstream_plot.png', width = 1000, height = 800, res = 150)
plot(boruta1, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta1$ImpHistory),function(i)
  boruta1$ImpHistory[is.finite(boruta1$ImpHistory[,i]),i])
names(lz) <- colnames(boruta1$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta1$ImpHistory), cex.axis = 0.45)
dev.off()



##====================================
## SO what best explains gpp/ER?

GPP_df_stream

summary(GPP_df_stream)

GPP_df_stream1<- GPP_df_stream%>%
  drop_na(spc_mean) #lake_par_int

GPP_df_stream2<- GPP_df_stream1%>%
  drop_na(lake_par_int)

GPP_df_stream3<- GPP_df_stream2%>%
  drop_na(flow_mean)


summary(GPP_df_stream2)

## Select data for all four sites - no stream variables:
summary(GPP_df_stream2[,c(7:26,31:39,43,44,45)])

str(ER_df2[,c(7:16,31:39,43,44)])

boruta_streamGPP <- Boruta(GPP_df_stream2$middle_GPP ~ ., data = GPP_df_stream2[,c(7:26,31:39,43,44,45)], doTrace = 2, maxRuns = 500)
print(boruta_streamGPP)

summary(ER_df_stream[,c(7:26,28:39,43:45)])

png('./SFS24_Analysis/figures/GPP_RF_stream_plot.png', width = 1000, height = 800, res = 150)
plot(boruta_streamGPP, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_streamGPP$ImpHistory),function(i)
  boruta_streamGPP$ImpHistory[is.finite(boruta_streamGPP$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_streamGPP$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_streamGPP$ImpHistory), cex.axis = 0.45)
dev.off()


















summary(GPP_df_stream)

GPP_df_stream1<- GPP_df_stream%>%
  drop_na(spc_mean) #lake_par_int

GPP_df_stream2<- GPP_df_stream1%>%
  drop_na(lake_par_int)

GPP_df_stream3<- GPP_df_stream2%>%
  drop_na(flow_mean)


summary(GPP_df_stream2)

## Select data for all four sites - no stream variables:
summary(GPP_df_stream2[,c(7:26,31:39,43,44,45)])

str(ER_df2[,c(7:16,31:39,43,44)])

boruta_streamGPP <- Boruta(GPP_df_stream2$middle_GPP ~ ., data = GPP_df_stream2[,c(7:26,31:39,43,44,45)], doTrace = 2, maxRuns = 500)
print(boruta_streamGPP)

summary(ER_df_stream[,c(7:26,28:39,43:45)])

png('./SFS24_Analysis/figures/GPP_RF_stream_plot.png', width = 1000, height = 800, res = 150)
plot(boruta_streamGPP, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_streamGPP$ImpHistory),function(i)
  boruta_streamGPP$ImpHistory[is.finite(boruta_streamGPP$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_streamGPP$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_streamGPP$ImpHistory), cex.axis = 0.45)
dev.off()





summary(ER_df_stream)

ER_df_stream1<- ER_df_stream%>%
  drop_na(spc_mean) #lake_par_int

ER_df_stream2<- ER_df_stream1%>%
  drop_na(lake_par_int)

ER_df_stream3<- ER_df_stream2%>%
  drop_na(flow_mean)


summary(ER_df_stream2)

## Select data for all four sites - no stream variables:
summary(ER_df_stream2[,c(7:26,31:39,43,44,45)])

str(ER_df2[,c(7:16,31:39,43,44)])

boruta_streamER <- Boruta(ER_df_stream2$middle_ER ~ ., data = ER_df_stream2[,c(7:26,31:39,43,44,45)], doTrace = 2, maxRuns = 500)
print(boruta_streamER)

summary(ER_df_stream[,c(7:26,28:39,43:45)])

png('./SFS24_Analysis/figures/ER_RF_stream_plot.png', width = 1000, height = 800, res = 150)
plot(boruta_streamGPP, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_streamER$ImpHistory),function(i)
  boruta_streamER$ImpHistory[is.finite(boruta_streamER$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_streamER$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_streamER$ImpHistory), cex.axis = 0.45)
dev.off()

GPP_df_stream$week <- week(GPP_df_stream$date)

## 
GPP_df_stream_week <- GPP_df_stream %>%
  group_by(site, shore, week, year, position) %>%
  summarise(middle_GPP=mean(middle_GPP, na.rm = TRUE),
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
            



## 

GPP_df_stream_week <- GPP_df_stream_week %>%
  mutate(lag_middle_GPP = lag(middle_GPP, 1))


GPP_df_stream_week <- GPP_df_stream_week %>%
  arrange(site, shore, week, year, position) %>%
  mutate(lag_middle_GPP = lag(middle_GPP),
         lag_middle_GPP = if_else(is.na(lag_middle_GPP), lag(middle_GPP), lag_middle_GPP, missing = middle_GPP))

GPP_df_stream$GPP <- log(GPP_df_stream$middle_GPP+1)
GPP_df_stream$GPP_lag <- log(GPP_df_stream$lag_middle_GPP+1)

hist(GPP_df_stream_week$middle_GPP)

GPP_df_stream_weeks<- GPP_df_stream_week%>%
  filter(week>21 & week<37)

summary(GPP_df_stream_week)
hist(GPP_df_stream_weeks$middle_GPP)
hist(log(GPP_df_stream_weeks$middle_GPP+1))

##==========================
## GPP and stream influence 
## NULL:
## pull out just BW and gb 

# Double check the correlation coeff: 


gpp_stream_mod <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            lag(middle_GPP, 2)+
                            # lag(middle_GPP, 3) +
                            #   lag(middle_GPP, 4) +
                            #   lag(middle_GPP, 5) +
                            #   lag(middle_GPP, 6) +
                            #   lag(middle_GPP, 7) +
                            (1|site), data = GPP_df_stream_week, family = lognormal)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals3)

## NULL:
gpp_stream_mod1 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                            #lag(middle_GPP, 2)+
                            # lag(middle_GPP, 3) +
                            #   lag(middle_GPP, 4) +
                            #   lag(middle_GPP, 5) +
                            #   lag(middle_GPP, 6) +
                            #   lag(middle_GPP, 7) +
                            (1|site), data = GPP_df_stream_week, family = lognormal)
summary(gpp_stream_mod1)
residuals3 <- residuals(gpp_stream_mod1)
hist(residuals3)

AIC(gpp_stream_mod, gpp_stream_mod1)

acf(residuals3, main="ACF of Residuals")
pacf(residuals3, main="PACF of Residuals")


lakedf_cor <- GPP_df_stream_week[, c("stream_temp", "tmean_C","lake_par_int","windsp_mean", "flow_sum", "spc_mean", "precip_sum")]
chart.Correlation(lakedf_cor, histogram=TRUE, pch=19)

## Predictors
gpp_stream_mod2 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                             # lag(middle_GPP, 3) + 
                             scale(stream_temp)  + 
                             scale(flow_sum)+
                             scale(precip_sum)+
                             (1|site) + (1|year), data = GPP_df_stream_week, family = lognormal)
summary(gpp_stream_mod2)
residuals3 <- residuals(gpp_stream_mod2)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod2)

gpp_stream_mod4 <- glmmTMB((middle_GPP) ~ lag(middle_GPP, 1) + 
                             # lag(middle_GPP, 3) + 
                             scale(stream_temp)  + 
                             scale(flow_sum)+
                             scale(precip_sum)+
                             scale(year) +
                             (1|site), data = GPP_df_stream_week, family = lognormal)
summary(gpp_stream_mod4)
residuals3 <- residuals(gpp_stream_mod4)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod4)

AIC(gpp_stream_mod, gpp_stream_mod1, gpp_stream_mod2, gpp_stream_mod3, gpp_stream_mod4)

###
gpp_stream_mod <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_stream_mod)
residuals3 <- residuals(gpp_stream_mod)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod)

gpp_stream_mod1 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                             scale(stream_temp)  + 
                             scale(flow_sum)+
                             scale(precip_sum)+
                             (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_stream_mod1)
residuals3 <- residuals(gpp_stream_mod1)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod1)

gpp_stream_mod2 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                          scale(stream_temp)  + 
                          scale(flow_sum)+
                          scale(precip_sum)+
                          scale(year)+
                          (1|site), data = GPP_df_stream_week)
summary(gpp_stream_mod2)
residuals3 <- residuals(gpp_stream_mod2)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod2)

gpp_stream_mod4 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                          scale(stream_temp)  + 
                          scale(flow_sum)+
                          scale(precip_sum)+
                          scale(year)+
                          scale(windsp_mean)+
                          (1|site), data = GPP_df_stream_week)
summary(gpp_stream_mod4)
residuals3 <- residuals(gpp_stream_mod4)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod4)
vif(gpp_stream_mod4)


gpp_stream_mod3 <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                          scale(stream_temp)  + 
                          scale(flow_sum)+
                          scale(precip_sum)+
                          scale(year)+
                          (1|site), data = GPP_df_stream_week)
summary(gpp_stream_mod3)
residuals3 <- residuals(gpp_stream_mod3)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod3)
vif(gpp_stream_mod3)

AIC(gpp_stream_mod, gpp_stream_mod1, gpp_stream_mod2, gpp_stream_mod3, gpp_stream_mod4)

##### not stream stuff

gpp_s_mod1 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                          scale(tmean_C)  + 
                          scale(flow_sum)+
                          scale(precip_sum)+
                          scale(windsp_mean)+
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod1)
residuals3 <- residuals(gpp_s_mod1)
hist(residuals3)
r.squaredGLMM(gpp_s_mod1)
vif(gpp_s_mod1)



gpp_s_mod1a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                     scale(tmean_C)  + 
                     scale(flow_sum)+
                     scale(precip_sum)+
                     scale(windsp_mean)+
                     (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod1a)
residuals3 <- residuals(gpp_s_mod1a)
hist(residuals3)
r.squaredGLMM(gpp_s_mod1a)
vif(gpp_s_mod1a)



gpp_s_mod1b <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                      scale(tmean_C)  + 
                      #scale(flow_sum)+
                      scale(precip_sum)+
                      scale(windsp_mean)+
                      (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod1b)
residuals3 <- residuals(gpp_s_mod1b)
hist(residuals3)
r.squaredGLMM(gpp_s_mod1b)
vif(gpp_s_mod1b)


gpp_s_mod2 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                     scale(lake_par_int)  + 
                     scale(flow_sum)+
                     scale(precip_sum)+
                     scale(windsp_mean)+
                     (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod2)
residuals3 <- residuals(gpp_s_mod2)
hist(residuals3)
r.squaredGLMM(gpp_s_mod2)
vif(gpp_s_mod1)


gpp_s_mod2a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                     scale(lake_par_int)  + 
                     scale(flow_sum)+
                     scale(precip_sum)+
                     scale(windsp_mean)+
                     (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod2a)
residuals3 <- residuals(gpp_s_mod2a)
hist(residuals3)
r.squaredGLMM(gpp_s_mod2a)
vif(gpp_s_mod1a)


gpp_s_mod2b <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                      scale(lake_par_int)  + 
                      #scale(flow_sum)+
                      scale(precip_sum)+
                      scale(windsp_mean)+
                      (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod2b)
residuals3 <- residuals(gpp_s_mod2b)
hist(residuals3)
r.squaredGLMM(gpp_s_mod2b)
vif(gpp_s_mod1b)



gpp_s_mod3 <- lmer(log(middle_GPP+1) ~ log(lag_middle_GPP+1) + 
                      scale(light_mean)  + 
                      scale(flow_sum)+
                      scale(precip_sum)+
                      scale(windsp_mean)+
                      (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod3)
residuals3 <- residuals(gpp_s_mod3)
hist(residuals3)
r.squaredGLMM(gpp_s_mod3)
vif(gpp_s_mod3)



gpp_s_mod3a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                     scale(light_mean)  + 
                     scale(flow_sum)+
                     scale(precip_sum)+
                     scale(windsp_mean)+
                     (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_s_mod3a)
residuals3 <- residuals(gpp_s_mod3a)
hist(residuals3a)
r.squaredGLMM(gpp_s_mod3a)
vif(gpp_s_mod3a)



gpp_stream_mod5a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                          scale(stream_temp)  + 
                          scale(flow_sum)+
                          scale(spc_mean)+
                           scale(precip_sum)+
                          scale(windsp_mean)+
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(gpp_stream_mod5a)
residuals3 <- residuals(gpp_stream_mod5a)
hist(residuals3)
r.squaredGLMM(gpp_stream_mod5a)
vif(gpp_stream_mod5a)



# Filter out missing values from the original data frame
GPP_df_stream_complete <- na.omit(GPP_df_stream_week)
# Extract fitted values for the complete observations
fitted_values <- predict(gpp_stream_mod4a, newdata = GPP_df_stream_complete)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(flow_sum = GPP_df_stream_complete$flow_sum, middle_GPP = fitted_values)

# Plot
G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = flow_sum, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = log(middle_GPP+1)), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_flowplot_a.png",sep=""), width=4,height=3,dpi=300)


# Extract fitted values for the complete observations
fitted_values <- predict(gpp_stream_mod4a, newdata = GPP_df_stream_complete)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(stream_temp = GPP_df_stream_complete$stream_temp, middle_GPP = fitted_values)

# Plot
G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = stream_temp, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = log(middle_GPP+1), x = stream_temp), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_streamtemp_a.png",sep=""), width=4,height=3,dpi=300)


# Extract fitted values for the complete observations
fitted_values <- predict(gpp_stream_mod4a, newdata = GPP_df_stream_complete)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(precip_sum = GPP_df_stream_complete$precip_sum, middle_GPP = fitted_values)

# Plot
G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = precip_sum, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = log(middle_GPP+1), x = precip_sum), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_precipevent_a.png",sep=""), width=4,height=3,dpi=300)


# Extract fitted values for the complete observations
fitted_values <- predict(gpp_stream_mod4a, newdata = GPP_df_stream_complete)
# Create a new data frame with spc_mean and fitted values
plot_data <- data.frame(windsp_mean = GPP_df_stream_complete$windsp_mean, middle_GPP = fitted_values)

# Plot
G_plot <- ggplot(data = GPP_df_stream_complete, aes(x = windsp_mean, y = log(middle_GPP+1), color = shore)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=position)) +
  geom_smooth(data = plot_data, aes(y = log(middle_GPP+1), x = windsp_mean), method = "lm", se = F, color = "grey50") +  # Linear trend for spc_mean and middle_ER from the model
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() 
G_plot

# ggsave(plot = G_plot, filename = paste("./SFS24_Analysis/figures/GPP_windsp_a.png",sep=""), width=4,height=3,dpi=300)


AIC(gpp_s_mod1, gpp_s_mod1a, gpp_s_mod2, gpp_s_mod2a, gpp_s_mod3, gpp_s_mod3a, gpp_stream_mod4a)

AIC(gpp_s_mod1a, gpp_s_mod2a, gpp_s_mod3a, gpp_stream_mod4a)

AIC(gpp_s_mod1b, gpp_s_mod2b)



ER_df2$week <- week(ER_df2$date)

## 
ER_df_week <- ER_df2 %>%
  group_by(site, shore, week, year, position) %>%
  summarise(middle_ER=mean(middle_ER, na.rm = TRUE),
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

ER_df_week$middle_ER<- c(ER_df_week$middle_ER*-1)
ER_df_week$middle_ER_lag <- lag(ER_df_week$middle_ER)


ER_stream_mod4a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                           scale(stream_temp)  + 
                           scale(flow_sum)+
                          # scale(spc_mean)+
                           scale(precip_sum)+
                           scale(windsp_mean)+
                           (1|site) + (1|year), data = GPP_df_stream_week)
summary(ER_stream_mod4a)
residuals3 <- residuals(ER_stream_mod4a)
hist(residuals3)
r.squaredGLMM(ER_stream_mod4a)
vif(ER_stream_mod4a)




ER_stream_mod2a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                          scale(light_mean)  + 
                          scale(flow_sum)+
                          #scale(spc_mean)+
                          scale(precip_sum)+
                          scale(windsp_mean)+
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(ER_stream_mod2a)
residuals3 <- residuals(ER_stream_mod2a)
hist(residuals3)
r.squaredGLMM(ER_stream_mod2a)
vif(ER_stream_mod2a)




ER_stream_mod3a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                          scale(tmean_C)  + 
                          scale(flow_sum)+
                          #scale(spc_mean)+
                          scale(precip_sum)+
                          scale(windsp_mean)+
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(ER_stream_mod3a)
residuals3 <- residuals(ER_stream_mod3a)
hist(residuals3)
r.squaredGLMM(ER_stream_mod3a)
vif(ER_stream_mod3a)


ER_stream_mod3a <- lmer(log(middle_GPP+1) ~ #log(lag_middle_GPP+1) + 
                          scale(tmean_C)  + 
                          scale(flow_sum)+
                          #scale(spc_mean)+
                          scale(precip_sum)+
                          scale(windsp_mean)+
                          (1|site) + (1|year), data = GPP_df_stream_week)
summary(ER_stream_mod3a)
residuals3 <- residuals(ER_stream_mod3a)
hist(residuals3)
r.squaredGLMM(ER_stream_mod3a)
vif(ER_stream_mod3a)

AIC(ER_stream_mod3a, ER_stream_mod2a, ER_stream_mod4a)


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

chem_datQ <- chem_dat %>%
  group_by(shore, location, date, substrate) %>%
  summarise(
    NO3_mgL_dl = mean(NO3_mgL_dl, na.rm=T),
    NH3_mgL_dl= mean(NH3_mgL_dl, na.rm=T),
    NH4_mgL_dl= mean(NH4_mgL_dl, na.rm=T),
    PO4_ugL_dl = mean(PO4_ugL_dl, na.rm=T),
    DOC_mgL_dl = mean(DOC_mgL_dl, na.rm=T))


summary(chem_datQ)

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
ggsave(plot = chem_grid, filename = paste("./SFS24_Analysis/figures/Lake_streamChem_plot.png",sep=""),width=5,height=5.05,dpi=300)



