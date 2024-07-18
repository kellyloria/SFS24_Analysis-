
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


## 
## create a new df for complete GPP obs 
GPP_df <- SFS_datQ1%>%
  dplyr::select( -middle_ER, -upper_ER, -lower_ER)

GPP_df<- GPP_df%>%
  drop_na(middle_GPP)

ER_df <- SFS_datQ1%>%
  dplyr::select(-middle_GPP, -upper_GPP, -lower_GPP)

ER_df<- ER_df%>%
  drop_na(middle_ER)


ER_df$ER <- c(ER_df$middle_ER *-1)
ER_df$ER_low <- c(ER_df$lower_ER *-1)
ER_df$ER_up <- c(ER_df$upper_ER *-1)
#GPP_df1<- na.omit(GPP_df)
summary(GPP_df)
summary(ER_df)

#===========================================
## GPP glms for stream influence
#============================================
# subset for shore BW and GB
GPP_df_stream <- GPP_df %>%
  filter(shore=="BW" |shore=="GB")

GPP_df_stream$middle_GPP_lag <- lag(GPP_df_stream$middle_GPP)
ER_df_stream <- ER_df %>%
  filter(shore=="BW" |shore=="GB")

ER_df_stream$middle_ER_lag <- lag(ER_df_stream$ER)

#============================================

GPP_df_BW <- GPP_df_stream%>%
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



## Issues with normality... 
GPP_df_GB <- GPP_df_stream%>%
  filter(shore=="GB") 
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

columns_to_scale <- c("logGPP", "logGPPlag", "loglake_tempC", "light_mean", "windsp_mean","tmean_C","log_streamflow")
# Scale the specific columns
GB_dat_scaled <- GPP_df_GB %>%
  mutate(across(all_of(columns_to_scale), scale))

GB_dat_scale <- as.data.frame(GB_dat_scaled)
# Display the structure of the scaled data frame
str(GB_dat_scaled)
names(GB_dat_scaled)

summary(GB_dat_scale)


##### ER 
hist(ER_df_stream$lake_tempC)
hist(ER_df_stream$middle_ER)
hist(ER_df_stream$par_int_3m)

ER_df_stream$logER <- log(ER_df_stream$ER+1)
hist(ER_df_stream$logER)

ER_df_stream$logERlag <- log(ER_df_stream$middle_ER_lag+1)
ER_df_stream$loglake_tempC <- log(ER_df_stream$lake_tempC)
hist(ER_df_stream$loglake_tempC)


hist(ER_df_stream$par_int_3m)
summary(ER_df_stream)

ER_df_Q <- ER_df_stream%>%
  select(site, shore, date, year, ER, middle_ER_lag, logER, logERlag, lake_tempC, loglake_tempC, Kd_fill, par_int_3m, tmean_C, light_mean, flow_mean, log_streamflow, log_ppmt, windsp_mean)

ER_df_BW <- ER_df_Q%>%
  filter(shore=="BW")

ER_df_GB <- ER_df_Q%>%
  filter(shore=="GB")

ER_df_23 <- ER_df_Q %>%
  filter(date> as.Date("2023-02-01"))
summary(ER_df_23)

ER_df_23<- drop_na(ER_df_23)


######


ER_df_GB <- ER_df_GB[-1,]

ER_df_BW <- ER_df_BW[-1,]





#####################
#######################

# https://rpubs.com/jebyrnes/brms_bayes_sem

library(piecewiseSEM) 
library(DiagrammeR)
library(brms)
library(gridExtra)
library(ggplot2)

names(GPP_df_BWQ)

###
BW_cor <- scale(GPP_df_BWQ[, c(6:16)])
chart.Correlation(BW_cor, histogram=TRUE, pch=19)

BW_dat <- GPP_df_BWQ[, c(6:16)]

# Load necessary library
library(dplyr)

# Assume BW_dat is your data frame
# Display the structure of the data frame
str(BW_dat)

# Columns to scale
columns_to_scale <- c("logGPP", "logGPPlag", "loglake_tempC", "light_mean", "windsp_mean","tmean_C","log_streamflow")
# Scale the specific columns
BW_dat_scaled <- GPP_df_BWQ %>%
  mutate(across(all_of(columns_to_scale), scale))

BW_dat_scale <- as.data.frame(BW_dat_scaled)
# Display the structure of the scaled data frame
str(BW_dat_scaled)
names(BW_dat_scaled)

flow_sub <- GPP_df_BWQ %>%
  filter(log_streamflow<0.102145)

hist(GPP_df_BWQ$log_streamflow)

# temp mod
plot(GPP_df_BWQ$loglake_tempC~GPP_df_BWQ$log_streamflow)
plot(flow_sub$lake_tempC~flow_sub$log_streamflow)

summary(lmer(loglake_tempC ~ log_streamflow + (1|site), data=BW_dat_scale))
summary(lmer(lake_tempC ~ log_streamflow + (1|site), data=flow_sub))


plot(GPP_df_BWQ$lake_tempC~GPP_df_BWQ$log_streamflow)


# temp mod
plot(GPP_df_BWQ$loglake_tempC~GPP_df_BWQ$windsp_mean)
summary(lmer(loglake_tempC ~ windsp_mean + (1|site), data=BW_dat_scale))

# plot(GPP_df_BWQ$loglake_tempC~GPP_df_BWQ$log_ppmt)
# summary(lmer(loglake_tempC ~ log_ppmt + (1|site), data=GPP_df_BWQ))

plot(GPP_df_BWQ$loglake_tempC~GPP_df_BWQ$tmean_C)
summary(lmer(loglake_tempC ~ tmean_C + (1|site), data=BW_dat_scale))



plot(GPP_df_BWQ$light_mean~GPP_df_BWQ$log_streamflow)
summary(lmer(light_mean ~ log_streamflow + (1|site), data=BW_dat_scale))


plot(flow_sub$light_mean~flow_sub$log_streamflow)
summary(lmer(light_mean ~ log_streamflow + (1|site), data=flow_sub))


# plot(GPP_df_BWQ$par_int_3m~GPP_df_BWQ$windsp_mean)
# summary(lmer(par_int_3m ~ windsp_mean + (1|site), data=GPP_df_BWQ))
# 
# plot(GPP_df_BWQ$par_int_3m~GPP_df_BWQ$log_ppmt)
# summary(lmer(par_int_3m ~ windsp_mean + (1|site), data=GPP_df_BWQ))





#########
names(GPP_df_BWQ)

k_fit_psem <- psem(
  lmer(logGPP ~ logGPPlag + light_mean + lake_tempC + log_streamflow +(1|site), data=flow_sub),
  lmer(lake_tempC ~ log_streamflow + tmean_C +  windsp_mean + (1|site), data=flow_sub)
)

plot(k_fit_psem)
summary(k_fit_psem)
fisherC(k_fit_psem)

k_fit_psem <- psem(
  lmer(logGPP ~ logGPPlag + par_int_3m + lake_tempC + log_streamflow +(1|site), data=flow_sub),
  lmer(lake_tempC ~ log_streamflow+ (1|site), data=flow_sub),
  lmer(par_int_3m ~ log_streamflow+ (1|site), data=flow_sub)
)

plot(k_fit_psem)
summary(k_fit_psem)
fisherC(k_fit_psem)



# start bayes 

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(gridExtra)
library(ggplot2)
library(grid)
library(brms)


GPP_mod <- bf(logGPP ~ logGPPlag + scale(par_int_3m) + scale(lake_tempC) + scale(log_streamflow) +(1|site)) #rich_mod
temp_mod <- bf(scale(lake_tempC) ~ scale(log_streamflow) + (1|site)) #cover_mod
light_mod <- bf(scale(par_int_3m) ~ scale(log_streamflow) + (1|site)) #cover_mod


k_fit_brms <- brm(GPP_mod +
                    temp_mod + 
                    light_mod +
                    set_rescor(FALSE), 
                  data=flow_sub,
                  iter=9000, warmup = 4500,
                  cores=4, chains = 3)

plot(k_fit_brms)
summary(k_fit_brms)



# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms)
posterior_summaries <- posterior_summary(k_fit_brms)

#  Extract the relevant estimates and CIs
estimates <- posterior_summaries[, "Estimate"]
ci_lower <- posterior_summaries[, "Q2.5"]
ci_upper <- posterior_summaries[, "Q97.5"]


# Create the DAG plot with DiagrammeR
dag_plot <- grViz("
  digraph DAG {
    graph [layout = dot, rankdir = TB]
    
    node [shape = ellipse]
    flow_mean -> lake_tempC [label = 'mean = -3.02, CI = [-4.68, -1.42]']
    flow_mean -> par_int_3m [label = 'mean = 761.45, CI = [643.61, 878.27]']
    flow_mean -> logGPP [label = 'mean = 0.01, CI = [0.00, 0.01]']
    lake_tempC -> logGPP [label = 'mean = 0.01, CI = [0.00, 0.01]']
    logGPPlag -> logGPP [label = 'mean = 0.82, CI = [0.78, 0.85]']
    par_int_3m -> logGPP [label = 'mean = -0.00, CI = [-0.00, -0.00]']
    
    logGPP [label = 'logGPP']
    lake_tempC [label = 'lake_tempC']
    logGPPlag [label = 'logGPPlag']
    par_int_3m [label = 'par_int_3m']
    flow_mean [label = 'flow_mean']
  }
")

# Convert DiagrammeR plot to an SVG, then to a grob
dag_svg <- export_svg(dag_plot)
dag_grob <- grid::rasterGrob(rsvg::rsvg(charToRaw(dag_svg)), interpolate = TRUE)

# Extract the relevant posteriors
posteriors <- posterior_samples(k_fit_brms)
post_logGPP_logGPPlag <- posteriors$b_logGPP_logGPPlag
post_logGPP_par_int_3m <- posteriors$b_logGPP_scalepar_int_3m
post_logGPP_lake_tempC <- posteriors$b_logGPP_scalelake_tempC
post_laketempC_flow_mean <- posteriors$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors$b_scaleparint3m_scalelog_streamflow
post_logGPP_flow_mean <- posteriors$b_logGPP_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#3283a8", linewidth = 1,  fill = "#666666ff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) +
    theme(plot.title = element_text(size = 10))
}




# Create histograms with estimates and CIs
hist_logGPP_logGPPlag <- create_histogram(post_logGPP_logGPPlag, estimates["b_logGPP_logGPPlag"], ci_lower["b_logGPP_logGPPlag"], ci_upper["b_logGPP_logGPPlag"], "logGPPlag -> logGPP")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logGPP_flow_mean <- create_histogram(post_logGPP_flow_mean, estimates["b_logGPP_scalelog_streamflow"], ci_lower["b_logGPP_scalelog_streamflow"], ci_upper["b_logGPP_scalelog_streamflow"], "flow_mean -> logGPP")


hist_logGPP_par_int_3m <- create_histogram(post_logGPP_par_int_3m, estimates["b_logGPP_scalepar_int_3m"], ci_lower["b_logGPP_scalepar_int_3m"], ci_upper["b_logGPP_scalepar_int_3m"], "par_int_3m -> logGPP")
hist_logGPP_lake_tempC <- create_histogram(post_logGPP_lake_tempC, estimates["b_logGPP_scalelake_tempC"], ci_lower["b_logGPP_scalelake_tempC"], ci_upper["b_logGPP_scalelake_tempC"], "lake_tempC -> logGPP")

# Arrange the DAG and histograms in a grid
sem_post_plt <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logGPP_logGPPlag),
  ggplotGrob(hist_logGPP_par_int_3m),
  ggplotGrob(hist_logGPP_lake_tempC),
  ggplotGrob(hist_logGPP_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt, filename = paste("./figures/NS24_BW_SEM_post_est3.png",sep=""),height=3,width=5.75,dpi=300)

############### GB
##############


#########
names(GPP_df_BWQ)

k_fit_psem <- psem(
  lmer(logGPP ~ logGPPlag + light_mean + lake_tempC + log_streamflow +(1|site), data=BW_dat_scale),
  lmer(lake_tempC ~ log_streamflow + tmean_C +  windsp_mean + (1|site), data=BW_dat_scale)
)

plot(k_fit_psem)
summary(k_fit_psem)
fisherC(k_fit_psem)

k_fit_psem <- psem(
  lmer(logGPP ~ logGPPlag + par_int_3m + lake_tempC + log_streamflow +(1|site), data=BW_dat_scale),
  lmer(lake_tempC ~ log_streamflow+ (1|site), data=BW_dat_scale),
  lmer(par_int_3m ~ log_streamflow+ (1|site), data=BW_dat_scale)
)

plot(k_fit_psem)
summary(k_fit_psem)
fisherC(k_fit_psem)



# start bayes 

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(gridExtra)
library(ggplot2)
library(grid)
library(brms)

GB_dat_flow <- GPP_df_GB%>%
  filter(log_streamflow<0.0029609)

plot(GB_dat_flow$log_streamflow~ GB_dat_flow$logGPPlag)


GPP_mod <- bf(logGPP ~ logGPPlag + scale(par_int_3m) + scale(lake_tempC) + scale(log_streamflow) +(1|site)) #rich_mod
temp_mod <- bf(scale(lake_tempC) ~ scale(log_streamflow) + (1|site)) #cover_mod
light_mod <- bf(scale(par_int_3m) ~ scale(log_streamflow) + (1|site)) #cover_mod


k_fit_brms_gb <- brm(GPP_mod +
                    temp_mod + 
                    light_mod +
                    set_rescor(FALSE), 
                  data=GB_dat_flow,
                  iter=9000, warmup = 4500,
                  cores=4, chains = 3)

plot(k_fit_brms_gb)
summary(k_fit_brms_gb)



# Extract posterior summaries
summary_k_fit_brms_gb <- summary(k_fit_brms_gb)
posterior_summaries_gb <- posterior_summary(k_fit_brms_gb)

#  Extract the relevant estimates and CIs
estimates_gb <- posterior_summaries_gb[, "Estimate"]
ci_lower_gb <- posterior_summaries_gb[, "Q2.5"]
ci_upper_gb <- posterior_summaries_gb[, "Q97.5"]


# Create the DAG plot with DiagrammeR
dag_plot <- grViz("
  digraph DAG {
    graph [layout = dot, rankdir = TB]
    
    node [shape = ellipse]
    flow_mean -> lake_tempC [label = 'mean = -3.02, CI = [-4.68, -1.42]']
    flow_mean -> par_int_3m [label = 'mean = 761.45, CI = [643.61, 878.27]']
    flow_mean -> logGPP [label = 'mean = 0.01, CI = [0.00, 0.01]']
    lake_tempC -> logGPP [label = 'mean = 0.01, CI = [0.00, 0.01]']
    logGPPlag -> logGPP [label = 'mean = 0.82, CI = [0.78, 0.85]']
    par_int_3m -> logGPP [label = 'mean = -0.00, CI = [-0.00, -0.00]']
    
    logGPP [label = 'logGPP']
    lake_tempC [label = 'lake_tempC']
    logGPPlag [label = 'logGPPlag']
    par_int_3m [label = 'par_int_3m']
    flow_mean [label = 'flow_mean']
  }
")

# Convert DiagrammeR plot to an SVG, then to a grob
dag_svg <- export_svg(dag_plot)
dag_grob <- grid::rasterGrob(rsvg::rsvg(charToRaw(dag_svg)), interpolate = TRUE)

# Extract the relevant posteriors
posteriors_gb <- posterior_samples(k_fit_brms_gb)
post_logGPP_logGPPlag <- posteriors_gb$b_logGPP_logGPPlag
post_logGPP_par_int_3m <- posteriors_gb$b_logGPP_scalepar_int_3m
post_logGPP_lake_tempC <- posteriors_gb$b_logGPP_scalelake_tempC
post_laketempC_flow_mean <- posteriors_gb$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors_gb$b_scaleparint3m_scalelog_streamflow
post_logGPP_flow_mean <- posteriors_gb$b_logGPP_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#a67c18", linewidth = 1,  fill = "#528cdeff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) +
    theme(plot.title = element_text(size = 10))
}




# Create histograms with estimates and CIs
hist_logGPP_logGPPlag <- create_histogram(post_logGPP_logGPPlag, estimates_gb["b_logGPP_logGPPlag"], ci_lower_gb["b_logGPP_logGPPlag"], ci_upper_gb["b_logGPP_logGPPlag"], "logGPPlag -> logGPP")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates_gb["b_scalelaketempC_scalelog_streamflow"], ci_lower_gb["b_scalelaketempC_scalelog_streamflow"], ci_upper_gb["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates_gb["b_scaleparint3m_scalelog_streamflow"], ci_lower_gb["b_scaleparint3m_scalelog_streamflow"], ci_upper_gb["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logGPP_flow_mean <- create_histogram(post_logGPP_flow_mean, estimates_gb["b_logGPP_scalelog_streamflow"], ci_lower_gb["b_logGPP_scalelog_streamflow"], ci_upper_gb["b_logGPP_scalelog_streamflow"], "flow_mean -> logGPP")

hist_logGPP_par_int_3m <- create_histogram(post_logGPP_par_int_3m, estimates_gb["b_logGPP_scalepar_int_3m"], ci_lower_gb["b_logGPP_scalepar_int_3m"], ci_upper_gb["b_logGPP_scalepar_int_3m"], "par_int_3m -> logGPP")
hist_logGPP_lake_tempC <- create_histogram(post_logGPP_lake_tempC, estimates_gb["b_logGPP_scalelake_tempC"], ci_lower_gb["b_logGPP_scalelake_tempC"], ci_upper_gb["b_logGPP_scalelake_tempC"], "lake_tempC -> logGPP")

# Arrange the DAG and histograms in a grid
sem_post_plt_gb <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logGPP_logGPPlag),
  ggplotGrob(hist_logGPP_par_int_3m),
  ggplotGrob(hist_logGPP_lake_tempC),
  ggplotGrob(hist_logGPP_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt_gb, filename = paste("./figures/NS24_GB_SEM_post_est3.png",sep=""),height=3,width=5.75,dpi=300)




# Function to create histogram plots with estimates and CIs for multiple datasets
create_facet_histogram <- function(data_list, estimates, ci_lowers, ci_uppers, titles, colors) {
  combined_data <- data.frame()
  
  for (i in seq_along(data_list)) {
    temp_data <- data.frame(
      value = data_list[[i]], 
      dataset = titles[i], 
      estimate = estimates[i], 
      ci_lower = ci_lowers[i], 
      ci_upper = ci_uppers[i],
      group = paste0("Estimate: ", estimates[i])
    )
    combined_data <- rbind(combined_data, temp_data)
  }
  
  ggplot(combined_data, aes(x = value, fill = dataset, color = dataset)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    geom_vline(aes(xintercept = estimate, color = dataset), size = 0.75) +
    geom_vline(aes(xintercept = ci_lower, color = dataset), linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper, color = dataset), linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle("Faceted Histogram by Estimates") +
    theme(plot.title = element_text(size = 10)) +
    facet_wrap(~group, scales = "free")
}


# Data lists
data_list <- list(post_logGPP_logGPPlag, post_logGPP_logGPPlag, post_laketempC_flow_mean, post_laketempC_flow_mean, post_parint3m_flow_mean, post_parint3m_flow_mean, post_logGPP_flow_mean, post_logGPP_flow_mean, post_logGPP_par_int_3m, post_logGPP_par_int_3m, post_logGPP_lake_tempC, post_logGPP_lake_tempC)
estimates <- c(estimates["b_logGPP_logGPPlag"], estimates_gb["b_logGPP_logGPPlag"], estimates["b_scalelaketempC_scalelog_streamflow"], estimates_gb["b_scalelaketempC_scalelog_streamflow"], estimates["b_scaleparint3m_scalelog_streamflow"], estimates_gb["b_scaleparint3m_scalelog_streamflow"], estimates["b_logGPP_scalelog_streamflow"], estimates_gb["b_logGPP_scalelog_streamflow"], estimates["b_logGPP_scalepar_int_3m"], estimates_gb["b_logGPP_scalepar_int_3m"], estimates["b_logGPP_scalelake_tempC"], estimates_gb["b_logGPP_scalelake_tempC"])
ci_lowers <- c(ci_lower["b_logGPP_logGPPlag"], ci_lower_gb["b_logGPP_logGPPlag"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_lower_gb["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_lower_gb["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_logGPP_scalelog_streamflow"], ci_lower_gb["b_logGPP_scalelog_streamflow"], ci_lower["b_logGPP_scalepar_int_3m"], ci_lower_gb["b_logGPP_scalepar_int_3m"], ci_lower["b_logGPP_scalelake_tempC"], ci_lower_gb["b_logGPP_scalelake_tempC"])
ci_uppers <- c(ci_upper["b_logGPP_logGPPlag"], ci_upper_gb["b_logGPP_logGPPlag"], ci_upper["b_scalelaketempC_scalelog_streamflow"], ci_upper_gb["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], ci_upper_gb["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_logGPP_scalelog_streamflow"], ci_upper_gb["b_logGPP_scalelog_streamflow"], ci_upper["b_logGPP_scalepar_int_3m"], ci_upper_gb["b_logGPP_scalepar_int_3m"], ci_upper["b_logGPP_scalelake_tempC"], ci_upper_gb["b_logGPP_scalelake_tempC"])
titles <- c("logGPPlag -> logGPP", "logGPPlag -> logGPP", "flow_mean -> lake_tempC", "flow_mean -> lake_tempC", "flow_mean -> par_int_3m", "flow_mean -> par_int_3m", "flow_mean -> logGPP", "flow_mean -> logGPP", "par_int_3m -> logGPP", "par_int_3m -> logGPP", "lake_tempC -> logGPP", "lake_tempC -> logGPP")
colors <- c("#3283a8", "#a67c18", "#3283a8", "#a67c18", "#3283a8", "#a67c18", "#3283a8", "#a67c18", "#3283a8", "#a67c18", "#3283a8", "#a67c18")

# Create faceted histogram
faceted_histogram <- create_facet_histogram(data_list, estimates, ci_lowers, ci_uppers, titles, colors)

# Display the plot
print(faceted_histogram)



########
flow_sub_BW <- ER_df_BW %>%
  filter(log_streamflow<0.102145)


names(ER_df_BW)

ER_mod <- bf(logER ~ logERlag + scale(par_int_3m) + scale(lake_tempC) + scale(log_streamflow) +(1|site)) #rich_mod
temp_mod <- bf(scale(lake_tempC) ~ scale(log_streamflow) + (1|site)) #cover_mod
light_mod <- bf(scale(par_int_3m) ~ scale(log_streamflow) + (1|site)) #cover_mod


k_fit_brms_BWER <- brm(ER_mod +
                    temp_mod + 
                    light_mod +
                    set_rescor(FALSE), 
                  data=flow_sub_BW,
                  iter=9000, warmup = 4500,
                  cores=4, chains = 3)

plot(k_fit_brms_BWER)
summary(k_fit_brms_BWER)



# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms_BWER)
posterior_summaries_BWER <- posterior_summary(k_fit_brms_BWER)

#  Extract the relevant estimates and CIs
estimates <- posterior_summaries_BWER[, "Estimate"]
ci_lower <- posterior_summaries_BWER[, "Q2.5"]
ci_upper <- posterior_summaries_BWER[, "Q97.5"]

# Extract the relevant posteriors
posteriors_BWER <- posterior_samples(k_fit_brms_BWER)
post_logER_logERlag <- posteriors_BWER$b_logER_logERlag
post_logER_par_int_3m <- posteriors_BWER$b_logER_scalepar_int_3m
post_logER_lake_tempC <- posteriors_BWER$b_logER_scalelake_tempC
post_laketempC_flow_mean <- posteriors_BWER$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors_BWER$b_scaleparint3m_scalelog_streamflow
post_logER_flow_mean <- posteriors_BWER$b_logER_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#3283a8", linewidth = 1,  fill = "#528cdeff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) + ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



###
library(ggplot2)
library(scales)  # For number_format function

# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  # Round the numbers to have only two digits after the decimal point
  estimate_formatted <- round(estimate, 2)
  ci_lower_formatted <- round(ci_lower, 2)
  ci_upper_formatted <- round(ci_upper, 2)
  
  ggplot(data.frame(value = data), aes(x = value)) +
    # geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#3283a8", linewidth = 1, fill = "#528cdeff") +  # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff
    geom_vline(aes(xintercept = estimate_formatted), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower_formatted), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper_formatted), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    scale_x_continuous(labels = number_format(accuracy = 0.01)) +  
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +  
    ggtitle(title) +  ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



## Create histograms with estimates and CIs
# hist_logER_logERlag <- create_histogram(post_logER_logERlag, estimates["b_logER_logERlag"], ci_lower["b_logER_logERlag"], ci_upper["b_logER_logERlag"], "logERlag -> logER")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logER_flow_mean <- create_histogram(post_logER_flow_mean, estimates["b_logER_scalelog_streamflow"], ci_lower["b_logER_scalelog_streamflow"], ci_upper["b_logER_scalelog_streamflow"], "flow_mean -> logER")


hist_logER_par_int_3m <- create_histogram(post_logER_par_int_3m, estimates["b_logER_scalepar_int_3m"], ci_lower["b_logER_scalepar_int_3m"], ci_upper["b_logER_scalepar_int_3m"], "par_int_3m -> logER")
hist_logER_lake_tempC <- create_histogram(post_logER_lake_tempC, estimates["b_logER_scalelake_tempC"], ci_lower["b_logER_scalelake_tempC"], ci_upper["b_logER_scalelake_tempC"], "lake_tempC -> logER")

# Arrange the DAG and histograms in a grid
sem_post_plt_BWER <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logER_logERlag),
  ggplotGrob(hist_logER_par_int_3m),
  ggplotGrob(hist_logER_lake_tempC),
  ggplotGrob(hist_logER_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt_BWER, filename = paste("./figures/NS24_BW_ER_SEM_post_est3.png",sep=""),height=3.25,width=6.5,dpi=300)


########
######


########

names(ER_df_BW)

GB_dat_flow_E <- ER_df_GB%>%
  filter(log_streamflow<0.0029609)


ER_mod <- bf(logER ~ logERlag + scale(par_int_3m) + scale(lake_tempC) + scale(log_streamflow) +(1|site)) #rich_mod
temp_mod <- bf(scale(lake_tempC) ~ scale(log_streamflow) + (1|site)) #cover_mod
light_mod <- bf(scale(par_int_3m) ~ scale(log_streamflow) + (1|site)) #cover_mod


k_fit_brms_BWER <- brm(ER_mod +
                         temp_mod + 
                         light_mod +
                         set_rescor(FALSE), 
                       data=GB_dat_flow_E,
                       iter=9000, warmup = 4500,
                       cores=4, chains = 3)

plot(k_fit_brms_BWER)
summary(k_fit_brms_BWER)



# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms_BWER)
posterior_summaries_BWER <- posterior_summary(k_fit_brms_BWER)

#  Extract the relevant estimates and CIs
estimates <- posterior_summaries_BWER[, "Estimate"]
ci_lower <- posterior_summaries_BWER[, "Q2.5"]
ci_upper <- posterior_summaries_BWER[, "Q97.5"]

# Extract the relevant posteriors
posteriors_BWER <- posterior_samples(k_fit_brms_BWER)
post_logER_logERlag <- posteriors_BWER$b_logER_logERlag
post_logER_par_int_3m <- posteriors_BWER$b_logER_scalepar_int_3m
post_logER_lake_tempC <- posteriors_BWER$b_logER_scalelake_tempC
post_laketempC_flow_mean <- posteriors_BWER$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors_BWER$b_scaleparint3m_scalelog_streamflow
post_logER_flow_mean <- posteriors_BWER$b_logER_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#3283a8", linewidth = 1,  fill = "#666666ff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) + ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



###
library(ggplot2)
library(scales)  # For number_format function

# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  # Round the numbers to have only two digits after the decimal point
  estimate_formatted <- round(estimate, 2)
  ci_lower_formatted <- round(ci_lower, 2)
  ci_upper_formatted <- round(ci_upper, 2)
  
  ggplot(data.frame(value = data), aes(x = value)) +
    # geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#3283a8", linewidth = 1, fill = "#0b5394ff") +  # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff
    geom_vline(aes(xintercept = estimate_formatted), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower_formatted), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper_formatted), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    scale_x_continuous(labels = number_format(accuracy = 0.01)) +  
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +  
    ggtitle(title) +  ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



## Create histograms with estimates and CIs
# hist_logER_logERlag <- create_histogram(post_logER_logERlag, estimates["b_logER_logERlag"], ci_lower["b_logER_logERlag"], ci_upper["b_logER_logERlag"], "logERlag -> logER")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logER_flow_mean <- create_histogram(post_logER_flow_mean, estimates["b_logER_scalelog_streamflow"], ci_lower["b_logER_scalelog_streamflow"], ci_upper["b_logER_scalelog_streamflow"], "flow_mean -> logER")


hist_logER_par_int_3m <- create_histogram(post_logER_par_int_3m, estimates["b_logER_scalepar_int_3m"], ci_lower["b_logER_scalepar_int_3m"], ci_upper["b_logER_scalepar_int_3m"], "par_int_3m -> logER")
hist_logER_lake_tempC <- create_histogram(post_logER_lake_tempC, estimates["b_logER_scalelake_tempC"], ci_lower["b_logER_scalelake_tempC"], ci_upper["b_logER_scalelake_tempC"], "lake_tempC -> logER")

# Arrange the DAG and histograms in a grid
sem_post_plt_BWER <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logER_logERlag),
  ggplotGrob(hist_logER_par_int_3m),
  ggplotGrob(hist_logER_lake_tempC),
  ggplotGrob(hist_logER_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt_BWER, filename = paste("./figures/NS24_BW_ER_SEM_post_est2.png",sep=""),height=3.25,width=6.5,dpi=300)


########
##### GB


k_fit_brms_GBER <- brm(ER_mod +
                         temp_mod + 
                         light_mod +
                         set_rescor(FALSE), 
                       data=GB_dat_flow_E,
                       iter=9000, warmup = 4500,
                       cores=4, chains = 3)

plot(k_fit_brms_GBER)
summary(k_fit_brms_GBER)



# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms_GBER)
posterior_summaries_GBER <- posterior_summary(k_fit_brms_GBER)

#  Extract the relevant estimates and CIs
estimates <- posterior_summaries_GBER[, "Estimate"]
ci_lower <- posterior_summaries_GBER[, "Q2.5"]
ci_upper <- posterior_summaries_GBER[, "Q97.5"]

# Extract the relevant posteriors
posteriors_GBER <- posterior_samples(k_fit_brms_GBER)
post_logER_logERlag <- posteriors_GBER$b_logER_logERlag
post_logER_par_int_3m <- posteriors_GBER$b_logER_scalepar_int_3m
post_logER_lake_tempC <- posteriors_GBER$b_logER_scalelake_tempC
post_laketempC_flow_mean <- posteriors_GBER$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors_GBER$b_scaleparint3m_scalelog_streamflow
post_logER_flow_mean <- posteriors_GBER$b_logER_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#a67c18", linewidth = 1,  fill = "#666666ff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) + ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



## Create histograms with estimates and CIs
# hist_logER_logERlag <- create_histogram(post_logER_logERlag, estimates["b_logER_logERlag"], ci_lower["b_logER_logERlag"], ci_upper["b_logER_logERlag"], "logERlag -> logER")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logER_flow_mean <- create_histogram(post_logER_flow_mean, estimates["b_logER_scalelog_streamflow"], ci_lower["b_logER_scalelog_streamflow"], ci_upper["b_logER_scalelog_streamflow"], "flow_mean -> logER")


hist_logER_par_int_3m <- create_histogram(post_logER_par_int_3m, estimates["b_logER_scalepar_int_3m"], ci_lower["b_logER_scalepar_int_3m"], ci_upper["b_logER_scalepar_int_3m"], "par_int_3m -> logER")
hist_logER_lake_tempC <- create_histogram(post_logER_lake_tempC, estimates["b_logER_scalelake_tempC"], ci_lower["b_logER_scalelake_tempC"], ci_upper["b_logER_scalelake_tempC"], "lake_tempC -> logER")

# Arrange the DAG and histograms in a grid
sem_post_plt_GBER <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logER_logERlag),
  ggplotGrob(hist_logER_par_int_3m),
  ggplotGrob(hist_logER_lake_tempC),
  ggplotGrob(hist_logER_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt_GBER, filename = paste("./figures/NS24_GB_ER_SEM_post_est3.png",sep=""),height=3.25,width=6.5,dpi=300)


########## shore analysis 
names(ER_df_23)

ER_mod <- bf(logER ~ logERlag + (light_mean) + (lake_tempC) + (windsp_mean) + (log_ppmt) +(1|shore/site)) 

k_fit_brms_ER23 <- brm(ER_mod, 
                       data=ER_df_23,
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)

plot(k_fit_brms_ER23)
summary(k_fit_brms_ER23)

# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms_ER23)
posterior_summaries_ER <- posterior_summary(k_fit_brms_ER23)



ER_moda <- bf(logER ~ logERlag + (par_int_3m) + (lake_tempC) + (windsp_mean) +(1|shore/site)) 

k_fit_brms_ER23a <- brm(ER_moda, 
                       data=ER_df_23,
                       iter=10000, warmup = 5000,
                       control = list(adapt_delta = 0.95, max_treedepth = 15), 
                       cores=4, chains = 3)

plot(k_fit_brms_ER23a)
summary(k_fit_brms_ER23a)


ER_modb <- bf(logER ~ logERlag + (par_int_3m) + (lake_tempC) + (windsp_mean) +(1|shore)) 

k_fit_brms_ER23b <- brm(ER_modb, 
                        data=ER_df_23,
                        iter=10000, warmup = 5000,
                        control = list(adapt_delta = 0.95, max_treedepth = 15), 
                        cores=4, chains = 3)

plot(k_fit_brms_ER23b)
summary(k_fit_brms_ER23b)


# Extract posterior summaries
summary_k_fit_brms <- summary(k_fit_brms_ER23)
posterior_summaries_GB23 <- posterior_summary(k_fit_brms_ER23)

#  Extract the relevant estimates and CIs
estimates <- posterior_summaries_GB23[, "Estimate"]
ci_lower <- posterior_summaries_GB23[, "Q2.5"]
ci_upper <- posterior_summaries_GB23[, "Q97.5"]

# Extract the relevant posteriors
posteriors_ER23 <- posterior_samples(k_fit_brms_ER23)
post_logER_logERlag <- posteriors_ER23$b_logER_logERlag
post_logER_par_int_3m <- posteriors_ER23$b_logER_scalepar_int_3m
post_logER_lake_tempC <- posteriors_ER23$b_logER_scalelake_tempC
post_laketempC_flow_mean <- posteriors_ER23$b_scalelaketempC_scalelog_streamflow
post_parint3m_flow_mean <- posteriors_ER23$b_scaleparint3m_scalelog_streamflow
post_logER_flow_mean <- posteriors_ER23$b_logER_scalelog_streamflow



# Function to create histogram plots with estimates and CIs
create_histogram <- function(data, estimate, ci_lower, ci_upper, title) {
  ggplot(data.frame(value = data), aes(x = value)) +
    #geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
    geom_density(color = "#a67c18", linewidth = 1,  fill = "#0b5394ff") + # lake= #0b5394ff , stream= #528cdeff, gpp= #666666ff"
    geom_vline(aes(xintercept = estimate), color = "black", size = 0.75) +
    geom_vline(aes(xintercept = ci_lower), color = "black", linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = ci_upper), color = "black", linetype = "dotted", size = 0.75) +
    theme_minimal() +
    ggtitle(title) + ylab(NULL)+
    theme(plot.title = element_text(size = 10))
}



## Create histograms with estimates and CIs
# hist_logER_logERlag <- create_histogram(post_logER_logERlag, estimates["b_logER_logERlag"], ci_lower["b_logER_logERlag"], ci_upper["b_logER_logERlag"], "logERlag -> logER")

hist_laketempC_flow_mean <- create_histogram(post_laketempC_flow_mean, estimates["b_scalelaketempC_scalelog_streamflow"], ci_lower["b_scalelaketempC_scalelog_streamflow"], ci_upper["b_scalelaketempC_scalelog_streamflow"], "flow_mean -> lake_tempC")
hist_parint3m_flow_mean <- create_histogram(post_parint3m_flow_mean, estimates["b_scaleparint3m_scalelog_streamflow"], ci_lower["b_scaleparint3m_scalelog_streamflow"], ci_upper["b_scaleparint3m_scalelog_streamflow"], "flow_mean -> par_int_3m")
hist_logER_flow_mean <- create_histogram(post_logER_flow_mean, estimates["b_logER_scalelog_streamflow"], ci_lower["b_logER_scalelog_streamflow"], ci_upper["b_logER_scalelog_streamflow"], "flow_mean -> logER")


hist_logER_par_int_3m <- create_histogram(post_logER_par_int_3m, estimates["b_logER_scalepar_int_3m"], ci_lower["b_logER_scalepar_int_3m"], ci_upper["b_logER_scalepar_int_3m"], "par_int_3m -> logER")
hist_logER_lake_tempC <- create_histogram(post_logER_lake_tempC, estimates["b_logER_scalelake_tempC"], ci_lower["b_logER_scalelake_tempC"], ci_upper["b_logER_scalelake_tempC"], "lake_tempC -> logER")

# Arrange the DAG and histograms in a grid
sem_post_plt_GBER <-grid.arrange(
  #dag_grob,
  ggplotGrob(hist_logER_logERlag),
  ggplotGrob(hist_logER_par_int_3m),
  ggplotGrob(hist_logER_lake_tempC),
  ggplotGrob(hist_logER_flow_mean),
  ggplotGrob(hist_laketempC_flow_mean),
  ggplotGrob(hist_parint3m_flow_mean),
  nrow = 2,
  ncol = 3
)

# ggsave(plot = sem_post_plt_GBER, filename = paste("./figures/NS24_GB_ER_SEM_post_est2.png",sep=""),height=3.25,width=6.5,dpi=300)

#####

GPP_df$middle_GPP_lag <- lag(GPP_df$middle_GPP)


ER_df_stream$middle_ER_lag <- lag(ER_df_stream$ER)

#============================================

## Issues with normality... 
hist(GPP_df_BW$lake_tempC)
hist(GPP_df_BW$middle_GPP)
hist(GPP_df_BW$par_int_3m)

GPP_df$logGPP <- log(GPP_df$middle_GPP+1)
hist(GPP_df$logGPP)

GPP_df$logGPPlag <- log(GPP_df$middle_GPP_lag+1)

GPP_df_23 <- GPP_df%>%
  select(site, shore, date, middle_GPP, middle_GPP_lag, logGPP, logGPPlag, 
         lake_tempC, Kd_fill, par_int_3m, tmean_C, light_mean, log_ppmt, 
         windsp_mean, precip_bi, real_NS_depth) %>%
  filter(date > as.Date("2023-02-01"))

summary(GPP_df_23)

GPP_df_BWQ<- drop_na(GPP_df_23)

# # Load necessary libraries
# library(brms)
# library(ggplot2)
# library(dplyr)
# library(reshape2)
# 
# # Extract posterior samples
# posterior_samples <- posterior_samples(k_fit_brms_GPP23a)
# 
# # Assuming the posterior samples are not directly linked to shore levels,
# # we need to replicate the levels across the posterior samples length.
# 
# # Create a sequence of shore levels to match the posterior samples length
# shore_levels <- rep(c("BW", "GB", "SH", "SS"), length.out = nrow(posterior_samples))
# 
# # Prepare the posterior samples data frame with shore levels
# posterior_samples_df <- posterior_samples %>%
#   select(b_Intercept, b_logGPPlag, b_par_int_3m, b_lake_tempC, b_precip_bi) %>%
#   mutate(shore = shore_levels)
# 
# # Melt the data for easier plotting with ggplot2
# posterior_samples_melted <- melt(posterior_samples_df, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')
# 
# # Calculate overall mean and credible intervals for each coefficient
# ci_summary_overall <- posterior_samples_melted %>%
#   group_by(Coefficient) %>%
#   summarize(
#     mean = mean(Estimate),
#     l_95 = quantile(Estimate, 0.025),
#     u_95 = quantile(Estimate, 0.975)
#   )
# 
# # Plot histograms with credible intervals and mean estimates
# ggplot(posterior_samples_melted, aes(x = Estimate, fill = shore, color = shore)) +
#   facet_wrap(~ Coefficient, scales = 'free') +
#   geom_density(alpha = 0.25, linewidth = 0.05) +
#   scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
#   scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
#   theme_minimal() +
#   labs(title = 'Posterior Distributions of Coefficients by Site',
#        x = 'Estimate',
#        y = 'Density') +
#   geom_vline(data = ci_summary_overall, aes(xintercept = mean), linetype = "solid", size = 0.75, color = "black") +
#   geom_vline(data = ci_summary_overall, aes(xintercept = l_95), linetype = "dashed", size = 0.5, color = "black") +
#   geom_vline(data = ci_summary_overall, aes(xintercept = u_95), linetype = "dashed", size = 0.5, color = "black")
# 

tempm1 <- lmer(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + precip_bi +(1|shore), data=GPP_df_23) 
summary(tempm1)
vif(tempm1)


tempm2 <- lmer(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + precip_bi +(1|shore/site), data=GPP_df_23) 
summary(tempm2)
vif(tempm2)

tempm3 <- lmer(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + log_ppmt +(1|shore), data=GPP_df_23) 
summary(tempm3)
vif(tempm3)

tempm4 <- lmer(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + log_ppmt +(1|shore/site), data=GPP_df_23) 
summary(tempm4)

tempm5 <- lmer(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + log_ppmt + real_NS_depth +(1|shore/site), data=GPP_df_23) 
summary(tempm5)
vif(tempm5)


GPP_modc <- bf(logGPP ~ logGPPlag + (par_int_3m) + (lake_tempC) + log_ppmt + real_NS_depth +(1|shore/site)) 

k_fit_brms_GPP23c <- brm(GPP_modc, 
                         data=GPP_df_23,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)

plot(k_fit_brms_GPP23c)
summary(k_fit_brms_GPP23c)

summary_k_fit_brms <- summary(k_fit_brms_GPP23c)
posterior_summaries_GB23c <- posterior_summary(k_fit_brms_GPP23c)


# Extract posterior samples
posterior_samplesc <- posterior_samples(k_fit_brms_GPP23c)


# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW", "GB", "SH", "SS"), length.out = nrow(posterior_samplesc))

# Prepare the posterior samples data frame with shore levels
posterior_samples_dfc <- posterior_samplesc %>%
  select(b_Intercept, b_logGPPlag, b_par_int_3m, b_lake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted <- melt(posterior_samples_df, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

# Calculate overall mean and credible intervals for each coefficient
ci_summary_overall <- posterior_samples_melted %>%
  group_by(Coefficient) %>%
  summarize(
    mean = mean(Estimate),
    sd = sd(Estimate),
    se = sd(Estimate)/sqrt(length(Estimate)),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
RE_semplot<-ggplot(posterior_samples_melted, aes(x = Estimate, fill = shore, color = shore)) +
  facet_wrap(~ Coefficient, scales = 'free') +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(title = 'Posterior Distributions of Coefficients by Site',
       x = 'Estimate',
       y = 'Density') +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean), linetype = "solid", size = 0.75, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = l_95), linetype = "dashed", size = 0.5, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = u_95), linetype = "dashed", size = 0.5, color = "black")


####### break this out by site
GPP_mod_site <- bf(logGPP ~ logGPPlag + log(par_int_3m) + log(lake_tempC) + log_ppmt + (1 | site)) 




###
##


k_fit_psem <- psem(logGPP ~ logGPPlag + log(par_int_3m) + log(lake_tempC) + log_ppmt +(1|site), data=BW_df) 

plot(k_fit_psem)
summary(k_fit_psem)
fisherC(k_fit_psem)



BW_df <- GPP_df_23%>%
  filter(shore=="BW" & date > as.Date("2023-02-01"))

k_fit_brms_GPP_BW <- brm(GPP_mod_site, 
                         data=BW_df,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)



# Extract posterior samples
posterior_samples_BW <- posterior_samples(k_fit_brms_GPP_BW)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))

# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq <- posterior_samples_BW %>%
  select(b_Intercept, b_logGPPlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BWq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


## GB
GB_df <- GPP_df_23%>%
  filter(shore=="GB" & date > as.Date("2023-02-01"))

k_fit_brms_GPP_GB <- brm(GPP_mod_site, 
                         data=GB_df,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)


summary_k_fit_GB <- summary(k_fit_brms_GPP_GB)
posterior_summaries_GB <- posterior_summary(k_fit_brms_GPP_GB)


# Extract posterior samples
posterior_samples_GB <- posterior_samples(k_fit_brms_GPP_GB)


# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))

# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_Intercept, b_logGPPlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

## SS
SS_df <- GPP_df_23%>%
  filter(shore=="SS" & date > as.Date("2023-02-01"))

k_fit_brms_GPP_SS <- brm(GPP_mod_site, 
                         data=SS_df,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)


summary_k_fit_SS <- summary(k_fit_brms_GPP_SS)
posterior_summaries_SS <- posterior_summary(k_fit_brms_GPP_SS)

# Extract posterior samples
posterior_samples_SS <- posterior_samples(k_fit_brms_GPP_SS)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))

# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_Intercept, b_logGPPlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


## SH
SH_df <- GPP_df_23%>%
  filter(shore=="SH" & date > as.Date("2023-02-01"))

k_fit_brms_GPP_SH <- brm(GPP_mod_site, 
                         data=SH_df,
                         iter=10000, warmup = 5000,
                         control = list(adapt_delta = 0.95, max_treedepth = 15), 
                         cores=4, chains = 3)


summary_k_fit_SH <- summary(k_fit_brms_GPP_SH)
posterior_summaries_SH <- posterior_summary(k_fit_brms_GPP_SH)

# Extract posterior samples
posterior_samples_SH <- posterior_samples(k_fit_brms_GPP_SH)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SH"), length.out = nrow(posterior_samples_SH))

# Prepare the posterior samples data frame with shore levels
posterior_samples_SHq <- posterior_samples_SH %>%
  select(b_Intercept, b_logGPPlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SH <- melt(posterior_samples_SHq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


## All 


SEM_df <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW, posterior_samples_melted_SS, posterior_samples_melted_SH)

ci_summary_overall <- SEM_df %>%
  group_by(Coefficient) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates
site_semplot <-ggplot(SEM_df, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 1, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
       x = NULL,
       y = 'Density') +
  facet_wrap(Coefficient~., ncol=5, scales = 'free') 
  
ggsave(plot = site_semplot, filename = paste("./figures/NS24_metab_all_site_semplot.png",sep=""),width=5.5,height=3.75,dpi=300)

  
site_semplot

##############
##############
## ER 
####### break this out by site

#============================================

ER_df$ER_lag <- lag(ER_df$ER)
ER_df$logER <- log(ER_df$ER+1)
ER_df$logERlag <- log(ER_df$ER_lag+1)

ER_df_23 <- ER_df%>%
  select(site, shore, date, ER, ER_lag, logER, logERlag, 
         lake_tempC, Kd_fill, par_int_3m, tmean_C, light_mean, log_ppmt, 
         windsp_mean, precip_bi, real_NS_depth) %>%
  filter(date > as.Date("2023-02-01"))

summary(ER_df_23)

ER_df_23q<- drop_na(ER_df_23)


ER_mod_site <- bf(logER ~ logERlag + log(par_int_3m) + log(lake_tempC) + log_ppmt  + (1 | site)) 

BW_df_ER <- ER_df_23q%>%
  filter(shore=="BW")

# k_fit_brms_ER_BW <- brm(ER_mod_site, 
#                          data=BW_df_ER,
#                          iter=10000, warmup = 5000,
#                          control = list(adapt_delta = 0.95, max_treedepth = 15), 
#                          cores=4, chains = 3)

k_fit_brms_ER_BW
summary_k_fit_BW <- summary(k_fit_brms_ER_BW)
posterior_summaries_BW <- posterior_summary(k_fit_brms_ER_BW)

# Extract posterior samples
posterior_samples_BW <- posterior_samples(k_fit_brms_ER_BW)


# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("BW"), length.out = nrow(posterior_samples_BW))

# Prepare the posterior samples data frame with shore levels
posterior_samples_BWq <- posterior_samples_BW %>%
  select(b_logERlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_BW <- melt(posterior_samples_BWq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")



## GB

GB_df_ER <- ER_df_23q%>%
  filter(shore=="GB" & date > as.Date("2023-02-01"))

# k_fit_brms_GPP_GB_ER <- brm(ER_mod_site, 
#                          data=GB_df_ER,
#                          iter=10000, warmup = 5000,
#                          control = list(adapt_delta = 0.95, max_treedepth = 15), 
#                          cores=4, chains = 3)

k_fit_brms_GPP_GB_ER

summary_k_fit_GB <- summary(k_fit_brms_GPP_GB_ER)
posterior_summaries_GB <- posterior_summary(k_fit_brms_GPP_GB_ER)

# Extract posterior samples
posterior_samples_GB <- posterior_samples(k_fit_brms_GPP_GB_ER)


# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("GB"), length.out = nrow(posterior_samples_GB))

# Prepare the posterior samples data frame with shore levels
posterior_samples_GBq <- posterior_samples_GB %>%
  select(b_logERlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_GB <- melt(posterior_samples_GBq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")

## SS
SS_df_ER <- ER_df_23q%>%
  filter(shore=="SS" & date > as.Date("2023-02-01"))

# k_fit_brms_ER_SS <- brm(ER_mod_site, 
#                          data=SS_df_ER,
#                          iter=10000, warmup = 5000,
#                          control = list(adapt_delta = 0.95, max_treedepth = 15), 
#                          cores=4, chains = 3)


k_fit_brms_ER_SS

summary_k_fit_SS <- summary(k_fit_brms_ER_SS)
posterior_summaries_SS <- posterior_summary(k_fit_brms_ER_SS)

# Extract posterior samples
posterior_samples_SS <- posterior_samples(k_fit_brms_ER_SS)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SS"), length.out = nrow(posterior_samples_SS))

# Prepare the posterior samples data frame with shore levels
posterior_samples_SSq <- posterior_samples_SS %>%
  select(b_logERlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SS <- melt(posterior_samples_SSq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


## SH

SH_df_ER <- ER_df_23q%>%
  filter(shore=="SH" & date > as.Date("2023-02-01"))

# k_fit_brms_ER_SH <- brm(ER_mod_site,
#                          data=SH_df_ER,
#                          iter=10000, warmup = 5000,
#                          control = list(adapt_delta = 0.95, max_treedepth = 15),
#                          cores=4, chains = 3)

k_fit_brms_ER_SH
summary_k_fit_SH <- summary(k_fit_brms_ER_SH)
posterior_summaries_SH <- posterior_summary(k_fit_brms_ER_SH)

# Extract posterior samples
posterior_samples_SH <- posterior_samples(k_fit_brms_ER_SH)

# Create a sequence of shore levels to match the posterior samples length
shore_levels <- rep(c("SH"), length.out = nrow(posterior_samples_SH))

# Prepare the posterior samples data frame with shore levels
posterior_samples_SHq <- posterior_samples_SH %>%
  select(b_logERlag, b_logpar_int_3m, b_loglake_tempC, b_log_ppmt) %>%
  mutate(shore = shore_levels)

# Melt the data for easier plotting with ggplot2
posterior_samples_melted_SH <- melt(posterior_samples_SHq, id.vars = 'shore', variable.name = 'Coefficient', value.name = 'Estimate')%>%
  filter(Coefficient!="b_Intercept")


## All 


SEM_df_ER <- rbind(posterior_samples_melted_GB, posterior_samples_melted_BW, posterior_samples_melted_SS, posterior_samples_melted_SH)

ci_summary_overall <- SEM_df_ER %>%
  group_by(Coefficient) %>%
  summarize(
    mean = mean(Estimate),
    l_95 = quantile(Estimate, 0.025),
    u_95 = quantile(Estimate, 0.975)
  )

# Plot histograms with credible intervals and mean estimates


# Plot histograms with credible intervals and mean estimates
site_semplot_ER <-ggplot(SEM_df_ER, aes(x = Estimate, fill = shore, color = shore)) +
  geom_density(alpha = 0.75, linewidth = 0.25) +
  geom_vline(xintercept=0, linetype = "solid", size = 1, color = "black") +
  geom_vline(data = ci_summary_overall, aes(xintercept = mean, color=shore), linetype = "solid", size = 0.5) +
  scale_fill_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(name = "Site", values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_minimal() +
  labs(
    x = 'Estimate',
    y = 'Density') +
  facet_wrap(Coefficient~., ncol=5, scales = 'free') 

# ggsave(plot = site_semplot, filename = paste("./figures/NS24_metab_all_site_semplot.png",sep=""),width=5.5,height=3.75,dpi=300)

site_semplot
## total grid:
SEM_grid <- ggarrange(site_semplot, site_semplot_ER,
                      ncol = 1, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")


# ggsave(plot = SEM_grid, filename = paste("./figures/NS24_metab_SEM_23.png",sep=""),width=7.5,height=4.5,dpi=300)
######
# summary stuff at the end:


