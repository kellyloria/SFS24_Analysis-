
library(tidyverse)
library(lubridate)
# plotting packages:
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(gridExtra)
library(ggpubr)
library(stringr)

##############
SFS_datQ <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/SFS24_data_T.rds")
summary(SFS_datQ)
str(SFS_datQ)


se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}


SFS_datQ1 <- SFS_datQ %>%
  mutate(
    position = case_when(
      site %in% c("BWNS1","GBNS1", "SSNS1", "SHNS1") ~ "north",
      site %in% c("BWNS2","GBNS2", "SSNS2", "SHNS2" ) ~ "center",
      site %in% c("BWNS3","GBNS3", "SHNS3") ~ "south",
      TRUE ~ as.character(site)
    )
  )


names(SFS_datQ)
SFS_datQ1 <-SFS_datQ%>%
  filter(middle_GPP <= 30 & middle_ER >= -30)


SFS_datQ_agg <- SFS_datQ1%>%
  group_by(shore, date)%>%
  summarise(
    GPP_m = mean(middle_GPP, na.rm = TRUE),
    GPP_sd = sd(middle_GPP, na.rm = T),
    GPP_se =  GPP_sd/sqrt(length(middle_GPP)),
    GPP_m_low = mean(lower_GPP, na.rm = TRUE),
    GPP_m_up = mean(upper_GPP, na.rm = TRUE),
    ER_m = mean(middle_ER, na.rm = TRUE),
    ER_sd = sd(middle_ER, na.rm = TRUE),
    ER_se =  se(middle_ER),
    ER_m_low = mean(lower_ER, na.rm = TRUE),
    ER_m_up = mean(upper_ER, na.rm = TRUE),
  )

SFS_datQ_agg$year <- year(SFS_datQ_agg$date)
SFS_datQ_agg$yday <- yday(SFS_datQ_agg$date)
SFS_datQ_agg$week <- week(SFS_datQ_agg$date)


SFS_datQ_agg21<- SFS_datQ_agg%>%
  filter(year==2021)

SFS_datQ_agg22<- SFS_datQ_agg%>%
  filter(year==2022)

SFS_datQ_agg23<- SFS_datQ_agg%>%
  filter(year==2023)



#####################
### weekly plots ###
#####################
# Create two-week intervals
SFS_datQ_agg$weekyr <- as.numeric(paste0(SFS_datQ_agg$year, sprintf("%02d", SFS_datQ_agg$week)))
aggregate_weeks <- function(weekyr) {
  ceiling(weekyr / 4)
}

# Apply the function to create the weekyr3 variable
SFS_datQ_agg$weekyr3 <- aggregate_weeks(SFS_datQ_agg$weekyr)
# Back calculate the date for each label in the x-axis based on weekyr3
starting_date <- as.Date("2021-06-12")
ending_date <- as.Date("2023-09-06")
# Create a sequence of dates occurring every 3 weeks
week_dates <- seq(starting_date, ending_date, by = "4 weeks")
range(week_dates)
range(SFS_datQ_agg$date)
formatted_labels <- format(week_dates, "%b-%y")

# Plot with aggregated intervals
G_plot <- ggplot(data = SFS_datQ_agg, aes(x = as.factor(weekyr3), y = GPP_m, fill = shore, color = shore)) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + ylim(0,20)+
  geom_hline(yintercept = 0, size = 0.3, color = "gray25") +
  labs(y=expression(GPP~mmol~O[2]~m^-3~d^-1), x = NULL) + 
  scale_x_discrete(labels =formatted_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
G_plot


E_plot <- ggplot(data = SFS_datQ_agg, aes(x = as.factor(weekyr3), y = ER_m, fill = shore, color = shore)) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + ylim(-20,0)+
  geom_hline(yintercept = 0, size = 0.3, color = "gray25") +
  scale_x_discrete(labels =formatted_labels) +
  labs(y=expression(ER~mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
E_plot


NEP_plot <- ggplot(data = SFS_datQ_agg, aes(x = as.factor(weekyr3), y = c(GPP_m + ER_m), fill = shore, color = shore)) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + ylim(-20,20)+
  geom_hline(yintercept = 0, size = 0.3, color = "gray25") +
  scale_x_discrete(labels =formatted_labels) +
  labs(y=expression(NEP~mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
NEP_plot

metab23a <- ggarrange(G_plot,
                      E_plot,
                      NEP_plot,
                      ncol = 1, nrow = 3,
                      common.legend = TRUE, 
                      labels=c("A", "B", "C"),
                      legend = "bottom")

# ggsave(plot = metab23a, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/Figures/NS24_metab_all_RF.png",sep=""),width=9,height=8,dpi=300)




###########
names(SFS_datQ1)
SFS_datQ1$DOY <- yday(SFS_datQ1$date)
SFS_datQ1$year <- year(SFS_datQ1$date)

#   se <- sd(dat)/sqrt(length(dat))
tmean_Cse=sd(SFS_datQ1$tmean_C)/ sqrt(length(SFS_datQ1$tmean_C))

SFS_env <- SFS_datQ1 %>%
  group_by(shore, DOY, year) %>%
  summarise(
    middle_GPP_m = mean(middle_GPP, na.rm = TRUE),
    middle_GPP_se = sd(middle_GPP_m, na.rm = TRUE) / sqrt(n()), # SE for middle_GPP
    lower_GPP_m = mean(lower_GPP, na.rm = TRUE),
    lower_GPP_se = sd(lower_GPP, na.rm = TRUE) / sqrt(n()), # SE for lower_GPP
    upper_GPP_m = mean(upper_GPP, na.rm = TRUE),
    upper_GPP_se = sd(upper_GPP, na.rm = TRUE) / sqrt(n()), # SE for upper_GPP
    tmean_Cm = mean(tmean_C, na.rm = TRUE),
    tmean_Cm_se = sd(tmean_C, na.rm = TRUE) / sqrt(n()), # SE for tmean_Cm
    light_mean = mean(light_mean, na.rm = TRUE),
    light_se = sd(light_mean, na.rm = TRUE) / sqrt(n()), # SE for light_mean
    stream_temp = mean(stream_temp, na.rm = TRUE),
    stream_temp_se = sd(stream_temp, na.rm = TRUE) / sqrt(n()), # SE for stream_temp
    spc_mean = mean(spc_mean, na.rm = TRUE),
    spc_mean_se = sd(spc_mean, na.rm = TRUE) / sqrt(n()), # SE for spc_mean
    lake_do = mean(lake_do, na.rm = TRUE),
    lake_do_se = sd(lake_do, na.rm = TRUE) / sqrt(n()), # SE for lake_do
    lake_wtemp = mean(lake_wtemp, na.rm = TRUE),
    lake_wtemp_se = sd(lake_wtemp, na.rm = TRUE) / sqrt(n()), # SE for lake_wtemp
    lake_par = mean(lake_par, na.rm = TRUE),
    lake_par_se = sd(lake_par, na.rm = TRUE) / sqrt(n()), # SE for lake_par
    lake_par_int = mean(lake_par_int, na.rm = TRUE),
    lake_par_int_se = sd(lake_par_int, na.rm = TRUE) / sqrt(n()), # SE for lake_par_int
    windsp_mean = mean(windsp_mean, na.rm = TRUE),
    windsp_mean_se = sd(windsp_mean, na.rm = TRUE) / sqrt(n()), # SE for windsp_mean
    flow_sum = sum(flow_mean, na.rm = TRUE),
    ppt_mean = mean(ppt_mm, na.rm = TRUE),
    flow_mean = mean(flow_mean, na.rm = TRUE),
    flow_mean_se = sd(flow_mean, na.rm = TRUE) / sqrt(n()) # SE for flow_mean
  )


names(SFS_env)

plot_ws <- ggplot(data = SFS_env, aes(x =DOY, y = windsp_mean, color = shore)) +
  geom_point(alpha=0.5)+
  #geom_boxplot(width = 0.75, alpha = 0.5) +  
  #scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +
  labs(x=expression(Day~of~year), y=expression(Mean~daily~wind~sp~m~s^-1)) 
plot




plot_at <- ggplot(data = SFS_env, aes(x =DOY, y = tmean_Cm, color = shore)) +
  geom_point(alpha=0.5)+
  #geom_boxplot(width = 0.75, alpha = 0.5) +  
  #scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +
  labs(x=expression(Day~of~year), y=expression(Mean~daily~air~temp~C)) 
plot


plot_par <- ggplot(data = SFS_env, aes(x =DOY, y = light_mean, color = shore)) +
  geom_point(alpha=0.5)+
  #geom_boxplot(width = 0.75, alpha = 0.5) +  
  #scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +
  labs(x=expression(Day~of~year), y=expression(Mean~daily~PAR~umol~m^-2~s^-1)) 
plot


SFS_env$precip_bi <- ifelse(SFS_env$ppt_mean > 0, 1, 0)


plot_ppt <- ggplot(data = SFS_env, aes(x = DOY, y = log(ppt_mean+1), fill = shore)) +
  geom_bar(stat = "identity", position = "dodge", width = 2) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() +
  labs(x = expression(Day ~ of ~ year), y = expression(log(Mean~daily~precip+1)~mm))


cliamte <- ggarrange(plot_at, 
                     plot_par,
                     plot_ws,
                     plot_ppt,
                     ncol = 1, nrow = 4,
                     common.legend = TRUE, 
                     labels=c("A", "B", "C", "D"),
                     #x = 0.95, y = 0.95,
                     legend = "bottom")


# ggsave(plot = cliamte, filename = paste("./figures/NS24_climate_F.png",sep=""),width=6,height=12,dpi=300)

GB <-(SFS_env%>%filter(shore=="GB")%>%select(windsp_mean))
summary(GB)

BW <-(SFS_env%>%filter(shore=="BW")%>%select(windsp_mean))
summary(BW)

SS <-(SFS_env%>%filter(shore=="SS")%>%select(windsp_mean))
summary(SS)

SH <-(SFS_env%>%filter(shore=="SH")%>%select(windsp_mean))
summary(SH)

