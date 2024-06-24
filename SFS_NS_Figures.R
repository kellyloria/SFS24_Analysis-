#' Script for creating figure for SFS 2024 conference 
#' 
#' @description nearshore metabolism output visualization 
#' @param lake in file name of .cvs of "\LittoralMetabModeling\FinalInputs\"
#' 
#' @return Returns .r file for running metabolism model 
#' @export 
#' 

##==============================================================================
## Created  01/31/2024 by KAL
#===============================================================================

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
# 
# ##============================================
# ## read in lake metabolism estimates
# #============================================
# # * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/
# ## temporary local path:
# ## setwd("/Users/kellyloria/Documents/UNR/MSMmetab/")
# getwd()
# 
# # Function to add a column based on the file name to a data frame
# add_site_column <- function(file_name) {
#   df <- read.csv(file_name)
#   # Extract information from the file name:
#   site_info <- str_match(file_name, "([A-Za-z0-9_]+)__daily_full.csv")
#   # Check if there is a match and extract the relevant part from the match
#   if (!is.null(site_info) && length(site_info) == 2) {
#     site_name <- site_info[2]
#   } else {
#     site_name <- "UnknownSite"
#   }
#   # Add a new column to the data frame with the extracted site name
#   df$site <- site_name
#   return(df)
# }
# 
# # Specify the folder path:
# folder_path <- "./plotdata_24"
# # Get a list of CSV files in the folder
# csv_files <- list.files(folder_path, pattern = "__daily_full.csv$", full.names = TRUE)
# # Apply the function to each CSV file and store the results in a list
# processed_data_list <- lapply(csv_files, add_site_column)
# # Combine the list of data frames into a single data frame
# LM_data <- do.call(rbind, processed_data_list)
# 
# ## subset data:
# LM_data <-LM_data %>% 
#   subset(name=="ER" | name=="GPP"|name=="NEP") %>%
#   subset(middle <= 30 & middle >= -30)
# 
# # add in column for date:
# LM_data$origin <- as.Date(paste0(LM_data$year, "-01-01"),) 
# LM_data$date <-as.Date(LM_data$yday, origin = LM_data$origin) 
# 
# 
# # quick plot of all data in long form: 
# p2 <- ggplot(data = LM_data %>% drop_na(year),aes(date, middle, color = name))+
#   geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = name),
#               linetype = 0, alpha = 0.2)+
#   geom_line()+ geom_point(size= 3, alpha = 0.6)+
#   # geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
#   #            color = "#4c4d4c") +
#   #geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
#   scale_color_manual(values = c("#982649", "#003F91", "#333333")) +
#   scale_fill_manual(values = c("#982649","#003F91","black")) +
#   scale_x_date(date_breaks = "3 month", date_labels = "%b-%y")+ #new
#   ylim(-35, 35) +
#   labs(y=expression(mmol~O[2]~m^-3~d^-1), x=NULL) + 
#   theme_classic() + 
#   theme(axis.text=element_text(size=18),
#         axis.title=element_text(size=20),
#         legend.position = 'bottom', 
#         legend.direction = "horizontal") +
#   facet_grid((site~.))
# p2
# 
# # ggsave(plot = p2, filename = paste("./SFS24_figures/24NS_all.png",sep=""),width=12.5,height=13,dpi=300)
# 
# 
# ## total grid:
# # all_grid <- ggarrange(tempgrid,
# #                       Vel_grid,
# #                       PS_grid,
# #                       NDVIgrid,
# #                       ncol = 1, nrow = 4,
# #                       widths = c(1,1, 0.7),
# #                       common.legend = TRUE, 
# #                       legend = "bottom")
# # 
# 
# 
# ##===========================================
# ## read data aggregated data for the project:
# #============================================
# # * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/

getwd()
SFS_datQ <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/SFS24_data_T.rds")
summary(SFS_datQ)
str(SFS_datQ)

##===========================================
## Begin plotting  
#============================================

## Best fit glm 
##  log(GPP+1) = β0 + β1 *scale(lag_GPP)+β3*scale(light_mean)+ β4*scale(windsp_cv)+ β6*scale(spc_mean)+b0i +ϵi 

se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}


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





GPP_plot <- ggplot(data = SFS_datQ_agg, aes(date, GPP_m, color = shore))+
  geom_ribbon(aes(ymin = GPP_m-GPP_sd, ymax = GPP_m+GPP_sd, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%y")+ #new
  ylim(-0, 35) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP", mmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right', 
        legend.direction = "vertical") + facet_grid(year~.)

# ggsave(plot = GPP_plot, filename = paste("./figures/NS24_GPP_all.png",sep=""),width=12,height=5,dpi=300)


GPP_plot <- ggplot(data = SFS_datQ_agg21, aes(yday, GPP_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = GPP_m_low, ymax = GPP_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(0, 32) +xlim(0,360)+
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP", mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)

# ggsave(plot = GPP_plot, filename = paste("./figures/NS24_GPP_all.png",sep=""),width=12,height=4,dpi=300)


ER_plot <- ggplot(data = SFS_datQ_agg21, aes(yday, ER_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = ER_m_low, ymax = ER_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(-32,0) + xlim(0,360)+
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER", mmol~O[2]~m^-3~d^-1), x = "Day of year") + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)


metab21 <- ggarrange(GPP_plot,
                     ER_plot,
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, 
                    legend = "bottom")

# ggsave(plot = metab21, filename = paste("./figures/NS24_metab21.png",sep=""),width=14,height=7,dpi=300)


SFS_datQ_agg$ER_m_t <- (SFS_datQ_agg$ER_m *-1)

GPP_ER_plot <- ggplot(data = SFS_datQ_agg, aes(ER_m_t, GPP_m, color = shore))+
  geom_point(size= 3, alpha = 0.6) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression(GPP~mmol~O[2]~m^-3~d^-1), x= expression(ER~mmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right', 
        legend.direction = "vertical") + facet_grid(shore~.)


# Add the R-squared value and 1:1 line to each facet
GPP_ER_plot <- GPP_ER_plot +
  stat_cor(method = "pearson", label.x = 0.8, label.y = 24, size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")

# ggsave(plot = GPP_ER_plot, filename = paste("./figures/GPP_ER_plot_24.png",sep=""),width=4.25,height=8,dpi=300)



GPP_plot <- ggplot(data = SFS_datQ_agg22, aes(yday, GPP_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = GPP_m_low, ymax = GPP_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(0, 32) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP", mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)

# ggsave(plot = GPP_plot, filename = paste("./figures/NS24_GPP_all.png",sep=""),width=12,height=4,dpi=300)


ER_plot <- ggplot(data = SFS_datQ_agg22, aes(yday, ER_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = ER_m_low, ymax = ER_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(-32,0) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER", mmol~O[2]~m^-3~d^-1), x = "Day of year") + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)


metab22 <- ggarrange(GPP_plot,
                     ER_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     legend = "bottom")

# ggsave(plot = metab22, filename = paste("./figures/NS24_metab22.png",sep=""),width=14,height=7,dpi=300)




GPP_plot <- ggplot(data = SFS_datQ_agg23, aes(yday, GPP_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = GPP_m_low, ymax = GPP_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(0, 32) +  xlim(0, 360) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP", mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)

ER_plot <- ggplot(data = SFS_datQ_agg23, aes(yday, ER_m, color = shore))+
  geom_line() +
  geom_ribbon(aes(ymin = ER_m_low, ymax = ER_m_up, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_point(size= 3, alpha = 0.6)+
  ylim(-32,0) +  xlim(0, 360) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER", mmol~O[2]~m^-3~d^-1), x = "Day of year") + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") + facet_grid(year~.)


metab23 <- ggarrange(GPP_plot,
                     ER_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     legend = "bottom")

# ggsave(plot = metab23, filename = paste("./figures/NS24_metab23.png",sep=""),width=14,height=7,dpi=300)






GPP_plot <- ggplot(data = SFS_datQ_agg, aes(yday, GPP_m, color = shore, shape=as.factor(year)))+
  geom_point(size= 3, alpha = 0.6)+
  ylim(0, 20) +  xlim(0, 360) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP", mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") 

ER_plot <- ggplot(data = SFS_datQ_agg, aes(yday, ER_m, color = shore, shape=as.factor(year)))+
  geom_point(size= 3, alpha = 0.6)+
  ylim(-20,0) +  xlim(0, 360) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER", mmol~O[2]~m^-3~d^-1), x = "Day of year") + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal")


metab23 <- ggarrange(GPP_plot,
                     ER_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     legend = "bottom")

# ggsave(plot = metab23, filename = paste("./figures/NS24_metab_all_shapeyear.png",sep=""),width=14,height=7,dpi=300)

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
  labs(y=expression(GPP~mmol~O[2]~m^-3~d^-1), x = NULL) + 
  scale_x_discrete(labels =formatted_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
G_plot


E_plot <- ggplot(data = SFS_datQ_agg, aes(x = as.factor(weekyr3), y = ER_m, fill = shore, color = shore)) +
  geom_boxplot(width = 0.75, alpha = 0.5) +  
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_color_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  theme_bw() + ylim(-20,0)+
  scale_x_discrete(labels =formatted_labels) +
  labs(y=expression(ER~mmol~O[2]~m^-3~d^-1), x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
E_plot

metab23a <- ggarrange(G_plot,
                     E_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     legend = "bottom")

# ggsave(plot = metab23a, filename = paste("./figures/NS24_metab_all_F.png",sep=""),width=9,height=5,dpi=300)



GGP_b4_plot <- ggplot(data = SFS_datQ1, aes(windsp_cv, log(middle_GPP+1), color = shore)) +
  #geom_pointrange(aes(ymin = log(lower_GPP+1), ymax = log(upper_GPP+1), fill = site.x), alpha = 0.2)+
  geom_point(size = 3, alpha = 0.6) +
  labs(y = "log(GPP+1)")+
  theme_bw() + 
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20), 
    legend.position = 'right', 
    legend.direction = "vertical"
  )



GGP_b6_plot <- ggplot(data = SFS_datQ1, aes(spc_mean, log(middle_GPP+1), color = site.x)) +
  geom_pointrange(aes(ymin = log(lower_GPP+1), ymax = log(upper_GPP+1), fill = site.x), alpha = 0.2)+
  geom_point(size = 3, alpha = 0.6) +
  labs(y = "log(GPP+1)")+
  theme_bw() + 
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20), 
    legend.position = 'right', 
    legend.direction = "vertical"
  )




GGP_b1_plot <- ggplot(data = SFS_datQ1, aes(light_cv, log(middle_GPP+1), color = site.x)) +
  geom_pointrange(aes(ymin = log(lower_GPP+1), ymax = log(upper_GPP+1), fill = site.x), alpha = 0.2)+
  geom_point(size = 3, alpha = 0.6) +
  labs(y = "log(GPP+1)")+
  theme_bw() + 
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20), 
    legend.position = 'right', 
    legend.direction = "vertical"
  )


#############################################################
############
############ Weather covaraites for the report 
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










#######
######
######
## JUNK CODE BELOW 



# LAKE NEP drivers

# streamflow?
# DOdat from 2022Figures 





sumNEPplot <- ggplot(data = sumNEPF %>% drop_na(year),aes(log(dischargeCFS), middle.NEP, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_pointrange(aes(ymin = lower, ymax = upper, fill = shore), alpha = 0.2)+
  #geom_line()+ 
  geom_point(size= 1, alpha = 0.4)+
  #geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
  # color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  geom_smooth(method="lm", se=F) +
  # scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  #ylim(-15, 15) +
  scale_y_continuous(breaks = seq(-15, 24, by = 3))+
  labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        #axis.text.x=element_text(angle=60, hjust=1), # new
        legend.position = 'right', 
        legend.direction = "vertical") 
#+ facet_grid((shore~.))

# ggsave(plot = sumNEPplot, filename = paste("./figures/NEP_flow_all.png",sep=""),width=4.25,height=4,dpi=300)



sumNEPplot <- ggplot(data = sumNEPFSPC %>% drop_na(year),aes(log(SPCmean), middle.NEP, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_pointrange(aes(ymin = lower, ymax = upper, 
                      fill = shore), alpha = 0.2)+
  #geom_line()+ 
  geom_point(size= 1, alpha = 0.4)+
  #geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
  # color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  geom_smooth(method="lm", se=F) +
  # scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  #ylim(-15, 15) +
  scale_y_continuous(breaks = seq(-15, 24, by = 3))+
  labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        #axis.text.x=element_text(angle=60, hjust=1), # new
        legend.position = 'right', 
        legend.direction = "vertical") 
#+ facet_grid((shore~.))

# ggsave(plot = sumNEPplot, filename = paste("./figures/NEP_SPC_all.png",sep=""),width=4.25,height=4,dpi=300)


mod1<- lmer(middle.NEP~ scale(SPCmean)+ shore+(1|site), data=sumNEPFSPC)
summary(mod1)


mod2<- lmer(middle.NEP~ scale(SPCmean) + scale(dischargeCFS)+ shore + (1|site), data=sumNEPFSPC)
summary(mod2)

hist(residuals(mod1))
r.squaredGLMM(mod1)
# 
# Soil_dat_l$date2<- Soil_dat_l$date - 20
# 
# sumNEPFSPC_A <- left_join(sumNEPFSPC, Soil_dat_l[,c("date2", "AFDMmgmL", "shore")],  by=c('date'='date2', 'shore'='shore'))
# summary(sumNEPFSPC_A)
# 
# sumNEPFSPC_A <- left_join(Soil_dat_l, sumNEPFSPC[,c("date", "middle.NEP", "shore")],  by=c('date2'='date', 'shore'='shore'))
# summary(sumNEPFSPC_A)
# 
# 
# sumNEPAplot <- ggplot(data = sumNEPFSPC %>% drop_na(year),aes((AFDMmgmL), middle.NEP, color = shore))+
#   geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
#   geom_pointrange(aes(ymin = lower, ymax = upper, 
#                       fill = shore), alpha = 0.2)+
#   #geom_line()+ 
#   geom_point(size= 1, alpha = 0.4)+
#   #geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
#   # color = "#4c4d4c") +
#   scale_color_manual(values=c("#a67d17", "#3283a8")) +
#   scale_fill_manual(values = c("#a67d17", "#3283a8")) +
#   geom_smooth(method="lm", se=F) +
#   # scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
#   #ylim(-15, 15) +
#   scale_y_continuous(breaks = seq(-15, 24, by = 3))+
#   labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
#   theme_classic() + 
#   theme(axis.text=element_text(size=11),
#         axis.title=element_text(size=11),
#         #axis.text.x=element_text(angle=60, hjust=1), # new
#         legend.position = 'right', 
#         legend.direction = "vertical") 
# 
# 




## STREAM PLOTS ##

BWL <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/plotDat/BWL_daily.csv")
BWL$date <- as.Date(BWL$date, origin="2021-01-01")
BWL$site <- "BWL"
BWL$shore <- "west"


GBL <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/plotDat/GBL_daily.csv")
GBL$date <- as.Date(GBL$date, origin="2021-01-01")
GBL$site <- "GBL"
GBL$shore <- "east"


Eplot_sp <- ggplot(GBL, aes(x = GPP_mean, y = ER_mean)) + ylab("Ecosystem respiration") + 
  xlab("Gross primary productivity") + #ylim(-25, 0) + xlim(0,5) + 
  geom_point(aes(x = GPP_mean, y = ER_mean), shape= 17, col = alpha(c("#a67d17"),0.5)) +
  geom_abline(intercept = 0, slope = -1, col = "grey50")+
  stat_density2d(aes(colour = ..level..)) +
  scale_colour_gradient(
    low = "#faefd4",
    high = "#a67d17",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )  +
  theme_classic()


Odd_plot <- ggplot(GBL, aes(x = K600_daily_mean, color = shore, fill = shore)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 15) +
  scale_fill_manual(values = alpha(c("#a67d17"), 0.2)) +
  scale_color_manual(values = alpha(c("#a67d17"), 0.9)) + theme_bw() + 
  geom_vline(data = GBL, aes(xintercept = mean(na.omit(K600_daily_mean))), linetype = "dashed") 

Gplot_sp <- ggplot(GBL, aes(x = K600_daily_mean, y = (ER_mean*-1))) +  
  geom_point(shape= 17, col = alpha(c("#a67d17"),0.5)) +
  geom_smooth(method ="lm", se=F)+  facet_grid(.~site)+
  theme_bw()

GB_lm <- lm(K600_daily_mean~(ER_mean*-1), data=GBL)
summary(GB_lm)


# Extract R-squared value
r_gsquared <- summary(GB_lm)$r.squared

GB_plot <- Gplot_sp + 
  geom_text(aes(x = max((GBL$K600_daily_mean), na.rm=T), y = min((GBL$ER_mean*-1), na.rm=T), 
                label = paste("R-squared =", round(r_gsquared, 2))),
            hjust = 1, vjust = 0, size = 4, col = "blue",
            parse = T, check_overlap = T, na.rm = T)



Wplot_sp <- ggplot(BWL, aes(x = GPP_mean, y = ER_mean)) + ylab("Ecosystem respiration") + 
  xlab("Gross primary productivity") + #ylim(-25, 0) + xlim(0,5) + 
  geom_point(aes(x = GPP_mean, y = ER_mean), shape= 17, col = alpha(c("#3283a8"),0.5)) +
  geom_abline(intercept = 0, slope = -1, col = "grey50")+
  stat_density2d(aes(colour = ..level..)) +
  scale_colour_gradient(
    low = "#3283a8",
    high = "#3283a8",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )  +
  theme_classic()


Odd_plot <- ggplot(BWL, aes(x = K600_daily_mean, color = shore, fill = shore)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 15) +
  scale_fill_manual(values = alpha(c("#3283a8"), 0.2)) +
  scale_color_manual(values = alpha(c("#3283a8"), 0.9)) + theme_bw() + 
  geom_vline(data = BWL, aes(xintercept = mean(na.omit(K600_daily_mean))), linetype = "dashed") 


Wplot_sp <- ggplot(BWL, aes(x = K600_daily_mean, y = (ER_mean*-1))) +  
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.5)) +
  geom_smooth(method ="lm", se=F) + facet_grid(.~site)+
  #geom_abline(intercept = 0, slope = 1, col = "grey50")  +
  theme_bw()

bw_lm <- lm(K600_daily_mean~(ER_mean*-1), data=BWL)
summary(bw_lm)


# Extract R-squared value
r_squared <- summary(bw_lm)$r.squared

BW_plot <- Wplot_sp + 
  geom_text(aes(x = max((BWL$K600_daily_mean), na.rm=T), y = min((BWL$ER_mean*-1), na.rm=T), 
                label = paste("R-squared =", round(r_squared, 2))),
            hjust = 1, vjust = 0, size = 4, col = "blue",
            parse = T, check_overlap = T, na.rm = T) 


## total grid:
k_grid <- ggarrange(GB_plot,
                      BW_plot,
                      ncol = 1, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")



ggsave(plot = k_grid, filename = paste("./figures/24_streamMetab_k_all.png",sep=""),width=5,height=7,dpi=300)



BWL$NEP<- BWL$GPP_daily_mean + BWL$ER_daily_mean
BWL$NEP_sd<- sd(BWL$GPP_daily_mean + BWL$ER_daily_mean)

GBL$NEP<- GBL$GPP_daily_mean + GBL$ER_daily_mean
GBL$NEP_sd<- sd(GBL$GPP_daily_mean + GBL$ER_daily_mean)


names(GBL)

streamNEP <- rbind(GBL, BWL)

StreamGPPplot <- ggplot(data = streamNEP, aes(date, GPP_daily_mean, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50") +
  geom_point(size= 1, alpha = 0.6) +
  geom_line() +
  geom_pointrange(aes(ymin =(GPP_2.5pct), 
                      ymax = (GPP_97.5pct)), alpha = 0.75) +
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
             color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y")+ #new
  scale_y_continuous(breaks = seq(-30, 25, by = 5))+
  labs(y=expression(GPPmmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),# new
        legend.position = 'right')

## total grid:
NEP_grid <- ggarrange(StreamGPPplot,
                      StreamERplot,
                      ncol = 1, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")

# save plot to figures folder
#  ggsave(plot = NEP_grid, filename = paste("./NEP_grid.png",sep=""),width=8,height=6.5,dpi=300)

###

###
### Gas exchange:
###

StreamERplot <- ggplot(data = streamNEP, aes(K600_daily_mean, (ER_daily_mean*-1), color = shore))+
  geom_point(size= 1, alpha = 0.6) +
  geom_pointrange(aes(ymin =(-1*ER_2.5pct), 
                      ymax = (-1*ER_97.5pct)), alpha = 0.75) +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  scale_y_continuous(breaks = seq(-30, 25, by = 5))+
  labs(y=expression(ERmmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right')+ facet_grid(.~shore)

###
Odd_plot <- ggplot(streamNEP, aes(x = K600_daily_mean, color = shore, fill = shore)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 15) +
  scale_fill_manual(values = alpha(c("#3283a8", "#a67d17"), 0.2)) +
  scale_color_manual(values = alpha(c("#3283a8", "#a67d17"), 0.9)) + theme_bw() + 
  geom_vline(data = streamNEP, aes(xintercept = mean(na.omit(K600_daily_mean))), linetype = "dashed") +
  facet_grid(.~shore)



sumNEP_SM <- left_join(sumNEPF, streamNEP[,c("date", "shore", "NEP", "GPP_daily_mean", "ER_daily_mean")],  by=c('date'='date', 'shore'='shore'))



sumNEPplot <- ggplot(data = sumNEP_SM %>% drop_na(year),aes(GPP_daily_mean, middle.NEP, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_pointrange(aes(ymin = lower, ymax = upper, 
                      fill = shore), alpha = 0.2)+
  #geom_line()+ 
  geom_point(size= 1, alpha = 0.4)+
  #geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
  # color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  geom_smooth(method="lm", se=F) +
  # scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  #ylim(-15, 15) +
  #scale_y_continuous(breaks = seq(-25, 25, by = 5))+
  #scale_x_continuous(breaks = seq(-30, 25, by = 5))+
  labs(y="Nearshore NEP", expression(mmol~O[2]~m^-3~d^-1),
       x="Stream NEP", expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11), # new
        legend.position = 'right', 
        legend.direction = "vertical") + facet_grid((shore~.))

ggsave(plot = sumNEPplot, filename = paste("./figures/LA_SM_GNEP_all.png",sep=""),width=5,height=7,dpi=300)
