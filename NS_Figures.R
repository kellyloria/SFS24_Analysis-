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

##============================================
## read in lake metabolism estimates
#============================================
# * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/
## temporary local path:
## setwd("/Users/kellyloria/Documents/UNR/MSMmetab/")
getwd()

# Function to add a column based on the file name to a data frame
add_site_column <- function(file_name) {
  df <- read.csv(file_name)
  # Extract information from the file name:
  site_info <- str_match(file_name, "([A-Za-z0-9_]+)__daily_full.csv")
  # Check if there is a match and extract the relevant part from the match
  if (!is.null(site_info) && length(site_info) == 2) {
    site_name <- site_info[2]
  } else {
    site_name <- "UnknownSite"
  }
  # Add a new column to the data frame with the extracted site name
  df$site <- site_name
  return(df)
}

# Specify the folder path:
folder_path <- "./plotdata_24"
# Get a list of CSV files in the folder
csv_files <- list.files(folder_path, pattern = "__daily_full.csv$", full.names = TRUE)
# Apply the function to each CSV file and store the results in a list
processed_data_list <- lapply(csv_files, add_site_column)
# Combine the list of data frames into a single data frame
LM_data <- do.call(rbind, processed_data_list)

## subset data:
LM_data <-LM_data %>% 
  subset(name=="ER" | name=="GPP"|name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)

# add in column for date:
LM_data$origin <- as.Date(paste0(LM_data$year, "-01-01"),) 
LM_data$date <-as.Date(LM_data$yday, origin = LM_data1$origin) 


# quick plot of all data in long form: 
p2 <- ggplot(data = LM_data1 %>% drop_na(year),aes(date, middle, color = name))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = name),
              linetype = 0, alpha = 0.2)+
  geom_line()+ geom_point(size= 3, alpha = 0.6)+
  # geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
  #            color = "#4c4d4c") +
  #geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
  scale_color_manual(values = c("#982649", "#003F91", "#333333")) +
  scale_fill_manual(values = c("#982649","#003F91","black")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%y")+ #new
  ylim(-35, 35) +
  labs(y=expression(mmol~O[2]~m^-3~d^-1), x=NULL) + 
  theme_classic() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'bottom', 
        legend.direction = "horizontal") +
  facet_grid((site~.))
p2

# ggsave(plot = p2, filename = paste("./SFS24_figures/24NS_all.png",sep=""),width=12.5,height=13,dpi=300)


## total grid:
# all_grid <- ggarrange(tempgrid,
#                       Vel_grid,
#                       PS_grid,
#                       NDVIgrid,
#                       ncol = 1, nrow = 4,
#                       widths = c(1,1, 0.7),
#                       common.legend = TRUE, 
#                       legend = "bottom")
# 


##===========================================
## read data aggregated data for the project:
#============================================
# * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/

getwd()
SFS_datQ <- read.csv("./SFS24_analysis_dat/SFS24_analysis_dat.csv") %>%
  mutate(date = as.Date(date))

str(SFS_datQ)
##===========================================
## Begin plotting  
#============================================

## Best fit glm 
##  log(GPP+1) = β0 + β1 *scale(lag_GPP)+β3*scale(light_mean)+ β4*scale(windsp_cv)+ β6*scale(spc_mean)+b0i +ϵi 

names(SFS_datQ)
SFS_datQ1 <-SFS_datQ%>%
  filter(middle_GPP <= 30 & middle_ER >= -30)

GGP_plot <- ggplot(data = SFS_datQ1, aes(date, middle_GPP, color = site.x))+
  geom_ribbon(aes(ymin = lower_GPP, ymax = upper_GPP, fill = site.x),
              linetype = 0, alpha = 0.2)+
  geom_line()+ geom_point(size= 3, alpha = 0.6)+
  #scale_color_manual(values=c("#a67d17", "#3283a8")) +
  #scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%y")+ #new
  #ylim(-15, 15) +
  #scale_y_continuous(breaks = seq(-15, 24, by = 3))+
  labs(y=expression("log(GPP+1)",mmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right', 
        legend.direction = "vertical")

# ggsave(plot = sumNEPplot, filename = paste("./figures/LA_NEP_all.png",sep=""),width=8.5,height=6.5,dpi=300)


GGP_b4_plot <- ggplot(data = SFS_datQ1, aes(windsp_cv, log(middle_GPP+1), color = site.x)) +
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




## STREAM PLOT

BWL <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/plotDat/BWL_daily.csv")
BWL$date <- as.Date(BWL$date, origin="2021-01-01")
BWL$site <- "BWL"
BWL$shore <- "west"


GBL <- read.csv("./GBL_daily.csv")
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


Wplot_sp <- ggplot(BWL, aes(x = GPP_mean, y = ER_mean)) + ylab("Ecosystem respiration") + 
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



BWL$NEP<- BWL$GPP_daily_mean + BWL$ER_daily_mean
BWL$NEP_sd<- sd(BWL$GPP_daily_mean + BWL$ER_daily_mean)

GBL$NEP<- GBL$GPP_daily_mean + GBL$ER_daily_mean

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




StreamERplot <- ggplot(data = streamNEP, aes(date, ER_daily_mean, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50") +
  geom_point(size= 1, alpha = 0.6) +
  geom_line() +
  geom_pointrange(aes(ymin =(ER_2.5pct), 
                      ymax = (ER_97.5pct)), alpha = 0.75) +
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
             color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y")+ #new
  scale_y_continuous(breaks = seq(-30, 25, by = 5))+
  labs(y=expression(ERmmol~O[2]~m^-3~d^-1)) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right')




## total grid:
NEP_grid <- ggarrange(StreamGPPplot,
                      StreamERplot,
                      ncol = 1, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")

# save plot to figures folder
#  ggsave(plot = NEP_grid, filename = paste("./NEP_grid.png",sep=""),width=8,height=6.5,dpi=300)



# ggsave(plot = StreamNEPplot, filename = paste("./figures/SM_NEP_all.png",sep=""),width=9.5,height=6.5,dpi=300)

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
