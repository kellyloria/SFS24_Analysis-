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
library(dataRetrieval)

##===========================================
## read in lake metabolism estimates
#============================================
# * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/

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

#### Convert data from long to wide
str(LM_data)
LM_data_wide <- LM_data %>%
  pivot_wider(
    id_cols = c(date, site),
    names_from = name,
    values_from = c(middle, upper, lower)
  )

##===========================================
## Bring in co variate data
#============================================
# precip, streamflow, stream temp, stream spc, chem, stream metab

# Define the site numbers for "GB" and "BW"
siteNo_GB <- "10336730"
siteNo_BW <- "10336660"

# Define the parameter codes for flow and stage
pCode_flow <- "00060"
pCode_stage <- "00065"

# Set the start date to "1980-01-01" and end date to today (current date)
start.date <- "2020-10-01"
end.date <- "2023-10-01"  # Use
## Download hourly flow data for both sites separately and combine
HRflow_data_GB <- readNWISuv(siteNumbers = siteNo_GB, parameterCd = pCode_flow, startDate = start.date, endDate = end.date) %>%
  mutate(site = "GBL", shore = "GB") %>% 
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", site, shore) %>%
  mutate(
    date=as.Date(datetime),
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/11.2)) 


HRflow_data_BW <- readNWISuv(siteNumbers = siteNo_BW, parameterCd = pCode_flow, startDate = start.date, endDate = end.date) %>%
  mutate(site = "BWL", shore = "BW") %>%
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", site, shore) %>%
  mutate(
    date=as.Date(datetime),
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/11.2)) 

# Combine data for both sites
HRflow_data <- rbind(HRflow_data_GB, HRflow_data_BW) 

HRflow_sum <- HRflow_data %>%
  group_by(site, date) %>%
  summarise(
    meanflow = mean(scale_Q, na.rm = TRUE),
    sdflow = sd(scale_Q, na.rm = TRUE),
    cvflow = sdflow / meanflow * 100  # CV as a percentage
  )

######### can edit
### 
PRISM20 <- read.csv("/Users/kellyloria/Downloads/PRISM_ppt_tmin_tmean_tmax_vpdmin_vpdmax_stable_4km_20201001_20210901.csv", skip=10)
PRISM21 <- read.csv("/Users/kellyloria/Downloads/PRISM_ppt_tmin_tmean_tmax_vpdmin_vpdmax_stable_4km_20211001_20220901.csv", skip=10)
PRISM22 <- read.csv("/Users/kellyloria/Downloads/PRISM_ppt_tmin_tmean_tmax_vpdmin_vpdmax_provisional_4km_20221001_20230901.csv", skip=10)

PRISMdata <- rbind(PRISM20,PRISM21, PRISM22)
str(PRISMdata)
PRISMdata <- PRISMdata %>%
  mutate(date=as.Date(Date, format="%Y-%m-%d"))

# need to find actual precpi for this 
unique(PRISMdata$Name)












########
########
#######
# BW 2021: all, 2022: NS1 and 2
BW1sum <- rbind(BWNS1, BWNS1_22)
BWs1NEP<- BW1sum%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
# long to wide

#library(reshape2)
BW1NEP <- reshape(data=BWs1NEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")


BW2sumNEP<- BWNS2%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
BW2NEP <- reshape(data=BW2sumNEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")

BW3sumNEP <- rbind(BWNS3, BWNS3_22)
BW3sumNEP<- BW3sumNEP%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
BW3NEP <- reshape(data=BW3sumNEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")


# BW 2021: all, 2022: NS1 and 2
BW1sum <- rbind(BWNS1, BWNS1_22)
BWs1NEP<- BW1sum%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
# long to wide

#library(reshape2)
BW1NEP <- reshape(data=BWs1NEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")


BW2sumNEP<- BWNS2%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
BW2NEP <- reshape(data=BW2sumNEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")


GB1sumNEP<- GBNS1%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
GB1NEP <- reshape(data=GB1sumNEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")

GBNS3_nep <- rbind(GBNS3, GBNS3_22)

GB3sumNEP<- GBNS3_nep%>%
  subset(name=="NEP") %>%
  subset(middle <= 40 & middle >= -40)
GB3NEP <- reshape(data=GB3sumNEP, idvar="date",
                  v.names = "middle",
                  timevar = "name",
                  direction="wide")

# long to wide



se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

BW_NEP <- rbind(BW1NEP, BW2NEP,BW3NEP)
BW_NEP$shore <-"west"

# subset for oct 01
BW_NEP21<- subset(BW_NEP, date>as.Date("2022-03-01"))

mean(BW_NEP21$middle.NEP)
se(BW_NEP21$middle.NEP)


GB_NEP <- rbind(GB1NEP, GB3NEP)
GB_NEP$shore <- "east"
# subset for oct 01
GB_NEP21<- subset(GB_NEP, date>as.Date("2021-03-01"))
mean(GB_NEP21$middle.NEP)
se(GB_NEP21$middle.NEP)





sumNEP <-rbind(BW_NEP,GB_NEP)


sumNEPplot <- ggplot(data = sumNEP %>% drop_na(year),aes(date, middle.NEP, color = shore))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = shore),
              linetype = 0, alpha = 0.2)+
  geom_line()+ geom_point(size= 3, alpha = 0.6)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
             color = "#4c4d4c") +
  scale_color_manual(values=c("#a67d17", "#3283a8")) +
  scale_fill_manual(values = c("#a67d17", "#3283a8")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  #ylim(-15, 15) +
  scale_y_continuous(breaks = seq(-15, 24, by = 3))+
  labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust=1), # new
        legend.position = 'right', 
        legend.direction = "vertical")

# ggsave(plot = sumNEPplot, filename = paste("./figures/LA_NEP_all.png",sep=""),width=8.5,height=6.5,dpi=300)



p2 <- ggplot(data = GB_NEP %>% drop_na(year),aes(date, middle.NEP, color = site))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site),
              linetype = 0, alpha = 0.2)+
  geom_line()+ geom_point(size= 3, alpha = 0.6)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
             color = "#4c4d4c") +
  #geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
  scale_color_manual(values = c("#982649", "#003F91", "#333333")) +
  scale_fill_manual(values = c("#982649","#003F91","black")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  ylim(-35, 35) +
  labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust=1), # new
        legend.position = 'bottom', 
        legend.direction = "horizontal") +
  facet_grid((site~.))
p2


max(BW_NEP$middle.NEP)





p2 <- ggplot(data = BW_NEP %>% drop_na(year),aes(date, middle.NEP, color = site))+
  geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = site),
              linetype = 0, alpha = 0.2)+
  geom_line()+ geom_point(size= 3, alpha = 0.6)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")),
             color = "#4c4d4c") +
  #geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
  scale_color_manual(values = c("#982649", "#003F91", "#333333")) +
  scale_fill_manual(values = c("#982649","#003F91","black")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y")+ #new
  ylim(-35, 35) +
  labs(y=expression(mmol~O[2]~m^-3~d^-1)) + 
  theme_classic() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60, hjust=1), # new
        legend.position = 'bottom', 
        legend.direction = "horizontal") +
  facet_grid((site~.))
p2

# LAKE NEP drivers

# streamflow?
# DOdat from 2022Figures 

## FLOW
# library(dataRetrieval)
## Add in flow data
siteNo <- "10336730"
pCode <- c("00060")
start.date <- "2021-06-01"
end.date <- "2022-10-01"

GBflow <- readNWISdata(siteNumbers = siteNo,
                       parameterCd = pCode,
                       startDate = start.date,
                       endDate = end.date)

GBflow.ts <- GBflow %>% select("dateTime", "X_00060_00003") %>% 
  dplyr::rename(datetime = "dateTime", dischargeCFS = "X_00060_00003") %>%
  select("datetime", "dischargeCFS")

GBflow.ts$Site<- "GB"
GBflow.ts$shore<-"east"


siteNo <- "10336660"
pCode <- c("00060")
start.date <- "2021-06-01"
end.date <- "2022-10-01"


flow <- readNWISdata(siteNumbers = siteNo,
                     parameterCd = pCode,
                     startDate = start.date,
                     endDate = end.date)

BWflow.ts <- flow %>% select("dateTime", "X_00060_00003") %>% 
  dplyr::rename(datetime = "dateTime", dischargeCFS = "X_00060_00003") %>%
  select("datetime", "dischargeCFS")

BWflow.ts$Site<- "BW"
BWflow.ts$shore<-"west"

flow<- rbind(BWflow.ts, GBflow.ts)
str(flow)
sumNEPF <- left_join(sumNEP, flow,  by=c('date'='datetime', 'shore'='shore'))

library(lmerTest)
library(lme4)
library(MuMIn)
hist(sumNEPF$dischargeCFS)
hist(log(sumNEPF$dischargeCFS +1))
hist((sumNEPF$middle.NEP))


mod1<- glm(middle.NEP~ scale(dischargeCFS)+ (shore), data=sumNEPF)
summary(mod1)

mod2<- lmer(middle.NEP~ scale(dischargeCFS)+ (shore) + (1|site), data=sumNEPF)
summary(mod2)


mod3<- lmer(middle.NEP~ scale(dischargeCFS)+ (1|site), data=sumNEPF)
summary(mod3)

hist(residuals(mod2))
r.squaredGLMM(mod2)





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




# stream SPC
sumNEPF <- left_join(sumNEP, flow,  by=c('date'='datetime', 'shore'='shore'))



GBL22_SPC<- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/CleanDat/22_GBL_SPC.csv")
GBL22_SPC$solar.time <- as.POSIXct(GBL22_SPC$datetime,
                                   format = "%Y-%m-%d %H:%M:%S",
                                   tz = "UTC")
GBL22_SPC$site<- "GBL"
GBL22_SPC$shore<- "east"

BWL22_SPC<- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/CleanDat/22_BWL_SPC.csv")
BWL22_SPC$solar.time <- as.POSIXct(BWL22_SPC$datetime,
                                   format = "%Y-%m-%d %H:%M:%S",
                                   tz = "UTC")
BWL22_SPC$site<- "BWL"
BWL22_SPC$shore<- "west"


SPC_dat<- rbind(BWL22_SPC, GBL22_SPC)
SPC_datD <- SPC_dat %>% 
  mutate(date= as.Date(datetime)) %>%
  group_by(date, site, shore) %>% 
  summarise(SPCmean = mean(SPC,na.rm = T),
            SPCstd = sd(SPC, na.rm = T))

sumNEPFSPC <- left_join(sumNEPF, SPC_datD[,c("date", "shore", "SPCmean", "SPCstd")],  by=c('date'='date', 'shore'='shore'))



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
