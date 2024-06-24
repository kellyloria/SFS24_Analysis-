#' Script for aggregating data for SFS 2024 conference 
#' 
#' @description creating dataframe for nearshore metabolism output  
#' @param 
#' 
#' @return  
#' @export 
#' 

##==============================================================================
## Created  02/01/2024 by KAL
## modified 05/08/2024
#===============================================================================

library(tidyverse)
library(lubridate)

# plotting packages:
library(ggplot2)
library(reshape2)

# stats:
se <- function(dat){
  se <- sd(dat)/sqrt(length(dat))
  return(se)}

##===========================================
## read in lake metabolism estimates
#============================================
# * consider deleting: /Users/kellyloria/Documents/LittoralMetabModeling/plotDat/


add_site_column <- function(file_name) {
  df <- read.csv(file_name)
  # Extract the file name from the full file path
  file_name_only <- basename(file_name)
  # Extract information from the file name:
  site_info <- str_match(file_name_only, "^(.*?)__daily_5ms_fourthlake_Offset_")
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
folder_path <- "/Users/kellyloria/Documents/LittoralMetabModeling/24_plotdata"
# Get a list of CSV files in the folder
csv_files <- list.files(folder_path, pattern = ".csv$", full.names = TRUE)
# Apply the function to each CSV file and store the results in a list
processed_data_list <- lapply(csv_files, add_site_column)
# Combine the list of data frames into a single data frame
LM_data <- do.call(rbind, processed_data_list)

unique(LM_data$site)

## subset data:
LM_data <-LM_data %>% 
  subset(name=="ER" | name=="GPP"|name=="NEP") %>%
  subset(middle <= 50 & middle >= -50)

# add in column for date:
LM_data$origin <- as.Date(paste0(LM_data$year, "-01-01"),) 
LM_data$date <-as.Date(LM_data$yday, origin = LM_data$origin) 


# get rid of NA's
LM_data1 <- LM_data %>%
  dplyr::group_by(date, site, name) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

summary(LM_data)

LM_data_na <- LM_data %>% drop_na(date)
summary(LM_data_na)

#### Convert data from long to wide
str(LM_data)
LM_data_wide <- LM_data_na %>%
  pivot_wider(
    id_cols = c(site, date),
    names_from = name,
    values_from = c(middle, upper, lower)) %>%
  mutate(
    shore = case_when(
      site %in% c("BWNS1","BWNS2", "BWNS3") ~ "BW",
      site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH",
      site %in% c("SSNS1","SSNS2", "SSNS3") ~ "SS",
      site %in% c("GBNS1","GBNS2","GBNS3") ~ "GB",
      TRUE ~ as.character(site)))

##===========================================
lake_DO_temp <- read_rds("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NS_miniDOT/24_NS_flited_dailyDO.rds") %>%
  dplyr::rename(site="Site", lake_tempC = "Temperature_deg_C", lake_DO = "Dissolved_O_mg_L") %>%
  filter(site=="BWNS1" | site=="BWNS2" |site=="BWNS3" |
           site=="GBNS1" | site=="GBNS2" |site=="GBNS3"|
           site=="SHNS1" | site=="SHNS2" |site=="SHNS3"|
           site=="SSNS1" | site=="SSNS2" |site=="SSNS3")

unique(lake_DO_temp$site)

lake_SPC <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/SPC/Lake_SPC_dat.rds")%>%
  group_by(date, site)%>%
  summarise(lake_SPC = mean(SPC_full_uscm, na.rm =T),
            lake_SPC_cv= sd(SPC_full_uscm, na.rm =T)/mean(SPC_full_uscm, na.rm =T),
            lake_SPCtemp_cv=sd(wtr_C, na.rm=T)/mean(wtr_C, na.rm=T),
            lake_SPCtemp =mean(wtr_C, na.rm=T)) 

unique(lake_SPC$site)

# Depth:
depth_dat <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/RBR\ profiles/24NS_depth_dat.rds") %>%
  dplyr::rename(real_NS_depth = sensor_depth, site = Site) %>%
  group_by(date, site) %>%
  summarize(
    real_NS_depth = mean(real_NS_depth, na.rm=T))

unique(depth_dat$site)


# Benthic light: 
PAR_dat <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/benthic_light/PAR_calc_dat.rds") %>%
  dplyr::select(shore, date, Kd_fill, in_par, sensor_depth, par_int_3m) %>%
  group_by(date, shore) %>%
  summarize(
    Kd_fill = mean(Kd_fill, na.rm=T),  #
    in_par = mean(in_par, na.rm=T),    # 
    par_int_3m = mean(par_int_3m, na.rm=T))

unique(PAR_dat$shore)




lake_df <- LM_data_wide %>%
  left_join(lake_DO_temp, by = c("site", "date")) %>%
  left_join(depth_dat, by = c("site", "date")) 

lake_df1 <- lake_df %>% 
  left_join(lake_SPC, by = c("site", "date")) 


lake_df2 <- lake_df1 %>% 
  left_join(PAR_dat, by = c("shore", "date")) 
  

unique(lake_df$shore)
unique(lake_df2$site)

summary(lake_df2)

##===========================================
## Bring in stream covariate data
#============================================
# precip, streamflow, stream temp, stream spc, chem, stream metab
library(dataRetrieval)
# Streamflow
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
  mutate(site = "GBL", shore = "GB", location = "stream") %>% 
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", site, shore) %>%
  mutate(
    date=as.Date(datetime),
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/11.2)) 


HRflow_data_BW <- readNWISuv(siteNumbers = siteNo_BW, parameterCd = pCode_flow, startDate = start.date, endDate = end.date) %>%
  mutate(site = "BWL", shore = "BW", location = "stream") %>%
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", site, shore) %>%
  mutate(
    date=as.Date(datetime),
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/11.2)) 

# Combine data for both sites
HRflow_data <- rbind(HRflow_data_GB, HRflow_data_BW) 

HRflow_sum <- HRflow_data %>%
  group_by(site, shore, date) %>%
  summarise(
    flow_mean = mean(scale_Q, na.rm = TRUE),
    flow_sd = sd(scale_Q, na.rm = TRUE),
    flow_cv = flow_sd / flow_mean  # CV 
  )

#============================================
# stream metab water quality:
BWL_dat <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/FinalInputs/24_BWL_modelInputs.rds") %>%
  mutate(date= as.Date(solar.time))%>%
  group_by(date)%>%
  summarise(stream_DO = mean(DO.obs, na.rm =T),
            stream_temp =mean(temp.water, na.rm=T),
            stream_temp_cv =(sd(temp.water, na.rm=T)/mean(temp.water, na.rm=T))) %>%
  mutate(site = "BWL", shore = "BW") 

GBL_dat <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/FinalInputs/24_GBL_modelInputs.rds") %>%
  mutate(date= as.Date(solar.time))%>%
  group_by(date)%>%
  summarise(stream_DO = mean(DO.obs, na.rm =T),
            stream_temp =mean(temp.water, na.rm=T),
            stream_temp_cv =(sd(temp.water, na.rm=T)/mean(temp.water, na.rm=T))) %>%
  mutate(site = "GBL", shore = "GB") 
  
stream_wq <- rbind(BWL_dat, GBL_dat)


# bring in SPC data
BWL_spc<- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWL_SPC.rds")%>%
  mutate(date= as.Date(datetime))%>% group_by(date)%>%
  summarise(stream_SPC = mean(SPC, na.rm =T),
            stream_SPC_cv = (sd(SPC, na.rm =T)/mean(SPC, na.rm =T)),
            Stream_temp_SPC =mean(wtr, na.rm=T)) %>%
  mutate(site = "BWL", shore = "BW") 


GBL_spc<- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_GBL_SPC.rds") %>%
  mutate(date= as.Date(datetime))%>% group_by(date)%>%
  summarise(stream_SPC = mean(SPC, na.rm =T),
            stream_SPC_cv = (sd(SPC, na.rm =T)/mean(SPC, na.rm =T)),
            Stream_temp_SPC =mean(wtr, na.rm=T)) %>%
  mutate(site = "GBL", shore = "GB") 

stream_spc <- rbind(GBL_spc, BWL_spc)

stream_df <- HRflow_sum %>%
  left_join(stream_spc, by = c("site", "shore", "date")) %>% 
  left_join(stream_wq, by = c("site", "shore", "date"))
  

##===========================================
## Lake covariates join stream data
#============================================
## Merge climate datasets and change some grouping variables 
PRISM <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/PRISM_precip/PRISM_Precip_dat.csv") %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  filter(site %in% c("BWNS1", "BWNS2", "BWNS3", "SSNS1", "SSNS2", "SSNS3", "SHNS1", "SHNS2", "SHNS3", "GBNS1", "GBNS2", "GBNS3")) %>%
  mutate(
    shore = case_when(
      site %in% c("BWNS1","BWNS2", "BWNS3", "BWL") ~ "BW",
      site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH",
      site %in% c("SSNS1","SSNS2", "SSNS3") ~ "SS",
      site %in% c("GBNS1","GBNS2","GBNS3", "GBL") ~ "GB",
      TRUE ~ as.character(site)
    )
  )
unique(PRISM$site)

## NLDAS data
light_dat <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/processed_light/nearshore_NLDAS_light.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site, date) %>%
  summarise(
    light_mean = mean(light, na.rm = TRUE),
    light_sd = sd(light, na.rm = TRUE),
    light_cv = light_sd / light_mean)

unique(light_dat$site)

wind_dat <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/processed_windsp/nearshore_NLDAS_windsp.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site, date) %>%
  summarise(
    windsp_mean = mean(windsp_ms, na.rm = TRUE),
    windsp_sd = sd(windsp_ms, na.rm = TRUE),
    windsp_cv = windsp_sd / windsp_mean)


unique(wind_dat$site)


##====
## Merge climate datasets and change some grouping variables 
## final data selection: ##
clim_dat <- light_dat %>%
  left_join(wind_dat, by = c("site", "date")) %>%
  mutate(
    shore = case_when( # create broad variable to lineup climate and DO dat
      site == "BWNS2" ~ "BW", # dat is ~4km resolution so called it from NLDAS based on center miniDOT in each array
      site == "SHNS2" ~ "SH",
      site == "SSNS2" ~ "SS",
      site == "GBNS2" ~ "GB",
      TRUE ~ as.character(site)))

clim_dat1 <- clim_dat[,c(-1)]
clim_sum1 <- PRISM %>%
  left_join(clim_dat1, by = c("shore", "date")) 

unique(clim_sum1$site)


#############################
#### merge back with lake! 
# lake only dataframe
NS_df <- lake_df2 %>%
  left_join(clim_sum1, by = c("site", "date", "shore")) %>%
  left_join(stream_df[,c(-1)], by = c("shore", "date")) 


names(NS_df)

# trim data 
SFS_datQ <- NS_df %>%
  dplyr::select(-Latitude, -Longitude, -X, -elevation_m, -miniDO_sat)
names(SFS_datQ)

hist(SFS_datQ$middle_GPP)
hist(log(SFS_datQ$middle_GPP+1))
hist(SFS_datQ$middle_ER)

summary(SFS_datQ)

SFS_datQ1 <- SFS_datQ %>%
  mutate(
    position = case_when(
      site %in% c("BWNS1","GBNS1", "SSNS1", "SHNS1") ~ "north",
      site %in% c("BWNS2","GBNS2", "SSNS2", "SHNS2" ) ~ "center",
      site %in% c("BWNS3","GBNS3", "SHNS3") ~ "south",
      TRUE ~ as.character(site)
    )
  )


SFS_datQ2 <- SFS_datQ1%>%
  filter(shore=="BW" | shore=="GB")



SFS_datQ1$precip_bi <- ifelse(SFS_datQ1$ppt_mm > 0, 1, 0)
SFS_datQ1$log_ppmt <- log(SFS_datQ1$ppt_mm+1)
SFS_datQ1$log_streamflow <- (log(SFS_datQ1$flow_mean)+1)

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
  mutate(year=year(date), 
           week= week(date))%>%
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



GPP_plot <- ggplot(data = SFS_week_dat%>%filter(middle_GPP<30 & shore=="BW"|shore=="GB"), aes(week, middle_GPP, color = shore, shape=position))+
  geom_ribbon(aes(ymin = lower_GPP, ymax = upper_GPP, fill = shore),
              linetype = 0, alpha = 0.1)+
  #geom_point(size= 2.5, alpha = 0.7)+
  geom_line(alpha = 0.7, size=1)+
  #scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+ #new
  ylim(0, 40) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP"~mmol~O[2]~m^-3~d^-1), x=NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16)) + facet_grid(.~year)





GPP_plot <- ggplot(data = SFS_datQ2%>%filter(middle_GPP<30 & shore=="BW"|shore=="GB"), aes(week, middle_GPP, color = shore, shape=position))+
  geom_ribbon(aes(ymin = lower_GPP, ymax = upper_GPP, fill = shore),
              linetype = 0, alpha = 0.1)+
  #geom_point(size= 2.5, alpha = 0.7)+
  #scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+ #new
  ylim(0, 40) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP"~mmol~O[2]~m^-3~d^-1), x=NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16)) + facet_grid(.~year)



ER_plot <- ggplot(data = SFS_week_dat%>%filter(middle_ER>-30 & shore=="BW"|shore=="GB"), aes(week, middle_ER, color = shore, shape=position))+
  geom_ribbon(aes(ymin = lower_ER, ymax = upper_ER, fill = shore),
              linetype = 0, alpha = 0.1)+
  geom_line(alpha = 0.7, size=1)+
 # geom_point(size= 2.5, alpha = 0.7)+
  # scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+ #new
  ylim(-40, 0) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER"~mmol~O[2]~m^-3~d^-1), x=NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))+facet_grid(.~year)


# library(ggpubr)


metab21 <- ggarrange(GPP_plot,
                     ER_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     labels = c("A", "B"),
                     label.x = c(0.07,0.07),
                     label.y = c(0.95,0.95),
                     legend = "bottom")

# ggsave(plot = metab21, filename = paste("./figures/NS24_metab_pre_v2.png",sep=""),width=14.5,height=6.5,dpi=300)






## save dat: 
# write.csv (x = SFS_datQ, file = "./SFS24_analysis_dat/SFS24_analysis_dat_v2.csv", row.names = TRUE)

# saveRDS(SFS_datQ, file = "./SFS24_analysis_dat/SFS24_analysis_dat.rds")
library(dplyr)
library(lubridate)

# Assuming SFS_week_dat has columns 'year' and 'week'

# Step 1: Convert week and year to date (start of the week)


SFS_week_dat <- SFS_week_dat %>%
  mutate(date = as.Date(paste(year, week, 1, sep = "-"), "%Y-%U-%u"))

# Step 2: Create a time grouping variable
SFS_week_dat <- SFS_week_dat %>%
  arrange(date) %>%
  mutate(time_group = row_number())

# Step 3: Filter the data based on a date range
start_date <- as.Date("2021-06-01")
end_date <- as.Date("2023-09-14")

filtered_data <- SFS_week_dat %>%
  filter(date > start_date & date < end_date)

# Print the filtered data
print(filtered_data)

SFS_week_dat_df <- as.data.frame(SFS_week_dat)
str(SFS_week_dat_df)



# Create the sequence of dates
start_date <- as.Date("2021-06-01")
end_date <- as.Date("2023-09-14")
date_sequence <- seq.Date(from = start_date, to = end_date, by = "day")

# Extract unique sites from your dataframe
unique_sites <- unique(SFS_week_dat_df$site)

# Create a dataframe with all combinations of sites and dates
site_date_df <- expand.grid(site = unique_sites, date = date_sequence) %>%
  mutate(
    shore = case_when(
      site %in% c("BWNS1","BWNS2", "BWNS3") ~ "BW",
      site %in% c("SHNS1", "SHNS2", "SHNS3") ~ "SH",
      site %in% c("SSNS1","SSNS2", "SSNS3") ~ "SS",
      site %in% c("GBNS1","GBNS2","GBNS3") ~ "GB",
      TRUE ~ as.character(site)))%>%
  select(date, shore)





# Merge the new dataframe with the original dataframe
SFS_week_dat_full <- merge(site_date_df, SFS_week_dat_df, by = c("date", "shore"), all.x = TRUE)

NA_dates <- SFS_week_dat_full %>%
  filter(is.na(site))


GPP_plot <- ggplot(data = SFS_week_dat_full%>%filter(middle_GPP<25), aes(date, middle_GPP, color = shore, shape=position))+
  geom_ribbon(aes(ymin = lower_GPP, ymax = upper_GPP, fill = shore),
              linetype = 0, alpha = 0.1)+
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+ #new
  ylim(0, 40) +
  geom_point(size= 1.75, alpha = 0.7)+
  geom_line(alpha = 0.7, size=1)+
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("GPP"~(mmol~O[2]~m^-3~d^-1)), x=NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16)) #+ facet_grid(.~year)



ER_plot <- ggplot(data = SFS_week_dat%>%filter(middle_ER>-25), aes(date, middle_ER, color = shore, shape=position))+
  geom_ribbon(aes(ymin = lower_ER, ymax = upper_ER, fill = shore),
              linetype = 0, alpha = 0.1)+
  geom_line(alpha = 0.7, size=1)+
  geom_point(size= 1.75, alpha = 0.7)+
   scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+ #new
  ylim(-40, 0) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  labs(y=expression("ER"~(mmol~O[2]~m^-3~d^-1)), x=NULL) + 
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))#+facet_grid(.~year)




metab21 <- ggarrange(GPP_plot,
                     ER_plot,
                     ncol = 1, nrow = 2,
                     common.legend = TRUE, 
                     labels = c("a", "b"),
                     label.x = c(0.07,0.07),
                     label.y = c(0.95,0.95),
                     legend = "bottom")

# ggsave(plot = metab21, filename = paste("./figures/NS24_metab_pre_v5.png",sep=""),width=14.5,height=7.25,dpi=300)



library(ggplot2)

# Existing dataframe for non-NA rows
non_NA_dates <- SFS_week_dat_full %>%
  filter(!is.na(site))

# Create the plot
GPP_plot <- ggplot() +
  # Plot for non-NA data
  geom_ribbon(data = non_NA_dates %>% filter(shore == "BW" | shore == "GB"), 
              aes(x = date, ymin = lower_GPP, ymax = upper_GPP, fill = shore),
              linetype = 0, alpha = 0.1) +
  geom_line(data = non_NA_dates %>% filter(shore == "BW" | shore == "GB"), 
            aes(x = date, y = middle_GPP, color = shore, shape = position), 
            alpha = 0.7, size = 1) +
  
  # Plot for NA data
  geom_line(data = NA_dates %>% filter(shore == "BW" | shore == "GB"), 
            aes(x = date, y = middle_GPP), 
            color = "white", size = 1, alpha = 0.7) +
  
  # Axes and scales
  scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
  ylim(0, 40) +
  scale_colour_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  scale_fill_manual(values = c(SS = "#136F63", BW = "#3283a8", GB = "#a67d17", SH = "#c76640")) +
  
  # Labels and theme
  labs(y = expression("GPP" ~ mmol ~ O[2] ~ m^-3 ~ d^-1), x = NULL) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# Print the plot
print(GPP_plot)




