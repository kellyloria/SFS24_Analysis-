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
  subset(middle <= 30 & middle >= -30)

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
## Bring in stream covariate data
#============================================
# precip, streamflow, stream temp, stream spc, chem, stream metab

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
    flow_mean = mean(scale_Q, na.rm = TRUE),
    flow_sd = sd(scale_Q, na.rm = TRUE),
    flow_cv = flow_sd / flow_mean  # CV 
  )

# stream SPC
BWL_SPC <- read.csv("./23_CleanDat/23b_BWL_SPC.csv") %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles"),
         date=as.Date(datetime),
         site="GBL") 

GBL_SPC <- read.csv("./23_CleanDat/23b_GBL_SPC.csv") %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles"),
         date=as.Date(datetime),
         site="GBL") 

SPC_data <- rbind(GBL_SPC, BWL_SPC) 

SPC_sum <- SPC_data %>%
  group_by(site, date) %>%
  summarise(
    spc_mean = mean(SPC_full_uscm, na.rm = TRUE),
    spc_sd = sd(SPC_full_uscm, na.rm = TRUE),
    spc_cv = spc_sd / spc_mean  # coefficient of variation
  )

# stream metab - later

## Full join stream data
streamDat <- HRflow_sum %>%
  full_join(SPC_sum) %>%
  mutate(
    shore = case_when(
      site %in% c("BWL") ~ "BW",
      site %in% c("GBL") ~ "GB",
      TRUE ~ as.character(site)
    ))

##===========================================
## Lake covariates join stream data
#============================================
## Merge climate datasets and change some grouping variables 
PRISM <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/PRISM_precip/PRISM_Precip_dat.csv") %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  filter(site %in% c("BWNS1", "BWNS2", "BWNS3", "BWL", "SSNS1", "SSNS2", "SSNS3", "GBL", "SHNS1", "SHNS2", "SHNS3", "GBNS1", "GBNS2", "GBNS3")) %>%
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
  dplyr::select(site, datetime, light)

wind_dat <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/processed_windsp/nearshore_NLDAS_windsp.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  dplyr::select(site, datetime, windsp_ms)

##====
## Merge climate datasets and change some grouping variables 
## final data selection: ##
clim_dat <- light_dat %>%
  left_join(wind_dat, by = c("datetime", "site")) %>%
  filter(datetime > as.POSIXct("2021-01-01 00:00:00"))  %>%
  mutate(
    shore = case_when( # create broad variable to lineup climate and DO dat
      site == "BWNS2" ~ "BW", # dat is ~4km resolution so called it from NLDAS based on center miniDOT in each array
      site == "SHNS2" ~ "SH",
      site == "SSNS2" ~ "SS",
      site == "GBNS2" ~ "GB",
      TRUE ~ as.character(site)),
    date=as.Date(datetime, format="%Y-%m-%d"))

clim_sum <- clim_dat %>%
  group_by(site,shore, date) %>%
  summarise(
    light_mean = mean(light, na.rm = TRUE),
    light_sd = sd(light, na.rm = TRUE),
    light_cv = light_sd / light_mean,
    windsp_mean = mean(windsp_ms, na.rm = TRUE),
    windsp_sd = sd(windsp_ms, na.rm = TRUE),
    windsp_cv = windsp_sd / windsp_mean)

clim_sum2 <- PRISM %>%
  left_join(clim_sum, by = c("shore", "date")) %>%
  filter(date > as.Date("2021-01-01"))  

# lake only dataframe
SFS_dat_ns <- LM_data_wide %>%
  left_join(clim_sum2, by = c("site"="site.x", "date"))

# stream dat
SFS_dat <- SFS_dat_ns %>%
  left_join(streamDat, by = c("shore"="shore", "date"))

# trim data 
SFS_datQ <- SFS_dat %>%
  dplyr::select(site.x, shore, date, middle_GPP, middle_ER, middle_NEP,
                upper_GPP, upper_ER, upper_NEP, lower_GPP,
                lower_ER, lower_NEP,ppt_mm,tmin_C,tmean_C,
                tmax_C,vpdmin_hPa,vpdmax_hPa, light_mean,light_cv,windsp_cv,
                flow_mean,flow_sd,flow_cv, spc_mean, spc_sd, spc_cv
  )

## save dat: 
# write.csv (x = SFS_datQ, file = "./SFS24_analysis_dat/SFS24_analysis_dat.csv", row.names = TRUE)

