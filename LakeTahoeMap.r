## Urban streams map

library(tidyverse)
library(ggrepel)
library(ggmap)
library(ggsci)
library(scales)

library(dplyr)
library(ggplot2)
library(ggrepel)

#devtools::install_github("briatte/tidykml")
library(tidykml)
library(sf)

# get ploygons for watersheds from .kml map 
kml_data <- st_read("/Users/kellyloria/Desktop/TempExtra/Tahoe_Soils_and_Hydro_Data.kml")

glimpse(kml_data)
#filter for blackwood
blackwood_creek <- kml_data %>%
  filter(Name == "BLACKWOOD CREEK")

GB_creek <- kml_data %>%
  filter(Name == "GLENBROOK CREEK")

# Plot BLACKWOOD CREEK
ggplot() +
  geom_sf(data = blackwood_creek)

ggplot() +
  geom_sf(data = GB_creek)

# read in site cordinates
dat <- read.csv("/Users/kellyloria/Downloads/TahoeCords2.csv", header = T, sep = ",")

range(dat$long)
range(dat$lat)

# Define the bounding box coordinates
bbox <- c(left = -120.3, bottom = 38.87, right = -119.81, top = 39.25)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = bbox, zoom = 12, maptype = 'stamen_terrain') 
# need api run for the first time

# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

population_map <- ggmap(map) + 
  geom_point(data = dat, aes(x = long, y = lat, color = site), size = 2, shape = 21) +
  # geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
  #                  color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#3283a8", "#3283a8", "#3283a8", "#3283a8",
                                "#a67d17", "#a67d17", "#a67d17", "#a67d17"),
                    guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

# Extracting coordinates from blackwood_creek
coords <- st_coordinates(blackwood_creek)

# Creating a data frame from the extracted coordinates
blackwood_creek_df <- as.data.frame(coords)[,c(1,2)]
names(blackwood_creek_df) <- c("long", "lat")  # Rename columns
head(blackwood_creek_df)

# Extracting coordinates from blackwood_creek
coords <- st_coordinates(GB_creek)

# Creating a data frame from the extracted coordinates
GB_df <- as.data.frame(coords)[,c(1,2)]
names(GB_df) <- c("long", "lat")  # Rename columns
head(GB_df)

# Plotting the ggmap object with the Watershed polygons
newmap <- population_map +
  geom_polygon(data = blackwood_creek_df, aes(x = long, y = lat), fill = "transparent", color = "#3283a8")+
  geom_polygon(data = GB_df, aes(x = long, y = lat), fill = "transparent", color = "goldenrod")

# ggsave(plot = newmap, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_MAPbig.png",sep=""),width=8,height=6,dpi=300)




# Define the bounding box coordinates
wbbox <- c(left = -120.18, bottom = 39.102, right = -120.145, top = 39.1155)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = wbbox,  zoom = 15, maptype = 'stamen_terrain') 


# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

dat_w <-dat%>%
  filter(site=="BWL"|site=="BWNS1"| site=="BWNS2"|site=="BWNS3")

population_map <- ggmap(map) + 
  geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
   geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                    color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#3283a8", "#3283a8", "#3283a8", "#3283a8"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_BWzoom.png",sep=""),width=8,height=6,dpi=300)





# Define the bounding box coordinates
ebbox <- c(left = -119.949, bottom = 39.0847, right = -119.927, top = 39.0925)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = ebbox,  zoom = 15, maptype = 'stamen_terrain') 


# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

dat_w <-dat%>%
  filter(site=="GBL"|site=="GBNS1"| site=="GBNS2"|site=="GBNS3")

population_map <- ggmap(map) + 
  geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
  geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                   color = "black", fontface = "bold", size = 3.5) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#a67d17", "#a67d17", "#a67d17", "#a67d17"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_GBzoom.png",sep=""),width=8,height=6,dpi=300)
