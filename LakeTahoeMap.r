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
bbox <- c(left = -120.300, bottom = 38.900, right = -119.800, top = 39.30)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = bbox, zoom = 12, maptype = 'stamen_toner_lite') 
# need api run for the first time

#library(ggmap)
# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

population_map <- ggmap(map) + 
  geom_point(data = dat, aes(x = long, y = lat, color = site), size = 2, shape = 19) +
  # geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
  #                  color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#4697bd", "#3283a8", "#3283a8", "#3283a8","#4697bd",
                                "#b8902c", "#a67d17", "#a67d17", "#a67d17", "#b8902c",
                                "#c76640", "#c76640", "#c76640",
                                "#136F63", "#136F63", "#136F63"),
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
  geom_polygon(data = blackwood_creek_df, aes(x = long, y = lat), color = "#3283a8", fill = "#3283a8", linewidth = 1.1) #+
  #geom_polygon(data = GB_df, aes(x = long, y = lat), fill = "#a67d17", color = "#a67d17", linewidth = 1.1)

# ggsave(plot = newmap, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/23_MAP_big_BW_W_z13.png",sep=""),width=8,height=6,dpi=300)




# Define the bounding box coordinates
wbbox <- c(left = -120.179, bottom = 39.102, right = -120.145, top = 39.1441)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = wbbox,  zoom = 18, maptype = 'stamen_terrain') 

# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

dat_w <-dat%>%
  filter(site=="BWL"|site=="BWNS1"| site=="BWNS2"|site=="BWNS3"|site=="BWO"|
           site=="SSNS1"| site=="SSNS2"| site=="SSNS3")

population_map <- ggmap(map) + 
  geom_point(data = dat_w, aes(x = long, y = lat, color = site), size = 2) +
   geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                    color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#4697bd", "#3283a8", "#3283a8", "#3283a8","#4697bd",
                                "#136F63", "#136F63", "#136F63"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

## ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_BWzoom.png",sep=""),width=8,height=6,dpi=300)


# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/23_west_big_z13.png",sep=""),width=8,height=6,dpi=300)



# Define the bounding box coordinates
ebbox <- c(left = -119.952, bottom = 39.08, right = -119.93, top = 39.105)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = ebbox,  zoom = 18, maptype = 'stamen_terrain') 


# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map

dat_e <-dat%>%
  filter(site=="GBL"|site=="GBNS1"| site=="GBNS2"|site=="GBNS3"| site=="GBO"|
         site=="SHNS1"|site=="SHNS2"| site=="SHNS3")

population_map <- ggmap(map) + 
  geom_point(data = dat_e, aes(x = long, y = lat, color = site), size = 2) +
  geom_label_repel(data = dat, aes(x = long, y = lat, label = site), 
                   color = "black", fontface = "bold", size = 3.5) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() + 
  scale_color_manual(values = c("#b8902c", "#a67d17", "#a67d17", "#a67d17", "#b8902c",
                                "#c76640", "#c76640", "#c76640"),
                     guide = guide_legend(override.aes = list(label = ""))) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  ) +
  # Add scale bar
  scalebar(data = dat_e, dist = 1, location = "bottomleft", dist_unit = "km", transform = TRUE)


population_map


# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/Ncycle/MSM_ncycle/figures/23_GBzoom.png",sep=""),width=8,height=6,dpi=300)

# ggsave(plot = population_map, filename = paste("/Users/kellyloria/Documents/UNR/MSMmetab/SFS24_Analysis/figures/23_GBzoom.png",sep=""),width=8,height=6,dpi=300)
