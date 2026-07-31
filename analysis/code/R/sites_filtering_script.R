####################################################################
# This script is for filtering sites                               #
#                                                                  #                                                   
# Author: Mathias Bellat and Biel Soriano-Elias                    #
# Affiliation : Tubingen University                                #
# Creation date : 17/09/2025                                       #
# E-mail:  mathias.archaeology@gmail.com & bielsoel28@gmail.com    #
####################################################################


# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Set folder direction

# Clean up workspace
rm(list=setdiff(ls()))

# 0.2 Install packages =========================================================

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(sf, mapview, readr, terra)

# 0.3 Show session infos =======================================================

sessionInfo()

# 1 Prepare the data ###########################################################

# 1.1 Import files =============================================================

Rhine_Bronze <- read_delim("./analysis/raw_data/inputs/Rhine_Bronze.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
Rhine_Iron <- read_delim("./analysis/raw_data/inputs/Rhine_Iron.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

Lazio_Bronze <- read.csv("./analysis/raw_data/inputs/Lazio_Bronze.csv")
Lazio_Iron <- read.csv("./analysis/raw_data/inputs/Lazio_Iron.csv")

Shephelah_Bronze <- read_delim("./analysis/raw_data/inputs/Shephelah_Bronze.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
Shephelah_Iron <- read_delim("./analysis/raw_data/inputs/Shephelah_Iron.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

Kurdistan_Bronze <- read.csv("./analysis/raw_data/inputs/Kurdistan_Bronze.csv")
Kurdistan_Iron <- read.csv("./analysis/raw_data/inputs/Kurdistan_Iron.csv")

# 1.2 Filter the files =========================================================

# For Rhine data
data_df <- Rhine_Bronze[Rhine_Bronze$LATITUDE < 48.5,]
data_df <- data_df[data_df$CARAC_NAME == "Immobilier",]
data_df <- data_df[data_df$CARAC_LVL1 == "Habitat",]
Rhine_Bronze_sf <- st_as_sf(data_df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) 

data_df <- Rhine_Iron[Rhine_Iron$LATITUDE < 48.5,]
data_df <- data_df[data_df$CARAC_NAME == "Immobilier",]
data_df <- data_df[data_df$CARAC_LVL1 == "Habitat",]
Rhine_Iron_sf <- st_as_sf(data_df, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) 

# For Lazio data
Lazio_Bronze <- Lazio_Bronze[Lazio_Bronze$Longitude < 12.83,]
Lazio_Bronze_sf <- st_as_sf(Lazio_Bronze, coords = c("Longitude", "Latitude"), crs = 4326)

Lazio_Iron <- Lazio_Iron[Lazio_Iron$Longitude < 12.83,]
Lazio_Iron_sf <- st_as_sf(Lazio_Iron, coords = c("Longitude", "Latitude"), crs = 4326) 

# For Shephelah data
Shephelah_Bronze_sf <- st_as_sf(Shephelah_Bronze, coords = c("Longitude", "Latitude"), crs = 4326)
Shephelah_Iron <- Shephelah_Iron[Shephelah_Iron$Type == "Settlement",]
Shephelah_Iron <- Shephelah_Iron[Shephelah_Iron$Period == "Iron Age II",]
Shephelah_Iron <- Shephelah_Iron[grepl("Farm|hamlet|structure", Shephelah_Iron$Subtype, ignore.case = TRUE),]
Shephelah_Iron_sf <- st_as_sf(Shephelah_Iron, coords = c("Longitude", "Latitude"), crs = 4326) 

# For Kurdistan data
Kurdistan_Bronze <- Kurdistan_Bronze[Kurdistan_Bronze$zone != "C",]
Kurdistan_Bronze <- Kurdistan_Bronze[Kurdistan_Bronze$Latitude > 36.86,]

Kurdistan_Iron <- Kurdistan_Iron[Kurdistan_Iron$zone != "C",]
Kurdistan_Iron <- Kurdistan_Iron[Kurdistan_Iron$Latitude > 36.86,]

Kurdistan_Bronze_sf <- st_as_sf(Kurdistan_Bronze, coords = c("Longitude", "Latitude"), crs = 4326)
Kurdistan_Iron_sf <- st_as_sf(Kurdistan_Iron, coords = c("Longitude", "Latitude"), crs = 4326) 

# 1.3 Plot data ================================================================

mapview(Rhine_Bronze_sf, col.regions = "lightblue") + mapview(Rhine_Iron_sf, col.regions = "pink")

mapview(Lazio_Bronze_sf, col.regions = "lightblue") + mapview(Lazio_Iron_sf, col.regions = "pink")

mapview(Shephelah_Bronze_sf, col.regions = "lightblue") + mapview(Shephelah_Iron_sf, col.regions = "pink")

mapview(Kurdistan_Bronze_sf, col.regions = "lightblue") + mapview(Kurdistan_Iron_sf, col.regions = "pink")

# 2 Compute basic statistics ###################################################

# 2.1 Reuse percent ============================================================

reoccupation <- data.frame(Rhine = 0, Lazio = 0, Shephelah = 0, Kurdistan = 0)

reoccupation[1,1] <- round(((length(intersect(Rhine_Iron_sf$SITE_SOURCE_ID, Rhine_Bronze_sf$SITE_SOURCE_ID)))/length(Rhine_Bronze_sf$SITE_SOURCE_ID))*100, digit = 2)
reoccupation[1,2] <- round(((length(intersect(Lazio_Iron_sf$Id, Lazio_Bronze_sf$Id)))/length(Lazio_Bronze_sf$Id))*100, digit = 2)
reoccupation[1,3] <- round(((length(intersect(Shephelah_Iron_sf$Name, Shephelah_Bronze_sf$Name)))/length(Shephelah_Bronze_sf$Name))*100, digit = 2)
reoccupation[1,4] <- round(((length(intersect(Kurdistan_Iron_sf$Name, Kurdistan_Bronze_sf$Name)))/length(Kurdistan_Bronze_sf$Name))*100, digit = 2)

reoccupation[2,1] <- round(((length(intersect(Rhine_Iron_sf$SITE_SOURCE_ID, Rhine_Bronze_sf$SITE_SOURCE_ID)))/length(Rhine_Iron_sf$SITE_SOURCE_ID))*100, digit = 2)
reoccupation[2,2] <- round(((length(intersect(Lazio_Iron_sf$Id, Lazio_Bronze_sf$Id)))/length(Lazio_Iron_sf$Id))*100, digit = 2)
reoccupation[2,3] <- round(((length(intersect(Shephelah_Iron_sf$Name, Shephelah_Bronze_sf$Name)))/length(Shephelah_Iron_sf$Name))*100, digit = 2)
reoccupation[2,4] <- round(((length(intersect(Kurdistan_Iron_sf$Name, Kurdistan_Bronze_sf$Name)))/length(Kurdistan_Iron_sf$Name))*100, digit = 2)


reoccupation
write.table(reoccupation, "./analysis/derived_data/pre_process/Reoccupation_stats.txt")

# 2.2 Surface of the study area ================================================

# For Rhine
Rhine_Iron_sf$era <- "Iron"
Rhine_Bronze_sf$era <- "Bronze"
Rhine_all <- rbind(Rhine_Iron_sf, Rhine_Bronze_sf)
Rhine_all <- st_transform(Rhine_all, 32632)
bbox <- st_bbox(Rhine_all)
bbox_sf <- st_as_sfc(bbox)
Rhine_bbox <- st_buffer(bbox_sf, 1000)
mapview(Rhine_bbox) + mapview(Rhine_all,  col.regions = "pink")

# For Lazio
Lazio_Iron_sf$era <- "Iron"
Lazio_Bronze_sf$era <- "Bronze"
Lazio_all <- rbind(Lazio_Iron_sf, Lazio_Bronze_sf)
Lazio_all <- st_transform(Lazio_all, 32633)
bbox <- st_bbox(Lazio_all)
bbox_sf <- st_as_sfc(bbox)
Lazio_bbox <- st_buffer(bbox_sf, 1000)
mapview(Lazio_bbox) + mapview(Lazio_all,  col.regions = "pink")

# For Shephelah
Shephelah_Iron_sf$era <- "Iron"
Shephelah_Bronze_sf$era <- "Bronze"
Shephelah_all <- rbind(Shephelah_Iron_sf, Shephelah_Bronze_sf)
Shephelah_all <- st_transform(Shephelah_all, 32636)
bbox <- st_bbox(Shephelah_all)
bbox_sf <- st_as_sfc(bbox)
Shephelah_bbox <- st_buffer(bbox_sf, 1000)
mapview(Shephelah_bbox) + mapview(Shephelah_all,  col.regions = "pink")

# For Kurdistan
Kurdistan_Iron_sf$era <- "Iron"
Kurdistan_Bronze_sf$era <- "Bronze"
Kurdistan_all <- rbind(Kurdistan_Iron_sf, Kurdistan_Bronze_sf)
Kurdistan_all <- st_transform(Kurdistan_all, 32638)
bbox <- st_bbox(Kurdistan_all)
bbox_sf <- st_as_sfc(bbox)
Kurdistan_bbox <- st_buffer(bbox_sf, 1000)
mapview(Kurdistan_bbox) + mapview(Kurdistan_all,  col.regions = "pink")

surface <- data.frame(Rhine = 0, Lazio = 0, Shephelah = 0, Kurdistan = 0)

surface[1,1] <- round(sum(st_area(Rhine_bbox))/1000000, digit = 2)
surface[1,2] <- round(sum(st_area(Lazio_bbox))/1000000, digit = 2)
surface[1,3] <- round(sum(st_area(Shephelah_bbox))/1000000, digit = 2)
surface[1,4] <- round(sum(st_area(Kurdistan_bbox))/1000000, digit = 2)

surface[2,1] <- nrow(Rhine_Bronze_sf)
surface[2,2] <- nrow(Lazio_Bronze_sf)
surface[2,3] <- nrow(Shephelah_Bronze_sf)
surface[2,4] <- nrow(Kurdistan_Bronze_sf)

surface[3,1] <- round(nrow(Rhine_all[Rhine_all$era == "Bronze",])/(sum(st_area(Rhine_bbox))/1000000), digit = 2)
surface[3,2] <- round(nrow(Lazio_all[Lazio_all$era == "Bronze",])/(sum(st_area(Lazio_bbox))/1000000), digit = 2)
surface[3,3] <- round(nrow(Shephelah_all[Shephelah_all$era == "Bronze",])/(sum(st_area(Shephelah_bbox))/1000000), digit = 2)
surface[3,4] <- round(nrow(Kurdistan_all[Kurdistan_all$era == "Bronze",])/(sum(st_area(Kurdistan_bbox))/1000000), digit = 2)

surface[4,1] <- nrow(Rhine_Iron_sf)
surface[4,2] <- nrow(Lazio_Iron_sf)
surface[4,3] <- nrow(Shephelah_Iron_sf)
surface[4,4] <- nrow(Kurdistan_Iron_sf)

surface[5,1] <- round(nrow(Rhine_all[Rhine_all$era == "Iron",])/(sum(st_area(Rhine_bbox))/1000000), digit = 2)
surface[5,2] <- round(nrow(Lazio_all[Lazio_all$era == "Iron",])/(sum(st_area(Lazio_bbox))/1000000), digit = 2)
surface[5,3] <- round(nrow(Shephelah_all[Shephelah_all$era == "Iron",])/(sum(st_area(Shephelah_bbox))/1000000), digit = 2)
surface[5,4] <- round(nrow(Kurdistan_all[Kurdistan_all$era == "Iron",])/(sum(st_area(Kurdistan_bbox))/1000000), digit = 2)

surface[6,1] <- nrow(Rhine_all) - length(intersect(Rhine_Iron_sf$SITE_SOURCE_ID, Rhine_Bronze_sf$SITE_SOURCE_ID))
surface[6,2] <- nrow(Lazio_all) - length(intersect(Lazio_Iron_sf$Id, Lazio_Bronze_sf$Id))
surface[6,3] <- nrow(Shephelah_all) - length(intersect(Shephelah_Iron_sf$Name, Shephelah_Bronze_sf$Name))
surface[6,4] <- nrow(Kurdistan_all) - length(intersect(Kurdistan_Iron_sf$Name, Kurdistan_Bronze_sf$Name))

surface[7,1] <- round(surface[6,1]/(sum(st_area(Rhine_bbox))/1000000), digit = 2)
surface[7,2] <- round(surface[6,2]/(sum(st_area(Lazio_bbox))/1000000), digit = 2)
surface[7,3] <- round(surface[6,3]/(sum(st_area(Shephelah_bbox))/1000000), digit = 2)
surface[7,4] <- round(surface[6,4]/(sum(st_area(Kurdistan_bbox))/1000000), digit = 2)

row.names(surface) <- c("Surface of the AOI", "BA sites", "BA sites density (km2)", "IA sites", "IA sites density (km2)", "Total sites", "Total sites density (km2)")
surface
write.table(surface, "./analysis/derived_data/pre_process/Surface_stats.txt")

# 2.3 Export AOI ===============================================================

bbox_sf <- st_as_sfc(st_bbox(Rhine_bbox))
bbox_sf <- st_transform(bbox_sf, 4326)
bbox_vect <- vect(bbox_sf)
r <- rast(ext(bbox_vect), resolution = 50, crs = st_crs(bbox_sf)$wkt)
values(r) <- 1
writeRaster(r, file = "./analysis/raw_data/grid/Rhine_AOI.tif", overwrite = TRUE)

bbox_sf <- st_as_sfc(st_bbox(Lazio_bbox))
bbox_sf <- st_transform(bbox_sf, 4326)
bbox_vect <- vect(bbox_sf)
r <- rast(ext(bbox_vect), resolution = 50, crs = st_crs(bbox_sf)$wkt)
values(r) <- 1
writeRaster(r, file = "./analysis/raw_data/grid/Lazio_AOI.tif", overwrite = TRUE)

bbox_sf <- st_as_sfc(st_bbox(Shephelah_bbox))
bbox_sf <- st_transform(bbox_sf, 4326)
bbox_vect <- vect(bbox_sf)
r <- rast(ext(bbox_vect), resolution = 50, crs = st_crs(bbox_sf)$wkt)
values(r) <- 1
writeRaster(r, file = "./analysis/raw_data/grid/Shephelah_AOI.tif", overwrite = TRUE)

bbox_sf <- st_as_sfc(st_bbox(Kurdistan_bbox))
bbox_sf <- st_transform(bbox_sf, 4326)
bbox_vect <- vect(bbox_sf)
r <- rast(ext(bbox_vect), resolution = 50, crs = st_crs(bbox_sf)$wkt)
values(r) <- 1
writeRaster(r, file = "./analysis/raw_data/grid/Kurdistan_AOI.tif", overwrite = TRUE)

save(Rhine_all, Lazio_all, Shephelah_all, Kurdistan_all, file = "./analysis/derived_data/save/Sites_filter.RData")
