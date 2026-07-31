##########################################################################
# This script is for pre-processing the covariates                       #    
#                                                                        #                                              
# Author: Mathias Bellat & Biel Soriano Elias                            #
# Affiliation : Tübingen University & Autonomous University of Barcelona #
# Creation date : 21/10/2025                                             #
# E-mail:  mathias.archaeology@gmail.com & bielsoel28@gmail.com          #
##########################################################################


# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Set folder direction
setwd()

# Clean up workspace
rm(list = ls(all.names = TRUE))

# 0.2 Install packages =========================================================

install.packages("pacman")        #Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(terra, sf, googledrive, mapview, leastcostpath) # Specify required packages and download it if needed

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Import data sets ##########################################################

eras <- c("BA", "IA")
areas <- c("Rhine", "Lazio", "Shephelah", "Kurdistan")
projections <- list(Rhine = "EPSG:32632", Lazio = "EPSG:32633", Shephelah = "EPSG:32636", Kurdistan = "EPSG:32638")
load("./data/derived_data/save/Sites_filter.RData")

# 01.1 Import from Google Drive ================================================
drive_auth()
files <- drive_ls()

for (area in areas) {
  
  # Import LSWI 1613
  file <- files[files$name == paste0(area, "_Sentinel2_LSWI_1613.tif"),]
  drive_download(
    as_id(file$id),
    path = paste0("./data/raw_data/covariates_raw/",area, "_Sentinel2_LSWI_1613_raw.tif"), 
    overwrite = TRUE
  )
  
  # Import LSWI 2202
  file <- files[files$name == paste0(area, "_Sentinel2_LSWI_2202.tif"),]
  drive_download(
    as_id(file$id),
    path = paste0("./data/raw_data/covariates_raw/",area, "_Sentinel2_LSWI_2202_raw.tif"), 
    overwrite = TRUE
  )
  
  # Import NDWI
  file <- files[files$name == paste0(area, "_Sentinel2_NDWI.tif"),]
  drive_download(
    as_id(file$id),
    path = paste0("./data/raw_data/covariates_raw/",area, "_Sentinel2_NDWI_raw.tif"), 
    overwrite = TRUE
  )
  
  # Import TCTW
  file <- files[files$name == paste0(area, "_Sentinel2_TCTW.tif"),]
  drive_download(
    as_id(file$id),
    path = paste0("./data/raw_data/covariates_raw/",area, "_Sentinel2_TCTW_raw.tif"), 
    overwrite = TRUE
  )
}

# 1.2 Compute each bbox ========================================================

bbox_list <- list()

for (area in areas) {
  file <- paste0(area,"_all")
  sites <- get(file)
  bbox <- st_bbox(sites)
  bbox_sf <- st_as_sfc(bbox)
  bbox_list[[area]] <- st_buffer(bbox_sf, 750)
}

# 2 Crop the covariates rivers and lakes with AOIs #############################

# 02.1 Resize and reproject all layers =========================================

area_layers <- list()
layers_names <- c("BA_bio03", "BA_bio12", "BA_bio15", "BA_bio16", "IA_bio03", "IA_bio12", "IA_bio15", 
                  "IA_bio16", "Aridity", "Aspect", "River_cost", "River_cost_permanent", "River_cost_seasonal",
                  "MSRM", "LSWI.1613", "LSWI.2202", "NDWI", "TCTW", "Slope", "TPI", "TWI")

type <- as.factor(c("cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", 
                    "disc", "cont", "cont", "cont", "cont", "cont", "cont", "cont", 
                    "cont", "cont", "disc", "cont"))


for (area in areas) {
  DEM <- rast(paste0("./data/raw_data/covariates_raw/", area,"_DEM_raw.tif"))
  DEM <- project(DEM, projections[[area]])
  mask <- vect(bbox_list[[area]])
  DEM <- crop(DEM, mask)
  
  list <- list.files("./data/raw_data/covariates_raw/", paste0(area, ".*\\.tif$"), recursive = TRUE, full.names = TRUE)
  
  if (area == "Shephelah") {
    list <- list[!grepl("DEM|rivers_seasonal", list)] # Remove the DEM
    layers_names_updated <- layers_names[!grepl("River_cost_seasonal|River_cost_permanent", layers_names, ignore.case = TRUE)]
    type_updated <- type[-c(12,13)]
  } else {
    list <- list[!grepl("DEM", list)] # Remove the DEM
    layers_names_updated <- layers_names
    type_updated <- type
  }
  
  provisional_layers <- list()
   for (i in 1:length(type_updated)) {
     r <- rast(list[[i]])
     layer_type <- type_updated[i]
     
     if (layer_type == "cont") {
       r <- project(r, projections[[area]], method = "bilinear")
       r <- crop(r, DEM)
       r <- resample(r, DEM)
     } else if (layer_type == "disc") {
       r <- project(r, projections[[area]], method = "mode")
       r <- crop(r, DEM)
       r <- resample(r, DEM, method = "mode")
     }
     r <- mask(r, DEM)
     names(r) <- layers_names_updated[[i]]
     provisional_layers[[layers_names_updated[[i]]]] <- r
   }
  provisional_layers[["DEM"]] <- DEM
  r <- rast(provisional_layers)
  r$TWI[r$TWI > 10, ] <- NA     # Remove water area
  area_layers[[area]] <- r
}


for (area in areas) {
writeRaster(area_layers[[area]], paste0("./data/raw_data/covariates/",area, "_covariates.tif"), overwrite=TRUE)
}

