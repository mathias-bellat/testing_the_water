##########################################################################
# This script is for downloading the covariates from present and         #
# past climates accessible via an API and computing DTM derived          #         
# covariates                                                             # 
#                                                                        #   
# Author: Mathias Bellat & Biel Soriano Elias                            #
# Affiliation : Tübingen University & Autonomous University of Barcelona #
# Creation date : 1/06/2026                                             #
# E-mail: mathias.archaeology@gmail.com & bielsoel28@gmail.com           #
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
pacman::p_load(pastclim, terra, sf, httr, googledrive, geodata, cli, leastcostpath, dplyr) # Specify required packages and download it if needed

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Import data sets ##########################################################

# 01.1 Prepare all data format  ================================================

eras <- c("BA", "IA")
areas <- c("Rhine", "Lazio", "Shephelah", "Kurdistan")
dir.create("./analysis/data/raw_data/covariates_raw/CHELSA_PAST")

set_data_path(path_to_nc = "./analysis/data/raw_data/covariates_raw/CHELSA_PAST")
download_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi")

for (period in eras) {
  dir.create(paste0("./analysis/data/raw_data/covariates_raw/CHELSA_PAST/",period))
}

time_period <- list(Rhine = list(BA = c(-4100, -2800), IA = c(-2200, -2000)),      # Bronze  (2200 - 800) and la Tène C and D (260 - 26)
                    Lazio = list(BA = c(-4200, -3000), IA = c(-2900, -2600)),     # Bronze  (2300 - 1000) and Iron (1000 - 600)
                    Shephelah = list(BA = c(-4900, -4500), IA = c(-2800, -2500)),    # Bronze  (3000 - 2500) and Iron (900 - 500)
                    Kurdistan = list(BA = c(-4900, -4600), IA = c(-2800, -2300)) # Bronze  (3000 - 2600) and Iron (900 - 300)
)

projections <- list(Rhine = "EPSG:32632", Lazio = "EPSG:32633", Shephelah = "EPSG:32636", Kurdistan = "EPSG:32638")

# 01.2 Download the past CHELSA  ===============================================

# Get vars
get_vars_for_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi")
bio_variables <- c("bio01", "bio03", "bio12", "bio15", "bio16")

# 02 Download and compute past Climate CHELSA ##################################
# Two are need as the reprojection of the CHELSA with a smaller extend exclude 
# a part of the layers

for (area in areas) {
  grid <- rast(paste0("./analysis/data/raw_data/grid/",area,"_AOI.tif"))
  crs(grid) <- "EPSG:4326"
  ext_large <- ext(
    xmin(grid) - 0.5,
    xmax(grid) + 0.5,
    ymin(grid) - 0.5,
    ymax(grid) + 0.5
  )
  
  grid_large <- extend(grid, ext_large)

  for (variable in bio_variables) {
    cli_text("Running {.val {variable}} : ({length(bio_variables)} layers)")
    cli_progress_bar(
      format = "Running CHELSA {.val {period}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
      total = 2, clear = FALSE)
    for (period in eras) {
      x <- region_series(
        time_bp = list(min = min(time_period[[area]][[period]]), max = max(time_period[[area]][[period]])),
        bio_variables = variable,
        dataset = "CHELSA_trace21k_1.0_0.5m_vsi"
      )
      x <- crop(x, grid_large)
      x <- rast(x)
      x <- mean(x, na.rm = TRUE)
      writeRaster(x, paste0("./nalaysis/data/raw_data/covariates_raw/CHELSA_PAST/",period, "/",area, "_", variable, ".tif" ), overwrite=T)
      cli_progress_update()
    }
    cli_progress_done()
    cat("\n")
  }
}

# 03 DTM derived covariates ####################################################

# 03.1 Export the DEM from Google Drive ========================================

drive_auth()
files <- drive_ls()
for (area in areas) {
  
  # Import DEM
  file <- files[files$name == paste0(area, "_DEM_raw.tif"),]
  drive_download(
    as_id(file$id),
    path = paste0("./analysis/data/raw_data/covariates_raw/",area, "_DEM_raw.tif"), 
    overwrite = TRUE
  )
}

# 03.2 Import list of DEMs =====================================================

# Define your folder path
folder_path <- "analysis/data/raw_data/covariates_raw"

# List all .tif files in that folder ending with "_DEM_raw.tif"
dtm_list <- list.files(path = folder_path, 
                       pattern = "_DEM_raw\\.tif$", 
                       full.names = TRUE)
#Prepare the names
base_name <- sub("_DEM_raw\\.tif$", "", basename(dtm_list))

# 03.3 Slope  ==================================================================

for (f in 1:length(dtm_list)) {
  
  # Load the raster
  r <- rast(dtm_list[[f]])
  
  # Compute slope in degrees
  slope <- terrain(r, v = "slope", unit = "degrees")
  
  # Create output filename based on original
  out_name <- file.path(folder_path, paste0(base_name[f], "_slope_deg.tif"))
  
  # Save to file
  writeRaster(slope, out_name, overwrite = TRUE)
  
}

# 03.4 MSRM  ====================================================================

#Function to select planar crs for correct msrm computation
ensure_planar_crs <- function(rast_obj) {
  # 1️⃣ Check CRS
  if (is.na(crs(rast_obj))) {
    stop("Input raster has no CRS defined!")
  }
  
  # 2️⃣ Check if CRS is geographic (lon/lat)
  if (is.lonlat(rast_obj)) {
    # Get centroid of raster extent
    ext <- ext(rast_obj)
    centroid_x <- (ext[1] + ext[2]) / 2
    centroid_y <- (ext[3] + ext[4]) / 2
    
    # Determine UTM zone based on longitude
    utm_zone <- floor((centroid_x + 180) / 6) + 1
    
    # Use appropriate hemisphere EPSG
    if (centroid_y >= 0) {
      epsg <- 32600 + utm_zone  # Northern Hemisphere
    } else {
      epsg <- 32700 + utm_zone  # Southern Hemisphere
    }
    
    message("Geographic CRS detected. Reprojecting to planar CRS EPSG:", epsg)
    rast_obj <- project(rast_obj, paste0("EPSG:", epsg))
    
  } else {
    message("Planar CRS detected. No transformation needed.")
  }
}

#Function to msrm computation (Always forced to 1 fmin)
msrm <- function(dem, fmin = 5, fmax = 100, x = 1.6, outdir = tempdir()) {
  library(terra)
  
  # --- Compute resolution (m) ---
  rr <- mean(res(dem))
  if (fmin <= rr) {
    message("fmin smaller than pixel size; setting fmin = raster resolution")
    fmin <- rr
  }
  
  # --- Calculate i and n values ---
  i <- floor(((fmin - rr) / (2 * rr))^(1 / x))
  n <- ceiling(((fmax - rr) / (2 * rr))^(1 / x))
  
  # Ensure valid kernel indices
  if (is.nan(i) || i < 1) i <- 1
  if (is.nan(n) || n <= i) n <- i + 1
  
  message("Raster resolution: ", round(rr, 3), " m")
  message("i = ", i, ", n = ", n)
  
  reliefs <- c()
  prev <- NULL
  
  for (ndx in i:n) {
    rad <- round((ndx ^ x))
    if (rad < 1) rad <- 1
    
    if (rad > 50) {
      warning("Kernel radius > 50 pixels skipped to avoid memory issues (rad = ", rad, ").")
      next
    }
    
    k <- matrix(1, nrow = rad * 2 + 1, ncol = rad * 2 + 1)
    k <- k / sum(k)
    
    # Create unique temp file
    fpath <- tempfile(pattern = paste0("LP_", ndx, "_"), tmpdir = outdir, fileext = ".tif")
    f <- focal(dem, w = k, fun = sum, na.rm = TRUE, filename = fpath, overwrite = TRUE)
    
    if (!is.null(prev)) {
      rpath <- tempfile(pattern = paste0("RM_", ndx, "_"), tmpdir = outdir, fileext = ".tif")
      rj <- prev - f
      writeRaster(rj, rpath, overwrite = TRUE)
      reliefs <- c(reliefs, rpath)
    }
    prev <- f
  }
  
  # Combine and compute mean of relief models
  msrm_files <- rast(reliefs)
  msrm_raw <- mean(msrm_files)
  msrm <- round(msrm_raw * 1000) / 1000
  
  # Optional: cleanup intermediate files
  unlink(reliefs)
  
  return(msrm)
}

for (f in 1:length(dtm_list)) {
  
  # Load the raster
  r <- rast(dtm_list[[f]])
  orig_crs <- crs(r)
  r_planar <- ensure_planar_crs(r)
  
  #Compute MSRM
  msrm_result <- msrm(r_planar, fmin = 100, fmax = 1000, x = 1.5) #f in meters
  
  #Convert raster bakc to original crs
  msrm_result_back <- project(msrm_result, orig_crs)
  
  # Create output filename based on original
  out_name <- file.path(folder_path, paste0(base_name[f], "_MSRM.tif"))
  
  # Save to file
  writeRaster(msrm_result_back, out_name, overwrite = TRUE)
  
}


# 03.5 TPI  =====================================================================

## Load Slope rasters
# List all .tif files in that folder ending with "_AOI_cropped.tif"
slope_list <- list.files(path = folder_path, 
                         pattern = "_slope_deg\\.tif$", 
                         full.names = TRUE)

## Compute TPI
for (f in 1:length(dtm_list)) {
  
  # Load the rasters
  r <- rast(dtm_list[[f]]) #dtm
  
  SlopeDegrees <- rast(slope_list[[f]]) #slope
  
  # Compute TPI classified
  # Compute TPIBig and TPISmall
  TPIS <- focal(r, w=9, fun=\(x) x[41] - mean(x[-41], na.rm = TRUE)) #100 m
  
  TPIB <- focal(r, w=81, fun=\(x) x[3281] - mean(x[-3281], na.rm = TRUE)) #1000 m
  
  # Extract the values from the rasters as a numeric vectors
  slope <- values(SlopeDegrees)   
  tpi_values_B <- values(TPIB)
  tpi_values_S <- values(TPIS)
  
  # Calculate the mean and standard deviation
  mean_TPI_B <- mean(tpi_values_B)
  std_TPI_B  <- sd(tpi_values_B)
  
  mean_TPI_S <- mean(tpi_values_S)
  std_TPI_S  <- sd(tpi_values_S)
  
  # Normalize the TPI values
  TPI_B_normalized <- (tpi_values_B - mean_TPI_B) / std_TPI_B
  TPI_S_normalized <- (tpi_values_S - mean_TPI_S) / std_TPI_S
  
  TPI_B_normalized <- as.integer(round(TPI_B_normalized))
  TPI_S_normalized <- as.integer(round(TPI_S_normalized))
  
  # Initialize an empty vector for classification results
  classification <- vector("character", length = length(TPI_S_normalized))
  
  for (i in 1:length(TPI_S_normalized)) {
    TPI_SN <- TPI_S_normalized[i]   # Extract SN (TPIst) value
    TPI_LN <- TPI_B_normalized[i]   # Extract LN (TPIst) value
    slope_ind <- slope[i]               # Extract slope value
    
    # Classification logic based on TPI (SN and LN) and slope conditions
    if (TPI_SN <= -1) {
      if (TPI_LN <= -1) {
        classification[i] <- "Canyons, V-shaped valleys"
      } else if (-1 < TPI_LN && TPI_LN < 1) {
        classification[i] <- "Mid-slope drainage, shallow valley"
      } else if (TPI_LN >= 1) {
        classification[i] <- "Highland drainage, headwaters"
      }
    } else if (-1 < TPI_SN && TPI_SN < 1) {
      if (TPI_LN <= -1) {
        classification[i] <- "U-shaped valleys"
      } else if (-1 < TPI_LN && TPI_LN < 1) {
        # Check if slope is NA
        if (is.na(slope_ind)) {
          classification[i] <- "Undefined plain area"
        } else {
          if (slope_ind <= 5) {
            classification[i] <- "Plains"
          } else {
            classification[i] <- "Open slopes"
          }
        }
      } else if (TPI_LN >= 1) {
        classification[i] <- "Upper slopes"
      }
    } else if (TPI_SN >= 1) {
      if (TPI_LN <= -1) {
        classification[i] <- "Local ridges/Hills in valleys"
      } else if (-1 < TPI_LN && TPI_LN < 1) {
        classification[i] <- "Middle ridges, small hills in plains"
      } else if (TPI_LN >= 1) {
        classification[i] <- "Peaks, high ridges"
      }
    }
  }
  
  # Convert classification to a raster
  classification_raster <- TPIS
  values(classification_raster)[!is.na(values(classification_raster))] <- classification
  
  out_name <- file.path(folder_path, paste0(base_name[f], "_tpi_categ.tif"))
  
  # Save to file
  writeRaster(classification_raster, out_name, overwrite = TRUE)
  
}

# 03.6 Aspect  ==================================================================

#Load reclass matrix
rcl <- matrix(c(
  -1, 0,      1,   # Flat
  0,    22.5, 2,   # North
  22.5, 67.5, 3,   # NE
  67.5, 112.5, 4,  # E
  112.5, 157.5, 5, # SE
  157.5, 202.5, 6, # S
  202.5, 247.5, 7, # SW
  247.5, 292.5, 8, # W
  292.5, 337.5, 9, # NW
  337.5, 360, 2    # North again
), ncol=3, byrow=TRUE)

for (f in 1:length(dtm_list)) {
  
  # Load the raster
  r <- rast(dtm_list[[f]])
  
  r_asp <- terrain(r, v = "aspect", unit = "degrees")
  
  # Reclassify raster
  asp_class <- classify(r_asp, rcl, include.lowest=TRUE, right=FALSE)
  
  # Create output filename based on original
  out_name <- file.path(folder_path, paste0(base_name[f], "_aspect.tif"))
  
  # Save to file
  writeRaster(asp_class, out_name, overwrite = TRUE)
  
}

# 03.7 TWI  =====================================================================

for (f in 1:length(dtm_list)) {
  
  # Load the raster
  r <- rast(dtm_list[[f]])
  
  # Compute flow direction
  flow_dir <- terrain(r, v = "flowdir")
  
  # Compute flow accumulation (upslope area)
  flow_acc <- flowAccumulation(flow_dir)
  
  # Flow accumulation (number of cells)
  cell_area <- prod(res(r))
  flow_acc_area <- flow_acc * cell_area
  
  # Slope in radians
  slope_rad <- terrain(r, v = "slope", unit = "radians")
  
  # Compute TWI
  twi <- log((flow_acc_area + 1) / tan(slope_rad))
  
  twi[is.infinite(values(twi))] <- NA
  
  # Create output filename based on original
  out_name <- file.path(folder_path, paste0(base_name[f], "_TWI.tif"))
  
  # Save to file
  writeRaster(twi, out_name, overwrite = TRUE)
  
}

# 03.8 Cost distance from rivers ===============================================

# --- Function to ensure planar CRS function ---
ensure_planar_crs <- function(sf_obj) {
  current_crs <- st_crs(sf_obj)
  
  if (is.na(current_crs)) stop("Input object has no CRS defined!")
  
  if (current_crs$units_gdal == "degree") {
    # Pick an appropriate UTM zone based on the centroid of the data
    centroid <- st_centroid(st_union(sf_obj))
    lon <- st_coordinates(centroid)[1]
    utm_zone <- floor((lon + 180)/6) + 1
    epsg <- 32600 + utm_zone  # northern hemisphere
    message("Geographic CRS detected. Reprojecting to planar CRS EPSG:", epsg)
    sf_obj <- st_transform(sf_obj, epsg)
  } else {
    message("Planar CRS detected. No transformation needed.")
  }
  
  return(sf_obj)
}


# Crop the rivers

rivers_croped <- list()

for (area in areas) {
  
  # Ensure CRS match
  if (st_crs(Rivers) != st_crs(bbox_list[[area]])) {
    Rivers_transformed <- st_transform(Rivers, st_crs(bbox_list[[area]]))
  }
  
  # Spatial crop / intersection
  rivers_croped[[area]] <- st_crop(Rivers_transformed, bbox_list[[area]])
  
  # Rasterisation and distance cost
  rivers_rast <- rasterize(vect(rivers_croped[[area]]), area_layers[[area]]$DEM, field=1, background=NA)
  area_layers[[area]]$water.distance <- distance(rivers_rast)
  
  
  st_write(rivers_croped[[area]], "./analysis/data/raw_data/shapefiles/Rivers.gpkg", layer = area, overwrite = TRUE, append = FALSE)
  writeRaster(area_layers[[area]]$water.distance, paste0("./analysis/data/raw_data/covariates_raw/",area, "water_distance.tif"), overwrite=TRUE)
}


# Crop the lakes

rivers_croped <- list()

# Temporarily disable s2 if using projected CRS
sf_use_s2(FALSE)

# Repair invalid geometries
Lakes <- st_make_valid(Lakes)

# Optional: check which features were invalid
invalid_idx <- which(!st_is_valid(Lakes))
if(length(invalid_idx) > 0) {
  message("Fixed ", length(invalid_idx), " invalid lake geometries")
}

for (area in areas) {
  
  # Ensure CRS match
  if (st_crs(Lakes) != st_crs(bbox_list[[area]])) {
    lakes_transformed <- st_transform(Lakes, st_crs(bbox_list[[area]]))
  }
  
  # Spatial crop / intersection
  lakes_croped[[area]] <- st_crop(lakes_transformed, bbox_list[[area]])
  
  # Skip empty results
  if (nrow(lakes_croped[[area]]) == 0) {
    message("⚠️ No lakes inside bbox for ", region)
    next
  }
  st_write(lakes_croped[[area]], "./analysis/data/raw_data/shapefiles/Lakes.gpkg", layer = area, overwrite = TRUE, append = FALSE)
}

#Load the rivers
gpkg_path <- "analysis/data/raw_data/shapefiles/Rivers.gpkg"

#Load the lakes
gpkg_path_lakes <- "analysis/data/raw_data/shapefiles/Lakes.gpkg"

# Get layer names
layer_names <- st_layers(gpkg_path)$name
layer_names_lakes <- st_layers(gpkg_path_lakes)$name

# Read each layer into a list
gpkg_list <- lapply(layer_names, function(layer) {
  st_read(gpkg_path, layer = layer, quiet = TRUE)
})

gpkg_list_lakes <- lapply(layer_names_lakes, function(layer) {
  st_read(gpkg_path_lakes, layer = layer, quiet = TRUE)
})

# Optionally name the list elements
names(gpkg_list) <- layer_names
names(gpkg_list_lakes) <- layer_names_lakes

for (f in 1:length(dtm_list)) {
  
  if(f == 3) {
    # Load the raster
    r <- rast(dtm_list[[f]])
    r_lr <- aggregate(r, fact= 2, fun = mean)
    orig_crs <- st_crs(r) #Store original CRS
    
    #Load the rivers and change them to points
    rivers <- gpkg_list[[base_name[f]]]
    rivers_proj <- ensure_planar_crs(rivers)
    rivers_proj <- st_cast(rivers_proj, "LINESTRING") #Ensure is a linestring
    rivers_proj <-  rivers_proj[rivers_proj$ORD_FLOW %in% c(1,2,3,4,5,6,7,8), ]
    
    #Subset rivers from permanent and seasonal
    rivers_seasonal <- rivers_proj[rivers_proj$ORD_FLOW %in% c(5,6,7,8), ]
    
    #Sample lines
    samples <- st_line_sample(rivers_proj, density = 1 / 50, type = "regular")  #Samp le points every 50 m along each line
    samples_s <- st_line_sample(rivers_seasonal, density = 1 / 50, type = "regular")  #Samp le points every 50 m along each line
    
    #Convert all of the lines to point to calculate cost
    all_points <- st_sfc(crs = st_crs(rivers_proj))
    all_points_s <- st_sfc(crs = st_crs(rivers_proj))
    
    for (i in seq_along(samples)) {
      if (length(samples[[i]]) > 0) {
        # Wrap as sfc before casting — preserves CRS
        mp <- st_sfc(samples[[i]], crs = st_crs(rivers_proj))
        pts <- st_cast(mp, "POINT")
        all_points <- c(all_points, pts)
      }
    }
    
    for (i in seq_along(samples_s)) {
      if (length(samples_s[[i]]) > 0) {
        # Wrap as sfc before casting — preserves CRS
        mp_s <- st_sfc(samples_s[[i]], crs = st_crs(rivers_proj))
        pts_s <- st_cast(mp_s, "POINT")
        all_points_s <- c(all_points_s, pts_s)
      }
    }
    
    rivers_points <- st_sf(geometry = all_points)  # Convert to sf object
    rivers_points <- st_transform(rivers_points, orig_crs)
    
    rivers_points_s <- st_sf(geometry = all_points_s)  # Convert to sf object
    rivers_points_s <- st_transform(rivers_points_s, orig_crs)
    
    #Compute cost raster
    cost_raster <- create_slope_cs(r_lr, cost_function =  "tobler", neighbours = 16)
    
    # Initialize a raster to store the minimum accumulated cost values
    min_cc <- rasterise(cost_raster)
    min_cc_s <- rasterise(cost_raster)
    
    values(min_cc) <- Inf
    values(min_cc_s) <- Inf
    
    
    ## General rivers ############################################################
    
    for (i in 1:nrow(rivers_points)) {
      # Compute accumulated cost from each origin
      coords <- st_coordinates(rivers_points[i,])
      
      # Check if point is inside raster extent
      if (!all(coords[,1] >= xmin(min_cc) & coords[,1] <= xmax(min_cc) &
               coords[,2] >= ymin(min_cc) & coords[,2] <= ymax(min_cc))) {
        next  # skip this iteration
      }
      
      cc <- create_accum_cost(x = cost_raster, origins = rivers_points[i,], FUN = mean, rescale = FALSE)
      
      # Update min_cc with the minimum value between the existing and the new cc
      min_cc <- min(min_cc, cc, na.rm = TRUE)
      
      rm(cc)
      rm(coords)
      gc()
    }
    
    #Replace Inf with NA
    min_cc[values(min_cc) == Inf] <- NA
    
    #Resample to original resolution
    min_cc <- resample(min_cc, r, method = "bilinear")
    
    # Create output filename based on original
    out_name <- file.path(folder_path, paste0(base_name[f], "_cost_from_rivers.tif"))
    
    # Save to file
    writeRaster(min_cc, out_name, overwrite = TRUE)
    
    ## Seasonal rivers ###########################################################
    
    for (i in 1:nrow(rivers_points_s)) {
      # Compute accumulated cost from each origin
      coords <- st_coordinates(rivers_points_s[i,])
      
      # Check if point is inside raster extent
      if (!all(coords[,1] >= xmin(min_cc_s) & coords[,1] <= xmax(min_cc_s) &
               coords[,2] >= ymin(min_cc_s) & coords[,2] <= ymax(min_cc_s))) {
        next  # skip this iteration
      }
      
      cc <- create_accum_cost(x = cost_raster, origins = rivers_points[i,], FUN = mean, rescale = FALSE)
      
      # Update min_cc_s with the minimum value between the existing and the new cc
      min_cc_s <- min(min_cc_s, cc, na.rm = TRUE)
      
      rm(cc)
      rm(coords)
      gc()
    }
    
    #Replace Inf with NA
    min_cc_s[values(min_cc_s) == Inf] <- NA
    
    #Resample to original resolution
    min_cc_s <- resample(min_cc_s, r, method = "bilinear")
    
    # Create output filename based on original
    out_name <- file.path(folder_path, paste0(base_name[f], "_cost_from_rivers_seasonal.tif"))
    
    # Save to file
    writeRaster(min_cc_s, out_name, overwrite = TRUE)
    
    
  } else {
    
    # Load the raster
    r <- rast(dtm_list[[f]])
    r_lr <- aggregate(r, fact= 2, fun = mean)
    orig_crs <- st_crs(r) #Store original CRS
    
    #Load the rivers and change them to points
    rivers <- gpkg_list[[base_name[f]]]
    rivers_proj <- ensure_planar_crs(rivers)
    rivers_proj <- st_cast(rivers_proj, "LINESTRING") #Ensure is a linestring
    rivers_proj <-  rivers_proj[rivers_proj$ORD_FLOW %in% c(1,2,3,4,5,6,7,8), ]
    
    #Subset rivers from permanent and seasonal
    rivers_permanent <- rivers_proj[rivers_proj$ORD_FLOW %in% c(1,2,3,4), ]
    rivers_seasonal <- rivers_proj[rivers_proj$ORD_FLOW %in% c(5,6,7,8), ]
    
    #Sample lines
    samples <- st_line_sample(rivers_proj, density = 1 / 50, type = "regular")  #Samp le points every 50 m along each line
    samples_p <- st_line_sample(rivers_permanent, density = 1 / 50, type = "regular")  #Samp le points every 50 m along each line
    samples_s <- st_line_sample(rivers_seasonal, density = 1 / 50, type = "regular")  #Samp le points every 50 m along each line
    
    #Convert all of the lines to point to calculate cost
    all_points <- st_sfc(crs = st_crs(rivers_proj))
    all_points_p <- st_sfc(crs = st_crs(rivers_proj))
    all_points_s <- st_sfc(crs = st_crs(rivers_proj))
    
    for (i in seq_along(samples)) {
      if (length(samples[[i]]) > 0) {
        # Wrap as sfc before casting — preserves CRS
        mp <- st_sfc(samples[[i]], crs = st_crs(rivers_proj))
        pts <- st_cast(mp, "POINT")
        all_points <- c(all_points, pts)
      }
    }
    
    for (i in seq_along(samples_p)) {
      if (length(samples_p[[i]]) > 0) {
        # Wrap as sfc before casting — preserves CRS
        mp_p <- st_sfc(samples_p[[i]], crs = st_crs(rivers_proj))
        pts_p <- st_cast(mp, "POINT")
        all_points_p <- c(all_points_p, pts_p)
      }
    }
    for (i in seq_along(samples_s)) {
      if (length(samples_s[[i]]) > 0) {
        # Wrap as sfc before casting — preserves CRS
        mp_s <- st_sfc(samples_s[[i]], crs = st_crs(rivers_proj))
        pts_s <- st_cast(mp_s, "POINT")
        all_points_s <- c(all_points_s, pts_s)
      }
    }
    
    rivers_points <- st_sf(geometry = all_points)  # Convert to sf object
    rivers_points <- st_transform(rivers_points, orig_crs)
    
    rivers_points_p <- st_sf(geometry = all_points_p)  # Convert to sf object
    rivers_points_p <- st_transform(rivers_points_p, orig_crs)
    
    rivers_points_s <- st_sf(geometry = all_points_s)  # Convert to sf object
    rivers_points_s <- st_transform(rivers_points_s, orig_crs)
    
    #Sample Rhin lakes
    if (f == 4) {
      
      lakes <- gpkg_list_lakes[[1]]
      lakes_proj <- st_transform(lakes, st_crs(rivers_proj))
      lakes_proj <- st_cast(lakes_proj, "LINESTRING") #Ensure is a linestring
      sample_lakes <-  st_line_sample(lakes_proj, density = 1 / 50, type = "regular")
      all_points_lakes <- st_sfc(crs = st_crs(lakes_proj))
      
      for (i in seq_along(sample_lakes)) {
        if (length(sample_lakes[[i]]) > 0) {
          # Wrap as sfc before casting — preserves CRS
          mp_ll <- st_sfc(sample_lakes[[i]], crs = st_crs(rivers_proj))
          pts_ll <- st_cast(mp_ll, "POINT")
          all_points_lakes <- c(all_points_lakes, pts_ll)
        }
      }
      
      lakes_points <- st_sf(geometry = all_points_lakes)  # Convert to sf object
      lakes_points <- st_transform(lakes_points, orig_crs)
      
      #Unite lakes with permanent rivers and all rivers
      rivers_points <- st_union(rivers_points, lakes_points)
      rivers_points_p <- st_union(rivers_points_p, lakes_points)
      
    }
    
    #Compute cost raster
    cost_raster <- create_slope_cs(r_lr, cost_function =  "tobler", neighbours = 16)
    
    # Initialize a raster to store the minimum accumulated cost values
    min_cc <- rasterise(cost_raster)
    min_cc_p <- rasterise(cost_raster)
    min_cc_s <- rasterise(cost_raster)
    
    values(min_cc) <- Inf
    values(min_cc_p) <- Inf
    values(min_cc_s) <- Inf
    
    
    ## General rivers ############################################################
    
    for (i in 1:nrow(rivers_points)) {
      # Compute accumulated cost from each origin
      coords <- st_coordinates(rivers_points[i,])
      
      # Check if point is inside raster extent
      if (!all(coords[,1] >= xmin(min_cc) & coords[,1] <= xmax(min_cc) &
               coords[,2] >= ymin(min_cc) & coords[,2] <= ymax(min_cc))) {
        next  # skip this iteration
      }
      
      cc <- create_accum_cost(x = cost_raster, origins = rivers_points[i,], FUN = mean, rescale = FALSE)
      
      # Update min_cc with the minimum value between the existing and the new cc
      min_cc <- min(min_cc, cc, na.rm = TRUE)
      
      rm(cc)
      rm(coords)
      gc()
    }
    
    #Replace Inf with NA
    min_cc[values(min_cc) == Inf] <- NA
    
    #Resample to original resolution
    min_cc <- resample(min_cc, r, method = "bilinear")
    
    # Create output filename based on original
    out_name <- file.path(folder_path, paste0(base_name[f], "_cost_from_rivers.tif"))
    
    # Save to file
    writeRaster(min_cc, out_name, overwrite = TRUE)
    
    ## Permanent rivers ##########################################################
    
    for (i in 1:nrow(rivers_points_p)) {
      # Compute accumulated cost from each origin
      coords <- st_coordinates(rivers_points_p[i,])
      
      # Check if point is inside raster extent
      if (!all(coords[,1] >= xmin(min_cc_p) & coords[,1] <= xmax(min_cc_p) &
               coords[,2] >= ymin(min_cc_p) & coords[,2] <= ymax(min_cc_p))) {
        next  # skip this iteration
      }
      
      cc <- create_accum_cost(x = cost_raster, origins = rivers_points_p[i,], FUN = mean, rescale = FALSE)
      
      # Update min_cc_p with the minimum value between the existing and the new cc
      min_cc_p <- min(min_cc_p, cc, na.rm = TRUE)
      
      rm(cc)
      rm(coords)
      gc()
    }
    
    #Replace Inf with NA
    min_cc_p[values(min_cc_p) == Inf] <- NA
    
    #Resample to original resolution
    min_cc_p <- resample(min_cc_p, r, method = "bilinear")
    
    # Create output filename based on original
    out_name <- file.path(folder_path, paste0(base_name[f], "_cost_from_rivers_permanent.tif"))
    
    # Save to file
    writeRaster(min_cc_p, out_name, overwrite = TRUE)
    
    ## Seasonal rivers ###########################################################
    
    for (i in 1:nrow(rivers_points_s)) {
      # Compute accumulated cost from each origin
      coords <- st_coordinates(rivers_points_s[i,])
      
      # Check if point is inside raster extent
      if (!all(coords[,1] >= xmin(min_cc_s) & coords[,1] <= xmax(min_cc_s) &
               coords[,2] >= ymin(min_cc_s) & coords[,2] <= ymax(min_cc_s))) {
        next  # skip this iteration
      }
      
      cc <- create_accum_cost(x = cost_raster, origins = rivers_points[i,], FUN = mean, rescale = FALSE)
      
      # Update min_cc_s with the minimum value between the existing and the new cc
      min_cc_s <- min(min_cc_s, cc, na.rm = TRUE)
      
      rm(cc)
      rm(coords)
      gc()
    }
    
    #Replace Inf with NA
    min_cc_s[values(min_cc_s) == Inf] <- NA
    
    #Resample to original resolution
    min_cc_s <- resample(min_cc_s, r, method = "bilinear")
    
    # Create output filename based on original
    out_name <- file.path(folder_path, paste0(base_name[f], "_cost_from_rivers_seasonal.tif"))
    
    # Save to file
    writeRaster(min_cc_s, out_name, overwrite = TRUE)
    
    
  }
  
}


# 04 Raster stack of each AOI ##################################################

# Define output directory for stacked GeoTIFFs
output_dir <- "analysis/data/raw_data/covariates"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define your folder path where covaraites layers are 
folder_path <- "analysis/data/raw_data/covariates_raw"

raster_stacks <- list()

for (area in areas) {
  # List all raster files for this area
  rasters_for_area <- list.files(
    folder_path,
    pattern = paste0("^", area, ".*\\.tif$"),  # Match rasters starting with area name
    full.names = TRUE
  )
  
  if (length(rasters_for_area) == 0) {
    warning("No rasters found for area: ", area)
    next
  }
  
  # Read all rasters and stack them
  raster_stack <- rast(rasters_for_area)
  
  # Optionally name the layers after the files (without extensions)
  names(raster_stack) <- tools::file_path_sans_ext(basename(rasters_for_area))
  
  # Store in list
  raster_stacks[[area]] <- raster_stack
  
  # Define output file path
  out_path <- file.path(output_dir, paste0(area, "_stack.tif"))
  
  # Save stack as multi-band GeoTIFF
  writeRaster( raster_stacks[[area]], out_path, overwrite = TRUE)
}
