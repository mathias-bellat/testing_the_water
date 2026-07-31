##########################################################################
# This script is for archaeological models                               #
#                                                                        #                                                   
# Author: Mathias Bellat & Biel Soriano Elias                            #
# Affiliation : Tübingen University & Autonomous University of Barcelona #
# Creation date : 16/06/2026                                             #
# E-mail: mathias.archaeology@gmail.com & bielsoel28@gmail.com           #
##########################################################################

# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Clean up workspace
rm(list=ls())

# 0.2 Install packages =========================================================

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(sf, ggplot2, cli, mapview, readr, usdm, corrplot, doParallel, SDMtune, zeallot, dismo, rJava, patchwork, 
               terra, grid, gridExtra, blockCV, rasterVis, colorspace)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Import data sets ##########################################################
eras <- c("BA", "IA")
areas <- c("Rhine", "Lazio", "Shephelah", "Kurdistan")
projections <- list(Rhine = "EPSG:32632", Lazio = "EPSG:32633", Shephelah = "EPSG:32636", Kurdistan = "EPSG:32638")

cov_sites <- list()

for (area in areas) {
  for (era in eras) {
   cov <-  read.csv(paste0("./analysis/data/raw_data/inputs/", area,"_sites_cov_for_", era,".csv"))
   cov_sites[[area]][[era]] <- cov[,-1]
  }
}

sf <- list()

for (area in areas) {
  sf[[area]] <- st_read("./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area, "_covariates"))
}

# 01.1 Create the pseudo-absence points ========================================
PA_Maxent <- list()

for (area in areas) {
  for (era in eras) {
    sites <- st_as_sf(cov_sites[[area]][[era]], coords = c("X", "Y"), crs = projections[[area]])
    sites_buffer <- st_buffer(sites, dist = 500)
    buffer_union <- st_union(sites_buffer)
    mask <- st_difference(sf[[area]], buffer_union)
    
    set.seed(1070)
    pa.1 <- st_sample(mask, size = 10000, type = "hexagonal")
    pa.1_sf <- st_as_sf(pa.1)
    coords <- st_coordinates(pa.1_sf$x)
    PA.1 <- as.data.frame(coords)
    
    set.seed(2070)
    pa.2 <- st_sample(mask, size = 10000, type = "hexagonal")
    pa.2_sf <- st_as_sf(pa.2)
    coords <- st_coordinates(pa.2_sf$x)
    PA.2 <- as.data.frame(coords)
    
    PA_Maxent[[area]][[era]] <- list("PA.1" = PA.1,
                             "PA.2" = PA.2)
  }
}


# 01.2 Plot the absence and presences ==========================================

for (area in areas) {
  for (era in eras) {
    presence <- as.data.frame(cov_sites[[area]][[era]])
    colnames(presence) <- names(cov_sites[[area]][[era]])
    absence.1 <- as.data.frame(PA_Maxent[[area]][[era]]$PA.1)
    absence.2 <- as.data.frame(PA_Maxent[[area]][[era]]$PA.2)
    
    # Extract response variable and coordinates for training data
    
    presence_sf <- st_as_sf(presence, coords = c("X", "Y"), crs = projections[[area]]) 
    absence.1_sf <- st_as_sf(absence.1, coords = c("X", "Y"), crs = projections[[area]])
    absence.2_sf <- st_as_sf(absence.2, coords = c("X", "Y"), crs = projections[[area]]) 
    
    mapview(presence_sf, layer.name = paste0("Sites ", era), col.regions = "lightgreen", cex = 5) + 
    mapview(absence.1_sf, layer.name = paste0("Sites absence 1 ", era), col.regions = "red", cex = 2) +
    mapview(absence.2_sf, layer.name = paste0("Sites absence 2 ", era), col.regions = "orange", cex = 2)
  }
}


# 01.3 Extract pseudo-absence covariates =======================================

cov_pa <- list()

for (area in areas) {
  for (era in eras) {
    cov_name <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    PA <- list()
    
    for (i in 1:2) {
      pa_sf <- st_as_sf(PA_Maxent[[area]][[era]][[i]], coords = c("X", "Y"), crs = projections[[area]])
      pa_cov <- st_join(pa_sf, sf[[area]])
      df_cov <- st_drop_geometry(pa_cov)
      PA[[i]] <- df_cov[,c(1,10:length(df_cov))]
      cov_era <-  df_cov[cov_name]
      colnames(cov_era) <- c("Bio03", "Bio12", "Bio15", "Bio16")
      PA[[i]] <- cbind(PA[[i]], cov_era)
    }
    pa.1 <- cbind(PA_Maxent[[area]][[era]][[1]], PA[[1]])
    pa.2 <- cbind(PA_Maxent[[area]][[era]][[2]], PA[[2]])
    pa.1$pa <- 1
    pa.2$pa <- 2
    df_pa <- rbind(pa.1, pa.2)
    df_pa <- na.omit(df_pa)
    write_delim(df_pa, paste0("./analysis/data/raw_data/inputs/", area, "_PA_cov_",era,".csv"), delim = ";")
    
    cov_pa[[area]][[era]] <- list(
      "PA.1" = df_pa[df_pa$pa == 1, !(names(df_pa) %in% "pa")],
      "PA.2" = df_pa[df_pa$pa == 2, !(names(df_pa) %in% "pa")]
    )
  }
}


save(cov_sites, cov_pa, PA_Maxent, file = "./analysis/data/derived_data/save/PA_APM.RData")

rm(list=setdiff(ls(), c("eras", "areas", "projections")))

# 02 Feature selection #########################################################

# 02.1 Realise the VIF =========================================================

load("./analysis/data/derived_data/save/PA_APM.RData")

for (area in areas) {
  for (era in eras) {
    presence <- cov_sites[[area]][[era]]
    presence <- presence[,-c(1:3)]
    PA.1 <- cov_pa[[area]][[era]]$PA.1[,-c(1:3)]
    PA.2 <- cov_pa[[area]][[era]]$PA.2[,-c(1:3)]
    points <- rbind(presence, PA.1, PA.2)
    
    vif <- vifcor(points[,-c(2,12)], th=0.7)
    print(vif)
    cat(area, "excluded variables for", era,":", vif@excluded)
    vif_df <- as.data.frame(vif@results)
    vif_plot <- ggplot(vif_df, aes(x = reorder(Variables, VIF), y = VIF)) +
      geom_bar(stat = "identity", fill = "lightblue") +
      coord_flip() +
      theme_minimal() +
      labs(title = paste0(area," VIF values for ", era, " period"), x = "Variables", y = "VIF") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    vif_plot
    
    ggsave(paste0("./analysis/data/derived_data/pre_process/VIF/VIF_",area, "_", era,".png"), vif_plot, width = 12, height = 8)
    ggsave(paste0("./analysis/data/derived_data/pre_process/VIF/VIF_",area, "_", era,".pdf"), vif_plot, width = 12, height = 8)
    
    write.table(vif_df, paste0("./analysis/data/derived_data/pre_process/VIF/VIF_",area, "_", era,".txt"))
  }
}


# 02.2 Human validation =========================================================

# Leave this part empty for now

rm(list=setdiff(ls(), c("areas", "eras", "projections")))

# 03 Maxent modelling ##########################################################

# 03.1 Prepare the data ========================================================
load("./analysis/data/derived_data/save/PA_APM.RData")
seeds <- c(1070,2707)

predictors <- list()

for (area in areas) {
  predictors[[area]]  <- rast(paste0("./analysis/data/raw_data/covariates/",area,"_covariates.tif"))   
}

# 03.2 Format the data ========================================================

cli_progress_bar(
  format = "MaxEnt {.val {area}} for {.val {era}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = 8, 
  clear = FALSE
)

for (area in areas) {
  for (era in eras) {
  
  Maxent_models <- list("PA.1" = NA, 
                        "PA.2" = NA)

  #Extract predictors of area and select the period ones
  predictors_final <- predictors[[area]]
  
  if (era == "BA"){
    predictors_final <- predictors_final[[-c(5,6,7,8)]]
  } else {
    predictors_final <- predictors_final[[5:22]]
  }
  
  if (area == "Levant") {
    names(predictors_final) <- c(
      "Bio03", "Bio12", "Bio15", "Bio16", 
      "Aridity", "Aspect",
      "River_cost",
      "MSRM", "LSWI.1613", "LSWI.2202", "NDWI", "TCTW",
      "Slope", "TPI", "TWI", "DEM"
    )
    
  } else {names(predictors_final) <- c(
    "Bio03", "Bio12", "Bio15", "Bio16", 
    "Aridity", "Aspect",
    "River_cost", "River_cost_permanent", "River_cost_seasonal",
    "MSRM", "LSWI.1613", "LSWI.2202", "NDWI", "TCTW",
    "Slope", "TPI", "TWI", "DEM"
  )}
  
  # Erase predictors on base of VIF (this part has to be removed if NO VIF is performed)
  filter <- read.table(paste0("./analysis/data/derived_data/pre_process/VIF/VIF_",area,"_",era,".txt"))
  vars_to_keep <- unique(c(filter[,1], "TPI", "Aspect"))
  idx <- which(names(predictors_final) %in% vars_to_keep)
  predictors_final <- predictors_final[[idx]]
  
  tiff(paste0("./analysis/data/derived_data/pre_process/blocks/Blocks_",area, "_", era,".tiff"), width = 12*300, height = 8*300, res = 300) # Width and height in pixels
  sac <- cv_spatial_autocor(predictors_final, progress = FALSE) 
  dev.off()
  
  for (i in 1:2){
    nb <- seeds[i]
    ifelse(era == "IA", nb <- nb + 1000, NA)
    
    coords_all <- rbind(cov_sites[[area]][[era]][, c(1:2)], cov_pa[[area]][[era]][[i]][, c(1:2)])
    data_all   <- rbind(cov_sites[[area]][[era]][, -c(1:3)], cov_pa[[area]][[era]][[i]][, -c(1:3)])
    pa <- c(rep(1, nrow(cov_sites[[area]][[era]])), rep(0, nrow(cov_pa[[area]][[era]][[i]])))
    data_all <- data_all[vars_to_keep]
    
    #Create a swd object to train MaxEnt model
    data_swd <- new("SWD",
                    species = era,
                    coords  = coords_all,
                    data    = data_all,
                    pa      = as.numeric(pa))  
    
    #Create sf object to compute Blocks
    data_sf <- cbind(coords_all,data_all,pa)
    
    data_sf <- st_as_sf(
      data.frame(data_sf),
      coords = c("X", "Y"),  # replace with your actual coordinate column names
      crs = projections[[area]] # set correct CRS
    )
    
    # 03.3 Split the data into blocks ===========================================
  
    kfolds <- cv_spatial(
      x = data_sf,
      column = "pa",
      r = predictors_final,
      size = sac$range, 
      k = 4,
      selection = "random",
      iteration = 100,
      hexagon = FALSE,
      plot = FALSE,
      progress = FALSE
    )
    
    # Spatial blocks with sampling points
    gg2 <- cv_plot(kfolds, x = data_sf)
    ggsave(filename = paste0("./analysis/data/derived_data/pre_process/blocks/Folds_",area, "_", era,"_PA",i,".tiff"),
           plot = gg2, width = 8, height = 6, dpi = 300, units = "in")
    
    rm(data_sf)
    
    # 03.4 Train a first model model ===========================================
    
    c(train_swd,test_swd) %<-% trainValTest(data_swd, 
                                            test = 0.30, 
                                            only_presence = TRUE, 
                                            seed = nb)
    
    first_model <- train("Maxent", 
                         data = train_swd,
                         progress = FALSE)
    
    # 03.4.1 Tune with genetic algorithm =======================================
    
    h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"),
              reg = seq(0.2, 5, 0.2),
              iter = c(500,1000,2000))
    
    
    genetic <- optimizeModel(first_model, 
                             hypers = h, 
                             metric = "auc", 
                             test = test_swd,
                             pop = 20, 
                             gen = 5,
                             keep_best = 0.4,
                             keep_random = 0.2,
                             mutation_chance = 0.4,
                             interactive = FALSE,
                             progress = FALSE,
                             seed = nb)
    
    # 03.4.2 Select best parameters and train final model ============================
    
    final_model <- combineCV(genetic@models[[which.max(genetic@results$test_AUC)]])
    
    # Final model training
    final_model <- train("Maxent", 
                         data = data_swd, 
                         fc = final_model@model@fc , 
                         reg = final_model@model@reg,
                         folds = kfolds,
                         progress = FALSE)
    
    # Save the run
    Maxent_models[[i]] <- c("data"= data, 
                            "first_model" = first_model, 
                            "genetic" = genetic, 
                            "final_model" = final_model, 
                            "train" = train_swd, 
                            "test" = test_swd)

    

  }
  save(Maxent_models, file = paste0("./analysis/data/derived_data/save/",era,"_",area,"_MaxEnt.RData"))
  cli_progress_update()
  }
}
cli_progress_done()

# 03.5 Create a merge report ===================================================

# Load the functions
source("./analysis/code/R/report_figures/plotROC_kfold_script.R")
source("./analysis/code/R/report_figures/modelReportCV_script.R")
cl <- makeCluster(4)
registerDoParallel(cl)

for (area in areas) {
  for (era in eras) {
    load(paste0("./analysis/data/derived_data/save/",era,"_",area,"_MaxEnt.RData"))
  
    if (!dir.exists(paste0("./analysis/data/derived_data/model/",area))) {
      dir.create(paste0("./analysis/data/derived_data/model/",area))}
  
    for (i in 1:2) {
      modelReportCV(model = Maxent_models[[i]]$final_model,
                  folder = paste0("./analysis/data/derived_data/model/",area,"/",era,"_PA.",i),
                  test = Maxent_models[[i]]$test_swd,
                  type = "cloglog",
                  response_curves = TRUE,
                  only_presence = FALSE,
                  jk = FALSE,
                  clamp = TRUE,
                  permut = 4,
                  verbose = TRUE)
      }
  }
} 

stopCluster(cl)  

rm(list=setdiff(ls(), c("eras", "areas", "projections")))

# 03.6 Check both models =======================================================
source("./analysis/code/R/report_figures/plotROC_kfold_script.R")

cli_progress_bar(
  format = "MaxEnt figures {.val {area}} {.val {era}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = length(eras), 
  clear = FALSE
)

Figures <- list()

cl <- makeCluster(4)
registerDoParallel(cl)

for (area in areas) {
  for (era in eras) {
    load(paste0("./analysis/data/derived_data/save/",era,"_",area,"_MaxEnt.RData"))
    
    cat("Training AUC 1: ",area, era, SDMtune::auc(Maxent_models[[1]]$final_model), "\n")
    cat("Testing AUC 1: ",area, era, SDMtune::auc(Maxent_models[[1]]$final_model, test = Maxent_models[[1]]$test_swd), "\n")
    plot(Maxent_models[[1]]$genetic, title = paste0("My experiment PA1 for ",era), interactive = TRUE)
    
    cat("Training AUC 2: ",area, era ,SDMtune::auc(Maxent_models[[2]]$final_model), "\n")
    cat("Testing AUC 2: ",area, era, SDMtune::auc(Maxent_models[[2]]$final_model, test = Maxent_models[[2]]$test_swd), "\n")
    plot(Maxent_models[[2]]$genetic, title = paste0("My experiment PA2",era), interactive = TRUE)
    
    # Select your model
    ifelse(SDMtune::auc(Maxent_models[[1]]$final_model) > SDMtune::auc(Maxent_models[[2]]$final_model), x <- 1, x <- 2)
    model <- Maxent_models[[x]]$final_model
    test <- Maxent_models[[x]]$test
    
    # 03.7 Explore the results ===================================================
    # Change variables names for plots
    
    replace_names <- function(original_names) {
      
      replacements <- list(
        list(pattern = "River_cost_permanent", 
             replacement = "River distance \ncost permanent"),
        list(pattern = "River_cost_seasonal", 
             replacement = "River distance \ncost seasonal"),
        list(pattern = "River_cost", 
             replacement = "River distance \ncost"),
        list(pattern = "LSWI.1613", 
             replacement = "LSWI 1613"),
        list(pattern = "LSWI.2202", 
             replacement = "LSWI 2202")
      )
      
      new_names <- original_names
      
      for (rule in replacements) {
        matches <- grep(rule$pattern, new_names, fixed = TRUE)
        
        if (length(matches) > 0) {
          new_names[matches] <- rule$replacement
        }
      }
      return(new_names)
    }
  
    replacers <- replace_names(colnames(model@data@data))
    
    
    # Jacknife with all variables
    jK <- doJk(model,
               metric = "auc",
               with_only = TRUE,
               test = test)
    
    jKtrain <- plotJk(jK, type = c("train"))
    jKtrain$data$Variable <- as.factor(replacers)
    
    jKtest <- plotJk(jK, type = c("test"))
    jKtest$data$Variable <- as.factor(replacers)
    
    combined_plot <- jKtrain + jKtest +
      plot_layout(guides = "collect")
    
    combined_plot <- combined_plot +
      plot_annotation(title = paste0(area, " Jackknife Plot - AUC for ",era, " (Train vs Test)"),
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    print(combined_plot)
    
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/Jacknife",area, "_", era,".png"), combined_plot, width = 12, height = 5)
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/Jacknife",area, "_", era,".pdf"), combined_plot, width = 12, height = 5)
    
    # AUC for export
    ROC <- plotROC_kfold(model, test = test)
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/AUC_",area, "_", era,".png"), ROC , width = 6, height = 6)
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/AUC_",area, "_", era,".pdf"), ROC , width = 6, height = 6)
    
    # Var importance
    var <- varImp(model, permut = 10)
    
    for (i in seq_along(replacers)) {
      var[,1] <- gsub(colnames(model@data@data)[i], replacers[i], var[,1])
    }
    
    write.table(var, paste0("./analysis/data/derived_data/model/",area,"/Variables_importance_",area, "_", era,".txt"))
    
    var_plot <- plotVarImp(var, color = "#abd9e9")
    var_plot
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/Variables_importance_",area, "_", era,".png"), var_plot , width = 6, height = 6)
    ggsave(paste0("./analysis/data/derived_data/model/",area,"/Variables_importance_",area, "_", era,".pdf"), var_plot , width = 6, height = 6)
    
    Figures[[area]][[era]] <- list(jk = jK,
                                  roc = ROC,
                                  var_imp = var) 

  
  }
  cli_progress_update()
}
cli_progress_done()
stopCluster(cl)  

PA_selection <- data.frame(BA = c(1, 2, 2, 1), IA = c(2, 2, 1, 1), row.names = c("Rhine", "Lazio", "Shephelah", "Kurdistan"))
write.table(PA_selection,"./analysis/data/derived_data/model/PA_selection.txt")

save(Figures, file = "./analysis/data/derived_data/save/Figures_APM.RData")
rm(list=setdiff(ls(), c("areas", "eras", "projections", "PA_selection")))

# 04 Prediction model ##########################################################
# 04.1 Prepare the data ========================================================

pred <- list()
sf <- list()
for (area in areas) {
  sf[[area]] <- st_read("./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area, "_covariates"))
  df <- st_drop_geometry(sf[[area]])
  row.names(df) <- df[,1]
  pred[[area]] <- df[,-1]
}


cli_progress_bar(
  format = "Predictions {.val {area}} {.val {era}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = length(eras), 
  clear = FALSE
)

cl <- makeCluster(4)
registerDoParallel(cl)

# 04.2 Run the predictions =====================================================

ths_PA <- data.frame(NULL)

for (area in areas) {
  for (era in eras) {
    
    cov_name <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    cov_pred <- pred[[area]][,c(9:length(pred[[area]]))]
    cov_era <-  pred[[area]][cov_name]
    colnames(cov_era) <- c("Bio03", "Bio12", "Bio15", "Bio16")
    cov_pred <- cbind(cov_pred, cov_era)
    cov_pred <- na.omit(cov_pred)
    
    load(paste0("./analysis/data/derived_data/save/",era,"_",area,"_MaxEnt.RData"))
    nb <- PA_selection[area,era]
    model <- Maxent_models[[nb]]$final_model 
    cov_pred <- cov_pred[names(model@data@data)]
    
    maps <- predict(model,
                    data = cov_pred,
                    fun = c("mean", "sd", "min", "max", "median"),
                    type = "cloglog")
    
    maps_df <- as.data.frame(maps)

    ths <- list()
    for (i in 1:4) {
      ths[[i]] <- SDMtune::thresholds(model@models[[i]], 
                                      type = "cloglog")
    }
    
    id_col <- ths[[1]][, 1, drop = FALSE]
    
    numeric_matrices <- lapply(ths, function(df) {
      as.matrix(df[, -1]) |> apply(2, as.numeric)
    })
    
    numeric_mean <- Reduce("+", numeric_matrices) / length(numeric_matrices)
    ths_mean <- data.frame(id_col, numeric_mean, check.names = FALSE)
    ths_PA[area,era] <- ths_mean[3, 2]
    maps_df$ths <- maps_df$mean
    maps_df$ths <- ifelse(maps_df$ths> ths_PA[1,era], 1, 0)
    
    maps_df$id <- row.names(cov_pred)
    pred_era <- merge(sf[[area]], maps_df, by.x = 0, by.y ="id", all = TRUE)
    pred_era <- na.omit(pred_era)
    pred_era <- pred_era[,(ncol(pred_era)-6):ncol(pred_era)]
    
    st_write(pred_era, paste0("./analysis/data/derived_data/maps/",area,"_predictions_maps.gpkg"), layer = paste0("Prediction_",era),  append = FALSE)
  }
  cli_progress_update()
}
cli_progress_done()

stopCluster(cl)

write.table(ths_PA,  "./analysis/data/derived_data/model/Thresholds.txt")

rm(list=setdiff(ls(), c("areas", "eras", "projections", "PA_selection", "ths_PA")))
# 04.3 Create raster maps ======================================================


for (area in areas) {
  for (era in eras) {
    raster <- rast(paste0("./analysis/data/raw_data/covariates/",area,"_covariates.tif"))
    dem <- raster$DEM
    
    sf <- st_read(paste0("./analysis/data/derived_data/maps/",area,"_predictions_maps.gpkg"), layer = paste0("Prediction_",era))
    raster_mean  <- rasterize(vect(sf), dem, field = "mean", background=NA)
    raster_sd  <- rasterize(vect(sf), dem, field = "sd", background=NA)
    raster_binary  <- rasterize(vect(sf), dem, field = "ths", background=NA)
    
    # Plot mean map
    plotPred(raster_mean,
             lt = paste0("Site probality mean\nfor ", era, " period"),
             colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
    
    # Plot SD map
    plotPred(raster_sd,
             lt = paste0("Site probality sd\nfor ", era, " period"),
             colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
    
    plotPred(raster_binary,
             lt = paste0("Site binary values\nfor ", era, " period"),
             colorramp = c("#ffffbf", "darkgreen"))
    
    writeRaster(raster_mean, paste0("./analysis/data/derived_data/maps/",area,"_prediction_maps_",era,".tif"), overwrite=TRUE)
    writeRaster(raster_sd, paste0("./analysis/data/derived_data/maps/",area,"_prediction_maps_sd_",era,".tif"), overwrite=TRUE)
    writeRaster(raster_binary, paste0("./analysis/data/derived_data/maps/",area,"_prediction_maps_binary_",era,".tif"), overwrite=TRUE)
  }  
}
rm(list=setdiff(ls(), c("areas", "eras", "projections", "PA_selection", "ths_PA")))

# 04.4 Check model statistics ==================================================

# Function to evaluate metrics

eval_metrics <- function(values, map, model, test = NULL) {
  conf <- table(factor(values[["pred"]], levels = c(0, 1)), factor(values[["obs"]], levels = c(0, 1)))
  
  TP <- conf["1", "1"]
  TN <- conf["0", "0"]
  FP <- conf["1", "0"]
  FN <- conf["0", "1"]
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  accuracy    <- (TP + TN) / sum(conf)
  
  expected_acc <- (((TP + FP) * (TP + FN) + (FN + TN) * (FP + TN)) / sum(conf)^2)
  kappa <- (accuracy - expected_acc) / (1 - expected_acc)
  
  tss <- SDMtune::tss(model, test =  test)
  auc <- SDMtune::auc(model, test = test)
  
  
  # Kvame gain as: gain = 1 - (pm / ps)
  pm <- as.numeric(global(map, fun = function(x) sum(x == 1, na.rm = TRUE) / sum(!is.na(x))))
  ps <- sum(values$obs == 1 & values$pred == 1)/sum(values$obs == 1)
  
  KG = 1 - (pm/ps)
  
  results <- data.frame(
    Metric      = c("Accuracy", "Kappa", "Sensitivity", "Specificity", "TSS", "AUC", "KG"),
    Value       = c(accuracy, kappa, sensitivity, specificity, tss, auc, KG)
  )
  if (!is.null(test)){
    test <- "test"
  }
  
  if (is.null(test)){
    test <- "train"
  }
  write.table(results,  paste0("./analysis/data/derived_data/model/",area,"/",area,"_metrics_",test, "_" ,era,".txt"))
  
  return(results)
}

metrics <- list()

cli_progress_bar(
  format = "Metrics {.val {area}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = length(eras), 
  clear = FALSE
)

cl <- makeCluster(4)
registerDoParallel(cl)

for (area in areas) {
  for (era in eras) {
    load(paste0("./analysis/data/derived_data/save/",era,"_",area,"_MaxEnt.RData"))
    sf <- st_read(paste0("./analysis/data/derived_data/maps/",area,"_predictions_maps.gpkg"), layer = paste0("Prediction_",era))
    
    nb <- PA_selection[area,era]
    model <- Maxent_models[[nb]]$final_model
    test <- Maxent_models[[nb]]$test
    
    # For training
    train_df <- cbind(model@data@coords, model@data@pa)
    train_df$id <- as.numeric(row.names(train_df))
    
    train_sf <- st_as_sf(train_df, coords = c("X", "Y"), crs = projections[[area]]) 
    train_sf <- st_intersection(train_sf, sf["ths"])
    train_sf <- train_sf[,-2]
    train_df <- st_drop_geometry(train_sf)
    colnames(train_df) <- c("obs", "pred")
    
    # For test
    test_df <- cbind(test@coords, test@pa)
    test_df$id <- as.numeric(row.names(test_df))
    
    test_sf <- st_as_sf(test_df, coords = c("X", "Y"), crs = projections[[area]]) 
    test_sf <- st_intersection(test_sf, sf["ths"])
    test_sf <- test_sf[,-2]
    test_df <- st_drop_geometry(test_sf)
    colnames(test_df) <- c("obs", "pred")
    
    map <- rast(paste0("./analysis/data/derived_data/maps/",area,"_prediction_maps_binary_",era,".tif"))
    
    train_metrics <- eval_metrics(train_df, map, model)
    train_metrics
    metrics[[area]][[era]]["Train"] <- train_metrics
    
    test_metrics <- eval_metrics(test_df, map, model, test)
    test_metrics
    metrics[[area]][[era]]["Test"] <- test_metrics
    
    
  }
  cli_progress_update()
}
cli_progress_done()
stopCluster(cl)

save(metrics, ths_PA, file = ("./analysis/data/derived_data/save/Models_results.RData"))
rm(list=setdiff(ls(), c("areas", "eras", "projections", "PA_selection", "ths_PA")))

# 05 Final figures for publication #############################################

# 05.1 Import the layers for prediction ========================================

APM_list <- NULL

for (area in areas) {
  for (era in eras) {
    file <- list.files("./analysis/data/derived_data/maps/", 
                       pattern = paste0(area,"_prediction_maps_",era,".tif"), 
                       full.names = TRUE, 
                       ignore.case = TRUE)
    
    APM_list <- c(APM_list, file) 
  }
}

APM_maps <- lapply(APM_list, rast)


# 05.2 Run for the predictive values ===========================================

rasters <- list(Rhine = c(APM_maps[[1]],APM_maps[[2]]), Lazio = c(APM_maps[[3]],APM_maps[[4]]),
                Shephelah = c(APM_maps[[5]],APM_maps[[6]]), Kurdistan = c(APM_maps[[7]],APM_maps[[8]]))

# Map function
plots <- list()

for (area in areas) {
  
  for (i in 1:2) {
    
    era <- ifelse(i == 1, "BA", "IA")
    n <- ifelse(area == "Shephelah" | area == "Kurdistan", 0.15, 0.3)
    df_study <- as.data.frame(rasters[[area]][[i]], xy = TRUE)
    names(df_study)[3] <- "value"
    legend <- paste(area, era)
    
    x_range <- range(df_study$x, na.rm = TRUE)
    y_range <- range(df_study$y, na.rm = TRUE)
    
    plot <- ggplot() +
      geom_raster(data = df_study, aes(x = x, y = y, fill = value)) +
      scale_fill_viridis_c(name = NULL) + 
      coord_equal() +
      annotate("label",  
               x = x_range[1] + diff(x_range) * 0.05,  
               y = y_range[2] - diff(y_range) * 0.05,
               label = legend,
                 hjust = n, 
                 vjust = 0.5,
                 size = 1.5,
                 color = "black",
                 fontface = "bold",
                 fill = "grey80",        
                 alpha = 0.6,           
                 label.size = 0.1,       
                 label.padding = unit(0.3, "lines")) +  
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(0, 0, 0, 0)
      )
    plots[[area]][[era]] <- plot
  }
}

# Add the legend
p1_legend <- ggplot(df_study, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Probability\nof sites") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 4)
  )

g_full <- ggplotGrob(p1_legend)
legend_grob <- g_full$grobs[[which(sapply(g_full$grobs, function(x) x$name) == "guide-box")]]

# Create a grid 
for (area in areas) {
  grid <- arrangeGrob(plots[[area]][[1]], plots[[area]][[2]],
                      ncol = 2,
                      widths = unit(rep(1, 2), "null"),
                      padding = unit(0, "pt"))
  
  n <- ifelse(area == "Shephelah" | area == "Kurdistan", 400, 500)
  tiff(paste0("./analysis/data/derived_data/figures/",area,"_predictive_plot.tiff"), width = 1400, height = n, res = 300, bg = "transparent")
  #png(paste0("./analysis/data/derived_data/figures/",area,"_predictive_plot.png"), width = 1400, height = n, res = 300, bg = "transparent")
  
  grid.newpage()
  grid.draw(arrangeGrob(
    grid,
    legend_grob,
    ncol = 2,
    widths = unit.c(unit(10, "cm"), unit(1.75, "cm"))
  ))
  dev.off() 
}

# 05.3 Import the layers for binary ============================================

APM_list <- NULL

for (area in areas) {
  for (era in eras) {
    file <- list.files("./analysis/data/derived_data/maps/", 
                       pattern = paste0(area,"_prediction_maps_binary_",era,".tif"), 
                       full.names = TRUE, 
                       ignore.case = TRUE)
    
    APM_list <- c(APM_list, file) 
  }
}

APM_maps <- lapply(APM_list, rast)

# 05.4 Run for the binary values ===============================================

rasters <- list(Rhine = c(APM_maps[[1]],APM_maps[[2]]), Lazio = c(APM_maps[[3]],APM_maps[[4]]),
                Shephelah = c(APM_maps[[5]],APM_maps[[6]]), Kurdistan = c(APM_maps[[7]],APM_maps[[8]]))


# Map function
plots <- list()
plots_legend <- list()

# Map function
for (area in areas) {
  for (i in 1:2) {
    
    era <- ifelse(i == 1, "BA", "IA")
    n <- ifelse(area == "Shephelah" | area == "Kurdistan", 0.15, 0.3)
    df_study <- as.data.frame(rasters[[area]][[i]], xy = TRUE)
    names(df_study)[3] <- "value"
    legend <- paste(area, era)
    
    x_range <- range(df_study$x, na.rm = TRUE)
    y_range <- range(df_study$y, na.rm = TRUE)
    
    plot <- ggplot() +
      geom_raster(data = df_study, aes(x = x, y = y, fill = factor(value))) +
      scale_fill_manual(
        values = c("0" = "grey", "1" = "green4"),
        name = "Classe",
        labels = c("0" = "Absence", "1" = "Presence")
      ) +
      coord_equal() +
      annotate("label",  
               x = x_range[1] + diff(x_range) * 0.05,  
               y = y_range[2] - diff(y_range) * 0.05,
               label = legend,
               hjust = n, 
               vjust = 0.5,
               size = 1.5,
               color = "black",
               fontface = "bold",
               fill = "grey80",        
               alpha = 0.6,           
               label.size = 0.1,       
               label.padding = unit(0.3, "lines")) +  
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(0, 0, 0, 0))
    
    

    plots[[area]][[era]] <- plot
    
  }
  # Add the legend
  p1_legend <- ggplot(df_study, aes(x = x, y = y, fill = factor(value))) +
    geom_raster() +
    scale_fill_manual(
      values = c("0" = "grey", "1" = "green4"),
      name = "Site presence\nbased on threshold",
      labels = c("0" = "Absence", "1" = paste0("Presence \n threholds: \n",
                                               "EB = ", round(ths_PA[area,1], digit = 3), "\n", 
                                               "IA = ", round(ths_PA[area,2], digit = 3))))  +
    coord_equal() +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 4)
    )
  
  g_full <- ggplotGrob(p1_legend)
  legend_grob <- g_full$grobs[[which(sapply(g_full$grobs, function(x) x$name) == "guide-box")]]
  plots_legend[[area]] <- legend_grob
  
}


# Create a grid 
for (area in areas) {
  grid <- arrangeGrob(plots[[area]][[1]], plots[[area]][[2]],
    ncol = 2,
    widths = unit(rep(1, 2), "null"),
    padding = unit(0, "pt")
  )
  
  n <- ifelse(area == "Shephelah" | area == "Kurdistan", 400, 500)
  tiff(paste0("./analysis/data/derived_data/figures/",area,"_binary_plot.tiff"), width = 1400, height = n, res = 300, bg = "transparent")
  #png(paste0("./analysis/data/derived_data/figures/",area,"_binary_plot.png"), width = 1400, height = n, res = 300, bg = "transparent")
  
  grid.newpage()
  grid.draw(arrangeGrob(
    grid,
    plots_legend[[area]],
    ncol = 2,
    widths = unit.c(unit(10, "cm"), unit(1.75, "cm"))
  ))
  dev.off() 
}

# 06 de Martone index map ######################################################

predictors <- list()

for (area in areas) {
  predictors[[area]]  <- rast(paste0("./analysis/data/raw_data/covariates/",area,"_covariates.tif"))   
}

Aridity_index <- list()
levels <- c(15, 24, 30, 35, 40, 50, 60, 187)
labels <- c(
  "Not class \n(I < 15)",
  "Semi-arid \n(15 ≤ I ≤ 24)",
  "Moderately-arid \n(24 ≤ I ≤ 30)",
  "Slightly-arid \n(30 ≤ I ≤ 35)",
  "Moderately-humid \n(35 ≤ I ≤ 40)",
  "Humid \n(40 ≤ I ≤ 50)",
  "Very-humid \n(50 ≤ I ≤ 60)",
  "Excessively-humid \n(60 ≤ I ≤ 187)"
)

for (area in areas) {
  for (era in eras) {
    r <- rast(c(paste0("./analysis/data/raw_data/covariates_raw/CHELSA_PAST/",era,"/",area,"_bio01.tif"),
                paste0("./analysis/data/raw_data/covariates_raw/CHELSA_PAST/",era,"/",area,"_bio12.tif")))
    
    r <- project(r, projections[[area]], method = "bilinear")
    r <- crop(r, predictors[[area]]$DEM)
    Aridity_index[[area]][[era]] <- resample(r, predictors[[area]]$DEM)
    names(Aridity_index[[area]][[era]]) <- c("Ta", "P")
    Aridity_index[[area]][[era]]$IDM <- Aridity_index[[area]][[era]]$P / (Aridity_index[[area]][[era]]$Ta + 10)
    
    IM_values <- values(Aridity_index[[area]][[era]]$IDM)
    Aridity_index[[area]][[era]]$IMartone <- cut(IM_values, breaks = c(-Inf, levels), labels = labels, right = FALSE)
  }
}

df_list <- list()

for(region in names(Aridity_index)) {
  for(period in names(Aridity_index[[region]])) {
    r <- Aridity_index[[region]][[period]]$IMartone
    
    temp_df <- as.data.frame(r, xy = TRUE) %>%
      rename(IMartone = 3) %>%
      mutate(
        Region = region,
        Period = period
      )
    
    df_list[[paste(region, period, sep = "_")]] <- temp_df
  }
}

all_df <- bind_rows(df_list)

plot_list <- list()

all_classes <- sort(unique(all_df$IMartone))
n_classes <- length(all_classes)

color_palette <- brewer.pal(min(n_classes,9), "YlGnBu")
names(color_palette) <- all_classes

for(r in unique(all_df$Region)) {
  for(p in unique(all_df$Period)) {
    temp <- all_df %>% filter(Region == r, Period == p)
    
    p_plot <- ggplot(temp, aes(x = x, y = y, fill = IMartone)) +
      geom_raster() +
      scale_fill_manual(values = color_palette, drop = FALSE, na.value = "transparent") +
      coord_equal() +
      labs(title = paste(r, "-", p), fill = "Aridity") +
      theme_minimal()
    
    plot_list[[paste(r,p,sep="_")]] <- p_plot
  }
}

gg1 <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]]) 
gg2 <- (plot_list[[5]] | plot_list[[6]]) / (plot_list[[7]] | plot_list[[8]])     

ggsave(filename = "./analysis/data/derived_data/pre_process/Aridity_01.pdf", plot = gg1, width = 10, height = 10)
ggsave(filename = "./analysis/data/derived_data/pre_process/Aridity_02.pdf", plot = gg2, width = 10, height = 10)
