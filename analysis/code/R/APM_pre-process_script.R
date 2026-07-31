##########################################################################
# This script is for preparing the data before the APM process           # 
#                                                                        #   
# Author: Mathias Bellat & Biel Soriano Elias                            #
# Affiliation : Tübingen University & Autonomous University of Barcelona #
# Creation date : 15/06/2026                                             #
# E-mail: mathias.archaeology@gmail.com & bielsoel28@gmail.com           #
##########################################################################

# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Set folder direction

# Clean up workspace
rm(list=ls())

# 0.2 Install packages =========================================================

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(terra, cli, sf, doParallel, exactextractr, foreach, mapview, yardstick, dplyr, ggplot2, patchwork, 
               viridis, caret, MLmetrics, readr, corrplot, kableExtra)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Prepare the tiles and rasters #############################################

# 1.1 Import all predictors ====================================================

eras <- c("BA", "IA")
areas <- c("Rhine", "Lazio", "Shephelah", "Kurdistan")
projections <- list(Rhine = "EPSG:32632", Lazio = "EPSG:32633", Shephelah = "EPSG:32636", Kurdistan = "EPSG:32638")

predictors <- list()

for (area in areas) {
  predictors[[area]]  <- rast(paste0("./analysis/data/raw_data/covariates/",area,"_covariates.tif"))   
}

# 1.2 Convert into tiles =======================================================

hex_grid <- list()
for (area in areas) {
  hex_grid[[area]] <- st_make_grid(
    st_as_sf(as.polygons(ext(predictors[[area]]))),
    cellsize = 50,
    square = FALSE)
  
  
  hex_grid[[area]] <- st_sf(id = 1:length(hex_grid[[area]]), geometry = hex_grid[[area]], crs = projections[[area]])
  st_write(hex_grid[[area]], "./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area,"_grid_50"))
}

# Set the discrete and continuous raster
type <- as.factor(c("cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", 
                    "disc", "cont", "cont", "cont", "cont", "cont", "cont", "cont", "cont", 
                    "cont", "disc", "cont", "cont"))

cl <- makeCluster(4)
registerDoParallel(cl)

hex_sf <- list()

for (area in areas) {
  sf <- hex_grid[[area]]
  
  if (area == "Shephelah") {
    type_updated <- type[-c(12,13)]
  } else {
    type_updated <- type
  }
  
  cli_text("Running tessellation i = {.val {area}} ({.val {4}} layers)")
  cli_progress_bar(
    format = "Layer {.val {i}}: {.val {layer_name}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(type_updated), 
    clear = FALSE
  )
  
  r <- predictors[[area]]
  
  for(i in 1:length(type_updated)){
    layer_name <- names(r)[i]
    layer_type <- type_updated[i]
    if (layer_type == "cont") {
      sf[[layer_name]] <- exact_extract(r[[i]], sf, "mean")
    } else if (layer_type == "disc") {
      sf[[layer_name]] <- exact_extract(r[[i]], sf, "mode")
    }
    cli_progress_update()
  }
  hex_sf[[area]] <- sf
  cli_progress_done()
  cat("\n")
  
  st_write(sf, "./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area,"_covariates"), overwrite = TRUE, append = FALSE)
}
stopCluster(cl)


# 1.3 Centroid extract  ========================================================

cl <- makeCluster(4)
registerDoParallel(cl)

cent_sf <- list()

for (area in areas) {
  
  if (area == "Shephelah") {
    type_updated <- type[-c(12,13)]
  } else {
    type_updated <- type
  }
  
  cli_text("Running centroid i = {.val {area}} ({.val {4}} layers)")
  cli_progress_bar(
    format = "Centroid {.val {i}}: {.val {layer_name}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(type_updated), 
    clear = FALSE
  )
  
  sf <- st_centroid(hex_grid[[area]])
  sf_layer <- list()
  
  for(i in 1:length(type_updated)){
    layer_name <- names(predictors[[area]][[i]])
    layer_type <- type_updated[i]
    if (layer_type == "cont") {
      sf_layer[[layer_name]] <- extract(predictors[[area]][[i]], vect(sf), method = "bilinear")[,2]
    } else if (layer_type == "disc") {
      sf_layer[[layer_name]] <- extract(predictors[[area]][[i]], vect(sf), method = "near")[,2]
    }
    cli_progress_update()
  }
  cent_sf[[area]] <- sf_layer
  cli_progress_done()
  cat("\n")
  
}
stopCluster(cl)  


# 1.4 Validation of the transformation  ========================================


cl <- makeCluster(4)
registerDoParallel(cl)

fig_trans <- list()

for (area in areas) {
  hex_val <- hex_grid[[area]]
  metrics <- data.frame(Layer = "", Type = "", RMSE = 0, MAE = 0, RSQ = 0, F1 = 0, Accuracy = 0)
  summaries <- data.frame(Layer = "", Type = "", Min. = 0, `1st Qu.` = 0, Median = 0, Mean = 0, `3rd Qu.` = 0, Max. = 0)
  
  if (area == "Shephelah") {
    type_updated <- type[-c(12,13)]
  } else {
    type_updated <- type
  }
  
  cli_text("Running i = {.val {area}} ({.val {4}} layers)")
  cli_progress_bar(
    format = "Processing validation {.val {i}}: {.val {layer_name}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(type_updated), 
    clear = FALSE
  )
  
  for(i in 1:length(type_updated)){
    layer_name <- names(predictors[[area]][[i]])
    layer_type <- type_updated[i]
    if (layer_type == "cont") {
      comparison_table <- tibble(
        hex_mean = hex_sf[[area]][[layer_name]],
        cent_bilinear = cent_sf[[area]][[layer_name]]
      )
      comparison_table$id <- seq(1, nrow(comparison_table))
      comparison_table[comparison_table[,2] == "NaN", 2] <- NA
      comparison_data <- na.omit(comparison_table)
      values <- metric_set(rmse, mae, rsq)(comparison_data, hex_mean, cent_bilinear)
      metrics[i,] <- data.frame(Layer = layer_name, Type = layer_type, RMSE = round(values[[1,3]], digit = 3), MAE = round(values[[2,3]], digit = 3), 
                                RSQ = round(values[[3,3]], digit = 3), NA, NA)
      
      comparison_data$difference <- comparison_data$hex_mean - comparison_data$cent_bilinear
      comparison_data$abs_difference <- abs(comparison_data$difference)
      comparison_data$relative_diff_pct <- (comparison_data$difference / comparison_data$cent_bilinear) * 100
      
      set.seed(1070)
      plot_data <- comparison_data[sample(nrow(comparison_data), min(10000, nrow(comparison_data))),]
      
      gg1 <- ggplot(plot_data, aes(x = cent_bilinear, y = hex_mean)) +
        geom_point(alpha = 0.6) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_smooth(method = "lm", color = "blue") +
        labs(x = "Centroid extraction (bilinear)", 
             y = "Hexagon mean",
             title = paste0(area, " area correlation between raster and hexagon for ", layer_name)) +
        theme_classic()
      
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_correlation_for_",layer_name, ".png"), gg1, width = 6, height = 6)
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_correlation_for_",layer_name, ".pdf"), gg1, width = 6, height = 6)
      
      fig_trans[[area]][["hex-rast_cor"]][[layer_name]]  <- gg1
      
      r_df <- as.data.frame(predictors[[area]][[i]])
      r_df$Type <- "Raster"
      hex_df <- as.data.frame(comparison_data$hex_mean)
      hex_df$type <- "Hexagon"
      colnames(hex_df) <- c(layer_name, "Type")
      data_df <- rbind(r_df, hex_df)
      colnames(data_df) <- c("Value", "Type")
      
      gg2 <- ggplot(data_df, aes(x = Type, y = Value, fill = Type)) +
        geom_violin() +
        scale_fill_viridis(discrete = TRUE) +
        geom_boxplot(fill="white", width=.1, size=0.1) +
        labs(x = "Type", 
             y = ylab(paste0(layer_name)),
             title = paste0(area, " area correlation between \nraster and hexagon for ", layer_name)) +
        theme_minimal() +
        theme(legend.position="none")
      
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_distribution_for_",layer_name, ".png"), gg2, width = 5, height = 6)
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_distribution_for_",layer_name, ".pdf"), gg2, width = 5, height = 6)
      
      fig_trans[[area]][["hex-rast_dis"]][[layer_name]]  <- gg2
      
      x <- summary(comparison_data[[1]])
      n <- nrow(summaries)+1
      summaries[n,] <- data.frame(Layer = layer_name, Type = "Hexagon", Min. = round(x[[1]], digit = 3), `1st Qu.` = round(x[[2]], digit = 3), 
                                  Median = round(x[[3]], digit = 3), Mean = round(x[[4]], digit = 3), `3rd Qu.` = round(x[[5]], digit = 3), Max. = round(x[[6]], digit = 3))
      
      x <- summary(r_df[,1])
      n <- nrow(summaries)+1
      summaries[n,] <- data.frame(Layer = layer_name, Type = "Raster", Min. = round(x[[1]], digit = 3), `1st Qu.` = round(x[[2]], digit = 3), 
                                  Median = round(x[[3]], digit = 3), Mean = round(x[[4]], digit = 3), `3rd Qu.` = round(x[[5]], digit = 3), Max. = round(x[[6]], digit = 3))
      
      colnames(comparison_data) <- c("hex_mean", paste0(layer_name,"_centroid"), "id", paste0(layer_name,"_difference"), 
                                     paste0(layer_name,"_abs_diff"), paste0(layer_name,"_realtive_diff"))
      x <- merge(comparison_table[,-c(1:2)], comparison_data[,-1], by = "id", all = TRUE)
      hex_val <- cbind(hex_val, x[,c(2:5)]) 
      
    } else if (layer_type == "disc") {
      comparison_table <- tibble(
        hex_mean = hex_sf[[area]][[layer_name]],
        cent_bilinear = cent_sf[[area]][[layer_name]]
      ) 
      comparison_table$id <- seq(1, nrow(comparison_table))
      comparison_data <- as.data.frame(comparison_table)
      comparison_data[comparison_data[,2] == "NaN", 2] <- NA
      comparison_data <- na.omit(comparison_data)
      comparison_data[,1] <- as.factor(comparison_data[,1])
      comparison_data[,2] <- as.factor(comparison_data[,2])
      
      F1 <- F1_Score(comparison_data[,1], comparison_data[,2])
      accuracy <- sum(comparison_data[,1] == comparison_data[,2]) / nrow(comparison_data)
      metrics[i,] <- data.frame(Layer = layer_name, Type = layer_type, NA , NA, NA, F1 = round(F1, digit = 3), Accuracy = round(accuracy, digit = 3))
      
      cm <- conf_mat(comparison_data, cent_bilinear, hex_mean)
      gg1 <- autoplot(cm, type = "heatmap") +
        scale_fill_gradient(low = "white", high = "steelblue")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste0(area, " area confusion Matrix for ", layer_name),
             x = "Centroïd",
             y = "Hexagon mean",
             fill = "Count")
      
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_matrix_for_",layer_name, ".png"), gg1, width = 6, height = 6)
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_matrix_for_",layer_name, ".pdf"), gg1, width = 6, height = 6)
      
      r_df <- as.data.frame(predictors[[area]][[i]])
      r_df$Type <- "Raster"
      hex_df <- as.data.frame(comparison_data$hex_mean)
      hex_df$type <- "Hexagon"
      colnames(hex_df) <- c(layer_name, "Type")
      data_df <- rbind(r_df, hex_df)
      colnames(data_df) <- c("Value", "Type")
      
      levels <- levels(factor(data_df[,1]))
      viridis_colors <- viridis(length(levels), option = "D")
      
      gg2 <- ggplot(data_df, aes(x = Type, fill = Value)) +
        geom_bar(position = "fill", width = 0.7) +
        scale_fill_manual(values = setNames(viridis_colors, levels), drop = FALSE) +
        ggtitle(paste0(area, " area distribution for ", layer_name)) +
        xlab("Type") +
        ylab("Proportion") +
        theme_minimal() +
        theme(legend.position = "right") 
      
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_distribution_class_for_",layer_name, ".png"), gg2, width = 5, height = 7)
      ggsave(paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_raster_distribution_class_for_",layer_name, ".pdf"), gg2, width = 5, height = 7)
      
      fig_trans[[area]][["hex-rast_dis_class"]][[layer_name]] <- gg1
      
      colnames(comparison_data) <- c("hex_mean", paste0(layer_name,"_centroid"), "id")
      name <- colnames(comparison_data[2])
      x <- merge(comparison_table[,-c(1:2)], comparison_data[,-1], by = "id", all = TRUE)
      hex_val <- cbind(hex_val, setNames(data.frame(x[[name]]), name))
    }
    cli_progress_update()
  }
  
  
  summaries <- summaries[-1,]
  write.table(summaries, paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_summary_stats_cells.txt"))
  write.table(metrics, paste0("./analysis/data/derived_data/pre_process/vectorisation/",area,"_hexagon_transformation_results.txt"))
  st_write(hex_val, "./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area, "_validation"), overwrite = TRUE, append=FALSE)
  cli_progress_done()
  cat("\n") 
  
}

save(fig_trans, file = "./analysis/data/derived_data/save/Figure_hex_trans.RData")

rm(list=setdiff(ls(), c("areas", "eras", "projections", "type", "predictors")))

# 02 Compare the data ##########################################################

# 2.1 Import files =============================================================

cov <- list()
sf <- list()
fig_stats <- list()

for (area in areas) {
  sf[[area]] <- st_read("./analysis/data/raw_data/covariates/Hex_grid.gpkg", layer = paste0(area, "_covariates"))
  df <- st_drop_geometry(sf[[area]])
  df <- df[,-1]
  cov[[area]] <- df
}


# 2.2 Compute statistics =======================================================


summary_numeric_df <- function(data) {
  num_data <- data[, sapply(data, is.numeric), drop = FALSE]
  
  stats <- t(sapply(num_data, function(x) {
    x_no_na <- x[!is.na(x)]
    c(NAs = sum(is.na(x)),
      Mean = mean(x_no_na, na.rm = TRUE),
      SD = sd(x_no_na, na.rm = TRUE),
      Min = min(x_no_na, na.rm = TRUE),
      Q1 = quantile(x_no_na, 0.25, na.rm = TRUE),
      Median = median(x_no_na, na.rm = TRUE),
      Q3 = quantile(x_no_na, 0.75, na.rm = TRUE),
      Max = max(x_no_na, na.rm = TRUE)
    )
  }))
  
  df_out <- as.data.frame(stats)
  
  df_out[] <- lapply(df_out, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  df_out
}

df_summary <- list()

for (area in areas) {
  df_summary[[area]] <- summary_numeric_df(cov[[area]])
  write.table(df_summary[[area]], paste0("./analysis/data/derived_data/pre_process/statistics/zones/",area,"_area_stats.txt"))
}

# 2.3 Compute test and density plot =============================================

df_wilcox <- list()

# For covariates present in whole area
names <- colnames(cov[["Rhine"]])
names <- names[-c(12:13)]

for (variable in names) {
  
  area_a <- cov[["Rhine"]][[variable]]
  area_b <- cov[["Lazio"]][[variable]]
  area_c <- cov[["Shephelah"]][[variable]]
  area_d <- cov[["Kurdistan"]][[variable]]
  n_sample <- 20000 
  set.seed(2025)
  area_a <- sample(na.omit(area_a), n_sample)
  area_b <- sample(na.omit(area_b), n_sample)
  area_c <- sample(na.omit(area_c), n_sample)
  area_d <- sample(na.omit(area_d), n_sample)
  
  df <- data.frame(
    valeur = c(area_a, area_b, area_c, area_d),
    zone = factor(rep(c("Rhine", "Lazio", "Shephelah", "Kurdistan"), each = n_sample))
  )
  
  if (variable %in% c("TPI", "Aspect")) {
    
    gg1 <- ggplot(df, aes(x = valeur, fill = zone)) +
      geom_histogram(bins = 30, alpha = 0.8, color = "black") +
      facet_wrap(~zone, ncol = 1, scales = "free_y") +
      labs(
        title = paste0(variable, " distribution per area"),
        x = "Values",
        y = "Density"
      )
    theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
  } else {
    # Density plot for other variables
    gg1 <- ggplot(df, aes(x = valeur, fill = zone)) +
      geom_density(alpha = 0.3) +
      labs(
        title = paste0(variable, " distribution per area"),
        x = "Values",
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
      )
  }
  
  
  ggsave(paste0("./analysis/data/derived_data/pre_process/statistics/zones/",variable, "_distribution_density.png"), gg1, width = 8, height = 6)
  ggsave(paste0("./analysis/data/derived_data/pre_process/statistics/zones/",variable, "_distribution_density.pdf"), gg1, width = 8, height = 6)
  
  fig_stats[["zone"]][[variable]] <- gg1
  
  test_a_b <- wilcox.test(area_a, area_b)
  test_a_c <- wilcox.test(area_a, area_c)
  test_a_d <- wilcox.test(area_a, area_d)
  test_b_c <- wilcox.test(area_b, area_c)
  test_b_d <- wilcox.test(area_b, area_d)
  test_c_d <- wilcox.test(area_c, area_d)
    
  df_wilcox[[variable]] <- data.frame(
    Area1 = c("Rhine", "Rhine", "Rhine", "Lazio", "Lazio", "Shephelah"),
    Area2 = c("Lazio", "Shephelah", "Kurdistan", "Shephelah", "Kurdistan", "Kurdsitan"),
    variable = variable,
    p_value = c(test_a_b$p.value, test_a_c$p.value, test_a_d$p.value, test_b_c$p.value, test_b_d$p.value,
                test_c_d$p.value),
    result = c(ifelse(test_a_b$p.value < 0.05, "Significant", "Not significant"),
                ifelse(test_a_c$p.value < 0.05, "Significant", "Not significant"),
                ifelse(test_a_d$p.value < 0.05, "Significant", "Not significant"),
                ifelse( test_b_c$p.value < 0.05, "Significant", "Not significant"),
                ifelse(test_b_d$p.value < 0.05, "Significant", "Not significant"),
                ifelse(test_c_d$p.value < 0.05, "Significant", "Not significant"))
  )

}

# For river cost (absent from Shephelah area)
names <- c("River_cost_permanent","River_cost_seasonal")

for (variable in names) {
  area_a <- cov[["Rhine"]][[variable]]
  area_b <- cov[["Lazio"]][[variable]]
  area_c <- cov[["Kurdistan"]][[variable]]
  n_sample <- 20000 
  set.seed(2025)
  area_a <- sample(na.omit(area_a), n_sample)
  area_b <- sample(na.omit(area_b), n_sample)
  area_c <- sample(na.omit(area_c), n_sample)
  
  df <- data.frame(
    valeur = c(area_a, area_b, area_c),
    zone = factor(rep(c("Rhine", "Lazio", "Kurdistan"), each = n_sample))
  )
  
  gg1 <-  ggplot(df, aes(x = valeur, fill = zone)) +
    geom_density(alpha = 0.3) +
    labs(title = paste0(variable," distribution per area"), x = "Value", y = "Density") +
    theme_minimal()
  
  ggsave(paste0("./analysis/data/derived_data/pre_process/statistics/zones/",variable, "_distribution_density.png"), gg1, width = 8, height = 6)
  ggsave(paste0("./analysis/data/derived_data/pre_process/statistics/zones/",variable, "_distribution_density.pdf"), gg1, width = 8, height = 6)
  
  fig_stats[["zone"]][[variable]] <- gg1
  
  
  test_a_b <- wilcox.test(area_a, area_b)
  test_a_c <- wilcox.test(area_a, area_c)
  test_b_c <- wilcox.test(area_b, area_c)
  
  
  df_wilcox[[variable]] <- data.frame(
    Area1 = c("Rhine", "Rhine", "Lazio"),
    Area2 = c("Lazio", "Kurdistan", "Kurdistan"),
    variable = variable,
    p_value = c(test_a_b$p.value, test_a_c$p.value, test_b_c$p.value),
    result = c(ifelse(test_a_b$p.value < 0.05, "Significant", "Not significant"),
               ifelse(test_a_c$p.value < 0.05, "Significant", "Not significant"),
               ifelse( test_b_c$p.value < 0.05, "Significant", "Not significant"))
  ) 

}

test <- do.call(rbind, df_wilcox)
test$p_value  <- formatC(test$p_value, format = "f", digits = 3)

write.table(test, "./analysis/data/derived_data/pre_process/statistics/zones/Wilcox_test.txt")

test$result <- ifelse(test$result == "Significant",
                            cell_spec(test$result, "html", background = "lightgreen"),
                            as.character(test$result))

kable(test, "html", escape = FALSE) %>% 
  kable_styling(full_width = TRUE)  

# 03 Compute the graphical overview of the background ##########################

# 9 per page in pdf
for (area in areas) {
  
  # Set a reasonable PDF size (9x6 inches)
  pdf(file.path("data/derived_data/pre_process",
                paste0(area,"_graphical_overview.pdf")),
      width = 10, height = 10,  # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk")   # Color model )
  
  # Get number of predictor layers
  n_layers <- dim(predictors[[area]])[3]
  
  # Loop through layers 16 at a time
  for (start_idx in seq(1, n_layers, by = 9)) {
    
    # Define the plotting layout
    par(mfrow = c(3,3), mar = c(2, 2, 3, 1))
    
    # Determine which layers go on this page
    end_idx <- min(start_idx + 8, n_layers)
    
    for (i in start_idx:end_idx) {
      # Plot the layer
      plot(predictors[[area]][[i]],
           main = names(predictors[[area]])[i],
           axes = FALSE, box = FALSE)
    }
  }
  
  dev.off()
}

# All in one in png
for (area in areas) {
  
  # Get number of predictor layers
  n_layers <- dim(predictors[[area]])[3]
  
  # Compute layout size automatically (closest square grid)
  n_cols <- ceiling(sqrt(n_layers))
  n_rows <- ceiling(n_layers / n_cols)
  
  # Define output file path
  png(file.path("data/derived_data/pre_process",
                paste0(area,"_graphical_overview.png")),
      width = 300 * n_cols, height = 300 * n_rows, res = 100)
  
  # Layout and margins
  par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 1))
  
  # Loop through all layers and plot
  for (i in 1:n_layers) {
    plot(predictors[[area]][[i]],
         main = names(predictors[[area]])[i],
         axes = FALSE, box = FALSE)
  }
  
  dev.off()
}

# 04 Extract variables #########################################################

load("./analysis/data/derived_data/save/Sites_filter.RData")

sites <- list("Rhine" = Rhine_all[,c(1,4,22,23)], "Lazio" = Lazio_all[,c(1,2,12,13)], 
              "Shephelah" = Shephelah_all[,c(2,26:28)], "Kurdistan" = Kurdistan_all[,c(2,3,13,14)])

# 4.1 Extract variables for positive sites =====================================

sites_with_cov <- list()

# Remove NA and join site
for (area in areas) {
  sf[[area]] <- na.omit(sf[[area]])
  sites_with_cov[[area]] <- st_join(sites[[area]], sf[[area]])
  sites_with_cov[[area]]$era[sites_with_cov[[area]]$era == "Bronze"] <- "BA"
  sites_with_cov[[area]]$era[sites_with_cov[[area]]$era == "Iron"] <- "IA"
  print(sum(is.na(sites_with_cov[[area]]) == TRUE))
}

# Absence in the Lazio dataset is for the toponymes

cov_sites <- list()

for (area in areas) {
  
  for (era in eras) {
    cov_name <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    df_cov <- st_drop_geometry(sites_with_cov[[area]])
    coords <- st_coordinates(sites_with_cov[[area]])
    df_cov$X <- coords[, "X"]
    df_cov$Y <- coords[, "Y"]
    df_cov <- df_cov[df_cov$era == era,]
    cov_sites[[area]][[era]] <- df_cov[,c(2,15:length(df_cov)-2)]
    cov_era <-  df_cov[cov_name]
    colnames(cov_era) <- c("Bio03", "Bio12", "Bio15", "Bio16")
    cov_sites[[area]][[era]] <- cbind(df_cov[,c(length(df_cov)-1,length(df_cov))], cov_sites[[area]][[era]], cov_era)
    write.csv(cov_sites[[area]][[era]], paste0("./analysis/data/raw_data/inputs/", area, "_sites_cov_for_", era,".csv"))
  }
}

# 4.2 Compute basic statistics =================================================

summary_numeric_df <- function(data) {
  num_data <- data[, sapply(data, is.numeric), drop = FALSE]
  
  stats <- t(sapply(num_data, function(x) {
    x_no_na <- x[!is.na(x)]
    c(NAs = sum(is.na(x)),
      Mean = mean(x_no_na, na.rm = TRUE),
      SD = sd(x_no_na, na.rm = TRUE),
      Min = min(x_no_na, na.rm = TRUE),
      Q1 = quantile(x_no_na, 0.25, na.rm = TRUE),
      Median = median(x_no_na, na.rm = TRUE),
      Q3 = quantile(x_no_na, 0.75, na.rm = TRUE),
      Max = max(x_no_na, na.rm = TRUE)
    )
  }))
  
  df_out <- as.data.frame(stats)
  
  df_out[] <- lapply(df_out, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  df_out
}

df_summary <- list()

for (area in areas) {
  for (era in eras) {
    df_summary[[area]][[era]] <- summary_numeric_df(cov_sites[[area]][[era]][,-c(1:2)])
    write.table(df_summary[[area]][[era]], paste0("./analysis/data/derived_data/pre_process/statistics/sites/",area,"_sites_stats_for_",era,".txt"))
  }
  
}


# 4.3 Compute Wilcox test BA vs. IA ============================================

df_wilcox <- list()

for (area in areas) {
  site_BA <- cov_sites[[area]][["BA"]][-1]
  site_IA <- cov_sites[[area]][["IA"]][-1]
  
  names <- colnames(site_BA[,-c(1:2)])

  for (variable in names) {
    
    cov_BA <- as.numeric(site_BA[[variable]])
    cov_IA <- as.numeric(site_IA[[variable]])
    
    test <- wilcox.test(cov_BA, cov_IA)
    
    df_wilcox[[area]][[variable]] <- data.frame(
      Area = area,
      variable = variable,
      p_value = c(test$p.value),
      result = c(ifelse(test$p.value < 0.05, "Significant", "Not significant"))
    )
  }
}  


Rhine_test <- do.call(rbind, df_wilcox[["Rhine"]])
Lazio_test <- do.call(rbind, df_wilcox[["Lazio"]])
Shephelah_test <- do.call(rbind, df_wilcox[["Shephelah"]])
Kurdistan_test <- do.call(rbind, df_wilcox[["Kurdistan"]])

Final_test <- rbind(Rhine_test, Lazio_test, Shephelah_test, Kurdistan_test)
row.names(Final_test) <- c(1:nrow(Final_test))
Final_test$p_value  <- formatC(Final_test$p_value, format = "f", digits = 3)

write.table(Final_test, "./analysis/data/derived_data/pre_process/statistics/sites/Wilcox_test_sites_vs_sites.txt")


Final_test$result <- ifelse(Final_test$result == "Significant",
                    cell_spec(Final_test$result, "html", background = "lightgreen"),
                    as.character(Final_test$result))

kable(Final_test, "html", escape = FALSE) %>% 
  kable_styling(full_width = TRUE)


# 4.4 Compute Wilcox test sites vs. background =================================

df_wilcox <- list()

for (area in areas) {

  for (era in eras) {
    names <- colnames(cov_sites[[area]][[era]][,-c(1:3)])
    sites_era <- cov_sites[[area]][[era]][,-c(1:3)]
    
    cov_name <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    cov_era <- cov[[area]][cov_name]
    colnames(cov_era) <- c("Bio03", "Bio12", "Bio15", "Bio16")
    cov_era <- cbind(cov[[area]][,c(9:length(cov[[area]]))], cov_era)
    
    for (variable in names) {
      n_sample <- 20000 
      set.seed(2025)
      cov_background <- as.numeric(cov_era[[variable]])
      cov_background <- sample(na.omit(cov_background), n_sample)
      
      cov_sites_era <- as.numeric(sites_era[[variable]])
      
      test <- wilcox.test(cov_sites_era, cov_background)
      
      df_wilcox[[area]][[era]][[variable]] <- data.frame(
        Area = area,
        Era = era,
        variable = variable,
        p_value = c(test$p.value),
        result = c(ifelse(test$p.value < 0.05, "Significant", "Not significant"))
      )
    }
  }
  
}  


Rhine_BA_test <- do.call(rbind, df_wilcox[["Rhine"]][["BA"]])
Rhine_IA_test <- do.call(rbind, df_wilcox[["Rhine"]][["IA"]])
Lazio_BA_test <- do.call(rbind, df_wilcox[["Lazio"]][["BA"]])
Lazio_IA_test <- do.call(rbind, df_wilcox[["Lazio"]][["IA"]])
Shephelah_BA_test <- do.call(rbind, df_wilcox[["Shephelah"]][["BA"]])
Shephelah_IA_test <- do.call(rbind, df_wilcox[["Shephelah"]][["IA"]])
Kurdistan_BA_test <- do.call(rbind, df_wilcox[["Kurdistan"]][["BA"]])
Kurdistan_IA_test <- do.call(rbind, df_wilcox[["Kurdistan"]][["IA"]])

Final_test <- rbind(Rhine_BA_test, Lazio_BA_test, Shephelah_BA_test, Kurdistan_BA_test,
                    Rhine_IA_test, Lazio_IA_test, Shephelah_IA_test, Kurdistan_IA_test)
row.names(Final_test) <- c(1:nrow(Final_test))
Final_test$p_value  <- formatC(Final_test$p_value, format = "f", digits = 3)

write.table(Final_test, "./analysis/data/derived_data/pre_process/statistics/sites/Wilcox_test_sites_vs_background.txt")

Final_test$result <- ifelse(Final_test$result == "Significant",
                            cell_spec(Final_test$result, "html", background = "lightgreen"),
                            as.character(Final_test$result))

kable(Final_test, "html", escape = FALSE) %>% 
  kable_styling(full_width = TRUE)  

# 05 Visualisations ############################################################
# 5.1 Density plot for sites and background combined ===========================

#Extract names of the variables
names <- colnames(cov[["Rhine"]])
cov_sites <- list()

#Extract background and sites (BA and IA) values
for(area in areas){
  
  #Extract sites points
  sf[[area]] <- na.omit(sf[[area]])
  sites_with_cov[[area]] <- st_join(sites[[area]], sf[[area]])
  sites_with_cov[[area]]$era[sites_with_cov[[area]]$era == "Bronze"] <- "BA"
  sites_with_cov[[area]]$era[sites_with_cov[[area]]$era == "Iron"] <- "IA"
  print(sum(is.na(sites_with_cov[[area]]) == TRUE))
  
  #Extract all values of sites
  
  for (era in eras) {
    cov_name <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    df_cov <- st_drop_geometry(sites_with_cov[[area]])
    df_cov <- df_cov[df_cov$era == era,]
    cov_sites[[area]][[era]] <- df_cov[,c(2,13:length(df_cov))]
    cov_era <-  df_cov[cov_name]
    colnames(cov_era) <- c(paste0(era,"_bio03"), paste0(era,"_bio12"), paste0(era,"_bio15"), paste0(era,"_bio16"))
    cov_sites[[area]][[era]] <- cbind(cov_sites[[area]][[era]], cov_era)
  }
  
  # Divide values in BA and IA
  BA_vals <- cov_sites[[area]]$BA
  IA_vals <- cov_sites[[area]]$IA
  
  #Extract variables values for background and sites 
  for (variable in names) {
    
    if(area == "Shephelah" && (variable == "River_cost_permanent" || variable == "River_cost_seasonal")){
      next
    }
    
    background_variable <- cov[[area]][[variable]]
    n_sample <- 20000
    set.seed(2025)
    background_values <- sample(na.omit(background_variable), n_sample)
    
    if (variable %in% names(IA_vals) && variable %in% names(BA_vals)) {
      
      IA_var_val <- IA_vals[[variable]]
      BA_var_val <- BA_vals[[variable]]
      
    } else if (variable %in% names(IA_vals)) {
      
      IA_var_val <- IA_vals[[variable]]
      BA_var_val <- NULL
      
    } else if (variable %in% names(BA_vals)) {
      
      IA_var_val <- NULL
      BA_var_val <- BA_vals[[variable]]
      
    } 
    
    #Create df for plotting
    df <- data.frame(
      value = c(background_values, IA_var_val, BA_var_val),
      zone = factor(c(
        rep("Background", length(background_values)),
        rep("IA", length(IA_var_val)),
        rep("BA", length(BA_var_val))
      ), levels = c("Background", "IA", "BA"))
    )
    
    #Plot results
    if (variable %in% c("TPI", "Aspect")) {
      
      gg1 <- ggplot(df, aes(x = value, fill = zone)) +
        geom_histogram(bins = 30, alpha = 0.8, color = "black") +
        facet_wrap(~zone, ncol = 1, scales = "free_y") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none"
        )
      
    } else {
      # Density plot for other variables
      gg1 <- ggplot(df, aes(x = value, fill = zone)) +
        geom_density(alpha = 0.3) +
        labs(
          title = paste0(variable, " distribution in ", area),
          x = "Values",
          y = "Density"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()
        )
    }
    
    # save plot
    ggsave(filename = file.path("data/derived_data/pre_process/statistics/sites",
          paste0(area, "_", variable,"_bg_vs_sites.png")),
          plot = gg1, width = 8, height = 5, dpi = 300)
    
    ggsave(filename = file.path("data/derived_data/pre_process/statistics/sites",
                                 paste0(area, "_", variable,"_bg_vs_sites.pdf")),
            plot = gg1, width = 8, height = 5)
    
    fig_stats[[area]][[variable]] <- gg1
    
  }
}

save(fig_stats, cov_sites, file = "./analysis/data/derived_data/save/Figure_APM_pre-process_stats.RData")

# 5.2 Correlation for Maxent (positive only) ===================================

for (area in areas) {
 
  if (area == "Shephelah") {
    
    for (era in eras) {  
      pdf(paste0("./analysis/data/derived_data/pre_process/",area, "_corrplot_for_", era, ".pdf"),    # File name
          width = 10, height = 10,  # Width and height in inches
          bg = "white",          # Background color
          colormodel = "cmyk")   # Color model 
      
      # Correlation of the data (removed discrete variables)
      corrplot(cor(cov_sites[[area]][[era]][,-c(1,3,11)]),  method = "color", col = viridis(200), 
               type = "upper", 
               addCoef.col = "black", # Add coefficient of correlation
               tl.col = "black", tl.srt = 45, # Text label color and rotation
               number.cex = 0.7, # Size of the text labels
               cl.cex = 0.7) # Size of the color legend text
      
      dev.off()
    }
  } else {
    
    for (era in eras) {  
      pdf(paste0("./analysis/data/derived_data/pre_process/",area, "_corrplot_for_", era, ".pdf"),    # File name
          width = 10, height = 10,  # Width and height in inches
          bg = "white",          # Background color
          colormodel = "cmyk")   # Color model 
      
      # Correlation of the data (removed discrete variables)
      corrplot(cor(cov_sites[[area]][[era]][,-c(1,3, 13)]),  method = "color", col = viridis(200), 
               type = "upper", 
               addCoef.col = "black", # Add coefficient of correlation
               tl.col = "black", tl.srt = 45, # Text label color and rotation
               number.cex = 0.7, # Size of the text labels
               cl.cex = 0.7) # Size of the color legend text
      
      dev.off()
    }
  }
}
