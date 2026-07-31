
plotROC_kfold <- function (model, test = NULL, show_folds = NULL) 
{
  if (!requireNamespace("plotROC", quietly = TRUE)) {
    cli::cli_abort("Please install package {.pkg plotROC} to use this function", 
                   call = NULL)
  }
  
  if (!inherits(model, "SDMmodelCV")) {
    cli::cli_abort("This function is designed for k-fold cross-validation models (SDMmodelCV objects)", 
                   call = NULL)
  }
  
  if (inherits(model@models[[1]]@model, "Maxent")) {
    type <- "raw"
  } 
  else {
    type <- "link"
  }
  
  df_list <- list()
  auc <- numeric()
  all_points <- c()
  all_preds <- c()
  labels_cv <- c()
  
  for (i in seq_along(model@models)) {
    fold_model <- model@models[[i]]
    
    df_list[[i]] <- data.frame(set = paste0("K",i), pa = fold_model@data@pa, pred = predict(fold_model, 
                                                                                            data = fold_model@data, type = type), stringsAsFactors = FALSE)
    
    all_preds <- c(all_preds, df_list[[i]]$pred)
    all_points <- c(all_points, df_list[[i]]$pa)
    auc[[i]] <- auc(fold_model)
    labels_cv[i] <- paste0("K-fold ", i, " ",  round(auc[[i]], 4))
  }  
  
  df <- data.frame(
    set = rep("Train mean", length(all_points)),
    stringsAsFactors = FALSE
  )
  df <- as.data.frame(cbind(df, all_points, all_preds))
  colnames(df) <- c("set", "pa", "pred") 
  
  labels <- paste0("Mean train ",round(SDMtune::auc(model),4)) 
  
  if (!is.null(test)) {
    df_test <- data.frame(set = "Test", pa = test@pa, pred = predict(model, 
                                                                     data = test, type = type), stringsAsFactors = FALSE)
    df_test_list <- list()
    auc_test_list <- numeric()
    all_points_test <- c()
    all_preds_test <- c()
    labels_cv_test <- c()
    
    for (i in seq_along(model@models)) {
      fold_model <- model@models[[i]]
      df_test_list[[i]] <- data.frame(set = paste0("K",i, " test"), pa = test@pa, 
                                      pred = predict(fold_model, data = test, type = type), 
                                      stringsAsFactors = FALSE)
      all_preds_test <- c(all_preds_test, df_test_list[[i]]$pred)
      all_points_test <- c(all_points_test, df_test_list[[i]]$pa)
      auc_test_list[[i]] <- auc(fold_model)
    }
    
    df_test <- data.frame(
      set = rep("Test mean", length(all_points_test)),
      stringsAsFactors = FALSE
    )
    df_test <- as.data.frame(cbind(df_test, all_points_test, all_preds_test))
    colnames(df_test) <- c("set", "pa", "pred") 
    
    df <- rbind(df, df_test) 
    label_test <- paste0("Mean test ",round(SDMtune::auc(model, test = test),4))
    labels <- c(label_test, labels)
  }
  
  if (!is.null(show_folds)) {
    df_folds <- do.call(rbind, df_list) 
    df <- rbind(df, df_folds) 
    labels <- c(labels_cv, labels)
  }  
  
  my_plot <- ggplot(df, aes(m = .data$pred, d = .data$pa, group = .data$set)) + 
    plotROC::geom_roc(n.cuts = 0, aes(color = .data$set), size = 0.5) + 
    ggplot2::scale_color_discrete(name = "AUC",labels = labels) + 
    ggplot2::geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "grey", linetype = 2) + 
    ggplot2::labs(x = "False Positive Rate", y = "True Positive Rate") + 
    ggplot2::coord_fixed() + ggplot2::theme_minimal() + 
    ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))
  if (!is.null(test)) {
    my_plot <- my_plot + ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE))
  }
#  plot(my_plot)
  sd_auc <- sd(auc)
  sd_auc_test <- sd(SDMtune::auc(model, test = test))
  
  
  # Afficher les statistiques dans la console
  cat("Cross-Validation ROC Statistics:\n")
  cat("Mean AUC train:", round(SDMtune::auc(model),4), "\n")
  cat("Standard deviation AUC train:", round(sd_auc, 4), "\n")
  cat("Min AUC train:", round(min(auc), 4), "\n")
  cat("Max AUC train:", round(max(auc), 4), "\n")
  if (!is.null(test)) {
    cat("Mean AUC test:", round(SDMtune::auc(model, test = test),4), "\n")
    cat("Standard deviation AUC test:", round(sd_auc_test, 4), "\n")
    cat("Min AUC test:", round(min(auc_test_list), 4), "\n")
    cat("Max AUC test:", round(max(auc_test_list), 4), "\n")
  }  
  return(my_plot)
  
}
