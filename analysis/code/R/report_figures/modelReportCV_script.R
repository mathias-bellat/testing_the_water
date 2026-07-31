modelReportCV <- function (model, folder, test = NULL, type = NULL, response_curves = FALSE, 
                           only_presence = FALSE, jk = FALSE, env = NULL, clamp = TRUE, 
                           permut = 10, verbose = TRUE) 
{
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    cli::cli_abort("Please install package {.pkg kableExtra} to use this function", 
                   call = NULL)
  }
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    cli::cli_abort("Please install package {.pkg htmltools} to use this function", 
                   call = NULL)
  }
  if (file.exists(file.path(getwd(), folder))) {
    msg <- cli::cli_text(cli::col_red(cli::symbol$fancy_question_mark), 
                         " The folder {.file {folder}} already exists,", 
                         " do you want to overwrite it?")
    continue <- utils::menu(choices = c("Yes", "No"), title = msg)
  }
  else {
    continue <- 1
  }
  if (continue == 1) {
    template <- file.path(getwd(), "./code/R/modelReportCV.Rmd")
    folder <- file.path(getwd(), folder)
    plot_folder <- file.path(folder, "plots")
    dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
    species <- gsub(" ", "_", tolower(model@data@species))
    output_file <- paste0(species, ".html")
    if (verbose) 
      cli::cli_text("\f", cli::rule(
        left = paste("Model Report - class:", class(model)[1]),
        right = cli::style_italic(model@data@species),
        line_col = "#4bc0c0", col = "#f58410", width = 80))
    rmarkdown::render(template, output_file = output_file, 
                      output_dir = folder, params = list(model = model, 
                                                         type = type, test = test, folder = folder, plot_folder = plot_folder, 
                                                         env = env, jk = jk, response_curves = response_curves, 
                                                         only_presence = only_presence, clamp = clamp, 
                                                         permut = permut, verbose = verbose), quiet = TRUE)
    utils::browseURL(file.path(folder, output_file))
  }
}