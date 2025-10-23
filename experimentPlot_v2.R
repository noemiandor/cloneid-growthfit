## ------------------------------------------------------------------
##  experimentPlot_v2.R
## ------------------------------------------------------------------

## ------------------------------------------------------------------
##  Helper functions
## ------------------------------------------------------------------
# get_script_path <- function() {
#   if (requireNamespace("rstudioapi", quietly = TRUE) &&
#       rstudioapi::isAvailable() &&
#       !is.null(rstudioapi::getSourceEditorContext()$path)) {
#     return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
#   } else {
#     args <- commandArgs(trailingOnly = FALSE)
#     path <- sub("--file=", "", args[grep("--file=", args)])
#     return(dirname(normalizePath(path)))
#   }
# }

get_script_path <- function() {
  
  ## a) If we’re knitting, use knitr’s input
  if (requireNamespace("knitr", quietly = TRUE)) {
    f <- knitr::current_input()
    if (!is.null(f) && nzchar(f)) {
      return(dirname(normalizePath(f)))
    }
  }
  
  ## b) If interactive in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    f <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(f) && nzchar(f)) {
      return(dirname(normalizePath(f)))
    }
  }
  
  ## c) If run by Rscript --file=...
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && nzchar(f)) {
    return(dirname(normalizePath(f)))
  }
  
  ## d) Fallback: current working directory
  getwd()
}

load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File '%s' does not exist.", filepath))
  }
  lines <- tryCatch(readLines(filepath, warn = FALSE),
                    error = function(e) stop("Error: Failed to read file '", filepath, "'. ", e$message))
  if (length(lines) == 0) {
    stop(sprintf("Error: File '%s' is empty.", filepath))
  }
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines))) {
    stop("Error: All lines must be in the format KEY=value, with valid variable names.")
  }
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  keys   <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  vars   <- setNames(values, keys)
  
  required_keys <- c("HOST", "DBNAME", "USER", "PASSWORD")
  missing_keys  <- setdiff(required_keys, names(vars))
  if (length(missing_keys) > 0) {
    stop("Error: Missing required keys: ", paste(missing_keys, collapse = ", "))
  }
  return(vars)
}

## ------------------------------------------------------------------
##  Main wrapper function
## ------------------------------------------------------------------
#' Run the experiment-plot pipeline.
#'
#' @param param_json_path Character. Path to the JSON parameter file.
#' @param db_creds_path   Character. Path to the db_creds.txt file.
#' @return Invisibly returns the main data frame produced (df), but most
#'         users will just care about the side-effects (plots and files).
#'
experimentPlot_V2 <- function(param_json_path, db_creds_path) {
  
  ## ----- libraries -----
  suppressPackageStartupMessages({
    library(RMySQL)
    library(jsonlite)
    library(ggplot2)
    library(stringr)
    library(dplyr)
    library(patchwork)
    library(rlang)
  })
  
  ## ----- set working directory to script’s location -----
  path <- get_script_path()
  setwd(path)
  
  ## ----- load config & DB credentials -----
  config  <- fromJSON(param_json_path,  flatten = FALSE, simplifyVector = FALSE)
  dbvars  <- load_db_vars(db_creds_path)
  
  ## ----- connect to the database -----
  db <- dbConnect(
    MySQL(),
    host     = dbvars["HOST"],
    user     = dbvars["USER"],
    password = dbvars["PASSWORD"],
    dbname   = dbvars["DBNAME"]
  )
  
  ## ----- original analysis pipeline (verbatim) -----
  # Retrieve Passaging & Media data
  x <- dbGetQuery(db, config$database$queries$passaging_query)
  m <- dbGetQuery(db, config$database$queries$media_query)
  
  filters <- config$filters[[1]]$ranges  # Read filter ranges
  
  ## (all original helper functions and processing code remain unchanged)
  recover_lineage <- function(row, data) {
    current_id <- row$last
    lineage <- c(current_id)
    
    while (!is.na(current_id) && current_id != row$first) {
      parent <- data$passaged_from_id1[which(data$id == current_id)]
      if (length(parent) == 0 || is.na(parent)) {
        stop(paste("Could not trace from", row$last, "to", row$first))
      }
      lineage <- c(parent, lineage)
      current_id <- parent
    }
    
    if (current_id != row$first) {
      stop(paste("Could not find", row$first, "starting from", row$last))
    }
    
    seeding_lineage <- lineage[data$event[match(lineage, data$id)] == "seeding"]
    
    filtered_x <- do.call(rbind, lapply(seq_along(seeding_lineage), function(i) {
      dfi <- data[data$id %in% seeding_lineage[i] | data$passaged_from_id1 %in% seeding_lineage[i], ]
      dfi$adjPass <- stringr::str_pad(i, width = 2)
      dfi$date <- as.Date(dfi$date)
      dfi$num_date <- as.numeric(as.Date(dfi$date))
      dfi$num_date <- dfi$num_date - min(dfi$num_date)
      dfi$intercept <- NaN
      dfi$g <- NaN
      if (nrow(dfi) < 2) return(dfi)
      fit <- lm(log(pmax(1, dfi$correctedCount)) ~ dfi$num_date)
      dfi$intercept <- exp(coef(fit)[1])
      dfi$g <- coef(fit)[2]
      dfi
    }))
    
    filtered_x$label_value    <- row$label
    filtered_x$sublabel_value <- row$label2
    return(filtered_x)
  }
  
  ## -- Process each filter entry
  db_ids <- lapply(seq_along(filters), function(i) recover_lineage(filters[[i]], x))
  df     <- do.call(rbind, db_ids)
  
  ## -- Group by seeding events
  group_by_seeding <- function(data) {
    seeds  <- data$id[data$event == "seeding"]
    groups <- lapply(seeds, function(si) {
      subset            <- data[data$id == si | data$passaged_from_id1 == si, ]
      subset$passage_id <- si
      subset
    })
    groups <- do.call(rbind, groups)
    groups <- groups[!is.na(groups$date), ]
    return(groups)
  }
  df <- group_by_seeding(df)
  
  ## -- Merge media information
  media_df <- m
  colnames(media_df)[colnames(media_df) == "id"] <- "media"
  df <- base::merge(df, media_df)
  
  ## -- Manual filters
  manual_filters <- config$manual_filters
  for (fi in manual_filters) {
    df <- df[!df[, fi$key] %in% fi$values, ]
  }
  
  ## -- Fit growth curves for each seed group
  df <- do.call(rbind, lapply(split(df, df$passage_id), function(fi) {
    if (nrow(fi) > 1) {
      fit <- lm(log(pmax(1, fi$correctedCount)) ~ fi$num_date)
      fi$intercept <- exp(coef(fit)[1])
      fi$g         <- coef(fit)[2]
      fi$pred      <- exp(predict(fit, newdata = fi))
    } else {
      fi$intercept <- NA
      fi$g         <- NA
      fi$pred      <- NA
    }
    fi
  }))
  
  ## -- Determine color variable
  color_var <- NaN
  if (config$plotting$color_by %in% colnames(df)) {
    color_var <- config$plotting$color_by
  }
  df[, color_var] <- as.character(df[, color_var])
  
  ## -- Convenience axis labels
  x_lab       <- if (!is.null(config$plotting$axis_labels$x_label))
    config$plotting$axis_labels$x_label else "Days"
  y_lab_cell  <- if (!is.null(config$plotting$axis_labels$y_label_cellcount))
    config$plotting$axis_labels$y_label_cellcount else "Cell Count"
  
  ## ------------------------------
  ##  Plotting
  ## ------------------------------
  
  #— enforce numeric order on your color_var across all facets —#
  raw_vals <- unique(df[[color_var]])
  nums     <- suppressWarnings(as.numeric(raw_vals))
  
  if (!any(is.na(nums))) {
    # sort by numeric value
    ordered_levels <- as.character(sort(nums))
  } else {
    # fallback to alphabetical
    ordered_levels <- sort(raw_vals)
  }
  
  # now turn the column into a factor with those levels:
  df[[color_var]] <- factor(df[[color_var]], levels = ordered_levels)
  
  if (config$plotting$make_growth_curves) {
    facet_row_var <- config$plotting$facet_grid$rows
    facet_col_var <- config$plotting$facet_grid$cols
    scales        <- config$plotting$facet_grid$scales
    ncol_wrap     <- 12
    
    plots   <- list()
    heights <- c()
    row_levels <- unique(df[[facet_row_var]])
    
    for (row_val in row_levels) {
      df_sub <- df %>% dplyr::filter(.data[[facet_row_var]] == row_val)
      n_facets       <- length(unique(df_sub[[facet_col_var]]))
      n_rows_needed  <- ceiling(n_facets / ncol_wrap)
      
      p <- ggplot(df_sub, aes(x = num_date,
                              y = correctedCount,
                              color = !!sym(color_var))) +
        geom_point() +
        geom_line(aes(y = pred, group = passage_id)) +
        scale_color_viridis_d() +
        facet_wrap(vars(!!sym(facet_col_var)), ncol = ncol_wrap,
                   scales = scales) +
        ggtitle(as.character(row_val)) +
        labs(x = x_lab, y = y_lab_cell) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots[[length(plots) + 1]] <- p
      heights <- c(heights, n_rows_needed)
    }
    
    final_plot <- wrap_plots(plots, ncol = 1, heights = heights)
    
    if (!is.null(config$plotting$titles$passage_vs_growth)) {
      final_plot <- final_plot +
        plot_annotation(title = "Growth Curves")
    }
    
    ggsave(config$plotting$output_path,
           plot   = final_plot,
           width  = 20, height = 20, units = "in")
    
    print(final_plot)          # <-- this line makes it appear in-line
    
  }
  
  ## -- Generation bar-plots ------------------------------------------
  zz <- df[!duplicated(df$passage_id), ]
  
  if (config$plotting$make_generation_barplot) {
    
    ## sort rows by sublabel_value, then by adjusted passage -------------
    zz <- zz %>% 
      mutate(passage_num = suppressWarnings(
        as.numeric(gsub("\\D", "", passage))  # "P06" → 6
      )) %>% 
      arrange(sublabel_value, adjPass)
    
    p2 <- ggplot(
      zz,
      aes(x     = sublabel_value,
          y     = g,
          group = adjPass,      # keeps bars separate
          fill  = !!sym(color_var))
    ) +
      facet_grid(rows = vars(label_value)) +
      geom_col(position = "dodge", colour = "black") +
      scale_y_continuous(expression(growth~rate~(day^{-1}))) +
      scale_fill_viridis_d(drop = FALSE) +
      labs(fill = config$plotting$color_by)
    
    print(p2)
  }
  
  ## -- Passage-vs-Growth plot
  if (config$plotting$make_passage_vs_growth) {
    p3 <- ggplot(df, aes(x = passage, y = g, color = !!sym(color_var))) +
      geom_col() +
      scale_color_viridis_d() +
      labs(x = if (!is.null(config$plotting$axis_labels$passage_label))
        config$plotting$axis_labels$passage_label else "Passage",
        y = if (!is.null(config$plotting$axis_labels$y_label_growth))
          config$plotting$axis_labels$y_label_growth else "Growth Rate",
        title = if (!is.null(config$plotting$titles$passage_vs_growth))
          config$plotting$titles$passage_vs_growth else
            "Passage vs Growth")
    print(p3)
  }
  
  ## ----- disconnect if requested -----
  if (config$database_settings$disconnect_on_completion) {
    dbDisconnect(db)
  }
  
  invisible(df)   # Return something useful without changing behaviour
}
