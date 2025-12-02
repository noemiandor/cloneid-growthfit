## ------------------------------------------------------------------
##  analyze_fitness_dynamics.R
## ------------------------------------------------------------------

get_script_path <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    f <- knitr::current_input()
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && nzchar(f)) return(dirname(normalizePath(f)))
  getwd()
}

load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) stop(sprintf("Error: File '%s' does not exist.", filepath))
  lines <- tryCatch(readLines(filepath, warn = FALSE), error = function(e) stop("Read failed"))
  valid_lines <- lines[grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines)]
  kv_pairs <- strsplit(valid_lines, "=", fixed = TRUE)
  keys   <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  setNames(values, keys)
}

analyze_fitness <- function(param_json_path, db_creds_path, output_csv = "data/fitness_dynamics_metrics.csv", diagnostic_pdf = "figures/fitness_diagnostics.pdf", recovery_plot_pdf = "figures/fitness_recovery_trends.pdf") {
  
  suppressPackageStartupMessages({
    library(RMySQL)
    library(jsonlite)
    library(dplyr)
    library(stringr)
    library(strucchange) 
    library(tidyr)
    library(ggplot2)
  })
  
  config  <- fromJSON(param_json_path, flatten = FALSE, simplifyVector = FALSE)
  dbvars  <- load_db_vars(db_creds_path)
  
  db <- dbConnect(MySQL(), host=dbvars["HOST"], user=dbvars["USER"], password=dbvars["PASSWORD"], dbname=dbvars["DBNAME"])
  
  raw_passaging <- dbGetQuery(db, config$database$queries$passaging_query)
  raw_media     <- dbGetQuery(db, config$database$queries$media_query)
  
  filters <- config$filters[[1]]$ranges
  
  recover_lineage <- function(row, data) {
    current_id <- row$last
    lineage <- c(current_id)
    while (!is.na(current_id) && current_id != row$first) {
      parent <- data$passaged_from_id1[which(data$id == current_id)]
      if (length(parent) == 0) stop(paste("Broken lineage:", row$last))
      lineage <- c(parent, lineage)
      current_id <- parent
    }
    seeding_lineage <- lineage[data$event[match(lineage, data$id)] == "seeding"]
    
    do.call(rbind, lapply(seq_along(seeding_lineage), function(i) {
      dfi <- data[data$id %in% seeding_lineage[i] | data$passaged_from_id1 %in% seeding_lineage[i], ]
      dfi$adjPass <- i
      dfi$date <- as.Date(dfi$date)
      dfi$num_date <- as.numeric(dfi$date) # Keep absolute date numeric first
      
      if (nrow(dfi) < 2) {
        dfi$g <- NA 
        dfi$g_r2 <- NA
      } else {
        # Fit per-passage growth relative to the start of THAT passage
        passage_time <- dfi$num_date - min(dfi$num_date)
        fit <- lm(log(pmax(1, dfi$correctedCount)) ~ passage_time)
        dfi$g <- coef(fit)[2]
        dfi$g_r2 <- summary(fit)$r.squared 
      }
      dfi$label_value <- row$label
      return(dfi)
    }))
  }
  
  df_list <- lapply(seq_along(filters), function(i) recover_lineage(filters[[i]], raw_passaging))
  df <- do.call(rbind, df_list)
  
  # Process Passages and Define Epochs
  df_passages <- df %>%
    filter(event == "seeding") %>%
    left_join(raw_media, by = c("media" = "id")) %>%
    arrange(label_value, date) %>%
    group_by(label_value) %>% 
    mutate(
      prev_media = lag(media),
      is_new_epoch = (media != prev_media) | is.na(prev_media),
      epoch_id = cumsum(is_new_epoch)
    ) %>%
    group_by(label_value, epoch_id) %>% 
    mutate(
      epoch_start_date = min(num_date),
      days_in_epoch = num_date - epoch_start_date,
      plot_group = paste(label_value, "Epoch", epoch_id, sep=" - ")
    ) %>%
    ungroup()
  
  lineage_metrics <- df_passages %>%
    group_by(label_value, epoch_id) %>%
    summarise(
      media_id = first(media),
      n_passages = n(),
      
      # Metrics for Deficit Calculation
      g_start = first(g),
      # We capture the LAST growth rate of the current epoch to serve as the baseline for the NEXT epoch
      g_end_of_epoch = last(g),
      
      # Calculate R_rec
      R_rec = if(n() > 2) coef(lm(g ~ days_in_epoch))[2] else NA,
      
      n_breakpoints = if(n() > 4) {
        tryCatch({
          bp <- breakpoints(g ~ 1, breaks = 2) 
          if(is.na(bp$breakpoints[1])) 0 else length(bp$breakpoints)
        }, error = function(e) 0)
      } else { 0 }
    ) %>%
    ungroup() %>%
    arrange(label_value, epoch_id) %>%
    group_by(label_value) %>%
    mutate(
      # UPDATED LOGIC: Use the last growth rate from the PREVIOUS epoch as baseline
      # This captures the state of adaptation just before the transition
      baseline_g = lag(g_end_of_epoch), 
      
      # Calculate Deficit without artificial NA replacement
      delta_F_init = (g_start - baseline_g) / baseline_g
    ) %>%
    ungroup()
  
  write.csv(lineage_metrics, output_csv, row.names = FALSE)
  
  # --- DIAGNOSTIC PLOTS ---
  pdf(diagnostic_pdf, width=10, height=8)
  
  # Updated P1 to reflect new baseline logic
  p1 <- ggplot(lineage_metrics, aes(x=baseline_g, y=delta_F_init)) +
    geom_point(alpha=0.6) +
    labs(title="Stability Check: Baseline Growth (Last of Prev Epoch) vs. Deficit",
         x="Baseline Growth Rate (g_end_prev)", y="Calculated Deficit (Delta F)") +
    theme_minimal()
  print(p1)
  
  p2 <- ggplot(df_passages, aes(x=num_date, y=g, color=as.factor(epoch_id))) +
    geom_line(alpha=0.5) + geom_point() +
    facet_wrap(~label_value, scales="free") +
    labs(title="Growth Rate Trajectories by Lineage",
         x="Absolute Date (Numeric)", y="Growth Rate (g)", color="Epoch") +
    theme_minimal()
  print(p2)
  
  p3 <- ggplot(df_passages, aes(x=g_r2)) +
    geom_histogram(bins=20, fill="steelblue", color="white") +
    labs(title="Quality of Per-Passage Growth Fits",
         x="R-squared of Exponential Fit", y="Count") +
    theme_minimal()
  print(p3)
  
  dev.off()
  
  # --- RECOVERY TREND PLOT (PAGINATED: MAX 9 PER PAGE) ---
  
  pdf(recovery_plot_pdf, width=12, height=10) # Good size for 3x3
  
  # 1. Get unique plot groups (Facet panels)
  unique_groups <- unique(df_passages$plot_group)
  n_groups <- length(unique_groups)
  
  # 2. Define page parameters
  plots_per_page <- 9
  n_pages <- ceiling(n_groups / plots_per_page)
  
  cat(sprintf("\nGenerating %d pages of plots (%d total panels)...\n", n_pages, n_groups))
  
  for (i in 1:n_pages) {
    # Calculate start/end indices for this page
    start_idx <- (i - 1) * plots_per_page + 1
    end_idx <- min(i * plots_per_page, n_groups)
    
    # Subset data for ONLY these groups
    page_groups <- unique_groups[start_idx:end_idx]
    df_page <- df_passages %>% filter(plot_group %in% page_groups)
    
    # Create the plot for this page
    p_page <- ggplot(df_page, aes(x=days_in_epoch, y=g, color=as.factor(epoch_id))) +
      geom_point(alpha=0.7, size=2) + 
      geom_smooth(method = "lm", se = FALSE, size = 1) +
      
      # Use the pre-calculated unique identifier for faceting
      facet_wrap(~plot_group, scales="free", ncol=3, nrow=3) +
      
      labs(title=sprintf("Fitness Recovery Trends per Epoch (Page %d of %d)", i, n_pages),
           subtitle="3x3 Grid - Separated by Lineage and Epoch ID",
           x="Days within Epoch", 
           y="Growth Rate (g)",
           color="Epoch ID") +
      theme_bw() +
      theme(legend.position="bottom")
    
    print(p_page)
  }
  
  dev.off()
  cat(sprintf("\nRecovery trend plot saved to: %s\n", recovery_plot_pdf))
  
  
  # --- REPORT ---
  cat("\n=== MANUSCRIPT REPORT ===\n")
  stress_epochs <- lineage_metrics %>% filter(epoch_id > 1, delta_F_init < 0)
  
  cat(sprintf("\n--- Initial Fitness Deficit ---\n"))
  if(nrow(stress_epochs) > 0) {
    mean_val <- mean(stress_epochs$delta_F_init, na.rm = TRUE) * 100
    min_val  <- min(stress_epochs$delta_F_init, na.rm = TRUE) * 100
    max_val  <- max(stress_epochs$delta_F_init, na.rm = TRUE) * 100
    cat(sprintf("Mean deficit: %.1f%%\n", mean_val))
    cat(sprintf("Range: %.1f%% – %.1f%%\n", min_val, max_val))
  } else {
    cat("N/A (No valid stress epochs detected)\n")
  }
  
  recovery_epochs <- lineage_metrics %>% filter(R_rec > 0)
  cat(sprintf("\n--- Rate of Recovery ---\n"))
  if(nrow(recovery_epochs) > 0) {
    max_rec <- max(recovery_epochs$R_rec, na.rm = TRUE)
    mean_rec <- mean(recovery_epochs$R_rec, na.rm = TRUE)
    cat(sprintf("Max Recovery: %.4f day^-1\n", max_rec))
    cat(sprintf("Mean Recovery: %.4f day^-1\n", mean_rec))
  } else {
    cat("N/A (No positive recovery)\n")
  }
  
  stepwise_lineages <- lineage_metrics %>% 
    group_by(label_value) %>%
    summarise(has_plateaus = any(n_breakpoints > 0)) %>%
    filter(has_plateaus == TRUE)
  
  cat(sprintf("\n--- Stepwise Plateaus ---\n"))
  cat(sprintf("Lineages with plateaus: %d out of %d\n", nrow(stepwise_lineages), length(unique(lineage_metrics$label_value))))
  
  if (config$database_settings$disconnect_on_completion) {
    dbDisconnect(db)
  }
}