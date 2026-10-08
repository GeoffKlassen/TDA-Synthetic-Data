library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(tidyselect)
library(htmlwidgets)
library(TDAstats)
library(patchwork)

# Step 1: Generate synthetic data ####
synthetic_3D_3000n <- generate_synthetic_dataset(n = 3000,
                                                 participant_count = 1,
                                                 d = 3,
                                                 Sigma_hi_corr = Sigma_3D_hi_corr,
                                                 Sigma_med_corr = Sigma_3D_med_corr,
                                                 Sigma_no_corr = Sigma_3D_no_corr)

synthetic_3D_5000n <- generate_synthetic_dataset(n = 5000,
                                                participant_count = 1,
                                                d = 3,
                                                Sigma_hi_corr = Sigma_3D_hi_corr,
                                                Sigma_med_corr = Sigma_3D_med_corr,
                                                Sigma_no_corr = Sigma_3D_no_corr,
                                                standardize_columns = FALSE)

create_ph_and_UMAP_de_plots_v2 <- function(df,
                                           tau_vector = c(1),
                                           E_vector = c(3),
                                           FD_vector = c(FALSE)) {
  # Initialize a list for collecting the Persistent Homology plots.
  # These plots will be tabulated and displayed in a grid for comparison.
  ph_plot_list <- list()
  
  participant_count <- length(unique(df$Participant_ID))
  n_per_participant <- length(unique(df$Response_Date))
  # Step 2: Generate Persistent Homology (PH) plots on the synthetic data ####
  ## 2a: extract column groups from the synthetic dataset ####
  column_groups <- extract_column_groups(df)
  
  all_column_names <- names(df)

  for (group in column_groups) {
    lap <- Sys.time()
    cat(paste0("Group: ", group, "\n"))
    
    # Fetch the subset of the dataframe that matches the current column 'group'
    columns_in_group <- all_column_names[grepl(paste0(group, "_[0-9]+$"), all_column_names)]
    
    cat(paste(columns_in_group, sep = ", "))
    cat("\n")
    
    title_text_for_ph_raw_plot <- paste0("Persistent Homology\nRaw Data\n", compress_names(c_names = columns_in_group))
    
    cat("Making Persistent Homology Plot for raw data...\n")
    # Create a PH plot of the raw data for the current column group
    ph_raw_plot <- make_persistent_homology_plot(dataset = df %>% select(all_of(columns_in_group)),
                                  title_text = title_text_for_ph_raw_plot,
                                  create_png = FALSE)
    print(ph_raw_plot)
    
    lap <- lap_timer(lap)
    for (tau in tau_vector) {
      for (E in E_vector) {
        for (FD in FD_vector) {
          cat("Making Delay Embedding dataframe...\n")
          # Create delay embedding (DE) dataframe
          de_df <- create_delay_embedding_df(df = df,
                                             id_col = "Participant_ID",
                                             date_col = "Response_Date",
                                             de_cols = columns_in_group,
                                             tau = tau,
                                             E = E)
          lap <- lap_timer(lap)
          
          title_text_for_ph_de_plot <- paste0("PH on DE Data\n",
                                              "Participants=", participant_count, ", ",
                                              "n per participant=", n_per_participant, "\n",
                                             compress_names(c_names = columns_in_group), "\n",
                                             "tau=", tau, ", E=", E, "\n",
                                             "FD=", FD)
          
          cat("Making Persistent Homology Plot for Delay Embedded Data...\n")
          ph_de_plot <- make_persistent_homology_plot(dataset = de_df %>% select(starts_with(group)),
                                             title_text = title_text_for_ph_de_plot,
                                             create_png = TRUE)
          print(ph_de_plot)
          ph_plot_list <- append(ph_plot_list, ph_de_plot)
          
          lap <- lap_timer(lap)
          
          cat("Making UMAP from Delay Embedded dataframe...\n")
          umapped_de_df <- create_UMAP(df = de_df,
                                       cols = names(de_df %>% select(starts_with(group))),
                                       n_dimensions = 3)
          
          lap <- lap_timer(lap)
          
          title_text_for_umap_de_plot <- paste0("UMAP on Delay Embedding Data\n",
                                                "Participant Count=", participant_count, ", ",
                                                "n_rows per participant=", n_per_participant, "\n",
                                                compress_names(c_names = columns_in_group), "\n",
                                                "tau=", tau, ", E=", E, "\n",
                                                "First Differences =", FD)
          cat("Making UMAP Plot...\n")
          umap_de_plot <- plot_ly(umapped_de_df$umap_output,
                               x = ~x,
                               y = ~y,
                               z = ~z,
                               type = "scatter3d",
                               mode = "markers",
                               marker = list(color = ~df$Response_Date,
                                             colorscale = create_colour_gradient(c("orange", "lightgrey", "blue")),
                                             # reversescale = TRUE,
                                             showscale = TRUE,
                                             size = 5)
          ) %>% layout(
            title = title_text_for_umap_de_plot
          )
          lap <- lap_timer(lap)
          
          print(umap_de_plot)
          
          cat("Saving UMAP...\n")
          cols_as_title <- compress_names(columns_in_group)
          htmlwidgets::saveWidget(widget = umap_de_plot,
                                  file = paste0("UMAP on DE, ", cols_as_title, ", FD=", include_first_diff, ", tau=", tau, ", E=", E, ".html")
          )
          lap <- lap_timer(lap)
        } # for FD in FD_vector
      } # for E in E_vector
    } # for tau in tau_vector
   
  } # for group in column_groups
  
  plot_data_list <- map(ph_plot_list, ~ ggplot_build(.x)$data[[2]])
  
  all_x <- unlist(map(plot_data_list, "x"))
  all_y <- unlist(map(plot_data_list, "y"))
  
  global_range <- c(0, 1.01*max(all_x, all_y))
  
  
  
  minimal_plots <- map(
    ph_plot_list,
    ~ .x +
      coord_cartesian(
        xlim = global_range,
        ylim = global_range
      ) +
      theme(
        plot.title      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",        # ← remove legend
        panel.border    = element_rect(color = "grey60", fill = NA, linewidth = 0.3)
      )
  )
  
  
  col_labels <- c("High correlation\nρ=0.96", "Medium correlation\nρ=0.5", "No correlation\nρ=0")
  col_labels_plot <- wrap_plots(
    map(col_labels, ~
          ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = .x, size = 4) +
          theme_void() +
          theme(
            plot.margin = margin(0, 0, 0, 0)   # remove outer margins
            #panel.spacing = unit(0, "pt")       # remove internal spacing
          )
    ),
    ncol = 3
  ) + plot_layout(heights = unit(12, "mm"))
  
  row_labels <- c(" ", "Walk", "Spring\nκ=0.005", "Spring\nκ=0.01", "Spring\nκ=0.025", "Spring\nκ=0.05", "Spring\nκ=0.1", "Noise", "Markov mean\nswitching (p=0.05)")
  row_labels_plot <- wrap_plots(
    map(row_labels, ~ ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = .x, angle = 90, size = 4) +
          theme_void() +
          theme(plot.margin = margin(0, 0, 0, 0))),
    ncol = 1
  ) + plot_layout(widths = unit(12, "mm"))  # force narrow width
  
  # Combine column labels + grid first
  grid_with_col_labels <- (col_labels_plot / wrap_plots(minimal_plots, ncol = 3)) +
    plot_layout(heights = c(0.08, 1))  # adjust label height
  
  # Then add row labels to the left
  final_plot <- row_labels_plot | grid_with_col_labels +
    plot_layout(widths = c(0.05, 1))   # adjust label column width
  final_plot
} # function

## Step 2b: Output PH plots to PNG files ####
  
# Step 3: Generate Delay Embedding (DE) data on synthetic data ####

# Step 4: Generate UMAP plots on DE data ####
## Step 4b: Output PH plots to HTML files ####

create_ph_and_UMAP_de_plots_v2(df = synthetic_3D_500n)

create_ph_and_UMAP_de_plots_v2(df=synthetic_3D_3000n,
                               tau_vector = c(2),
                               E_vector = c(3))
library(plotly)
umap_iteration_plot <- plot_ly(as.data.frame(embeddings_by_epoch[[100]]),
                                  x = ~V1,
                                  y = ~V2,
                                  z = ~V3,
                                  type = "scatter3d",
                                  mode = "markers",
                                  marker = list(color = ~synthetic_3D_500n$Response_Date,
                                                colorscale = create_colour_gradient(c("orange", "lightgrey", "blue")),
                                                # reversescale = TRUE,
                                                showscale = TRUE,
                                                size = 5)
)
print(umap_iteration_plot)

library(reticulate)
reticulate::use_python("C:/Users/geoff/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)
reticulate::py_install('plotly')

plotly::save_image(umap_iteration_plot, file = "test export umap iteration plot.png") 
