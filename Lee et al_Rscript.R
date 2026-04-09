##### Manuscript Script #######

#### Environmental data (Figure 1, Table S1, Table S2) ####

folder_names <- c("a1", "a2", "e1", "e2")

# A function to load only CH₄ data from a specific directory
read_and_combine_CH4 <- function(folder_path, folder_name) {
  file_list <- list.files(folder_path, pattern = "\\.CSV$", full.names = TRUE, recursive = TRUE)
  
  if (length(file_list) == 0) {
    warning(paste("No file in the folder:", folder_path))
    return(NULL)
  }
  combined_data <- lapply(file_list, function(file) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    date <- sub(".*/(2024[0-9]{4})\\.CSV$", "\\1", file, ignore.case = TRUE)
    data$Date <- as.Date(date, format="%Y%m%d")
    
    # Select the specific column (ex. CH4)
    if ("CH4" %in% colnames(data)) {
      data <- data %>% dplyr::select(Date, CH4) %>% mutate(Source = folder_name)
    } else {
      warning(paste("No CH4 column in the file", file))
      return(NULL)
    }
    
    return(data)
  }) %>% bind_rows()
  
  return(combined_data)
}

CH4_data_list <- mapply(read_and_combine_CH4, folders, folder_names, SIMPLIFY = FALSE)

# Remove the outliers: CO2 above 2000 ppm
CH4_data_filtered <- CH4_data %>% filter(CH4 <= 2000)

# Graph
ggplot(CH4_summary_filtered, aes(x = Date, y = mean_CH4, color = Group, group = Group)) +
  geom_line(size = 1) +  
  geom_point(size = 3) +  
  # geom_errorbar(aes(ymin = mean_CH4 - sd_CH4, ymax = mean_CH4 + sd_CH4), width = 0.2) + 
  labs(title = "Daily Mean CO2 Concentration (Grouped: aCO2 vs eCO2)",
       x = "Date",
       y = "Mean CO2 Concentration (ppm)",
       color = "Group") +
  scale_color_manual(values = c("eCO2" = "#E69F00", 
                                "aCO2" = "#009E73")) + 
  scale_y_continuous(
    limits = c(0, 800),
    breaks = seq(0, 800, by =200),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),

    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
  )



#### Table S2
folders <- c("a1", "a2", "e1", "e2")

start_date_str <- "20241105"
end_date_str   <- "20241125"
start_date <- ymd(start_date_str)
end_date <- ymd(end_date_str)

process_folder_data <- function(folder_name) {
  full_path <- file.path(base_dir, folder_name)
  
  if (!dir.exists(full_path)) {
    warning(paste("No folder:", full_path))
    return(NULL)
  }
  
  file_list <- list.files(path = full_path, pattern = "\\.[cC][sS][vV]$", full.names = TRUE)
  data_list <- list()
  
  for (file in file_list) {
    filename <- basename(file)
    date_part <- str_extract(filename, "\\d{8}")
    
    if (!is.na(date_part)) {
      file_date <- ymd(date_part)
      
      if (file_date >= start_date & file_date <= end_date) {
        temp_data <- read_csv(file, show_col_types = FALSE)
        
        if (all(c("CH4", "CH5", "CH6") %in% names(temp_data))) {
          selected_data <- temp_data %>%
            dplyr::select(CH4, CH5, CH6) %>%
            mutate(Folder = folder_name, Date = file_date)
          
          data_list[[length(data_list) + 1]] <- selected_data
        }
      }
    }
  }
  
  if (length(data_list) > 0) {
    return(bind_rows(data_list))
  } else {
    return(NULL)
  }
}


all_data <- purrr::map_dfr(target_folders, process_folder_data)

if (!is.null(all_data) && nrow(all_data) > 0) {
  
  # Remove the outliers: temperature above 30
  all_data <- all_data %>% 
    filter(CH5 <= 30)
  
  stats_long <- all_data %>%
    pivot_longer(cols = c("CH4", "CH5", "CH6"), names_to = "Channel", values_to = "Value") %>%
    group_by(Date, Folder, Channel) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SE = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    # 2) "Mean ± SE" string
    mutate(
      Result_Str = paste0(sprintf("%.2f", Mean), " ± ", sprintf("%.2f", SE))
    )
  stats_long$Channel <- factor(stats_long$Channel, levels = c("CH5", "CH6", "CH4"))
  
  final_table <- stats_long %>%
    dplyr::select(Date, Folder, Channel, Result_Str) %>%
    arrange(Folder, Channel) %>% 
    pivot_wider(
      names_from = c(Folder, Channel),
      values_from = Result_Str,
      names_sep = "_"
    ) %>%
    arrange(Date) 
  
  desired_order <- c("Date")
  for (f in target_folders) {
    for (ch in c("CH5", "CH6", "CH4")) {
      col_name <- paste(f, ch, sep = "_")
      if (col_name %in% names(final_table)) {
        desired_order <- c(desired_order, col_name)
      }
    }
  }
  
  final_table <- final_table %>% dplyr::select(all_of(desired_order))
  
  output_filename <- "sensor_mean_se_matrix.xlsx"
  output_path <- file.path(base_dir, output_filename)
  write_xlsx(final_table, path = output_path)
  
}

#### Quality control of bioassay data ####
# Filtering for samples with successful M. persicae settlement
filter_condition <- function(df) (
  df %>% filter(
    (herbivory %in% c("ctrl")) |
      herbivory %in% c("spo")|
      (herbivory == "gpa" & gpa_num > 3) |
      (herbivory == "both" & gpa_num > 3 )
  )
)


#### Effect size (Figure 2, Table S3) ####
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

subset_data <- subset(data, !is.na(plant_freshmass))
subset_data <- subset(data, 
                      (herbivory == "ctrl") |
                        (herbivory == "spo" & !is.na(spo_mass)) |
                        (herbivory == "gpa" & gpa_num > 3) |
                        (herbivory == "both" & !is.na(spo_mass) & gpa_num > 3)
)
model <- lmerTest::lmer(log(plant_freshmass) ~ co2 * her_spo * her_gpa + exp+ (1|chamber_id),
                              data = subset_data)
model <- lmerTest::lmer(log(spo_mass) ~ co2 * her_gpa + exp + (1|chamber_id),
                              data = subset_data)
model <- lmerTest::lmer(log(gpa_num) ~ co2 * her_spo + exp + (1|chamber_id),
                              data = subset_data)



#### Regression between feature abundance and S. litura larval mass (Table S5) ####
library(tidyverse)
library(broom)

target_ids <- annot_unique_EB$id
sample_names <- mmo$metadata$sample_col

long_features <- mmo$feature_data %>%
  filter(id %in% target_ids) %>%
  select(id, all_of(sample_names)) %>% 
  pivot_longer(
    cols = -id, 
    names_to = "sample_col", 
    values_to = "abundance"
  )

# Integration with mmo$metadata
merged_data <- long_features %>%
  left_join(mmo$metadata %>% select(sample_col, spo_mass), by = "sample_col")

merged_data <- merged_data %>%
  mutate(spo_mass = as.numeric(spo_mass)) %>% 
  filter(!is.na(abundance) & !is.na(spo_mass))

regression_results <- merged_data %>%
  group_by(id) %>%
  nest() %>% 
  mutate(
    model = map(data, ~ lm(abundance ~ spo_mass, data = .)),
    tidied = map(model, tidy) 
  ) %>%
  unnest(tidied) %>%
  select(-data, -model)

final_results <- regression_results %>%
  filter(term == "spo_mass") %>%
  arrange(p.value)

# Merge annotation data to final_results based on matching IDs.
final_results_annot <- final_results %>%
  left_join(mmo$sirius_annot, by = "id")

 