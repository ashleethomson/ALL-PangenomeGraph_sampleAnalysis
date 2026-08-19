library(tidyverse)
library(dplyr)
library(ggplot2)



colours = c("#62109F", "#FE6244")


sample_ids <- c(
  "A2189-1026",
  "A2241-10841",
  "A2638-12600",
  "A4551-18835",                  
  "A5058-21220",
  "A5322-22937",
  "A5560-24589",
  "A5654-25234",
  "A6314-28788",
  "ADU-0705-201321",
  "ADU-0919-201785",
  "ADU-1262-202602",
  "ADU-1281-202610",
  "AYA-0846-201738",
  "AYA-0892-201728",
  "AYAII-0454-DIA1-BM",
  "CHI-0796-202173",
  "CHI-1383-202781")




get_variant_size_csv <- function(ref, alt) {
  size <- mapply(function(r, a) {
    # Deletion: Alt is "-"
    if (a == "-") {
      return(-nchar(r))
    }
    # Insertion: Ref is "-"
    if (r == "-") {
      return(nchar(a))
    }
    # Complex indel (e.g. GACGGG -> A, or delins)
    return(nchar(a) - nchar(r))
  }, ref, alt)
  
  return(as.integer(size))
}


load_unfiltered_tables <- function(sample_ids) {
  data_list <- list()
  
  for (sample_id in sample_ids) {
    for (type in c("Linear", "Graph")) {
      dir_path  <- file.path("Data/exclude_Homopolymer/")
      file_name <- paste0(sample_id, "_", type, "_exclude_HP.csv")
      file_path <- file.path(dir_path, file_name)
      
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      df$size   <- get_variant_size_csv(df$Ref, df$Alt)
      df$sample <- sample_id
      df$type   <- type
      
      data_list[[paste(sample_id, type, sep = "_")]] <- df
    }
  }
  
  bind_rows(data_list)
}



load_filtered_tables <- function(sample_ids) {
  data_list <- list()
  
  for (sample_id in sample_ids) {
    for (type in c("Linear", "Graph")) {
      dir_path  <- file.path("Data/exclude_Homopolymer/")
      file_name <- paste0(sample_id, "_", type, "_exclude_HP_filtered.csv")
      file_path <- file.path(dir_path, file_name)
      
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      df$size   <- get_variant_size_csv(df$Ref, df$Alt)
      df$sample <- sample_id
      df$type   <- type
      
      data_list[[paste(sample_id, type, sep = "_")]] <- df
    }
  }
  
  bind_rows(data_list)
}


filtered_data <- load_filtered_tables(sample_ids)


filtered_data_v2 <- filtered_data %>%
  select(sample, VAF, Gene, AA.Change, size, type)



# plot - VAF vs size
VAFvsINDELsize <- ggplot(filtered_data_v2, aes(x = size, y = VAF, colour = type, shape = type)) +
  geom_point(alpha = 0.7, size = 3) +
  #geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Graph" = "#62109F", "Linear" = "#FE6244")) +
  scale_shape_manual(values = c("Graph" = 16, "Linear" = 17)) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(c(0, 0.05))) +
  labs(x = "INDEL Size (bp)", y = "VAF", colour = "Type", shape  = "Type") +
  theme_linedraw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", family="Verdana"),
        axis.title.y = element_text(size = 18, face = "bold", family="Verdana"),
        axis.text.x = element_text(size = 13, face = "bold", family="Verdana"),
        axis.text.y = element_text(size = 13, face = "bold", family="Verdana"))



ggsave("VariantLength/VAFvsINDELsize_04.png", VAFvsINDELsize,
       width = 20, height = 10, units = "cm", dpi = 100)


