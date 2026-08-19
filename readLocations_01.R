library(tidyverse)
library(dplyr)
library(ggridges)
library(ggplot2)
library(tidyr)
library(forcats)
library(scales)




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


############## Get variant size ################

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



load_filtered_tables <- function(sample_ids) {
  data_list <- list()
  
  for (sample_id in sample_ids) {
    for (type in c("Linear", "Graph")) {
      dir_path  <- file.path("Data/exclude_Homopolymer/")
      file_name <- paste0(sample_id, "_", type, "_exclude_HP_filtered.csv")
      file_path <- file.path(dir_path, file_name)
      
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      df$size   <- get_variant_size_csv(df$Ref, df$Alt)
      df$ID <- sample_id
      df$type   <- type
      
      data_list[[paste(sample_id, type, sep = "_")]] <- df
    }
  }
  
  bind_rows(data_list)
}


filtered_data <- load_filtered_tables(sample_ids)


filtered_data_summary <- filtered_data %>%
  select(Chr, Start, End, ID, size) %>%
  unique()




filtered_data_summary <- filtered_data_summary %>%
  arrange(ID, Chr, Start) %>%
  group_by(ID, Chr) %>%
  mutate(
    new_group = Start - lag(End, default = first(Start) - 2) > 1,  # TRUE if NOT adjacent to previous row
    group_id  = cumsum(new_group)
  ) %>%
  group_by(ID, Chr, group_id) %>%
  mutate(size_c = sum(size, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-new_group, -group_id) %>%
  mutate(
    range = case_when(
      size_c >= 16                  ~ 8,
      size_c >= 11 & size_c <= 15   ~ 7,
      size_c >= 6  & size_c <= 10   ~ 6,
      size_c >= 0  & size_c <= 5    ~ 5,
      size_c >= -5  & size_c < 0    ~ 4,
      size_c >= -10 & size_c <= -6  ~ 3,
      size_c >= -15 & size_c <= -11 ~ 2,
      size_c <= -16                 ~ 1,
      TRUE ~ NA_real_
    )
  )


filtered_data_summary <- filtered_data_summary %>%
  mutate(range_label = factor(range, levels = 1:8,
                              labels = c("≤ -16", "-15 to -11", "-10 to -6", "-5 to 0",
                                         "0 to 5", "6 to 10", "11 to 15", "≥ 16")))





############## Get Position and CIGAR strings ##############

# Get all graph output files
graph_file <- list.files(
  "readLocs/Data",
  pattern = "_graphOutput_01\\.txt$",
  full.names = TRUE
)

# Function to process one graph/linear file pair
process_pair <- function(graph_file) {
  
  # Create matching linear filename
  linear_file <- str_replace(graph_file, "graphOutput", "linearOutput")
  
  # Extract metadata from filename
  file_base <- basename(graph_file) %>%
    str_remove("_graphOutput_01\\.txt")
  
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  sample_id <- parts[1]
  gene <- parts[2]
  
  # Everything after gene becomes location
  location <- paste(parts[3:length(parts)], collapse = "_")
  
  cols <- paste0("X", 1:9)
  
  # Read graph file
  Graph_df <- read_tsv(graph_file, show_col_types = FALSE, col_names = cols,) %>%
    #filter(X9 > 0, X5 == 60) %>%
    filter(X9 > 0) %>%
    select(1, 3, 4, 6) %>%
    rename("ReadName" = X1, "GChr" = X3, "GPos" = X4, "GCigar" = X6) %>%
    mutate(GPos = as.numeric(GPos), GChr = as.character(GChr)) %>%
    mutate(GChr = if_else(str_starts(GChr, "chr"),
                          GChr, paste0("chr", GChr))) %>%
    mutate(.after = 1, ID = sample_id, Gene = gene, Location = location) %>%
    unique()
  
  # Read linear file
  Linear_df <- read_tsv(linear_file, show_col_types = FALSE, col_names = cols,) %>%
    #filter(X9 > 0, X5 == 60) %>%
    filter(X9 > 0) %>%
    select(1, 3, 4, 6) %>%
    rename("ReadName" = X1, "LChr" = X3, "LPos" = X4, "LCigar" = X6) %>%
    mutate(LPos = as.numeric(LPos), LChr = as.character(LChr)) %>%
    mutate(LChr = if_else(str_starts(LChr, "chr"),
                          LChr, paste0("chr", LChr))) %>%
    mutate(.after = 1, ID = sample_id, Gene = gene, Location = location) %>%
    unique()
  
  # Combine and annotate
  combined_df <- Graph_df %>%
    full_join(Linear_df, by = c("ReadName", "ID", "Gene", "Location")) %>%
    mutate(Chr = str_extract(Location, "^[^_]+"),
           Start = as.numeric(str_extract(Location, "(?<=_)\\d+(?=-)")),
           End = as.numeric(str_extract(Location, "(?<=-)\\d+$"))) %>%
    mutate(in_both = !is.na(GPos) & !is.na(LPos),
           same_chr = GChr == Chr & LChr == Chr,
           position_match = abs(GPos - LPos) <= 5 & GChr == LChr,
           CIGAR_match = GCigar == LCigar & position_match == TRUE)
  
  return(combined_df)
}



############## Process all files and combine into one dataframe ############## 

all_combined_df <- map_dfr(graph_file, process_pair)

all_combined_filtered <- all_combined_df %>%
  filter(!(grepl("^([0-9]+[MN])+$", GCigar) & grepl("^([0-9]+[MN])+$", LCigar)))

all_combined_filtered <- all_combined_filtered %>%
  left_join(filtered_data_summary, by = c("Chr", "Start", "End", "ID"))



# Unmapped vs Mismatched
LPos_count <- all_combined_filtered %>%
  filter(!is.na(GPos)) %>%
  group_by(ID, Gene, size_c) %>%
  summarise(n_present = sum(!is.na(LPos)),
            n_missing = sum(is.na(LPos)),
            .groups = "drop")


same_chr_count <- all_combined_filtered %>%
  filter(!is.na(LPos), !is.na(GPos)) %>%
  group_by(ID, Gene, size_c) %>%
  summarise(n_diff = sum(!same_chr, na.rm = TRUE),
            n_same = sum(same_chr, na.rm = TRUE),
            .groups = "drop")






##### PossDiff plots
# Calculate bp distance of seq read from Variant location
DiffChr <- all_combined_df %>%
  filter(same_chr == FALSE)


SameChr <- all_combined_df %>%
  filter(same_chr == TRUE | is.na(same_chr))


same_pos_count <- all_combined_filtered %>%
  filter(!is.na(LPos), !is.na(GPos)) %>%
  filter(same_chr == TRUE | is.na(same_chr)) %>%
  group_by(ID, Gene, size_c) %>%
  summarise(n_diff = sum(!position_match, na.rm = TRUE),
            n_same = sum(position_match, na.rm = TRUE),
            .groups = "drop")



PosDiff_v2 <- SameChr %>%
  filter(position_match == FALSE) %>%
  mutate(GDist = GPos - Start, LDist = LPos - Start)



PosDiff_summary <- PosDiff_v2 %>%
  group_by(Gene, LDist) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(sign = ifelse(LDist < 0, "Downstream", "Upstream"),
         LDist = abs(LDist)) %>%
  filter(LDist >= 5)



######## Plots
PosDiff_Dist_plot <- ggplot(PosDiff_summary, aes(x = LDist, y = Gene,
                            size = count, colour = Gene,
                            shape = sign)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 10),
                        name = "Count",
                        breaks = c(1, 25, 50, 75, 100)) +
  scale_color_manual(values = palette_10, guide = "none") +
  scale_x_log10(labels = scales::comma, n.breaks = 10) +
  labs(x = "Distance of read from variant location (log10)", y = "") +
  theme_linedraw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15, face = "bold", family="Verdana"),
        axis.title.y = element_text(size = 15, face = "bold", family="Verdana"),
        axis.text.x = element_text(size = 12, family="Verdana"),
        axis.text.y = element_text(size = 12, face = "italic", family="Verdana"))



ggsave('plots/PosDiff_Dist_05.png', PosDiff_Dist_plot,
       width = 20, height = 10, units = "cm", dpi = 100)


