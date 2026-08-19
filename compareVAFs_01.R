library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)
library(GGally)
library(ggrepel)
library(paletteer)


colours = c("#62109F", "#FE6244")
palette_gen <- colorRampPalette(colours)

palette_7 <- palette_gen(7)
palette_11 <- palette_gen(11)


info <- read_tsv("sampleInfo.tsv")

samples <- c(
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




#### Functions
Graph_samples <- function(s) {
  
  # Read data
  Graph <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP_filtered.csv")) %>%
    select(Gene, AA.Change, VAF) %>%
    rename(Graph_VAF = VAF) %>%
    mutate(.before = 1, ID = s) 
}


Linear_samples <- function(s) {
  # Read data
  Linear <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP_filtered.csv")) %>%
    select(Gene, AA.Change, VAF) %>%
    rename(Linear_VAF = VAF) %>%
    mutate(.before = 1, ID = s) 
}


process_samples <- function(s) {
  
  # Read data
  Graph <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP_filtered.csv")) %>%
    select(Gene, AA.Change, VAF) %>%
    rename(Graph_VAF = VAF) %>%
    mutate(.before = 1, ID = s) 

  # Read data
  Linear <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP_filtered.csv")) %>%
    select(Gene, AA.Change, VAF) %>%
    rename(Linear_VAF = VAF) %>%
    mutate(.before = 1, ID = s) 
  
  
  VAFs <- full_join(by = c("ID", "Gene", "AA.Change"), Graph, Linear)
}






##### Process samples
all_VAFs <- map_dfr(samples, process_samples)

Graph_VAFs <- map_dfr(samples, Graph_samples)
Linear_VAFs <- map_dfr(samples, Linear_samples)


Graph_VAFs <- Graph_VAFs %>%
  group_by(ID, Gene, AA.Change) %>%
  slice_max(order_by = Graph_VAF, n = 1, with_ties = FALSE) %>%
  ungroup()

Linear_VAFs <- Linear_VAFs %>%
  group_by(ID, Gene, AA.Change) %>%
  slice_max(order_by = Linear_VAF, n = 1, with_ties = FALSE) %>%
  ungroup()


VAFs <- full_join(by = c("ID", "Gene", "AA.Change"), Graph_VAFs, Linear_VAFs)


VAFs$highlight <- ifelse(rowSums(is.na(VAFs[, 4:5])) > 0, "NoMatch", "Match")


VAFs_02 <- full_join(by = "ID", VAFs, info)



VAFs_long <- pivot_longer(VAFs_02, cols = c(Graph_VAF, Linear_VAF),
                         names_to = "Source",
                         values_to = "Value") %>%
  mutate(Xpos = ifelse(highlight == "NoMatch", 1, 2),
         nudge_dir = ifelse(Xpos == 1, -0.3, 0.3))



VAFs_long_02 <- pivot_longer(VAFs_02, cols = c(Graph_VAF, Linear_VAF),
                          names_to = "Source",
                          values_to = "Value")


####### Plots

## Plots by Age

VAFplot_CHI <- ggparcoord(VAFs_02 %>% filter(AgeGroup %in% c("CHI")), columns = c(5,4),
           groupColumn = "ID",
           scale = "globalminmax",
           boxplot = FALSE,
           mapping = ggplot2::aes(linewidth = 0.5),
           showPoints = FALSE) +
  geom_point(data = VAFs_long_02 %>% filter(AgeGroup %in% c("CHI")),
             aes(x = Source, y = Value, group = ID, fill = highlight, size = highlight),
             shape = 21,
             colour = "lightgrey",
             stroke = 0.1,
             inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(rep("darkgrey", n_distinct(VAFs_02$ID)), unique(VAFs_02$ID))) +
  scale_fill_manual(values = c("Match" = "darkgrey", "NoMatch" = "red")) +
  scale_size_manual(values = c("Match" = 1, "NoMatch" = 3)) +
  scale_linewidth_identity() +
  labs(title = "CHI",
       y = "Variant Allele Frequency (VAF)", 
       x = "") +
  scale_x_discrete(labels = c("Linear_VAF" = "Linear","Graph_VAF" = "Graph")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0.1, 1.0)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 18, face = "bold", family="Verdana"),
        axis.title.x = element_text(size = 18, face = "bold", family="Verdana"),
        axis.title.y = element_text(size = 18, face = "bold", family="Verdana"),
        axis.text.x = element_text(size = 13, face = "bold", family="Verdana"),
        axis.text.y = element_text(size = 13, face = "bold", family="Verdana"))



VAFplot_AYA_ADU <- ggparcoord(VAFs_02 %>% filter(AgeGroup %in% c("AYA", "ADU")), columns = c(5,4),
           groupColumn = "ID",
           scale = "globalminmax",
           boxplot = FALSE,
           mapping = ggplot2::aes(linewidth = 0.5),
           showPoints = FALSE) +
  geom_point(data = VAFs_long_02 %>% filter(AgeGroup %in% c("AYA", "ADU")),
             aes(x = Source, y = Value, group = ID, fill = highlight, size = highlight),
             shape = 21,
             colour = "lightgrey",
             stroke = 0.1,
             inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(rep("darkgrey", n_distinct(VAFs_02$ID)), unique(VAFs_02$ID))) +
  scale_fill_manual(values = c("Match" = "darkgrey", "NoMatch" = "red")) +
  scale_size_manual(values = c("Match" = 1, "NoMatch" = 3)) +
  scale_linewidth_identity() +
  labs(title = "AYA & ADU",
       y = "Variant Allele Frequency (VAF)", 
       x = "") +
  scale_x_discrete(labels = c("Linear_VAF" = "Linear","Graph_VAF" = "Graph")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0.1, 1.0)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 18, face = "bold", family="Verdana"),
        axis.title.x = element_text(size = 18, face = "bold", family="Verdana"),
        axis.title.y = element_text(size = 18, face = "bold", family="Verdana"),
        axis.text.x = element_text(size = 13, face = "bold", family="Verdana"),
        axis.text.y = element_text(size = 13, face = "bold", family="Verdana"))



ggsave('plots/VAFs_CHI_05.png', VAFplot_CHI,
       width = 15, height = 20, units = "cm", dpi = 100)

ggsave('plots/VAFs_AYA_ADU_05.png', VAFplot_AYA_ADU,
       width = 15, height = 20, units = "cm", dpi = 100)










