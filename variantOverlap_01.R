
library(VennDiagram)
library(dplyr)
library(tidyverse)
library(here)
library(purrr)



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
process_sample <- function(s) {
  
  # Read data
  Graph <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP.csv")) %>%
    select(AAChange.refGene)
    ##unique()
  
  Linear <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP.csv")) %>%
    select(AAChange.refGene)
    ##unique()
  
  # Overlap
  both <- intersect(Graph, Linear)
  only_graph  <- setdiff(Graph, Linear)
  only_linear <- setdiff(Linear, Graph)
  
  # Return counts
  data.frame(
    Sample = s,
    graph_total = length(Graph$AAChange.refGene),
    only_graph  = length(only_graph$AAChange.refGene),
    linear_total = length(Linear$AAChange.refGene),
    only_linear = length(only_linear$AAChange.refGene),
    both = length(both$AAChange.refGene)
  )
}

process_sample_KG <- function(s) {
  
  Graph <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP_filtered.csv")) %>%
    select(AA.Change)
    ##unique()
  
  Linear <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP_filtered.csv")) %>%
    select(AA.Change)
    ##unique()
  
  both <- intersect(Graph, Linear)
  only_graph  <- setdiff(Graph, Linear)
  only_linear <- setdiff(Linear, Graph)
  
  data.frame(
    Sample = s,
    graph_total = length(Graph$AA.Change),
    only_graph  = length(only_graph$AA.Change),
    linear_total = length(Linear$AA.Change),
    only_linear = length(only_linear$AA.Change),
    both = length(both$AA.Change)
  )
}


KG_counts <- read_delim("KG_counts.txt")


##### Process
all_counts <- map_dfr(samples, process_sample)

all_counts_KG <- map_dfr(samples, process_sample_KG)


###### Calculate values
venn_totals <- all_counts %>%
  summarise(
    graph_total = sum(graph_total),
    only_graph  = sum(only_graph),
    linear_total = sum(linear_total),
    only_linear = sum(only_linear),
    both = sum(both)
  )

venn_totals_KG <- KG_counts %>%
  summarise(
    graph_total = sum(graph_total),
    only_graph  = sum(only_graph),
    linear_total = sum(linear_total),
    only_linear = sum(only_linear),
    both = sum(both)
    )



###### Draw Venn diagram with classic Venn layout (not Euler)
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = venn_totals$graph_total,
  area2 = venn_totals$linear_total,
  cross.area = venn_totals$both,
  category = c("Graph", "Linear"),
  fill = c("#62109F", "#FE6244"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2, # Inside font size
  cat.just = list(c(0, -3), c(1.5, -3)),
  lwd = 3,  # circle line width
  fontface = rep("bold"),
  main = "Venn Diagram of Overlap",
  scaled = FALSE,
  euler.d = FALSE # This forces classic Venn layout
)




grid.newpage()
venn.plot_KG <- draw.pairwise.venn(
  area1 = venn_totals_KG$graph_total,
  area2 = venn_totals_KG$linear_total,
  cross.area = venn_totals_KG$both,
  category = c("Graph", "Linear"),
  fill = c("#62109F", "#FE6244"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2, # Inside font size
  cat.just = list(c(0, -3), c(1.5, -3)),
  lwd = 3,  # circle line width
  fontface = rep("bold"),
  main = "Venn Diagram of Overlap",
  scaled = FALSE,
  euler.d = FALSE # <--- This forces classic Venn layout
)
grid.newpage()




plot_file <- paste0("plots/Venn/overlapCounts_all_07.jpg")
plot_file_KG <- paste0("plots/Venn/overlapKeyGenes_all_06.jpg")


ggsave(plot_file, venn.plot, width = 15, height = 15, units = "cm", dpi = 100)
ggsave(plot_file_KG, venn.plot_KG, width = 15, height = 15, units = "cm", dpi = 100)

