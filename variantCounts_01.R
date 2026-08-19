

library(ggplot2)
library(dplyr)
library(here)
library(stringr)
library(tidyverse)


######### Graph code #######

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



process_sample <- function(s) {
  
  # Read data
  Graph_Counts <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP.csv"))
  Linear_Counts <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP.csv"))
  
  Graph_Counts_KG <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Graph_exclude_HP_filtered.csv"))
  Linear_Counts_KG <- read_delim(paste0("Data/exclude_Homopolymer/", s,"_Linear_exclude_HP_filtered.csv"))
  
  counts <- data.frame(
    Source = c("Graph_Counts", "Linear_Counts"),
    Counts = c(nrow(Graph_Counts), nrow(Linear_Counts)),
    KeyGenes = c(nrow(Graph_Counts_KG), nrow(Linear_Counts_KG)),
    ID = s)
}



allCounts <- map_dfr(samples, process_sample)


countsTotal <- allCounts %>%
  group_by(Source) %>%
  summarise(
    Counts = sum(Counts),
    KeyGenes = sum(KeyGenes),
    .groups = "drop"
  )


counts_long <- countsTotal %>%
  pivot_longer(cols = c(Counts, KeyGenes), names_to = "Category", values_to = "Value") %>%
  mutate(Source = str_remove(Source, "_Counts$"))

facetLabels <- list(
  'Counts'="Total Counts",
  'KeyGenes'="Counts of Key Genes"
)

facetLabels <- c("Total Counts", "Counts of Key Genes")
names(facetLabels) <- c("Counts", "KeyGenes")



# Plot: facet by Category, each with its own Y-axis
plot_counts <- ggplot(counts_long, aes(x = Source, y = Value, fill = Source)) +
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Category, scales = "free_y", labeller = labeller(Category = facetLabels)) +
  labs(y = "Number of INDELs") +
  scale_fill_manual(values = c("Graph" = "#62109F", "Linear" = "#FE6244"),
                    name = "Source") +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  geom_text(aes(label = Value), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 4) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )




ggsave("plots/variantCounts_all_07.jpg", plot_counts, width = 15, height = 18, units = "cm", dpi = 100)




