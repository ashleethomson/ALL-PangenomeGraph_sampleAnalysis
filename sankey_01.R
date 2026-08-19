# Library
library(dplyr)
library(grid)
library(tidyverse)
library(PantaRhei)




nodes <- tribble(
  ~ID,     ~label,                 ~x,    ~y,      ~label_pos,
  "total", "Total",                "0",   "0",     "left",
  "fo",    "In Final\nOutput",      "1.5", "total", "above",
  "nfo",   "Not in\nFinal Output", "1.5", "-1",    "below",
  "rg",    "In raw\nGATK",          "3.5", "nfo",   "above",
  "nrg",   "Not in\nraw GATK",     "3.5", "-2.25", "below",
  "p",     "PASS",                 "5.5", "-0.3",  "right",
  "qd",    "QD",                   "5.5", "rg",    "right",
  "fs",    "FS",                   "5.5", "-1.5",  "right",
  "b",     "In BAM",               "5.5", "nrg",   "right",
  "nb",    "Not in\nBAM",           "5.5", "-3",    "right"
)



flows_G <- tribble(
  ~from, ~to, ~ substance, ~quantity,
  "total", "fo", "Graph", 96,
  "total", "nfo", "", 5,
  "nfo", "rg", "", 4,
  "rg", "p", "", 3,
  "rg", "qd", "", 1,
  "nfo", "nrg", "", 1,
  "nrg", "b", "", 1
)



flows_L <- tribble(
  ~from, ~to, ~ substance, ~quantity,
  "total", "fo", "Linear", 63,
  "total", "nfo", "", 38,
  "nfo", "rg", "", 11,
  "rg", "p", "", 7,
  "rg", "qd", "", 2,
  "rg", "fs", "", 2,
  "nfo", "nrg", "", 27,
  "nrg", "b", "", 4,
  "nrg", "nb", "", 23
)

## aesthetics
colours_G <- tribble(
  ~substance, ~color,
  "Graph",    "#62109F"
)

colours_L <- tribble(
  ~substance, ~color,
  "Linear",    "#FE6244"
)


ns <- list(type = "arrow",              # node type
           length = 0.3,
           gp = gpar(fill = "#FFDEB9",   # node colour
                     col = "white",        # node border colour
                     lwd = 4),             # border width around node
           mag_pos = "inside",
           mag_fmt = "%1.0f")



sankey(nodes, flows_L, colours_L, node_style = ns,
       legend = TRUE, max_width = 0.1, rmin = 0.3,
       page_margin=c(0.1, 0.3, 0.2, 0.1))

sankey(nodes, flows_G, colours_G, node_style = ns,
       legend = TRUE, max_width = 0.1, rmin = 0.3,
       page_margin=c(0.1, 0.3, 0.2, 0.1))


pdf("linearSankey_01.pdf",
    height = 9, width = 15) # Set up PDF device
sankey(nodes, flows_L, colours_L, node_style = ns, legend = TRUE)           # plot diagram
dev.off() 


pdf("rgraphSankey_01.pdf",
    height = 4, width = 15) # Set up PDF device
sankey(nodes, flows_G, colours_G, node_style = ns, legend = TRUE)          # plot diagram
dev.off() 




