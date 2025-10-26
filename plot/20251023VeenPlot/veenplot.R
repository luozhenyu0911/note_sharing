
# install.packages("ggVennDiagram")
library(ggVennDiagram)

dataf = read.table('../data/seven_dataset.csv', sep = ',', header = 1)

colnames(dataf)
x <- list(`A (5218)`=dataf$A[dataf$A != ""],
          `B (2976)`=dataf$B[dataf$B != ""],
          `C (10850)`=dataf$C[dataf$C != ""],
          `D (14132)`=dataf$D[dataf$D != ""],
          `E (40451)`=dataf$E[dataf$E != ""],
          `F (19321)`=dataf$F[dataf$F != ""],
          `G (7,266)`=dataf$G[dataf$G != ""])

library(ggplot2)
ggVennDiagram(x, label = "count",label_alpha=0,label_geom='label',edge_size=1.5,
              set_size = 6,
              set_color = c("#1f77b4", "#ff7f0e", "#2ca02c", 
                            "#d62728", "#9467bd", "#8c564b", 
                            "#e377c2"))+#scale_fill_distiller(palette = "PiYG")
  scale_fill_gradientn(colours = c("white")) + theme(legend.position = '')






