library(ggVennDiagram)

setwd('D:\\03_mtDNA\\bigcsMT\\analysis\\mtDNA_BIGCS_II_4x100\\output\\outputfile')

dataf = read.table('seven_dataset.csv', sep = ',', header = 1)
tail(dataf$BIGCS_II[dataf$BIGCS_II != ""])
dataf$BIGCS_II[dataf$BIGCS_II != ""]

colnames(dataf)
x <- list(`BIGCS II (5218)`=dataf$BIGCS_II[dataf$BIGCS_II != ""],
          `1KGP (2976)`=dataf$X1000G[dataf$X1000G != ""],
          `GnomAD (10850)`=dataf$Gnomad[dataf$Gnomad != ""],
          `Helix (14132)`=dataf$Helix[dataf$Helix != ""],
          `HmtVar (40451)`=dataf$Hmtvar[dataf$Hmtvar != ""],
          `MitoMap (19321)`=dataf$Mitomap[dataf$Mitomap != ""],
          `NyuWa (7,266)`=dataf$NyuWa[dataf$NyuWa != ""])

library(ggplot2)
ggVennDiagram(x, label = "count",label_alpha=0,label_geom='label',edge_size=1.5,
              set_size = 6,
              set_color = c("#1f77b4", "#ff7f0e", "#2ca02c", 
                            "#d62728", "#9467bd", "#8c564b", 
                            "#e377c2"))+#scale_fill_distiller(palette = "PiYG")
  scale_fill_gradientn(colours = c("white")) + theme(legend.position = '')

?ggVennDiagram






