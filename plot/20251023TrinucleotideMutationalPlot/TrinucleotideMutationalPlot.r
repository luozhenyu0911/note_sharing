library(forcats)
library(ggplot2)
library(ggpubr)
library(spgs)
library(stringr)

paper_theme <- theme(plot.title = element_text(size = 8),
                     axis.title.x = element_text(size = 8), 
                     axis.text.x  = element_text(size = 8), 
                     axis.title.y = element_text(size = 8),
                     axis.text.y  = element_text(size = 8),
                     axis.ticks.length = unit(0.15, "cm"),
                     panel.background = element_blank(), 
                     plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
                     axis.line = element_line(colour = "black"),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 8), 
                     strip.text.x = element_text(size = 8),
                     strip.text.y = element_text(size = 8))

# Figure 1a - plot the mutational signature estimated by the model, across the reference mtDNA (excluding OriB-OriH)

file <- read.delim(file = './mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# convert to pyrimidine for plotting
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
file$pyr_mut <- ifelse(file$mut == "G>T","C>A",
                     ifelse(file$mut == "G>A", "C>T",
                            ifelse(file$mut == "G>C", "C>G",
                                   ifelse(file$mut == "A>T", "T>A",
                                          ifelse(file$mut == "A>C", "T>G",
                                                 ifelse(file$mut == "A>G", "T>C", file$mut))))))
file$pyr_tri <- ifelse(file$REF == "G" | file$REF == "A", 
                       reverseComplement(file$trinucleotide, case = "upper"), 
                       as.character(file$trinucleotide))
file$strand <- factor(ifelse(file$REF == "G" | file$REF == "A", "Heavy", "Light"), 
                      levels = c("Light", "Heavy"), 
                      labels = c("Reference / Light", "Reverse complement / Heavy"))

# exclude OriB-OriH m.191-16197
for_plot <- unique(file[file$POS > 191 & file$POS < 16197, 
                        c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")]) 
  
ggplot(data = for_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Pastel1", name = "Strand") +
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide", y = "Likelihood") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0)) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25)
