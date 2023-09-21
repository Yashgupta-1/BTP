library(BiocManager)
library(rtracklayer)
library(dplyr)
library(readr)
library(tidyr)
library(ggtranscript)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggpubr)

human_genome <- import("C:\\Users\\Yash\\Desktop\\BTP- Yash Gupa_Prof. Ishaan\\First_Task_BTP\\data_files\\gencode.v44.annotation.gtf")

mouse_genome <- import("C:\\Users\\Yash\\Desktop\\BTP- Yash Gupa_Prof. Ishaan\\First_Task_BTP\\data_files\\gencode.vM33.annotation.gtf")

human_df <- as.data.frame(human_genome)
mouse_df <- as.data.frame(mouse_genome)

human_genes <- unique(human_genome$gene_name)
mouse_genes <- unique(mouse_genome$gene_name)

common_genes <- intersect(human_genes, mouse_genes)

human_C9orf72 <- unique(human_genome[human_genome$gene_name == "C9orf72",]$gene_id)
mouse_C9orf72 <- unique(mouse_genome[mouse_genome$gene_name == "C9orf72",]$gene_id)
cat("Human C9orf72 gene_id:", human_C9orf72, "\n")
cat("Mouse C9orf72 gene_id:", mouse_C9orf72, "\n")



human_gene <- subset(human_df, gene_id == 'ENSG00000147894.17' & type == 'exon')
mouse_gene <- subset(mouse_df, gene_id == 'ENSMUSG00000028300.15' & type == 'exon')



human_plot <- human_gene %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
  geom_range(aes(fill = 'Human Exon'), color = 'darkblue', height = 0.3) +  # Reduced height for exons
  geom_intron(data = to_intron(human_gene, 'transcript_id'), aes(strand = strand), color = 'darkblue', arrow = NULL) +
  labs(title = 'Human C9orf72', y = NULL, x = NULL) +
  theme_minimal() +
  theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
  scale_fill_manual(values = c('Human Exon' = 'darkblue'))

# Plot mouse gene
mouse_plot <- mouse_gene %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
  geom_range(aes(fill = 'Mouse Exon'), color = 'red', height = 0.3) +  # Reduced height for exons
  geom_intron(data = to_intron(mouse_gene, 'transcript_id'), aes(strand = strand), color = 'red', arrow = NULL) +
  labs(title = 'Mouse C9orf72', y = NULL, x = NULL) +
  theme_minimal() +
  theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
  scale_fill_manual(values = c('Mouse Exon' = 'red'))




# Combine the plots
combined_plot <- grid.arrange(human_plot, human_plot1, mouse_plot, mouse_plot1, ncol = 1)

print(combined_plot)


###################################################################3

plot_gene_from_gtf <- function(gtf_path_1, gtf_path_2, gene_name, species_1 = "Species 1", species_2 = "Species 2") {
  # Load required libraries
  library(rtracklayer)
  library(dplyr)
  library(ggtranscript)
  library(ggplot2)
  library(gridExtra)
  
  # Import GTF files
  genome_1 <- import(gtf_path_1)
  genome_2 <- import(gtf_path_2)
  
  # Convert to data frames
  df_1 <- as.data.frame(genome_1)
  df_2 <- as.data.frame(genome_2)
  
  # Get gene_id for the specified gene_name
  gene_id_1 <- unique(genome_1[genome_1$gene_name == gene_name,]$gene_id)
  gene_id_2 <- unique(genome_2[genome_2$gene_name == gene_name,]$gene_id)
  
  # Filter for exons of the specified genes
  gene_1 <- subset(df_1, gene_id == gene_id_1 & type == 'exon')
  gene_2 <- subset(df_2, gene_id == gene_id_2 & type == 'exon')
  
  # Plot gene from genome 1
  plot_1 <- gene_1 %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = paste(species_1, 'Exon')), color = 'darkblue', height = 0.3) +
    geom_intron(data = to_intron(gene_1, 'transcript_id'), aes(strand = strand), color = 'darkblue', arrow = NULL) +
    labs(title = paste(species_1, gene_name), y = NULL, x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
    scale_fill_manual(values = c(paste(species_1, 'Exon') = 'darkblue'))
  
  # Plot gene from genome 2
  plot_2 <- gene_2 %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = paste(species_2, 'Exon')), color = 'red', height = 0.3) +
    geom_intron(data = to_intron(gene_2, 'transcript_id'), aes(strand = strand), color = 'red', arrow = NULL) +
    labs(title = paste(species_2, gene_name), y = NULL, x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
    scale_fill_manual(values = c(paste(species_2, 'Exon') = 'red'))
  
  # Combine the plots
  combined_plot <- grid.arrange(plot_1, plot_2, ncol = 1)
  
  return(combined_plot)
}