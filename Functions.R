#Plotting transcripts from two different sources.

plot_gene_from_gtf <- function(gtf_path_1, gtf_path_2, gene_name, species_1 = "Species 1", species_2 = "Species 2") {
  
  library(rtracklayer)
  library(dplyr)
  library(ggtranscript)
  library(ggplot2)
  library(gridExtra)
  
  
  genome_1 <- import(gtf_path_1)
  genome_2 <- import(gtf_path_2)
  
  
  df_1 <- as.data.frame(genome_1)
  df_2 <- as.data.frame(genome_2)
  
  
  gene_id_1 <- unique(genome_1[genome_1$gene_name == gene_name,]$gene_id)
  gene_id_2 <- unique(genome_2[genome_2$gene_name == gene_name,]$gene_id)
  
  
  gene_1 <- subset(df_1, gene_id == gene_id_1 & type == 'exon')
  gene_2 <- subset(df_2, gene_id == gene_id_2 & type == 'exon')
  
  
  plot_1 <- gene_1 %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = paste(species_1, 'Exon')), color = 'darkblue', height = 0.3) +
    geom_intron(data = to_intron(gene_1, 'transcript_id'), aes(strand = strand), color = 'darkblue', arrow = NULL) +
    labs(title = paste(species_1, gene_name), y = NULL, x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
    scale_fill_manual(values = c(paste(species_1, 'Exon') = 'darkblue'))
  
  
  plot_2 <- gene_2 %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = paste(species_2, 'Exon')), color = 'red', height = 0.3) +
    geom_intron(data = to_intron(gene_2, 'transcript_id'), aes(strand = strand), color = 'red', arrow = NULL) +
    labs(title = paste(species_2, gene_name), y = NULL, x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none', plot.title = element_text(hjust = 1), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm')) +
    scale_fill_manual(values = c(paste(species_2, 'Exon') = 'red'))
  
  
  combined_plot <- grid.arrange(plot_1, plot_2, ncol = 1)
  
  return(combined_plot)
}


#######Plotting transcripts along with their counts

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

plot_transcript_samples <- function(gtf_file_path, count_matrix, transcript_id) {
  
  gtf_file <- import(gtf_file_path)
  
  
  gtf_df <- as.data.frame(gtf_file)
  
  selected_transcript <- gtf_df %>% filter(transcript_id == transcript_id & type == "exon")
  
  all_colors <- c("exon" = "green", "intron" = "blue")
  
 
  plots_list <- list()
  

  for (i in 1:nrow(count_matrix)) {
    count <- count_matrix$count[i]
    
    
    replicated_transcript <- lapply(1:count, function(j) {
      df <- selected_transcript
      df$y_position <- j
      return(df)
    }) %>% bind_rows()
    
    
    plot <- ggplot(replicated_transcript) +
      geom_range(aes(xstart = start, xend = end, 
                     y = interaction(transcript_id, y_position), 
                     fill = ifelse(type == "exon", "exon", "intron")), height = 1) +
      geom_segment(data = subset(replicated_transcript, type == "intron"),
                   aes(x = start, xend = end, y = interaction(transcript_id, y_position), 
                       yend = interaction(transcript_id, y_position)), color = "black") +
      scale_fill_manual(values = all_colors, name = "Type") +
      theme_minimal() +
      labs(title = paste("Sample:", count_matrix$sample_name[i])) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 1, vjust = -0.5)) +
      guides(fill = FALSE)
    
    
    plots_list[[i]] <- plot
  }
  
  
  combined_plot <- do.call(grid.arrange, c(plots_list, list(ncol = 1)))
  
  return(combined_plot)
}