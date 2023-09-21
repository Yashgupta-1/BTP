install.packages("dplyr")
install.packages("readr")
install.packages("rtracklayer")
install.packages("BiocManager")
BiocManager::install("rtracklayer")

library(BiocManager)
library(rtracklayer)
library(dplyr)
library(readr)
library(dplyr)
#The code starts here.


gencode_file1 <- import("C:\\Users\\Yash\\Desktop\\BTP- Yash Gupa_Prof. Ishaan\\First_Task_BTP\\data_files\\gencode.v44.annotation.gtf")
class(gencode_file1)
gencode_file1 #Just viewing the gRange.
colnames(mcols(gencode_file1))

Gen_Chr_df <- data.frame(chromosome = seqnames(gencode_file1), 
                      gene_id = mcols(gencode_file1)$gene_id)

gene_counts <- Gen_Chr_df %>%
  group_by(chromosome) %>%
  summarise(num_genes = n_distinct(gene_id))

Chr13 <- gencode_file1[seqnames(gencode_file1) == "chr13"]
Chr13
gene_ids <- mcols(Chr13)$gene_id
gene_names <- mcols(Chr13)$gene_name
gene_types <- mcols(Chr13)$gene_type
length(unique(gene_ids))==length(unique(gene_names)) #Just Curious
length(unique(gene_ids))==length(unique(gene_types)) #Just Curious
length(unique(gene_ids)) #Just Curious
length(unique(gene_names)) #Just Curious
length(unique(gene_types)) #Just Curious

Chr13_df <- data.frame(gene_identity = mcols(Chr13)$gene_name, transcript_id = mcols(Chr13)$transcript_id)
gene_expression <- Chr13_df %>% group_by(gene_identity) %>% summarisesc(num_transcripts = n_distinct(transcript_id))
good_gene_expression = gene_expression[gene_expression$num_transcripts>mean(gene_expression$num_transcripts), ]

library(ggplot2)
ggplot(good_gene_expression, aes(x = gene_identity, y = num_transcripts)) +
     geom_bar(stat = 'identity') +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2), plot.title = element_text(hjust = 0.5)) +
     labs(title = "Gene Expression in Chr1", 
                   x = "Gene", 
                   y = "Expression Level")

#Plotting the exon-intron boundaries of specific genes in a chromosome and coloring constitutive exons as red in color.


install.packages("devtools")
devtools::install_github("dzhang32/ggtranscript")
library(ggtranscript)
Chr13
Chromosome13_df <- as.data.frame(Chr13)
gene_data <- Chromosome13_df %>% filter(gene_name == "BRCA2")
exons <- gene_data %>% filter(type == "exon")
exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_type)
  ) +
  scale_fill_manual(values = c("protein_coding" = "blue", "nonsense_mediated_decay" = "green", "retained_intron" = "red")) +  
  geom_intron(
    data = to_intron(exons, "transcript_name"),
    aes(strand = strand)
  )
