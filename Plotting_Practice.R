library(BiocManager)
library(rtracklayer)
library(dplyr)
library(readr)
library(tidyr)
library(ggtranscript)
library(ggplot2)
library(scales)

gencode_file1 <- import("C:\\Users\\Yash\\Desktop\\BTP- Yash Gupa_Prof. Ishaan\\First_Task_BTP\\data_files\\gencode.v44.annotation.gtf")
gencode_df <- as.data.frame(gencode_file1)

gene_counts <- gencode_df %>%
  filter(type == "gene") %>%
  group_by(seqnames) %>%
  summarise(num_genes = n_distinct(gene_id)) %>%
  arrange(-num_genes)

transcript_counts <- gencode_df %>%
  filter(type == "transcript") %>%
  group_by(seqnames) %>%
  summarise(num_transcripts = n_distinct(transcript_id)) %>%
  arrange(-num_transcripts)

Chr1_data <- subset(gencode_file1, seqnames == "chr1")

transcript_data <- subset(Chr1_data, type == "transcript")
transcript_gene_df <- data.frame(
  transcript_id = mcols(transcript_data)$transcript_id,
  gene_id = mcols(transcript_data)$gene_id
)
transcript_gene_df <- unique(transcript_gene_df)
transcript_gene_df <- transcript_gene_df %>%
  distinct(transcript_id, .keep_all = TRUE)

#calculating half of the total entries in the dataframe.
half_n <- nrow(transcript_gene_df) %/% 2 

#random selection of half of the entries
transcript_df1 <- transcript_gene_df %>%
  sample_n(half_n)

#creating count matrix
count_matrix <- transcript_df1 %>%
  group_by(gene_id) %>%
  summarise(
    transcripts = paste(transcript_id, collapse = ","),
    num_transcripts = n()
  ) %>%
  ungroup()
count_matrix1 <- as.data.frame(count_matrix)

#matrix 2
transcript_df2 <- transcript_gene_df %>%
  sample_n(half_n)
count_matrix2 <- transcript_df2 %>%
  group_by(gene_id) %>%
  summarise(
    transcripts = paste(transcript_id, collapse = ","),
    num_transcripts = n()
  ) %>%
  ungroup()
count_matrix2 <- as.data.frame(count_matrix2)

#matrix3
transcript_df3 <- transcript_gene_df %>%
  sample_n(half_n)
count_matrix3 <- transcript_df3 %>%
  group_by(gene_id) %>%
  summarise(
    transcripts = paste(transcript_id, collapse = ","),
    num_transcripts = n()
  ) %>%
  ungroup()
count_matrix3 <- as.data.frame(count_matrix3)

#matrix4
transcript_df4 <- transcript_gene_df %>%
  sample_n(half_n)
count_matrix4 <- transcript_df4 %>%
  group_by(gene_id) %>%
  summarise(
    transcripts = paste(transcript_id, collapse = ","),
    num_transcripts = n()
  ) %>%
  ungroup()
count_matrix4 <- as.data.frame(count_matrix4)

#matrix5
transcript_df5 <- transcript_gene_df %>%
  sample_n(half_n)
count_matrix5 <- transcript_df5 %>%
  group_by(gene_id) %>%
  summarise(
    transcripts = paste(transcript_id, collapse = ","),
    num_transcripts = n()
  ) %>%
  ungroup()
count_matrix5 <- as.data.frame(count_matrix5)

#creating a bigger matrix

names(count_matrix1)[names(count_matrix1) == "transcripts"] <- "transcripts_matrix1"
names(count_matrix2)[names(count_matrix2) == "transcripts"] <- "transcripts_matrix2"
names(count_matrix3)[names(count_matrix3) == "transcripts"] <- "transcripts_matrix3"
names(count_matrix4)[names(count_matrix4) == "transcripts"] <- "transcripts_matrix4"
names(count_matrix5)[names(count_matrix5) == "transcripts"] <- "transcripts_matrix5"

final_matrix <- count_matrix1 %>%
  select(gene_id, transcripts_matrix1) %>%
  inner_join(select(count_matrix2, gene_id, transcripts_matrix2), by = "gene_id") %>%
  inner_join(select(count_matrix3, gene_id, transcripts_matrix3), by = "gene_id") %>%
  inner_join(select(count_matrix4, gene_id, transcripts_matrix4), by = "gene_id") %>%
  inner_join(select(count_matrix5, gene_id, transcripts_matrix5), by = "gene_id")

#Plotting

gene_of_interest <- "ENSG00000007341.19"

transcripts_of_interest <- final_matrix %>%
  filter(gene_id == gene_of_interest) %>%
  select(-gene_id) %>%
  gather(matrix_id, transcript_id, everything()) %>%
  separate_rows(transcript_id, sep = ",")

transcript_structures <- Chr1_data[Chr1_data$transcript_id %in% transcripts_of_interest$transcript_id, ]
transcript_structures_df <- as.data.frame(transcript_structures)

transcript_structures_df <- transcript_structures_df %>%
  left_join(transcripts_of_interest, by = "transcript_id")

transcript_structures_df$matrix_id <- gsub("transcripts_", "", transcript_structures_df$matrix_id)


all_colors <- c(
  matrix1 = "yellow", matrix1_intron = "lightyellow",
  matrix2 = "red", matrix2_intron = "pink",
  matrix3 = "green", matrix3_intron = "lightgreen",
  matrix4 = "blue", matrix4_intron = "lightblue",
  matrix5 = "purple", matrix5_intron = "lavender"
)

plot <- ggplot(transcript_structures_df) +
  geom_range(aes(xstart = start, xend = end, y = interaction(transcript_id, matrix_id), 
                 fill = ifelse(type == "exon", matrix_id, paste0(matrix_id, "_intron")))) +
  scale_fill_manual(values = all_colors, name = "Matrix ID") +
  theme_minimal() +
  labs(y = "Transcript ID") +
  theme(axis.text.y = element_blank())  

print(plot)

#to highlight differences

reference_transcript_id <- "ENST00000360743.8"

reference_exons <- transcript_structures_df %>%
  filter(transcript_id == reference_transcript_id & type == "exon")

other_transcripts_exons <- transcript_structures_df %>%
  filter(transcript_id != reference_transcript_id & type == "exon")

transcript_diffs <- to_diff(
  exons = other_transcripts_exons,
  ref_exons = reference_exons,
  group_var = "transcript_id"
)

plot <- ggplot() +
  geom_range(
    data = transcript_structures_df %>% filter(type == "exon"),
    aes(xstart = start, xend = end, y = transcript_id)
  ) +
  geom_range(
    data = transcript_diffs,
    aes(xstart = start, xend = end, y = transcript_id, fill = diff_type),
    alpha = 0.2
  ) +
  theme_minimal() +
  labs(y = "Transcript ID") +
  theme(axis.text.y = element_blank())

print(plot)

#to create a box around certain exons

common_exons_50 <- transcript_structures_df %>%
  filter(type == "exon") %>%
  group_by(start, end) %>%
  summarise(count = n_distinct(transcript_id)) %>%
  filter(count / length(unique(transcript_structures_df$transcript_id)) >= 0.50) %>%
  ungroup()

print(common_exons_50)

buffer <- 500 #to increase the width of rectangle

plot <- ggplot(transcript_structures_df) +
  geom_range(aes(xstart = start, xend = end, y = interaction(transcript_id, matrix_id), 
                 fill = ifelse(type == "exon", matrix_id, paste0(matrix_id, "_intron")))) +
  geom_rect(data = common_exons_50, 
            aes(xmin = start - buffer, xmax = end + buffer, ymin = -Inf, ymax = Inf), 
            fill = NA, color = "black", linetype = "solid") +
  scale_fill_manual(values = all_colors, name = "Matrix ID") +
  theme_minimal() +
  labs(y = "Transcript ID") +
  theme(axis.text.y = element_blank())

print(plot)