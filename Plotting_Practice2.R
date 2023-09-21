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

gencode_file1 <- import("C:\\Users\\Yash\\Desktop\\BTP- Yash Gupa_Prof. Ishaan\\First_Task_BTP\\data_files\\gencode.v44.annotation.gtf")
gencode_df <- as.data.frame(gencode_file1)

# Extract gene_id and transcript_id from the metadata columns of the GTF object
transcript_gene_df <- data.frame(
  transcript_id = mcols(gencode_file1)$transcript_id,
  gene_id = mcols(gencode_file1)$gene_id
)

# Remove rows with NA values (if any)
transcript_gene_df <- transcript_gene_df %>%
  filter(!is.na(transcript_id) & !is.na(gene_id))

# Group by gene_id and count the number of unique transcripts for each gene
count_matrix <- transcript_gene_df %>%
  group_by(gene_id) %>%
  summarise(transcript_count = n_distinct(transcript_id)) %>%
  ungroup()

selected_genes <- count_matrix %>%
  filter(transcript_count > 5 & transcript_count < 20)

random_genes <- selected_genes %>%
  sample_n(3)


subset_transcript_gene_df <- transcript_gene_df %>%
  filter(gene_id %in% random_genes$gene_id)
subset_transcript_gene_df <- subset_transcript_gene_df %>%
  distinct(transcript_id, .keep_all = TRUE)

unique_transcripts_df <- subset_transcript_gene_df %>%
  distinct(transcript_id, .keep_all = TRUE)

set.seed(123)  # Setting a seed for reproducibility

# Extract the transcript counts from random_genes
transcript_counts <- random_genes$transcript_count

# Initialize table_df with random values within the transcript counts
table_df <- data.frame(
  Column1 = sample(1:transcript_counts[1], 5, replace = TRUE),
  Column2 = sample(1:transcript_counts[2], 5, replace = TRUE),
  Column3 = sample(1:transcript_counts[3], 5, replace = TRUE)
)

# Randomly Select Rows from subset_transcript_gene_df
set.seed(123)  
new_dataframes <- list()
for (i in 1:nrow(table_df)) {
  current_df <- data.frame()
  for (j in 1:nrow(random_genes)) {
    subset_df <- subset_transcript_gene_df[subset_transcript_gene_df$gene_id == random_genes$gene_id[j], ]
    selected_transcripts <- subset_df %>%
      sample_n(table_df[i, j])
    current_df <- rbind(current_df, selected_transcripts)
  }
  new_dataframes[[i]] <- current_df
}

# Convert Dataframes to Count Matrices
count_matrices <- list()
for (i in 1:length(new_dataframes)) {
  count_df <- new_dataframes[[i]] %>%
    mutate(`No. of Molecules` = sample(1:20, n(), replace = TRUE))
  count_matrices[[i]] <- count_df
}

# Assign each count matrix to a separate variable in the environment
for (i in 1:length(count_matrices)) {
  assign(paste("count_matrix_", i, sep = ""), count_matrices[[i]])
}

#joining the matrices

joined_count_matrix <- count_matrix_1 %>%
  full_join(count_matrix_2, by = c("transcript_id", "gene_id"), suffix = c("_1", "_2")) %>%
  full_join(count_matrix_3, by = c("transcript_id", "gene_id")) %>%
  full_join(count_matrix_4, by = c("transcript_id", "gene_id"), suffix = c("_3", "_4")) %>%
  full_join(count_matrix_5, by = c("transcript_id", "gene_id"))

#subsetting the matrices

desired_transcript_ids <- c("ENST00000552515.5", "ENST00000397997.6")  # Replace with your desired transcript_ids

# Subset the joined_count_matrix based on the desired transcript_id(s)
subset_joined_matrix <- joined_count_matrix %>%
  filter(transcript_id %in% desired_transcript_ids)

renamed_subset_matrix <- subset_joined_matrix %>%
  rename(
    Sample_1 = `No. of Molecules_1`,
    Sample_2 = `No. of Molecules_2`,
    Sample_3 = `No. of Molecules_3`,
    Sample_4 = `No. of Molecules_4`,
    Sample_5 = `No. of Molecules`
  )

#plotting gene expression


# Reshape the data to long format
long_data <- renamed_subset_matrix %>%
  pivot_longer(cols = starts_with("Sample"), names_to = "Sample", values_to = "Expression")

# Plot the histograms
ggplot(long_data, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~transcript_id, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(title = "Expression of Transcripts by Sample", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#######################################################################################
selected_transcript <- gencode_df %>% filter(transcript_id == "ENST00000552515.5")
count <- 10
replicated_transcript <- lapply(1:count, function(i) {
  df <- selected_transcript
  df$y_position <- i
  return(df)
}) %>% bind_rows()
replicated_transcript <- replicated_transcript %>% filter(type != "CDS")

replicated_transcript <- replicated_transcript %>% filter(type != "CDS")

# Define colors (adjust as needed)
all_colors <- c("exon" = "green", "intron" = "blue")

# Plot
plot <- ggplot(replicated_transcript) +
  geom_range(aes(xstart = start, xend = end, 
                 y = interaction(transcript_id, y_position), 
                 fill = ifelse(type == "exon", "exon", "intron")), height = 1) +  # Adjust height here
  scale_fill_manual(values = all_colors, name = "Type") +
  theme_minimal() +
  labs(y = "Transcript ID") +
  theme(axis.text.y = element_blank())

print(plot)

##############################################################################

n_samples <- 5  # Adjust this as needed

# Create the count matrix
count_matrix <- data.frame(
  sample_name = paste0("Sample_", 1:n_samples),
  count = sample(1:7, n_samples, replace = TRUE)
)

# Print the count matrix
print(count_matrix)

selected_transcript <- gencode_df %>% filter(transcript_id == "ENST00000552515.5")

# Remove CDS from the selected transcript
selected_transcript <- selected_transcript %>% filter(type != "CDS")

# Define colors (adjust as needed)
all_colors <- c("exon" = "green", "intron" = "blue")

# Create a list to store the plots
plots_list <- list()

# Generate plots based on the counts in count_matrix
for (i in 1:nrow(count_matrix)) {
  count <- count_matrix$count[i]
  
  # Replicate the transcript data based on the count
  replicated_transcript <- lapply(1:count, function(j) {
    df <- selected_transcript
    df$y_position <- j
    return(df)
  }) %>% bind_rows()
  
  # Plot
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
  
  # Store the plot in the list
  plots_list[[i]] <- plot
}

# Print the plots
for (i in 1:length(plots_list)) {
  print(plots_list[[i]])
}

grid.arrange(plots_list[[1]], plots_list[[2]], plots_list[[3]], plots_list[[4]], plots_list[[5]], ncol = 1)


###################################

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
  # Import GTF file
  gtf_file <- import(gtf_file_path)
  
  # Convert to data frame
  gtf_df <- as.data.frame(gtf_file)
  
  # Filter based on transcript ID
  selected_transcript <- gtf_df %>% filter(transcript_id == transcript_id & type == "exon")
  
  # Define colors (adjust as needed)
  all_colors <- c("exon" = "green", "intron" = "blue")
  
  # Create a list to store the plots
  plots_list <- list()
  
  # Generate plots based on the counts in count_matrix
  for (i in 1:nrow(count_matrix)) {
    count <- count_matrix$count[i]
    
    # Replicate the transcript data based on the count
    replicated_transcript <- lapply(1:count, function(j) {
      df <- selected_transcript
      df$y_position <- j
      return(df)
    }) %>% bind_rows()
    
    # Plot
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
    
    # Store the plot in the list
    plots_list[[i]] <- plot
  }
  
  # Combine the plots
  combined_plot <- do.call(grid.arrange, c(plots_list, list(ncol = 1)))
  
  return(combined_plot)
}

#count_matrix should have a specific form and the name of the fields should be sample_1,2,3,4,5.