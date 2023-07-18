source("process_public_data.R")
source("generic_process_functions.R")
library(Rsubread)

mm10_gene_features_df <- create_feature_table(genome="mm10", #Genome assembly
                                             num_bins=10,   # Number of bins for promoter and gene body
                                             tss_bin_size=100,   #Size of TSS bin
                                             min_gene_length=200)   #Discard genes/miRNAs whose length is below this threshold.

#bam_paths should be a named vector. For ex : bam_paths <- c("untreated_1"="/path/to/bam/1","untreated_2"="/path/to/bam/2",treated"="/path/to/bam/2")
results <- find_feature_read_counts(bam_paths,
                                    mm10_gene_features_df, 
                                         cell_line = "HCT116",  #Name of cell line         
                                         lib="H3K27ac",   #Type of library                          
                                dmel_genome_fai_dt=NULL,  #Ignore this option          
                                         use_spike_in=F)

#This computes a TPM-type normalization for the read counts in each bin
normalized_read_counts_df <- normalize_read_counts( results$read_counts, 
                                                   results$mapping_info, 
                                                   sample_names=sample_names,  #Typically, sample_names = names(bam_paths)
                                                   dmel_features_df=NULL )

#Set bin = "all" to sum up TPM across all bins. Set bin = "none" to get TPM values for each bin. 
summarized_normalized_counts_df <- summarize_counts( normalized_read_counts_df, sample_names, 
                                                    bin="all" ) 

