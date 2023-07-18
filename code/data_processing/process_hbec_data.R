create_hbec_bam_paths <- function() { 
    hbec_data_path <- "/data/Vishaka_Varun/Analysis_COMBINED"
    data_anno_path <- file.path("/data/Vishaka_Varun/data_for_running_analyses/")

    samples_list <- list("Control_1"="Sample_1_IP_Untreated1",
                              "Control_2"="Sample_2_IP_Untreated2",
                         "Control_3"="Sample_3_IP_Untreated3",
                         "Treated_1"="Sample_4_IP_Pano_1",
                         "Treated_2"="Sample_5_IP_Pano_2",
                         "Treated_3"="Sample_6_IP_Pano_3",
                         "Control_Input"="Sample_7_INP_Untreated2",
                         "Treated_Input"="Sample_8_INP_Pano2")

    bam_file_paths <- c()
    for (sample_name in names(samples_list)) {
        bam_path <- file.path(hbec_data_path,samples_list[[sample_name]])
        for (file_path in list.files(bam_path,full.names=T)) {
            if (grepl(".sorted.markedup.bam$",file_path)) {
                bam_file_paths <- c(bam_file_paths,file_path)
            }
        }
    }
    names(bam_file_paths) <- names(samples_list)

    return(bam_file_paths)
}

load_hbec_data <- function(genome_info,gene_features_df) {
    bam_paths <- create_hbec_bam_paths()
    ret <- load_data("HBEC",bam_paths,genome_info,gene_features_df,lib="H3K27ac",use_spike_in=T)

    ret$summarized_normalized_counts_df <- ret$summarized %>% mutate(Control_Mean=(Control_1+Control_3)/2, Treated_Mean=(Treated_1+Treated_2+Treated_3)/3) %>% mutate(logFC=(Treated_Mean - Control_Mean - Treated_Input + Control_Input)) 
    return(ret)
}
