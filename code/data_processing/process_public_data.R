library(GEOquery)
source("generic_process_functions.R")

public_path <- "/data/Vishaka_Varun/public_data"

create_rh4_bam_paths <- function() {
    rh4_k27ac_path <- file.path(public_path,"output","RH4","aligned")
    rh4_rna_seq_path <- file.path(public_path,"output","RH4","salmon_out")
    sample_sheet_dt <- fread(file.path(public_path,"RH4_sample_sheet.tsv")) %>% dplyr::filter(time_point=="6h") %>%
    mutate(sample_name=paste(cell_line,treatment,sep="_")) %>% data.table

    rh4_k27ac_paths <- c()
    rh4_rna_seq_paths <- list()
    rh4_k27ac_sample_names <- c()
    rh4_rna_seq_sample_names <- c()
    print(sample_sheet_dt)
    for (idx in 1:nrow(sample_sheet_dt)) {
        lib <- sample_sheet_dt[idx,library]
        sample_name <- sample_sheet_dt[idx,sample_name]
        srr_id <- sample_sheet_dt[idx,srr]
        if (lib == "H3K27ac") {
            rh4_k27ac_sample_names <- c(rh4_k27ac_sample_names,sample_name)
            rh4_k27ac_paths <- c(rh4_k27ac_paths,file.path(rh4_k27ac_path,paste(srr_id,"dupmarked.sorted.bam",sep=".")))
        }  else {
            rh4_rna_seq_sample_names <- c(rh4_rna_seq_sample_names, sample_name)
            rh4_rna_seq_paths <- c(rh4_rna_seq_paths, file.path(rh4_rna_seq_path,srr_id,"quant.sf"))
        }
    }
    names(rh4_k27ac_paths) <- rh4_k27ac_sample_names
    names(rh4_rna_seq_paths) <- rh4_rna_seq_sample_names

    return(list("k27ac"=rh4_k27ac_paths,"rna-seq"=rh4_rna_seq_paths))
}

load_rh4_data <- function(genome_info,gene_features_df,bin) {
    rh4_bam_paths <- create_rh4_bam_paths()
    ret <- load_data("RH4",rh4_bam_paths[["k27ac"]],genome_info,gene_features_df,"H3K27ac",use_spike_in=T,bin=bin)

    ret$summarized_normalized_counts_df <- ret$summarized %>% mutate(logFC=RH4_treated - RH4_untreated)

    ret$`RNA-seq` <- run_tximport( unlist(rh4_bam_paths[["rna-seq"]]) )
    
    return(ret)
}

create_hct116_paths <- function() {
    chip_seq_path <- file.path(public_path,"output","HCT116","aligned")
    rna_seq_path <- file.path(public_path,"output","HCT116","salmon_out")
    sample_sheet_dt <- fread(file.path(public_path,"backup","HCT116_sample_sheet.tsv")) %>% mutate(sample_name=paste(cell_line,treatment,sep="_")) %>% data.table

    file_path_list <- list()
    file_info_list <- list()

    for (idx in 1:nrow(sample_sheet_dt)) {
        lib <- sample_sheet_dt[idx,library]
        sample_name <- sample_sheet_dt[idx,sample_name]
        srr_id <- sample_sheet_dt[idx,srr]
        if (lib == "RNA-seq") {
            file_path <- file.path(rna_seq_path,srr_id,"quant.sf")
        } else {
            file_path <- file.path(chip_seq_path,paste(srr_id,"dupmarked.sorted.bam",sep="."))
        }
        file_path_list[[lib]] <- c(file_path_list[[lib]],file_path)
        file_info_list[[lib]] <- c(file_info_list[[lib]],sample_name)
    }
        
    for (lib in unique(sample_sheet_dt$library)) {
        names(file_path_list[[lib]]) <- file_info_list[[lib]]
    }

    return(file_path_list)
}

load_hct116_data <- function(genome_info,gene_features_df,bin) {
    hct116_paths <- create_hct116_paths()
    to_return <- list()

    dosages <- c("4.68nM","9.37nM","18.75nM","30nM","37.5nM","75nM","150nM","300nM")
    for (lib in names(hct116_paths)) {
        if (lib == "RNA-seq")
            next

        print(lib)
        info <- load_data("HCT116",hct116_paths[[lib]],genome_info,gene_features_df,lib,use_spike_in=F,bin=bin)
        flush.console()
        to_return[[lib]] <- info
        for (dosage in dosages) {
            dosage_col <- paste("HCT116_treated",dosage,sep="_")
            if (!dosage_col %in% colnames(info$summarized))
                next
            logFC_col <- paste("logFC",dosage,sep="_")
            to_return[[lib]]$summarized[[logFC_col]] <- info$summarized[[dosage_col]] - info$summarized$HCT116_untreated
        }
    }
    
    to_return$`RNA-seq` <- run_tximport( unlist(hct116_paths[["RNA-seq"]]) )
    return(to_return)
}

create_l3.6_bam_paths <- function() {
    k27ac_path <- file.path(public_path,"output","L3.6","aligned")
    rna_seq_path <- file.path(public_path,"output","L3.6","salmon_out")
    sample_sheet_dt <- fread(file.path(public_path,"L3.6_sample_sheet.tsv")) 

    k27ac_paths <- c()
    rna_seq_paths <- list()
    k27ac_sample_names <- c()
    for (idx in 1:nrow(sample_sheet_dt)) {
        lib <- sample_sheet_dt[idx,library]
        sample_name <- paste(c("L3.6",sample_sheet_dt[idx,.(treatment,replicate)]),collapse="_")
        srr_id <- sample_sheet_dt[idx,srr]
        if (lib == "H3K27ac") {
            k27ac_sample_names <- c(k27ac_sample_names,sample_name)
            k27ac_paths <- c(k27ac_paths,file.path(k27ac_path,paste(srr_id,"dupmarked.sorted.bam",sep=".")))
        } 
    }
    names(k27ac_paths) <- k27ac_sample_names
    return(k27ac_paths)
}
