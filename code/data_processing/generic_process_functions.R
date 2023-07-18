library(Rsubread)

read_genome_file <- function(genome="hg38") {
    if (genome == "hg38") {
        print("Loading hg38 genome")
        genome_fai_dt <- fread(file.path("../data_for_running_analyses/hg38_dm6.fa.fai"))
        dm_genome_fai_dt <- genome_fai_dt[grepl("dm",V1),.(GeneID=V1,Chr=V1,Start=1,End=V2,Strand="*")]
        return(list("human"=genome_fai_dt,"dmel"=dm_genome_fai_dt))
    } else if (genome == "mm10") {
        print("Loading mm10 genome")
        flush.console()
        genome_fai_dt <- fread("/fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa.fai") 
        return(list("mouse"=genome_fai_dt))
    }

}

run_tximport <- function(file_paths) {
    human_gene_dt <- read_human_genome_annotation()
    tx_out <- tximport( file_paths,type="salmon",
                 tx2gene=human_gene_dt %>% dplyr::select(transcript_stable_id,gene_stable_id), ignoreTxVersion=T, importer=fread )

    return(tx_out)
}

read_genome_annotation <- function(genome="hg38") {
    if (genome == "hg38") {
        gene_transcript_info_dt <- fread("../data_for_running_analyses/ensembl_101_hg38_human_genes.tsv.gz")
        gene_info_dt <- gene_transcript_info_dt[chr != "MT" & !grepl("CTG",chr) & !grepl("PATCH",chr) & end - start > 100,] %>% mutate(chr=paste0("chr",chr),
    strand=fifelse(strand == -1,"-","+")) %>% unique 
    } else if (genome == "mm10") {
        gene_transcript_info_dt <- fread("../data_for_running_analyses/ensembl_102_GRCm38_mouse_genes.tsv.gz")
        gene_info_dt <- gene_transcript_info_dt[chr != "MT" & !grepl("CTG",chr) & !grepl("PATCH",chr) & !grepl("CHORI",chr) & !grepl("GL",chr) & !grepl("JH",chr) & end - start > 100,] %>% mutate(chr=paste0("chr",chr), strand=fifelse(strand == -1,"-","+")) %>% unique 
    }
    
    #gene_info_dt <- gene_transcript_info_dt[chr != "MT" & !grepl("CTG",chr) & !grepl("PATCH",chr) & end - start > 100, .(gene_name,chr=paste0("chr",chr),start,end,strand=fifelse(strand == -1,"-","+"))] %>% unique 

    return(gene_info_dt)
}

create_feature_table <- function(num_bins=10,tss_bin_size=100,min_gene_length=200,genome="hg38") {
    feature_table_file_name <- file.path(paste(genome,num_bins,"bins_features_table.tsv",sep="_"))
    if (!file.exists(feature_table_file_name)) {
        genome_info <- read_genome_file(genome)
        gene_info_dt <- read_genome_annotation(genome)
        if (genome == "hg38") {
            genome_fai_dt <- genome_info$human
        } else if (genome == "mm10") {
            genome_fai_dt <- genome_info$mouse
        }

        gene_info_dt <- gene_info_dt %>% dplyr::filter(end-start >= min_gene_length)
        gene_granges_info <- makeGRangesFromDataFrame(gene_info_dt,keep.extra.columns=T)
        tss_granges_info <- flank(gene_granges_info,width = tss_bin_size,both=T)
        tss_bins_df <- as.data.frame(tss_granges_info)

        gene_granges_bin_info <- tile(narrow(gene_granges_info,start=tss_bin_size),n=num_bins) %>% unlist
        gene_bins_df <- as.data.frame(gene_granges_bin_info)
        promoter_granges_info <- promoters(shift(gene_granges_info,shift=-tss_bin_size),downstream=0,upstream=2000)
        promoter_granges_bin_info <- promoter_granges_info %>% tile(.,n=num_bins) %>% unlist
        promoter_bins_df <- as.data.frame(promoter_granges_bin_info)
        gene_names <- gene_info_dt$gene_name

        bin_names <- apply(gene_info_dt,1,function(row){return(paste(row["gene_name"], 
                                                       row["chr"],
                                                       row["start"],
                                                       row["end"],
                                                       1:num_bins,sep="|") %>% gsub("\ *","",.));})                                   
        gene_bin_names <- paste("Gene",bin_names,sep="_")
        promoter_bin_names <- paste("Prom",bin_names,sep="_")
        tss_bin_names <- paste("TSS",bin_names[seq(1,length(bin_names),by=num_bins)],sep="_")
        gene_bins_df$bin <- gene_bin_names
        promoter_bins_df$bin <- promoter_bin_names
        tss_bins_df$bin <- tss_bin_names

        gene_features_df <- rbindlist(list(gene_bins_df,promoter_bins_df,
                                           tss_bins_df %>% dplyr::select(all_of(colnames(gene_bins_df))))) %>%
        dplyr::select(GeneID=bin,Chr=seqnames,Start=start,End=end,Strand=strand)

        fwrite(gene_features_df,feature_table_file_name,sep="\t")
    } else {
        gene_features_df <- fread(feature_table_file_name)
    }

    return(gene_features_df)
}

find_feature_read_counts <- function(bam_file_paths,gene_features_df,cell_line,lib="H3K27ac",dmel_genome_fai_dt=NULL,n_threads=6,use_spike_in=T){
    sample_names <- names(bam_file_paths)
    to_return <- list()
    
    features_file_name <- paste(cell_line,lib,"feature_table_10_bins.rds",sep="_")
    print(features_file_name)
    flush.console()

    if (file.exists(features_file_name)) {
        gene_feature_table <- readRDS(features_file_name)
    } else {
        print("Counting reads in features.")
        flush.console()
        gene_feature_table <- featureCounts(files=bam_file_paths,annot.ext=gene_features_df,ignoreDup=T, nthreads=n_threads,allowMultiOverlap=T)

        colnames(gene_feature_table$counts) <- sample_names
        saveRDS(gene_feature_table,features_file_name)
    }

    gene_read_counts_df <- as.data.frame(gene_feature_table$counts) %>% tibble::rownames_to_column("GeneID") %>% mutate_at(all_of(sample_names),list(as.numeric))
    to_return[["read_counts"]] <- gene_read_counts_df

    mapping_info_file_name <- paste(cell_line,lib,"mapping_information.tsv",sep="_")
    if (file.exists(mapping_info_file_name)) {
        mapping_info_dt <- fread(mapping_info_file_name)
    } else {
        print("Counting total number of mapped reads to the human genome.")
        flush.console()
        mapping_info_dt <- propmapped(bam_file_paths)
        rownames(mapping_info_dt) <- sample_names
        fwrite(mapping_info_dt, mapping_info_file_name, sep="\t")
    }
    to_return[["mapping_info"]] <- mapping_info_dt

    if (use_spike_in==T) {
        spike_in_features_file_name <- paste(cell_line,"spikein_table.rds",sep="_")
        if (file.exists(spike_in_features_file_name)) {
            dmel_feature_table <- readRDS(spike_in_features_file_name)
            colnames(dmel_feature_table$counts) <- sample_names
        } else {
            print("Counting number of mapped reads to the fly genome.")
            flush.console()
            dmel_feature_table <- featureCounts(files=bam_file_paths,annot.ext=dmel_genome_fai_dt,ignoreDup=T, nthreads=n_threads)
            colnames(dmel_feature_table$counts) <- sample_names
            saveRDS(dmel_feature_table,spike_in_features_file_name)
        }
        to_return[["dmel_read_counts"]] <- dmel_feature_table$counts

        spike_in_info_file_name <- paste(cell_line,"dmel_spikein_info.tsv",sep="_")
        if (file.exists(spike_in_info_file_name)) {
            spike_in_info_dt <- fread(spike_in_info_file_name)
        } else {
            spike_in_info_dt <- tibble::enframe(colSums(dmel_feature_table$counts),name="treatment",value="dmel_reads") %>% mutate(bam_file_name=dmel_feature_table$targets)
            fwrite(spike_in_info_dt, sep="\t" )
        }
        to_return[["dmel_info"]] <- spike_in_info_dt
    }

    return(to_return)
}

summarize_counts <- function(df_,sample_names,bin="all") {
    split_gene_id_mat <- str_split(df_$GeneID,"\\|") %>% unlist %>% matrix(.,nrow=5) %>% t
    colnames(split_gene_id_mat) <- c("gene_name","chr","start","end","bin")
    
    gene_id_df_ <- as.data.frame(split_gene_id_mat) %>% mutate(element_type=gsub("_.*","",gene_name) %>% gsub("_","",.), gene_name=gsub("^..*?_","",gene_name))
    rm(split_gene_id_mat)
    
    if (bin == "all") {
         summarized_df_ <- cbind(df_,gene_id_df_) %>% 
 group_by(gene_name,element_type,chr,start,end) %>% summarize_at(sample_names,sum) %>% ungroup
    } else if (bin %in% seq(1,10)) {
        bin_ = bin
        summarized_df_ <- cbind(df_,gene_id_df_) %>% dplyr::filter(bin == bin_)
    } else if (bin == "none") {
         summarized_df_ <- cbind(df_,gene_id_df_)       
    }
    
    return(summarized_df_)
}

normalize_read_counts <- function( gene_features_df, mapping_info_df, 
                                   sample_names, dmel_features_df=NULL, pseudocount=1 ) {
    normalized_gene_features_df <- copy(gene_features_df)
    mapping_info_df <- mapping_info_df %>% mutate(row.name=sample_names)
    chip_scale_factor <- 1e6
    if (!is.null(dmel_features_df)) {
        dmel_counts <- apply(dmel_features_df,2,sum)
        dmel_scale_factor <- 10 ** ceiling(log10(min(dmel_counts))) 
        scale_factor <- chip_scale_factor * dmel_scale_factor
    } else {
        scale_factor <- chip_scale_factor
    }

    for (sample_name in sample_names) {
        total_read_count <- mapping_info_df %>% dplyr::filter(row.name==sample_name) %>% pull(NumMapped)
        if (!is.null(dmel_features_df)) {
            dmel_read_count <- dmel_counts[sample_name]
            normalized_counts <- ((scale_factor*gene_features_df[[sample_name]]/total_read_count)/dmel_read_count)
        } else {
            normalized_counts <- (scale_factor*gene_features_df[[sample_name]]/total_read_count)
        }
        normalized_gene_features_df[[sample_name]] <- normalized_counts
    }
    
    return(normalized_gene_features_df)
}

load_data <- function(cell_line,bam_paths,genome_info,gene_features_df,lib,use_spike_in=F,bin="all") {
    #bam_paths <- create_rh4_bam_paths()
    sample_names <- names(bam_paths)#names(rh4_bam_paths)
    flush.console()

    results <- find_feature_read_counts(bam_paths,gene_features_df,
                                             cell_line = cell_line,
                                             lib,
                                    dmel_genome_fai_dt=dmel_fai_dt,
                                             use_spike_in=use_spike_in)

    if (use_spike_in==F) {
        dmel_fai_dt <- NULL
        dmel_features_df <- NULL
    } else {
        dmel_fai_dt=genome_info$dmel
        dmel_features_df <- results$dmel_read_counts 
    }
    normalized_read_counts_df <- normalize_read_counts( results$read_counts, results$mapping_info, sample_names=sample_names, dmel_features_df=dmel_features_df )

    summarized_normalized_counts_df <- summarize_counts( normalized_read_counts_df, sample_names, bin=bin )
    return(list("normalized"=normalized_read_counts_df,"summarized"=summarized_normalized_counts_df, "read_counts"=results$read_counts))
}

create_HCT116_gene_categories <- function(counts_mat){
    log_tpm_df <- as.data.frame(rlog(floor(counts_mat)))
    log_tpm_df$gene_id <- rownames(counts_mat)

    original_tpm_df <- as.data.frame(counts_mat)
    colnames(original_tpm_df) <- paste("original",colnames(counts_mat),sep="_")
    plot_df <- cbind( log_tpm_df, original_tpm_df )
    non_zero_genes <- plot_df %>% dplyr::filter(HCT116_untreated > 0) %>% pull(gene_id)

    log_tpm_df <- log_tpm_df %>% dplyr::filter(gene_id %in% non_zero_genes)
    num_genes <- 1000
    rna_seq_logFC_cols <- c()
    dosage_cols <- c()
    for (column in colnames(log_tpm_df)) {
        category_vec <- rep("unchanged",length=nrow(log_tpm_df))
        if (grepl("_treated",column)) {
            logFC <- log_tpm_df[[column]] - log_tpm_df$HCT116_untreated
            dosage <- str_match(column,"_[0-9].*?nM") %>% gsub("_","",.) %>% as.character
            dosage_cols <- c(dosage_cols,dosage)
            logFC_col <- paste("rna_seq_logFC",dosage,sep="_")
            rna_seq_logFC_cols <- c(rna_seq_logFC_cols,logFC_col)
            log_tpm_df[[logFC_col]] <- logFC
            
            up_thresh <- sort(logFC,decreasing=T)[num_genes]
            down_thresh <- sort(logFC)[num_genes]
            category_vec[logFC > up_thresh] <- "up"
            category_vec[logFC < down_thresh] <- "down"
            log_tpm_df[[paste("category",dosage,sep="_")]] <- category_vec
            log_tpm_df <- dplyr::rename(log_tpm_df,!!(paste0("rna_seq_",column)):=all_of(column))
        }
    }
    log_tpm_df <- dplyr::rename(log_tpm_df,rna_seq_HCT116_untreated=HCT116_untreated)

    return(log_tpm_df)
}

