# Directories and file paths
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_dir> <output_folder>")
}

input_dir <- args[1]
output_folder <- args[2]

cdna <- file.path(input_dir, "Homo_sapiens.GRCh38.106.cdna.all.clean.fa")
count_file <- file.path(input_dir, "transcript_expression.csv")

C_folder <- "simseq_sc/"
C_compiled_code <- "build/simseq_sc"
read_len <- 91
fragment_mean <- 400
fragment_std <- 100




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("R.utils", "GenomicFeatures", "Biostrings", "ShortRead")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  )
}



flat_file <- function(file_name,is_gz,read_no, output_file_name) {
  con <- if (is_gz) gzfile(file_name, open = "r") else file(file_name, open = "r")
  on.exit(close(con))
  
  umi_set <- character()  # Store unique UMIs for this function call
  
  generate_unique_umi <- function(umi_length) {
    repeat {
      umi <- paste0(sample(c("A", "T", "G", "C"), umi_length, replace = TRUE), collapse = "")
      if (!(umi %in% umi_set)) {
        umi_set <<- c(umi_set, umi)
        return(umi)
      }
    }
  }
  
  output <- character(2e7)  # preallocate large enough vector (adjust if needed)
  idx <- 1
  barcode_idx <- 1
  read_mat <- as.matrix(read_count[, barcode_idx, drop = FALSE])
  filter_matrix <- read_mat[rowSums(read_mat) > 0, , drop = FALSE]
  total_current_read <- sum(filter_matrix)
  seq_buffer <- character()

  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, ">")) {
      if (length(seq_buffer) > 0) {
        if(read_no == 1) {
          read <- idx / 2
          while(read > total_current_read) {
            barcode_idx <- barcode_idx + 1
            read_mat <- as.matrix(read_count[, barcode_idx, drop = FALSE])
            filter_matrix <- read_mat[rowSums(read_mat) > 0, , drop = FALSE]
            total_current_read <- total_current_read + sum(filter_matrix)
            umi_set <- character()
          }
          barcode <- colnames(read_count)[barcode_idx]
          umi <- generate_unique_umi(12)
          output[idx] <- paste0(barcode,umi)
        } else {
          output[idx] <- paste0(seq_buffer, collapse = "")
        }
        idx <- idx + 1
        seq_buffer <- character()
      }
      output[idx] <- line
      idx <- idx + 1
    } else {
      seq_buffer <- paste0(seq_buffer, gsub("\\s", "", line))
    }
    
    if(read_no == 1) print(paste0("Apply barcode and UMI for read: ", idx))
    else if(read_no == 2) print(paste0("Formatting read ", idx))
  }
  
  if (length(seq_buffer) > 0) {
    if(read_no == 1) {
      barcode <- colnames(read_count)[barcode_idx]
      umi <- generate_unique_umi(12)
      output[idx] <- paste0(barcode,umi)
    } else {
      output[idx] <- paste0(seq_buffer, collapse = "")
    }
  }
  writeLines(output, output_file_name)
}


# Load transcript sequences
tx.all.fasta <- readDNAStringSet(cdna)
all_fasta_name <- names(tx.all.fasta)

# preprocessing
print("!!!Begin preprocessing:")
read_count <- read.csv(count_file, header = TRUE, row.names = 1)
# Remove rows where all values are 0
read_count <- read_count[rowSums(read_count != 0) > 0, ]
# Remove columns where all values are 0
read_count <- read_count[, colSums(read_count != 0) > 0]
print("!!!Preprocessing done!")

# read generating
print("!!!Begin simulating reads:")
dir.create(paste0(output_folder))
dir.create(paste0(output_folder,"/fasta"))

for (i in 1:ncol(read_count)) {
  readmat_input <- as.matrix(read_count[, i, drop = FALSE])
  selected_fasta_file <- paste0(output_folder,"selected_fasta.fa")
  read_count_file <- paste(output_folder,"/readcount.txt",sep="")

  non_zero_rows <- rowSums(readmat_input) > 0
  filtered_matrix <- readmat_input[non_zero_rows, , drop = FALSE]
  transcript_df <- data.frame(
    TranscriptID = paste0(rownames(filtered_matrix),"Cell",i),
    Count = as.integer(filtered_matrix[, 1])
  )
  write.table(
    transcript_df,
    file = read_count_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )

  transcriptID <- rownames(filtered_matrix)
  x <- strsplit(all_fasta_name, " ")
  x <- do.call(rbind, x)
  y <- x[,2]
  x <- x[,1]
  m <- match(transcriptID, x)
  selected_transcripts <- tx.all.fasta[ m ]
  tx_names <- paste0(x[ m ], "Cell", i, " ", y[ m ])
  names(selected_transcripts) <- tx_names
  writeXStringSet(selected_transcripts, selected_fasta_file)
  flat_file(selected_fasta_file,0,0, selected_fasta_file)

  gzip(selected_fasta_file)
  
  system(paste0("cat ",selected_fasta_file,".gz >> ",output_folder,"/fasta/final_fasta.fa.gz"))
  file.remove(paste0(selected_fasta_file,".gz"))

  print(paste0("Simulating read: Cell ",i, " done!"))
}

select_fasta_folder <- paste0(output_folder,"/fasta/")
read_count_file <- paste(output_folder,"/readcount.txt",sep="")
command <- paste("./",C_folder,C_compiled_code," ",select_fasta_folder, " -o ",output_folder," -b 0 -rc ",read_count_file," -r ",read_len," -fm ",fragment_mean," -fs ",fragment_std,sep = "")
system(command)
print("!!!Simulating read done!")


# apply UMI and barcode
print("!!!Begin apply UMI and barcode for each read:")
barcode <- colnames(read_count)[i]
# Read 1
read_1 <- paste(output_folder, "final_fasta.fa_sim_2.fasta.gz", sep = "")
output_read1_fasta <- paste(output_folder, "read1.fasta", sep = "")
flat_file(read_1, 1, 1, output_read1_fasta)

# Read 2
read_2 <- paste(output_folder, "final_fasta.fa_sim_1.fasta.gz", sep = "")
output_read2_fasta <- paste(output_folder, "read2.fasta", sep = "")
flat_file(read_2, 1, 2, output_read2_fasta)

print("!!!Apply UMI and barcode for each read done!")
