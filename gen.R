library(R.utils)
library("GenomicFeatures")
library("Biostrings")
library(ShortRead)

C_folder <- "simseq_sc/"
C_compiled_code <- "build/simseq_sc"
output_folder <- "output/"
read_len <- 91
fragment_mean <- 400
fragment_std <- 100

flat_file <- function(file_name,is_gz,read_no) {
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
    print(paste("Read:",idx,"; Barcode_idx:",barcode_idx))
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
  writeLines(output, file_name)
}

# Directories and file paths
input_dir = "input/"
cdna <- "/home/vunguyen/tin_sinh/input/fasta_folder/Homo_sapiens.GRCh38.106.cdna.all.clean.fa"
count_file <- paste(input_dir,"transcript_expression.csv",sep = "")

# Load transcript sequences
tx.all.fasta <- readDNAStringSet(cdna)
all_fasta_name <- names(tx.all.fasta)

# Read transcript expression count
read_count <- read.csv(count_file, header = TRUE, row.names = 1)
# Remove rows where all values are 0
read_count <- read_count[rowSums(read_count != 0) > 0, ]
# Remove columns where all values are 0
read_count <- read_count[, colSums(read_count != 0) > 0]

dir.create(paste0(output_folder,"/fasta"))

# Simulate reads and generate FASTQ files
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
  flat_file(selected_fasta_file,0,0)
  
  gzip(selected_fasta_file)
  
  system(paste0("cat ",selected_fasta_file,".gz >> ",output_folder,"/fasta/final_fasta.fa.gz"))
  file.remove(paste0(selected_fasta_file,".gz"))

  print(paste("Cell:",i))
}

select_fasta_folder <- paste0(output_folder,"/fasta/")
read_count_file <- paste(output_folder,"/readcount.txt",sep="")
command <- paste("./",C_folder,C_compiled_code," ",select_fasta_folder, " -o ",output_folder," -b 0 -rc ",read_count_file," -r ",read_len," -fm ",fragment_mean," -fs ",fragment_std,sep = "")
system(command)

input_fasta <- paste0(output_folder,"final_fasta.fa_sim_1.fasta.gz")
barcode <- colnames(read_count)[i]
  
  # Read 1
  read_1 <- paste(output_folder, "final_fasta.fa_sim_1.fasta.gz", sep = "")
  output_read1_fasta <- (paste(output_folder, "read1.fasta", sep = ""))
  flat_file(read_1, 1, 1)
  
  # Read 2
  read_2 <- paste(output_folder, "selected_fasta_sim_1.fasta.gz", sep = "")
  output_read2_fasta <- (paste(output_folder, "read2.fasta", sep = ""))
  flat_file(read_2, 1, 0)