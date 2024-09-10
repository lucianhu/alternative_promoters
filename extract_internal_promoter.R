# Load necessary libraries
library(proActiv)   # For promoter analysis
library(SummarizedExperiment)  # For working with SummarizedExperiment objects
library(biomaRt)    # For accessing BioMart databases

# Define the GTF file path for promoter annotation
gtf_file_hg19 <- '/home/lucianhu/rna_snakemake/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz'

# Check if the GTF file exists
if (!file.exists(gtf_file_hg19)) {
  stop("GTF file does not exist at the specified path.")
}

# Prepare promoter annotation using the GTF file for Homo sapiens
promoterAnnotation_gencode_v19 <- preparePromoterAnnotation(
  file = gtf_file_hg19, 
  species = 'Homo_sapiens'
)

# Define the TxDb file path for promoter annotation
sqlite_file_hg19 <- '/home/lucianhu/rna_snakemake/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz.sqlite'

# Check if the TxDb file exists
if (!file.exists(sqlite_file_hg19)) {
  stop("TxDb file does not exist at the specified path.")
}

# Load the TxDb object from the SQLite file
txdb_hg19 <- loadDb(sqlite_file_hg19)

# Prepare promoter annotation using the TxDb object for Homo sapiens
promoterAnnotation_gencode_v19 <- preparePromoterAnnotation(
  txdb = txdb_hg19, 
  species = 'Homo_sapiens'
)

# Convert the promoterAnnotation_gencode_v19 object to a data frame
test <- data.frame(promoterAnnotation_gencode_v19@promoterCoordinates)

# Check the names of the columns to confirm they exist
print(names(test))

# Proceed if 'seqnames' and 'internalPromoter' columns exist
if ("seqnames" %in% names(test) & "internalPromoter" %in% names(test)) {
  
  # Create a new column 'pos' combining 'seqnames' and 'start' with a colon separator
  test$pos <- paste0(test$seqnames, ":", test$start)
  
  # Filter the data to keep only rows where 'internalPromoter' is FALSE, NA, or 'seqnames' is NA
  test <- test[test$internalPromoter == TRUE | is.na(test$internalPromoter) | is.na(test$seqnames), ]
  
  # Select only the columns 'seqnames', 'start', and 'end' for further processing
  selected_data <- test[, c("seqnames", "start", "end")]
  
  # Convert 'seqnames' to a factor and set the order from 'chr1' to 'chrY'
  selected_data$seqnames <- factor(selected_data$seqnames, 
                                   levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
  
  # Sort the data first by 'seqnames' and then by 'start'
  sorted_data <- selected_data[order(selected_data$seqnames, selected_data$start), ]
  
  # Print the sorted data without row names and headers, and save to a BED file
  write.table(sorted_data, file = "/home/lucianhu/rna_snakemake/internal_promoters.v19.bed", 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
} else {
  # If columns are missing, print an error message
  stop("The required columns 'seqnames' or 'internalPromoter' are not present in the data.")
}
