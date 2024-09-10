# Extract promoter coordinates from the promoterAnnotation_gencode_v19 object
# and create a data frame called promote_full_table
promote_full_table <- data.frame(promoterAnnotation_gencode_v19@promoterCoordinates)

# Ensure 'seqnames' is a factor with the specified order
# The order is set to chromosomes 1 through 22, followed by 'chrX', 'chrY', and 'chrM'
promote_full_table$seqnames <- factor(promote_full_table$seqnames, 
                                      levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))

# Sort the data frame first by 'seqnames' and then by 'start' within each 'seqnames'
promote_full_table_sorted <- promote_full_table[order(promote_full_table$seqnames, promote_full_table$start), ]

# Convert any list columns in the data frame to character vectors
# Lists are converted to comma-separated strings for each element
promote_full_table_sorted[] <- lapply(promote_full_table_sorted, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(unlist(y), collapse = ", "))  # Convert list to comma-separated string
  } else {
    x
  }
})

# Save the sorted data frame to a TSV (Tab-Separated Values) file
# The file is saved to the specified path with tab as the separator, no row names, and no quotes around values
write.table(promote_full_table_sorted, file = "/home/lucianhu/rna_snakemake/promote_full_table_sorted.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
