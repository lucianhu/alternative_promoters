library(AnnotationDbi)
library(GenomicFeatures)

# PREPARE SQLITE FILE FROM ANNOTATION GTF FILE

# Specify the path to your GTF file for the hg19 genome version
gtf_file_hg19 <- '/path/to/gencode.v19.annotation.gtf.gz'

# Check if the hg19 GTF file exists at the specified path
file.exists(gtf_file_hg19)

# Create a TxDb object from the hg19 GTF file, which contains transcript annotations
txdb_hg19 <- txdbmaker::makeTxDbFromGFF(gtf_file_hg19, format = "gtf", dataSource="gencode", organism = "Homo sapiens")

# Specify the path to save the SQLite database file for hg19
sqlite_file_hg19 <- '/path/to/gencode.v19.annotation.gtf.gz.sqlite'

# Save the TxDb object as a SQLite database file for hg19
saveDb(txdb_hg19, file = sqlite_file_hg19)

# Check if the hg19 SQLite file was successfully created
file.exists(sqlite_file_hg19)
