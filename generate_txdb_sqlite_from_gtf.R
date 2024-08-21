library(AnnotationDbi)
library(GenomicFeatures)

# PREPARE SQLITE FILE FROM ANNOTATION GTF FILE

# Specify the path to your GTF file for the hg38 genome version
gtf_file_hg38 <- 'C:/Users/lucia/OneDrive/Desktop/phD/annotation/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz'

# Check if the hg38 GTF file exists at the specified path
file.exists(gtf_file_hg38)

# Specify the path to your GTF file for the hg19 genome version
gtf_file_hg19 <- 'C:/Users/lucia/OneDrive/Desktop/phD/annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz'

# Check if the hg19 GTF file exists at the specified path
file.exists(gtf_file_hg19)

# Create a TxDb object from the hg38 GTF file, which contains transcript annotations
txdb_hg38 <- txdbmaker::makeTxDbFromGFF(gtf_file_hg38, format = "gtf")

# Create a TxDb object from the hg19 GTF file, which contains transcript annotations
txdb_hg19 <- txdbmaker::makeTxDbFromGFF(gtf_file_hg19, format = "gtf")

# Specify the path to save the SQLite database file for hg38
sqlite_file_hg38 <- 'C:/Users/lucia/OneDrive/Desktop/phD/annotation/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.sqlite'

# Save the TxDb object as a SQLite database file for hg38
saveDb(txdb_hg38, file = sqlite_file_hg38)

# Check if the hg38 SQLite file was successfully created
file.exists(sqlite_file_hg38)

# Specify the path to save the SQLite database file for hg19
sqlite_file_hg19 <- 'C:/Users/lucia/OneDrive/Desktop/phD/annotation/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.sqlite'

# Save the TxDb object as a SQLite database file for hg19
saveDb(txdb_hg19, file = sqlite_file_hg19)

# Check if the hg19 SQLite file was successfully created
file.exists(sqlite_file_hg19)
