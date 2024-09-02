# Install CRAN packages
install.packages(c("gridExtra", "Rtsne", "ggplot2", "tidyr", "SummarizedExperiment"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(c("biomaRt", "AnnotationDbi", "GenomicFeatures", "proActiv"))
