args <- commandArgs(TRUE)
mode <- if(length(args)) args[1] else "check"
genome <- if(length(args)>1) args[2] else "hg38"
txdb <- paste0("TxDb.Hsapiens.UCSC.", genome, ".knownGene")
cran <- c("jsonlite", "ggplot2", "statmod")
bioc <- c("csaw", "edgeR", "ChIPseeker", "rtracklayer", "GenomicFeatures",
          "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "enrichplot", "Gviz", txdb)
if(mode=="install") {
  if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
  missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
  missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
}
status <- vapply(c(cran,bioc),requireNamespace,logical(1),quietly=TRUE)
print(status)
if(!all(status)) quit(status=1)
