# Optional ID/description mapping only: does not replace genomic coordinates.
args<-commandArgs(TRUE)
if(length(args)!=2)stop('Usage: Rscript biomart_mapping.R config.json ENSEMBL_RELEASE_NUMBER')
suppressPackageStartupMessages({library(jsonlite);library(biomaRt)})
cp<-normalizePath(args[1]);cfg<-fromJSON(cp)
resolve_path <- function(base, path, mustWork=FALSE) {
  if(length(path)!=1L || is.na(path) || !nzchar(path)) stop("Expected one nonempty path")
  candidate <- if(startsWith(path,"/")) path else file.path(normalizePath(base,mustWork=TRUE),path)
  normalizePath(candidate,mustWork=mustWork)
}
root<-resolve_path(dirname(cp),cfg$project_root,mustWork=TRUE);out<-resolve_path(root,cfg$output_dir)
d<-read.delim(file.path(out,'annotation','consensus.nearest_gene.tsv'),stringsAsFactors=FALSE)
ids<-unique(as.character(na.omit(d$geneId)))
version<-as.integer(args[2]);if(is.na(version))stop('A numeric Ensembl release is required for reproducibility')
mart<-useEnsembl(biomart='genes',dataset='hsapiens_gene_ensembl',version=version)
res<-getBM(attributes=c('entrezgene_id','ensembl_gene_id','external_gene_name','description'),filters='entrezgene_id',values=ids,mart=mart)
# Keep one-to-many mappings; do not silently select one Ensembl gene per Entrez gene.
write.table(res,file.path(out,'annotation','BioMart.gene_mapping.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('Ensembl release:',version),'ID mapping only; UCSC coordinates remain unchanged.',capture.output(sessionInfo())),file.path(out,'logs','BioMart.sessionInfo.txt'))
