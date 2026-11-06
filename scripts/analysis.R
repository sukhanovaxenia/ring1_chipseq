args <- commandArgs(TRUE)
if(length(args)!=2) stop("Usage: Rscript analysis.R /path/config.json counts|annotate|enrich|loci")
suppressPackageStartupMessages({library(jsonlite); library(ggplot2); library(GenomicRanges); library(SummarizedExperiment)})
configfile <- normalizePath(args[1]); stage <- args[2]
cfg <- fromJSON(configfile)
resolve_path <- function(base, path, mustWork=FALSE) {
  if(length(path)!=1L || is.na(path) || !nzchar(path)) stop("Expected one nonempty path")
  candidate <- if(startsWith(path,"/")) path else file.path(normalizePath(base,mustWork=TRUE),path)
  normalizePath(candidate,mustWork=mustWork)
}
root <- resolve_path(dirname(configfile),cfg$project_root,mustWork=TRUE)
out <- resolve_path(root,cfg$output_dir)
fp <- function(...) file.path(out,...)
for(d in c("tables","figures","beds","annotation","loci","logs")) dir.create(fp(d),recursive=TRUE,showWarnings=FALSE)
file.copy(configfile,fp("logs",paste0(stage,".config.json")),overwrite=TRUE)
samples <- cfg$samples
save_tsv <- function(x,path) write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
save_plot <- function(plot,name,w=8,h=5) ggsave(fp("figures",paste0(name,".pdf")),plot,width=w,height=h)
base_pdf <- function(name,fun,w=7,h=6) {pdf(fp("figures",paste0(name,".pdf")),width=w,height=h);on.exit(dev.off());fun()}
read_bed <- function(path) {
  if(!file.exists(path) || file.info(path)$size==0) return(GRanges())
  d <- read.delim(path,header=FALSE,comment.char="#",stringsAsFactors=FALSE)
  gr <- GRanges(d[[1]],IRanges(d[[2]]+1L,d[[3]]))
  if(ncol(d)>=6) strand(gr) <- ifelse(d[[6]]==".","*",d[[6]])
  if(ncol(d)>=4) names(gr) <- d[[4]]
  gr
}
write_bed <- function(gr,path,named=FALSE) {
  d <- data.frame(chr=as.character(seqnames(gr)),start=start(gr)-1L,end=end(gr))
  if(named){d$name<-names(gr); d$score<-0; d$strand<-as.character(strand(gr));d$strand[d$strand=="*"]<-"."}
  write.table(d,path,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
coords <- function(gr) data.frame(region_id=names(gr),chr=as.character(seqnames(gr)),start0=start(gr)-1L,end=end(gr),width=width(gr))
regions <- read_bed(fp("beds","consensus.bed"))
names(regions) <- sprintf("region_%07d",seq_along(regions))
colors <- c(BF="#2166ac",HS="#b2182b")

if(stage=="counts") {
  suppressPackageStartupMessages({library(csaw);library(edgeR)})
  if(!isTRUE(cfg$biological_replicates_confirmed) && file.exists(fp("tables","differential.HS_vs_BF.tsv"))) stop("Previous differential results exist. Use a new output_dir for a descriptive-only rerun.")
  libs <- read.delim(fp("tables","libraries.tsv"),stringsAsFactors=FALSE)
  chip <- libs[match(samples$sample,libs$sample),]
  control <- libs[match(paste0(samples$sample,"_Input"),libs$sample),]
  chroms <- read.delim(fp("references","chrom.sizes"),header=FALSE)$V1
  # Prepared BAMs have already been filtered identically for MAPQ, flags, chromosomes and blacklist.
  param <- readParam(pe="none",minq=NA,dedup=FALSE,restrict=chroms)
  ext <- if(is.null(cfg$fragment_length)) NA_integer_ else as.integer(cfg$fragment_length)
  rc <- regionCounts(chip$bam,regions,ext=ext,param=param)
  bg <- windowCounts(chip$bam,bin=TRUE,width=cfg$background_bin,ext=ext,param=param)
  stopifnot(identical(as.numeric(rc$totals),as.numeric(bg$totals)))
  f <- normFactors(bg,se.out=FALSE)
  if(any(!is.finite(f)) || any(f<=0) || any(rc$totals<=0)) stop("Invalid normalization factors/depths")
  counts <- assay(rc);colnames(counts)<-samples$sample;rownames(counts)<-names(regions)
  group <- factor(samples$condition,levels=c("BF","HS"))
  y <- DGEList(counts=counts,lib.size=rc$totals,group=group,norm.factors=f)
  design <- model.matrix(~group)
  keep <- filterByExpr(y,design=design)
  if(sum(keep)<10) stop("Too few testable consensus regions; inspect peaks and read counts.")
  logcpm <- cpm(y,log=TRUE,prior.count=2)
  cpm_linear <- cpm(y)
  norm <- data.frame(sample=samples$sample,bam=chip$bam,kind="ChIP",library_size=rc$totals,TMM_factor=f,bamCoverage_scale=1e6/(rc$totals*f))
  # Input is shown as independently depth-normalized QC, not mixed into ChIP TMM fitting.
  inp <- regionCounts(control$bam,regions,ext=ext,param=param)
  if(any(inp$totals<=0))stop("Empty Input library")
  norm <- rbind(norm,data.frame(sample=paste0(samples$sample,"_Input"),bam=control$bam,kind="Input",library_size=inp$totals,TMM_factor=1,bamCoverage_scale=1e6/inp$totals))
  save_tsv(norm,fp("tables","normFactors.txt"))
  save_tsv(cbind(coords(regions),as.data.frame(counts)),fp("tables","consensus.raw_counts.tsv"))
  save_tsv(cbind(coords(regions),as.data.frame(logcpm)),fp("tables","consensus.logCPM.tsv"))
  inputcpm <- sweep(assay(inp),2,inp$totals/1e6,"/")
  input_ratio <- log2((cpm_linear+1)/(inputcpm+1))
  colnames(input_ratio)<-samples$sample
  save_tsv(cbind(coords(regions),as.data.frame(input_ratio)),fp("tables","ChIP_Input.exploratory_log2ratio.tsv"))
  save_tsv(data.frame(region_id=names(regions),tested=keep),fp("tables","tested_regions.tsv"))
  bg_y <- DGEList(assay(bg),lib.size=bg$totals,norm.factors=f)
  bg_log <- cpm(bg_y,log=TRUE,prior.count=2)
  colnames(bg_log)<-samples$sample
  qc <- data.frame(sample=samples$sample,condition=group,eligible_reads=rc$totals,TMM_factor=f,tested_regions=sum(keep))
  save_tsv(qc,fp("tables","normalization_QC.tsv"))
  save_plot(ggplot(qc,aes(sample,eligible_reads/1e6,fill=condition))+geom_col()+scale_fill_manual(values=colors)+theme_bw()+labs(y="Eligible reads (millions)",x=NULL),"library_depth")
  frip<-read.delim(fp("tables","read_QC.tsv"))
  save_plot(ggplot(frip,aes(sample,FRiP_consensus_read_overlap,fill=kind))+geom_col()+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y="Read-overlap FRiP in shared consensus"),"FRiP",10,5)
  for(kind in c("consensus","background")) {
    mat<-if(kind=="consensus")logcpm[keep,,drop=FALSE] else bg_log
    cor_mat<-cor(mat,method="spearman")
    save_tsv(data.frame(sample=rownames(cor_mat),cor_mat),fp("tables",paste0(kind,".spearman.tsv")))
    d<-as.data.frame(as.table(cor_mat));names(d)<-c("sample1","sample2","rho")
    save_plot(ggplot(d,aes(sample1,sample2,fill=rho))+geom_tile()+geom_text(aes(label=sprintf("%.3f",rho)))+scale_fill_gradient(low="white",high="steelblue",limits=c(-1,1))+theme_minimal()+labs(x=NULL,y=NULL),paste0(kind,".correlation"))
    v<-apply(mat,1,var);mat<-mat[is.finite(v)&v>0,,drop=FALSE]
    if(nrow(mat)>5000) mat<-mat[order(apply(mat,1,var),decreasing=TRUE)[1:5000],,drop=FALSE]
    if(nrow(mat)<2){warning("Too few variable regions for PCA: ",kind);next}
    pc<-prcomp(t(mat),center=TRUE,scale.=FALSE)
    frac<-100*pc$sdev^2/sum(pc$sdev^2)
    pp<-data.frame(sample=samples$sample,condition=group,PC1=pc$x[,1],PC2=pc$x[,2])
    save_plot(ggplot(pp,aes(PC1,PC2,color=condition,label=sample))+geom_point(size=3)+geom_text(vjust=-0.8)+scale_color_manual(values=colors)+theme_bw()+labs(x=sprintf("PC1 %.1f%%",frac[1]),y=sprintf("PC2 %.1f%%",frac[2])),paste0(kind,".PCA"))
  }
  # All samples evaluated on the SAME regions. Boxplot dots/rows are regions, not biological replicates.
  long<-data.frame(region_id=rep(names(regions),times=4),sample=rep(samples$sample,each=length(regions)),value=as.vector(logcpm),condition=rep(as.character(group),each=length(regions)))
  save_plot(ggplot(long,aes(sample,value,fill=condition))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=colors)+theme_bw()+labs(x=NULL,y="log2 normalized CPM (prior count 2)",subtitle="Distribution over shared regions; no per-peak pseudoreplication test"),"signal.boxplot")
  ir<-data.frame(sample=rep(samples$sample,each=length(regions)),value=as.vector(input_ratio),condition=rep(as.character(group),each=length(regions)))
  save_plot(ggplot(ir,aes(sample,value,fill=condition))+geom_boxplot(outlier.shape=NA)+geom_hline(yintercept=0,linetype=2)+scale_fill_manual(values=colors)+theme_bw()+labs(x=NULL,y="log2[(ChIP TMM CPM + 1)/(Input CPM + 1)]",subtitle="Exploratory background check; no inferential test"),"ChIP_Input.exploratory_boxplot")
  for(g in c("BF","HS")){
    ix<-which(samples$condition==g)
    discordance<-cbind(coords(regions),rep1_logCPM=logcpm[,ix[1]],rep2_logCPM=logcpm[,ix[2]],absolute_logCPM_difference=abs(logcpm[,ix[1]]-logcpm[,ix[2]]))
    discordance<-discordance[order(discordance$absolute_logCPM_difference,decreasing=TRUE),]
    save_tsv(discordance,fp("tables",paste0(g,".replicate_discordance.tsv")))
  }
  scatter<-data.frame(BF1=logcpm[,1],BF2=logcpm[,2],HS1=logcpm[,3],HS2=logcpm[,4])
  save_plot(ggplot(scatter,aes(BF1,BF2))+geom_bin_2d(bins=60)+geom_abline(slope=1,intercept=0,color="red")+coord_equal()+theme_bw()+labs(x="BF rep1 logCPM",y="BF rep2 logCPM"),"BF.replicate_scatter")
  save_plot(ggplot(scatter,aes(HS1,HS2))+geom_bin_2d(bins=60)+geom_abline(slope=1,intercept=0,color="red")+coord_equal()+theme_bw()+labs(x="HS rep1 logCPM",y="HS rep2 logCPM"),"HS.replicate_scatter")
  objects<-list(regions=regions,y=y,keep=keep,logcpm=logcpm,samples=samples,design=design,input_ratio=input_ratio)
  if(isTRUE(cfg$biological_replicates_confirmed)) {
    if(file.exists(fp("logs","differential.NOT_RUN.txt")))unlink(fp("logs","differential.NOT_RUN.txt"))
    yf<-y[keep,,keep.lib.sizes=TRUE]
    yf<-estimateDisp(yf,design,robust=TRUE)
    fit<-glmQLFit(yf,design,robust=TRUE)
    ql<-glmQLFTest(fit,coef=2)
    # Thresholded test: evidence for |log2FC| greater than configured biological effect threshold.
    tr<-glmTreat(fit,coef=2,lfc=cfg$minimum_log2_fold_change)
    tt<-topTags(ql,n=Inf,sort.by="none")$table
    tt$FDR_TREAT<-topTags(tr,n=Inf,sort.by="none")$table$FDR
    res<-cbind(coords(regions[keep]),tt)
    res$direction<-ifelse(res$FDR_TREAT<cfg$fdr,ifelse(res$logFC>0,"HS_higher","BF_higher"),"not_significant")
    save_tsv(res,fp("tables","differential.HS_vs_BF.tsv"))
    for(label in c("HS_higher","BF_higher"))write_bed(regions[keep][res$direction==label],fp("beds",paste0(label,".bed")))
    base_pdf("dispersion",function(){plotBCV(yf)})
    base_pdf("QL_dispersion",function(){plotQLDisp(fit)})
    save_plot(ggplot(res,aes(logFC,-log10(pmax(FDR_TREAT,1e-300)),color=direction))+geom_point(alpha=.5,size=.8)+geom_hline(yintercept=-log10(cfg$fdr),linetype=2)+theme_bw()+labs(x="log2FC HS / BF",y="-log10 TREAT FDR"),"differential.volcano")
    save_plot(ggplot(res,aes(logCPM,logFC,color=direction))+geom_point(alpha=.4,size=.8)+theme_bw()+labs(y="log2FC HS / BF"),"differential.MA")
    objects$results<-res
  } else {
    writeLines("Differential inference skipped: confirm independent biological replicates in config.json. QC remains descriptive.",fp("logs","differential.NOT_RUN.txt"))
  }
  saveRDS(objects,fp("tables","analysis.rds"))
}

if(stage=="annotate") {
  suppressPackageStartupMessages({library(ChIPseeker);library(GenomicFeatures);library(AnnotationDbi);library(org.Hs.eg.db)})
  pkg<-paste0("TxDb.Hsapiens.UCSC.",cfg$genome,".knownGene")
  if(!requireNamespace(pkg,quietly=TRUE))stop("Install ",pkg)
  txdb<-getExportedValue(pkg,pkg)
  org<-org.Hs.eg.db
  sizes<-read.delim(fp("references","chrom.sizes"),header=FALSE,stringsAsFactors=FALSE)
  genome_gr<-GRanges(sizes$V1,IRanges(1,sizes$V2))
  trim_to_genome<-function(gr){strand(gr)<-"*";gr<-gr[as.character(seqnames(gr))%in%sizes$V1];intersect(reduce(gr),genome_gr,ignore.strand=TRUE)}
  genes_gr<-genes(txdb);genes_gr<-genes_gr[as.character(seqnames(genes_gr))%in%sizes$V1]
  # One strand-aware gene model per retained Entrez gene; not every alternative transcript TSS.
  write_bed(genes_gr,fp("beds","genes.bed"),named=TRUE)
  exons_gr<-trim_to_genome(exons(txdb))
  utr5<-trim_to_genome(unlist(fiveUTRsByTranscript(txdb),use.names=FALSE))
  utr3<-trim_to_genome(unlist(threeUTRsByTranscript(txdb),use.names=FALSE))
  tx<-transcripts(txdb)
  prom<-trim_to_genome(promoters(tx,upstream=cfg$promoter_flank,downstream=cfg$promoter_flank+1))
  bodies<-trim_to_genome(genes_gr)
  introns<-setdiff(bodies,exons_gr,ignore.strand=TRUE)
  # Mutually exclusive territory hierarchy for bp-based density; separate from ChIPseeker nearest annotation.
  territories<-list(Promoter=prom,UTR5=utr5,UTR3=utr3,Exon=exons_gr,Intron=introns)
  used<-GRanges()
  for(k in names(territories)){
    territories[[k]]<-setdiff(territories[[k]],used,ignore.strand=TRUE)
    used<-reduce(c(used,territories[[k]]),ignore.strand=TRUE)
  }
  territories$Intergenic<-setdiff(genome_gr,used,ignore.strand=TRUE)
  black<-read_bed(resolve_path(root,cfg$blacklist))
  territories<-lapply(territories,function(gr)setdiff(gr,black,ignore.strand=TRUE))
  saveRDS(list(txdb=txdb,genes=genes_gr,territories=territories),fp("annotation","reference.rds"))
  setfiles<-list.files(fp("beds"),pattern="\\.bed$",full.names=TRUE)
  setfiles<-setfiles[!basename(setfiles)%in%c("genes.bed")]
  annotations<-list();densities<-list();totals<-list()
  for(file in setfiles){
    label<-sub("\\.bed$","",basename(file));gr<-read_bed(file)
    if(!length(gr))next
    names(gr)<-if(label=="consensus")names(regions) else paste0(label,"_",seq_along(gr))
    mcols(gr)$region_id<-names(gr)
    ann<-annotatePeak(gr,TxDb=txdb,annoDb="org.Hs.eg.db",tssRegion=c(-cfg$promoter_flank,cfg$promoter_flank),verbose=FALSE)
    annotations[[label]]<-ann
    tab<-as.data.frame(ann)
    tab$start0<-tab$start-1L
    # ChIPseeker keeps region_id metadata; refuse a positional fallback after sorting.
    if(!"region_id"%in%names(tab))stop("ChIPseeker did not preserve region_id metadata")
    save_tsv(tab,fp("annotation",paste0(label,".nearest_gene.tsv")))
    hits<-findOverlaps(gr,genes_gr,ignore.strand=TRUE)
    overlaps<-data.frame(region_id=names(gr)[queryHits(hits)],geneId=names(genes_gr)[subjectHits(hits)])
    if(nrow(overlaps))overlaps$SYMBOL<-unname(mapIds(org,keys=overlaps$geneId,column="SYMBOL",keytype="ENTREZID",multiVals="first"))
    save_tsv(overlaps,fp("annotation",paste0(label,".gene_body_overlaps.tsv")))
    mid<-resize(gr,1,fix="center");union<-reduce(gr,ignore.strand=TRUE)
    dd<-do.call(rbind,lapply(names(territories),function(k){
      ter<-territories[[k]];area<-sum(width(ter));covered<-sum(width(intersect(union,ter,ignore.strand=TRUE)))
      centers<-sum(countOverlaps(mid,ter,ignore.strand=TRUE)>0)
      data.frame(group=label,category=k,territory_bp=area,peak_centers=centers,centers_per_Mb=if(area)centers/(area/1e6) else NA,covered_bp=covered,percent_territory_covered=if(area)100*covered/area else NA)
    }))
    densities[[label]]<-dd
    totals[[label]]<-data.frame(group=label,n_peaks=length(gr),median_width=median(width(gr)),union_bp=sum(width(union)),nearest_genes=length(unique(na.omit(tab$geneId))),overlapped_genes=length(unique(overlaps$geneId)))
  }
  saveRDS(annotations,fp("annotation","annotations.rds"))
  save_tsv(do.call(rbind,totals),fp("tables","peak_set_summary.tsv"))
  dens<-do.call(rbind,densities);save_tsv(dens,fp("tables","genomic_territory_density.tsv"))
  comparisons<-list(pooled=c("BF_pool","HS_pool"),bdgdiff=c("BF_enriched","HS_enriched","common"),replicates=samples$sample,differential=c("BF_higher","HS_higher"))
  for(n in names(comparisons)){
    selected<-intersect(comparisons[[n]],names(annotations));if(!length(selected))next
    save_plot(plotAnnoBar(annotations[selected]),paste0(n,".annotation"),10,5)
    save_plot(plotDistToTSS(annotations[selected]),paste0(n,".distance_TSS"),10,5)
    dd<-dens[dens$group%in%selected,]
    save_plot(ggplot(dd,aes(category,centers_per_Mb,fill=group))+geom_col(position="dodge")+theme_bw()+labs(x=NULL,y="Peak centers per Mb of eligible territory"),paste0(n,".category_density"),10,5)
    save_plot(ggplot(dd,aes(category,percent_territory_covered,fill=group))+geom_col(position="dodge")+theme_bw()+labs(x=NULL,y="Territory covered by peaks (%)"),paste0(n,".category_coverage"),10,5)
  }
  # Jaccard computed from occupied BASE PAIRS, not asymmetric counts of overlapping intervals.
  labs<-intersect(samples$sample,names(annotations));jac<-matrix(NA,length(labs),length(labs),dimnames=list(labs,labs))
  for(i in labs)for(j in labs){x<-reduce(read_bed(fp("beds",paste0(i,".bed"))));y<-reduce(read_bed(fp("beds",paste0(j,".bed"))));den<-sum(width(reduce(c(x,y))));jac[i,j]<-if(den)sum(width(intersect(x,y)))/den else NA}
  if(length(labs)){
    save_tsv(data.frame(sample=rownames(jac),jac),fp("tables","replicate.Jaccard_bp.tsv"))
    jd<-as.data.frame(as.table(jac));names(jd)<-c("sample1","sample2","Jaccard")
    save_plot(ggplot(jd,aes(sample1,sample2,fill=Jaccard))+geom_tile()+geom_text(aes(label=sprintf("%.3f",Jaccard)))+scale_fill_gradient(low="white",high="steelblue",limits=c(0,1))+theme_minimal()+labs(x=NULL,y=NULL),"replicate.Jaccard_bp")
  }
  if(file.exists(fp("tables","analysis.rds"))){
    a<-readRDS(fp("tables","analysis.rds"));tab<-as.data.frame(annotations$consensus)
    category<-sub(" \\(.*$","",tab$annotation)
    cat_by_id<-setNames(category,tab$region_id)
    d<-data.frame(region_id=rep(names(regions),4),sample=rep(samples$sample,each=length(regions)),value=as.vector(a$logcpm),category=rep(unname(cat_by_id[names(regions)]),4))
    save_plot(ggplot(d,aes(sample,value,fill=sample))+geom_boxplot(outlier.shape=NA)+facet_wrap(~category)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y="log2 normalized CPM",subtitle="Same regions per category; descriptive distributions"),"signal.by_annotation",13,8)
    # Integrated table preserves all tested regions and raw/statistical direction.
    if(!is.null(a$results)){
      ann_tab<-read.delim(fp("annotation","consensus.nearest_gene.tsv"),check.names=FALSE)
      joined<-merge(a$results,ann_tab,by="region_id",all.x=TRUE,suffixes=c("",".annotation"))
      joined<-joined[order(joined$FDR_TREAT),]
      save_tsv(joined,fp("tables","differential.annotated.tsv"))
    }
  }
}

if(stage=="enrich"){
  suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot)})
  file<-fp("tables","differential.annotated.tsv")
  if(!file.exists(file))stop("Run differential counts (with confirmed biological repeats) and annotate first")
  d<-read.delim(file,check.names=FALSE)
  universe<-unique(as.character(na.omit(d$geneId)))
  for(label in c("BF_higher","HS_higher")){
    target<-unique(as.character(na.omit(d$geneId[d$direction==label])))
    if(length(target)<5){writeLines("Fewer than 5 mapped target genes; GO not run.",fp("logs",paste0(label,".GO.skipped.txt")));next}
    go<-enrichGO(gene=target,universe=universe,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=TRUE)
    tab<-as.data.frame(go);save_tsv(tab,fp("tables",paste0(label,".GO_BP.tsv")))
    sig<-tab[!is.na(tab$p.adjust)&tab$p.adjust<cfg$fdr,,drop=FALSE]
    save_tsv(sig,fp("tables",paste0(label,".GO_BP.significant.tsv")))
    if(nrow(sig)){
      go@result<-go@result[go@result$ID%in%sig$ID,,drop=FALSE]
      save_plot(dotplot(go,showCategory=15)+ggtitle(paste(label,"GO BP; BH FDR <",cfg$fdr)),paste0(label,".GO_BP"),10,7)
    } else writeLines("No GO BP term passes the configured BH FDR threshold.",fp("logs",paste0(label,".GO.no_significant_terms.txt")))
  }
}

if(stage=="loci"){
  suppressPackageStartupMessages({library(Gviz);library(rtracklayer)})
  ref<-readRDS(fp("annotation","reference.rds"))
  candidates<-fp("loci","regions_to_plot.tsv")
  if(!file.exists(candidates)){
    diff<-fp("tables","differential.annotated.tsv")
    if(!file.exists(diff))stop("Create loci/regions_to_plot.tsv with chr,start0,end,label columns; see README")
    d<-read.delim(diff,check.names=FALSE)
    d<-d[d$direction!="not_significant",]
    d<-d[order(d$FDR_TREAT),]
    # Up to 3 loci per direction; explicit selection can replace this file.
    selected<-do.call(rbind,lapply(c("BF_higher","HS_higher"),function(k)head(d[d$direction==k,],3)))
    if(!nrow(selected))stop("No significant candidate loci; select descriptive loci manually.")
    z<-selected[,c("chr","start0","end")];z$label<-selected$region_id
    save_tsv(z,candidates)
  }
  selected<-read.delim(candidates,stringsAsFactors=FALSE)
  sizes<-read.delim(fp("references","chrom.sizes"),header=FALSE)
  sizes<-setNames(sizes$V2,sizes$V1)
  for(i in seq_len(nrow(selected))){
    z<-selected[i,]; chr<-z$chr
    if(!chr%in%names(sizes)||is.na(z$start0)||z$start0<0||z$end<=z$start0||z$end>sizes[[chr]])stop("Invalid locus coordinates: ",z$label)
    lo<-max(1,z$start0+1-10000);hi<-min(sizes[[chr]],z$end+10000)
    if(is.null(hi)||!is.finite(hi)||hi<=lo)stop("Invalid locus")
    which<-GRanges(chr,IRanges(lo,hi))
    signals<-lapply(samples$sample,function(s)import(BigWigFile(fp("tracks",paste0(s,".bw"))),which=which))
    ymax<-max(c(1,unlist(lapply(signals,function(x)x$score))),na.rm=TRUE)
    tracks<-list(GenomeAxisTrack(),GeneRegionTrack(ref$txdb,genome=cfg$genome,chromosome=chr,start=lo,end=hi,name="UCSC genes",showId=TRUE,geneSymbol=FALSE))
    for(j in seq_along(signals)){
      gr<-signals[[j]]
      if(!length(gr))gr<-GRanges(chr,IRanges(lo,hi),score=0)
      tracks[[length(tracks)+1]]<-DataTrack(range=gr,data=gr$score,genome=cfg$genome,chromosome=chr,name=samples$sample[j],type="histogram",ylim=c(0,ymax),col.histogram=colors[[samples$condition[j]]],fill.histogram=colors[[samples$condition[j]]])
    }
    # Inputs are useful for local background checks; label their separate CPM scale.
    for(j in seq_len(nrow(samples))){
      path<-fp("tracks",paste0(samples$sample[j],"_Input.bw"));if(!file.exists(path))next
      gr<-import(BigWigFile(path),which=which);if(!length(gr))gr<-GRanges(chr,IRanges(lo,hi),score=0)
      tracks[[length(tracks)+1]]<-DataTrack(range=gr,data=gr$score,genome=cfg$genome,chromosome=chr,name=paste(samples$sample[j],"Input CPM"),type="histogram")
    }
    external<-fp("references","external_tracks.tsv")
    if(file.exists(external)){
      extra<-read.delim(external,stringsAsFactors=FALSE)
      for(k in seq_len(nrow(extra))){
        if(extra$assembly[k]!=cfg$genome)stop("External track assembly mismatch: ",extra$label[k])
        path<-resolve_path(root,extra$path[k])
        bfile<-BigWigFile(path)
        ext_sizes<-seqlengths(seqinfo(bfile))
        if(!"chr1"%in%names(ext_sizes)||ext_sizes[["chr1"]]!=sizes[["chr1"]])stop("External bigWig reference/chr naming mismatch: ",extra$label[k])
        gr<-import(bfile,which=which)
        if(!length(gr))next
        tracks[[length(tracks)+1]]<-DataTrack(range=gr,data=gr$score,genome=cfg$genome,chromosome=chr,name=paste(extra$label[k],"external"),type="histogram")
      }
    }
    tracks[[length(tracks)+1]]<-AnnotationTrack(GRanges(chr,IRanges(z$start0+1,z$end)),genome=cfg$genome,chromosome=chr,name="Selected region")
    name<-gsub("[^A-Za-z0-9_.-]","_",z$label)
    pdf(fp("loci",paste0(name,".pdf")),width=12,height=max(11,length(tracks)*.8))
    plotTracks(tracks,from=lo,to=hi,chromosome=chr,main=paste(z$label,chr,lo,hi))
    dev.off()
  }
}
writeLines(capture.output(sessionInfo()),fp("logs",paste0(stage,".sessionInfo.txt")))
