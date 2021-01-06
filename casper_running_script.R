#!/usr/bin/env Rscript


suppressMessages(library(CaSpER))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomeGraphs))
suppressMessages(library("tximport"))
suppressMessages(library("tximportData"))
suppressMessages(library("DESeq2"))
suppressMessages(library("vsn"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))


data("hg19_cytoband")

## change the following directories
#==============================================================================
## INPUTS
#==============================================================================
studyName="Magill_etal"

## output directory
output_dir="/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/casper_output/Magill_etal/"
system(paste("mkdir -p",output_dir))

## Tumor Expression salmon output directory
salmon_dir="/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/salmon_output/Magill_etal/"

## BAFextract directory
bafextract_dir="/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/bafextract/Magill_etal/"

#==============================================================================
#==============================================================================

### Obtain the expression profiles of tumors

cmd = paste0("ls ",salmon_dir)
tmp=system(cmd,intern = TRUE)
fnames = str_replace(tmp,pattern = "_quant","")
meta=data.frame(Run=fnames)

suppressMessages(library("ensembldb"))
txdb <- makeTxDbFromGFF("/n/data1/bch/genetics/lee/eam63/projects/lncRNA_project/RNA_seq/reference/gencode.v32.chr_patch_hapl_scaff.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

runs <- paste0(meta$Run,"_quant")
files <- file.path(salmon_dir, runs,"quant.sf")
names(files) <- meta$Run
meta$runs <- runs
all(file.exists(files))
txi <- tximport(files, type="salmon", tx2gene=tx2gene,countsFromAbundance="lengthScaledTPM",ignoreAfterBar = TRUE,dropInfReps = TRUE)

tpm_counts <- txi$counts %>% data.frame(stringsAsFactors = FALSE) %>% tibble::rownames_to_column(var = "ensmbl_gene")
tpm_counts$ensmbl_gene = gsub("\\.[0-9]*$","",tpm_counts$ensmbl_gene)
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = tpm_counts$ensmbl_gene,
                                           columns = c("SYMBOL","ENTREZID","GENENAME"),
                                           keytype = "ENSEMBL",
                                          multiVals="first")
tpm_counts <- merge(annotations_orgDb[,c("ENSEMBL","SYMBOL")],tpm_counts, by.x="ENSEMBL",by.y="ensmbl_gene")


write.table(tpm_counts, file=paste0(output_dir,studyName,"_TPM.csv"), sep=",")

data_expr = tpm_counts[,-2]



## download control samples expr
## Note these samples were generated with a combination of samples to be used across tumors
## Might need to build your own control reference TPM values and bafextract outputs 
## following a protocol similar to how we did it for tumors
control_expr <- read.csv(file="/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/sra_files/normal_samples_TPM.csv",row.names=1, stringsAsFactors = FALSE)
control_expr <- control_expr[,-2]
control_names <- colnames(control_expr)[2:ncol(control_expr)]

data_expr <- merge(data_expr,control_expr,by="ENSEMBL")

data_expr <- data_expr[!duplicated(data_expr$ENSEMBL),]
rownames(data_expr) <- data_expr[,1]
data_expr <- data_expr[,-which(colnames(data_expr)=="ENSEMBL")]

## BAFextract output
loh_raw <- readBAFExtractOutput(path=bafextract_dir, sequencing.type="bulk",suffix = ".bcf")

## BAFextract control samples
loh_control <- readBAFExtractOutput(path="/n/scratch3/users/e/eam63/clinical_projects/meningioma/RNAseq/bafextract/normal_samples/", sequencing.type="bulk", suffix=".bcf")

loh_raw <- c(loh_raw, loh_control)

names(loh_raw) <- gsub(".bcf","",names(loh_raw))
loh_raw.name.mapping <- data.frame(loh.name = names(loh_raw), sample.name= colnames(data_expr))

samps_raw = colnames(data_expr)
## generate annotation data.frame
### This part might take a while secondary to connection problems

annotation <- generateAnnotation(id_type="ensembl_gene_id",genes=rownames(data_expr),ishg19=T,centromere,host="uswest.ensembl.org")

data_expr <- data_expr[match(annotation$Gene,rownames(data_expr)),]


## create CaSpER object
object <- CreateCasperObject(raw.data=data_expr, loh.name.mapping= loh_raw.name.mapping, sequencing.type="bulk",
                            cnv.scale=3, loh.scale=3, matrix.type="normalized", expr.cutoff=4.5,annotation=annotation,
                            method="iterative", loh=loh_raw, filter="median",control.sample.ids=control_names,cytoband=cytoband)


## run CaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

##  large scale event summary
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
#common <- intersect(order.sampleNames, intersect(rownames(finalChrMat), rownames(genoMat)))
#finalChrMat <- finalChrMat[match(common, rownames(finalChrMat)), ]
fname = paste0(output_dir,"casper_",studyName,"_largeScale.csv")
write.table(x = finalChrMat,file = fname,sep=",",row.names = TRUE)

## segment based summary
gamma <- 7
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh

loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
loh.final <- loh[loh$count>=gamma, ]

chrom_sizes <- read_delim("/n/data1/bch/genetics/lee/eam63/reference/hg19.chrom.sizes",col_names = c("chrom","length"), delim="\t")
genome_length = sum(chrom_sizes$length)

size_df <- segment.summary$all.summary.loss %>% group_by(ID) %>% summarise(genome_fraction=sum(width)/genome_length) 
fname = paste0(output_dir,"casper_CNV_genomeburden_",studyName,".csv")
write.table(x=size_df,file=fname,sep=",",row.names=FALSE)

all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
    IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

fname = paste0(output_dir,"casper_",studyName,"_genelevel.csv")
write.table(x = rna.matrix,file = fname,sep=",",row.names = TRUE)

cat("CASPER has finished running \n")
cat("Output files are in ",output_dir, "\n")