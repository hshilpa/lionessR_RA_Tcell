#*******************************************************************************
#---------------- DESEQ2 SCRIPT -----------------------------------
#*******************************************************************************

# STEPS ----
# 0. load libraries
# 1. set paths for input and output files 
# 2. upload countsmatrix and metadata (sample info)
# 3. create DESeqDataSet object
# 4. perform deseq in three steps:
#   a. estimate sizeFactors
#   b. estimate dispersions
#   c. perform Wald test
# 5. extract results
# 6. Save tables
#

#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
### Load Libraries ----
#------------------------------------------------------------------------------#


library("DESeq2")



#------------------------------------------------------------------------------#
### Set paths ----
#------------------------------------------------------------------------------#




# setwd to a directory in your folder
setwd(".")
outdir = "./RA_RNAseq/GSE118829/counts"

# set IP directory 
indir = "./RA_RNAseq/GSE118829/counts"

# dataset name
GSE = "GSE118829"

#------------------------------------------------------------------------------#
### Load tables ----
#------------------------------------------------------------------------------#

# load the subset for which we are doing deseq2
metadata = read.csv(
  "./annot_tables/GSE118829_sampleDetails.csv")
metadata = metadata[,c(
  "Run","cell_type","gender","disease_status","tissue")]

metadata$condition <- paste(metadata$cell_type,metadata$disease_status,sep="_")
metadata$condition <- gsub(" ","_",metadata$condition)

metadata$condition <- as.factor(metadata$condition)
metadata$gender <-  as.factor(metadata$gender)


#------------------------------------------------------------------------------#
### Load tables ----
#------------------------------------------------------------------------------#

# load the gene counts dataframe
genecountTable0 = read.csv("./RA_RNAseq/GSE118829/counts/2021_04_23_13_26_28_geneCounts_STAR.tsv",sep="\t")

# covert to matrix
genecountTable = as.matrix(
  genecountTable0[,(colnames(genecountTable0) %in% metadata$Run)])

rownames(genecountTable) <- genecountTable0$GeneID

#------------------------------------------------------------------------------#
### Construct Deseq2 object ----
#------------------------------------------------------------------------------#

dds <- DESeqDataSetFromMatrix(genecountTable,
                              metadata,
                              design = ~ 0 + condition + gender)

#------------------------------------------------------------------------------#
### DESeq ----
#------------------------------------------------------------------------------#

# a. estimate size factors
dds <- estimateSizeFactors(dds)

# filter lowcount rows

keep <- (rowSums(counts(dds, normalized=TRUE) > 0))
# there should be no 0 count rows 

#
dds <- dds[keep,]

# b estimate dispersions
dds <- estimateDispersions(dds)


# perform Walds test
dds <- nbinomWaldTest(dds)

#------------------------------------------------------------------------------#
### Extract results ----
#------------------------------------------------------------------------------#
celltypes <- unique(gsub(" ","_", metadata$cell_type))

for (i in 1: length(celltypes)){
  
  RAvhealthy <- results(dds,
                      contrast = c("condition",
                                   paste0(celltypes[i],"_RA_non_treatment"),
                                   paste0(celltypes[i],"_healthy")))
  RAvhealthy <- cbind(geneID = rownames(RAvhealthy),RAvhealthy)                      
  
  RAvhealthy<- merge(RAvhealthy,
                   genecountTable0[,c("GeneID","gene_name")],
                   by.x = "geneID", by.y = "GeneID",
                   all.x = T, all.y = F)
  write.csv(RAvhealthy,
          file.path(outdir,
                    paste0(GSE,"/logFC/RAvHealthy/",GSE,"_STAR_deseq2_genes_healthyvRA_",celltypes[i],".csv")),
          row.names = F)
  
  IFX <- results(dds,
                        contrast = c("condition",
                                     paste0(celltypes[i],"_RA_IFX_treatment"),
                                     paste0(celltypes[i],"_RA_non_treatment")))
  IFX <- cbind(geneID = rownames(IFX),IFX)                      
  
  IFX<- merge(IFX,
                     genecountTable0[,c("GeneID","gene_name")],
                     by.x = "geneID", by.y = "GeneID",
                     all.x = T, all.y = F)
  write.csv(IFX,
            file.path(outdir,
                      paste0(GSE,"/logFC/treatedvUntreated/",GSE,"_STAR_deseq2_genes_IFX_",celltypes[i],".csv")),
            row.names = F)  
  
  MTX <- results(dds,
                 contrast = c("condition",
                              paste0(celltypes[i],"_RA_MTX_treatment"),
                              paste0(celltypes[i],"_RA_non_treatment")))
  MTX <- cbind(geneID = rownames(MTX),MTX)                      
  
  MTX<- merge(MTX,
              genecountTable0[,c("GeneID","gene_name")],
              by.x = "geneID", by.y = "GeneID",
              all.x = T, all.y = F)
  write.csv(MTX,
            file.path(outdir,
                      paste0(GSE,"/logFC/treatedvUntreated/",GSE,"_STAR_deseq2_genes_MTX_",celltypes[i],".csv")),
            row.names = F)  
  
  TCZ <- results(dds,
                 contrast = c("condition",
                              paste0(celltypes[i],"_RA_TCZ_treatment"),
                              paste0(celltypes[i],"_RA_non_treatment")))
  TCZ <- cbind(geneID = rownames(TCZ),TCZ)                      
  
  TCZ<- merge(TCZ,
              genecountTable0[,c("GeneID","gene_name")],
              by.x = "geneID", by.y = "GeneID",
              all.x = T, all.y = F)
  write.csv(TCZ,
            file.path(outdir,
                      paste0(GSE,"/logFC/treatedvUntreated/",GSE,"_STAR_deseq2_genes_TCZ_",celltypes[i],".csv")),
            row.names = F)  
  
  
}


#------------------------------------------------------------------------------#
### SAVE ----
#------------------------------------------------------------------------------#
# normalised counts 
norm_counts <- as.data.frame(counts(dds, normalized=TRUE))
norm_counts$GeneID <- rownames(counts(dds))

write.csv(norm_counts,"GSE118829_STAR_deseq2_genes_countsmat_normalized.csv",row.names = F)

rlog_counts <- rlog(dds, blind = F)
rlog_counts$GeneID <- rownames(counts(dds))

write.csv(rlog_counts,"GSE118829_STAR_deseq2_genes_rlog_countsmat.csv",row.names = F)

