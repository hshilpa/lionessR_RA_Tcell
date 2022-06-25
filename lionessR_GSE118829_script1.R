## subset genes and perform lionessR
## perform limma between edge weights to compare means


# script for lionessR on the server
# steps
# inputs: 
# 1. maxMean table
# 2. covdesc file: paste source and celltype
# 3. key enzymes table
# 4. KEGGmet table for pathway IDs
# 5. KEGG pathway annotation for maxMean genes done in local machine using annotation hub
# 6. GTRD transcription factor IDs
# Step 1: subset from maxMean: glycolysis, keyenzymes, TFs
# Step 2: create summarized experiment
# Step 3: calculate correlation for RA and ctrl celltype, find difference
# Step 4: run lionessR
# Step 5: run limma on edge weights
# Step 6: save topTable
# Step 7: convert to edge table and save
# Step 9: Limma on genes
# from edge table:
# find edges with fdr <= 0.1
# find genes with fdr <= 0.1 and in edges with fdr <= 0.1

# steps -------------------------------------
# Step 0: inputs, libraries ----------
# setwd
setwd(".")
dir.create("withoutNetDiffgt05")
setwd("./withoutNetDiffgt05")
annot_dir <- "annot_tables"

# libraries
library(lionessR)
library(SummarizedExperiment)
library(reshape2)
library(igraph)
library(limma)
library(ggplot2)
library(DESeq2)

# inputs ---- 
## 1. maxMean table

maxmean0 <- read.csv("GSE118829_STAR_deseq2_genes_rlog_countsmat.csv",header=T)

## 2. covdesc file: paste source and celltype ----
covdesc = read.csv("../GSE118829_sampleDetails.csv",
                     header=T)
covdesc <- covdesc[,c("Run","cell_type",
                      "disease_status","gender",
                      "source_name","tissue")]
covdesc$condition <- paste(covdesc$disease_status,covdesc$cell_type,sep="_")
covdesc$condition <- gsub(" ","_",covdesc$condition)
covdesc$cell_type <- gsub(" ","_",covdesc$cell_type)

## remove synovial fluid cells because of low samples size -----
covdesc = covdesc[which(covdesc$tissue == "blood"),]

## arrange maxMean columns according to covdesc order -----
covdesc = covdesc[which(covdesc$Run %in% colnames(maxmean0)),]
col.order <- covdesc$Run
row.names(maxmean0) <- maxmean0$GeneID
maxmean0 <- maxmean0[,c(col.order)]


## 3. key enzymes table -----
keyenz_df  <- read.delim(file.path(annot_dir,"glycolysis_keyEnzymes.tsv"),
                         sep = "\t",header=T)
keyenz <- keyenz_df$ENSEMBL

## 4. KEGGmet table for pathway IDs ----
KEGGmet <- read.delim(file.path(annot_dir,"2021_09_03_metabolic_pathways.tsv"),
                      sep="\t",header=T)
KEGGmet$pathwayclass <- gsub("Metabolism; ","",KEGGmet$pathwayclass)
KEGGmet$name <- gsub(" - Homo sapiens \\(human\\)","",KEGGmet$name)

metpathIDs <- KEGGmet$pathways

## 5. KEGG pathway annotation for maxMean genes done in local machine using annotation hub
metpathannot <- read.delim(file.path(annot_dir,"GSE118829_rlogcountsmat_KEGGpathway_mapping.tsv"),
                           sep="\t",header=T)
glyGenes <- metpathannot$ENSEMBL[grep("hsa00010",metpathannot$metPaths)]

## 6. GTRD transcription factor IDs
GTRDTF <- read.delim(file.path(annot_dir,"2021_09_1_TFonly_annotation.tsv"),
                     sep = "\t", header = T)
TF <- GTRDTF$TF_ENSEMBL[!is.na(GTRDTF$TF_ENSEMBL)]

# Step 1: subset from maxMean: glycolysis, keyenzymes, TFs -------
all_sel <- unique(c(keyenz,glyGenes,TF))
maxMean <- maxmean0[rownames(maxmean0) %in% all_sel,]
rm(maxmean0)
rm(all_sel)
# Step 2: create summarized experiment----------
rowData <- DataFrame(row.names = rownames(maxMean), gene = rownames(maxMean))
colData <- DataFrame(row.names = factor(covdesc$Run), 
                     sample = as.character(covdesc$Run), 
                     condition = covdesc$condition,
                     sex = covdesc$gender)

dat <- SummarizedExperiment(assays = list(counts = as.matrix(maxMean)), 
                           colData = colData, rowData = rowData)

rm(rowData)
rm(colData)


# for each contrast ----
celltypes <- unique(covdesc$cell_type)
celltypes <- gsub(" ","_",celltypes)
contrastlist_RAvhealthy <- sapply(celltypes,
                                  function(x){
                                    paste0("groupRA_non_treatment_",x," - ","grouphealthy_",x)})

contrastlist_RAvIFX <- sapply(celltypes,
                       function(x){
                         paste0("groupRA_IFX_treatment_",x," - ","groupRA_non_treatment_",x)})

contrastlist_RAvMTX <- sapply(celltypes,
                              function(x){
                                paste0("groupRA_MTX_treatment_",x," - ","groupRA_non_treatment_",x)})

contrastlist_RAvTCZ <- sapply(celltypes,
                              function(x){
                                paste0("groupRA_TCZ_treatment_",x," - ","groupRA_non_treatment_",x)})


sapply(celltypes,dir.create)           
# Step 3: calculate correlation for RA and ctrl celltype, find difference ----

# condition specific networks
# netRA <- cor(t(assay(dat)[, grep("RA",dat$condition)]))
# netctrl  <- cor(t(assay(dat)[, grep("control",dat$condition)]))
# netdiff <- netRA-netctrl
# rm(netRA)
# rm(netctrl)

# Step 4: get edgelabels to be selected from the square matrix ----
edgeLabels0 = sapply(rownames(dat),
                     function(x){paste(x,
                                       rownames(dat),
                                       sep="_")})

edgelabels = melt(upper.tri(edgeLabels0))
edgelabels = edgelabels[which(edgelabels$value),]
values = edgeLabels0[which(upper.tri(edgeLabels0))]
edgelabels = cbind(edgelabels[,1:2],values)
tosel = edgelabels$values


rm(edgeLabels0)
rm(edgelabels)
rm(values)


# Step 5: run lionessR -----
cormat <- lioness(dat, netFun)

# Step 6: subset edges with diff > 0.5 ----
corsub <- assay(cormat[which(row.names(cormat) %in% tosel), ])
rm(cormat)
corsub_df <- data.frame(corsub)
corsub_df$edge <- rownames(corsub)
write.table(corsub_df,"GSE118829_edgeweight_table.tsv",
            sep="\t",row.names = F,col.names=T,quote=F)
rm(corsub_df)


# Step 7: run limma on edge weights ----
group <- factor(dat$condition)
sex <- factor(dat$sex)
design <- model.matrix(~0+group+sex)

for (i in 1:length(celltypes)) {

  ## RA untreated v healthy ----
  cont.matrix <- makeContrasts(contrasts = contrastlist_RAvhealthy[i],
                               levels = design)
  fit <- lmFit(corsub, design)

  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2e <- eBayes(fit2)
  toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")

  rm(fit)
  rm(fit2)
  rm(fit2e)
  

  # Step 8: save topTable 
  toptable$edges <- row.names(toptable)
  write.table(toptable,
              file.path(celltypes[i],
                        paste0(celltypes[i],"_edgeweightLIMMA_table_RAvhealthy.tsv")),
              sep="\t",row.names = F,col.names=T,quote=F)

 rm(toptable)
 


 

 
 ## RA v IFX ----
 cont.matrix <- makeContrasts(contrasts = contrastlist_RAvIFX[i],
                              levels = design)
 fit <- lmFit(corsub, design)
 
 fit2 <- contrasts.fit(fit, cont.matrix)
 fit2e <- eBayes(fit2)
 toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")
 
 rm(fit)
 rm(fit2)
 rm(fit2e)
 
 
 # Step 8: save topTable 
 toptable$edges <- row.names(toptable)
 write.table(toptable,
             file.path(celltypes[i],
                       paste0(celltypes[i],"_edgeweightLIMMA_table_IFXvUntreated.tsv")),
             sep="\t",row.names = F,col.names=T,quote=F)
 
 rm(toptable)
 
 ## RA v MTX ----
 cont.matrix <- makeContrasts(contrasts = contrastlist_RAvMTX[i],
                              levels = design)
 fit <- lmFit(corsub, design)
 
 fit2 <- contrasts.fit(fit, cont.matrix)
 fit2e <- eBayes(fit2)
 toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")
 
 rm(fit)
 rm(fit2)
 rm(fit2e)

 # Step 8: save topTable 
 toptable$edges <- row.names(toptable)
 write.table(toptable,
             file.path(celltypes[i],
                       paste0(celltypes[i],"_edgeweightLIMMA_table_MTXvUntreated.tsv")),
             sep="\t",row.names = F,col.names=T,quote=F)
 
 rm(toptable)
 
 
 
 
 
 ## RA v TCZ ----
 cont.matrix <- makeContrasts(contrasts = contrastlist_RAvTCZ[i],
                              levels = design)
 fit <- lmFit(corsub, design)
 
 fit2 <- contrasts.fit(fit, cont.matrix)
 fit2e <- eBayes(fit2)
 toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")
 
 rm(fit)
 rm(fit2)
 rm(fit2e)

 # Step 8: save topTable 
 toptable$edges <- row.names(toptable)
 write.table(toptable,
             file.path(celltypes[i],
                       paste0(celltypes[i],"_edgeweightLIMMA_table_TCZvUntreated.tsv")),
             sep="\t",row.names = F,col.names=T,quote=F)
 
 rm(toptable)
}
















