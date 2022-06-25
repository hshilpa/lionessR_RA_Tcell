## create node attribute table, annotate with differential expression of genes

# upload FDR GTRD tables
# annotate DEGs
# check for keyenzyme + glycolysis enzymes edges


setwd("GTRD_annotated_graphs")


# keyenzymes + reg ptn list ----
GlyRegPtn <- read.delim("annot_tables/glycolysis_regulatoryPtns.tsv",
                        sep="\t",header=T)
GlyKeyEnzymes <- read.delim("annot_tables/glycolysis_keyEnzymes.tsv",
                            sep="\t",header=T)

all_enzreg <- unique(rbind(GlyRegPtn,GlyKeyEnzymes))
all_enzreg$enzreg <- ifelse(all_enzreg$SYMBOL %in% GlyKeyEnzymes$SYMBOL,
                            "enzyme","regPtn")

# get file names ----
noa_files <- list.files(pattern = "_glycolysisFDR_GTRD_nodes.tsv",
                        full.names = T)
names(noa_files) <- gsub("\\.|/","",sapply(strsplit(noa_files,"_glycolysisFDR_GTRD_nodes.tsv"),"[[",1))

DEGRAvHC_files <- list.files(path="GSE118829/logFC/RAvHealthy",
                           full.names = T,
                           pattern = "GSE118829_STAR_deseq2_genes_healthyvRA_")


names(DEGRAvHC_files) <- gsub("GSE118829_STAR_deseq2_genes_",
                            "",
                            unlist(strsplit(DEGRAvHC_files,
                                            split="GSE118829/logFC/RAvHealthy/"))[c(F,T)])

names(DEGRAvHC_files) <- gsub("\\.csv","",names(DEGRAvHC_files))
treatedvUntreated_files <- list.files(path="GSE118829/logFC/treatedvUntreated",
                             full.names = T,
                             pattern = "GSE118829_STAR_deseq2_genes_")


names(treatedvUntreated_files) <- gsub("GSE118829_STAR_deseq2_genes_treatedvuntreated_",
                              "",
                              unlist(strsplit(treatedvUntreated_files,
                                              split="GSE118829/logFC/treatedvUntreated/"))[c(F,T)])
names(treatedvUntreated_files) <- gsub("\\.csv","",names(treatedvUntreated_files))

DEGall_files = c(DEGRAvHC_files,treatedvUntreated_files)

celltypes_noa = names(noa_files)
celltypes_DEG = names(DEGall_files)

# order the two manually ----
celltypes_DEG = c("MTX_Tcm_CD4",
                  "healthyvRA_CD4Tcm",
                  "TCZ_Tcm_CD4",
                  "IFX_Tcm_CD8",
                  "MTX_Tcm_CD8",
                  "healthyvRA_CD8Tcm",
                  "TCZ_Tcm_CD8",
                  "IFX_Tem_CD4",
                  "MTX_Tem_CD4",
                  "healthyvRA_CD4Tem",
                  "TCZ_Tem_CD4",
                  "IFX_Tem_CD8",
                  "MTX_Tem_CD8",
                  "healthyvRA_CD8Tem",
                  "TCZ_Tem_CD8",
                  "IFX_Temra_CD8",
                  "MTX_Temra_CD8",
                  "healthyvRA_CD8Temra",
                  "TCZ_Temra_CD8",
                  "IFX_Tn_CD4",
                  "MTX_Tn_CD4",
                  "healthyvRA_CD4Tn",
                  "TCZ_Tn_CD4",
                  "IFX_Tn_CD8",
                  "MTX_Tn_CD8",
                  "healthyvRA_CD8Tn",
                  "TCZ_Tn_CD8")
#DEGall_files = DEGall_files[celltypes_DEG]

# create function to annotate ----
annotate_DEGs <- function(celltype,celltype_DEG) {
  
  # upload files
  noa <- read.delim(noa_files[celltype],header=T,sep="\t")
  DEGs <- read.delim(DEGall_files[celltype_DEG],header=T,sep=",")
  DEGs <- DEGs[!is.na(DEGs$padj),]
  print(dim(DEGs))
  upGenes <- DEGs$geneID[which(DEGs$padj <= 0.1 & DEGs$log2FoldChange >= log2(1.5))]
  downGenes <- DEGs$geneID[which(DEGs$padj <= 0.1 & DEGs$log2FoldChange <= -log2(1.5))]
  print(paste("upgenes of ",celltype_DEG,"=",length(upGenes)))
  print(upGenes)
  print(downGenes)
  print(paste("downGenes of ",celltype_DEG,"=",length(downGenes)))
  
  # annotate
  noa$DE <- ifelse(noa$ENSEMBL %in% upGenes,"UP",
                   ifelse(noa$ENSEMBL %in% downGenes,"DOWN","-"))
  
  noa$enzreg <- all_enzreg$enzreg[match(noa$nodes,all_enzreg$SYMBOL)]
  
  # save file
  write.table(noa,
              paste0(celltype,"_glycolysis_DE_FDRGTRD_nodes.tsv"),
              sep="\t",row.names = F,quote=F)
  
  
}

# call function ----
for (i in 1: length(celltypes_DEG)) {
  annotate_DEGs(celltype = celltypes_noa[i],
                celltype_DEG = celltypes_DEG[i])
  
}
