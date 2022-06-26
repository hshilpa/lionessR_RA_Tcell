

# identify edges with FDR <= 0.1
# find TFs and targets from GTRD in them
# subset TF-target and TF-TF edges
# assign direction to the edges
# identify DEGs in these.

# subsetting from FDR tables
# get list of genes in FDR table
# find edges in GTRD with one or both nodes in this list
# save table
# create noa for these nodes

setwd(".")
dir.create("GTRD_annotated_graphs")
# upload GTRD network


GTRD <- read.delim("annot_tables/20210901_TF_target_annotated.tsv",
                   header=T,sep="\t")
GTRDTF0 <- read.delim("annot_tables/2021_09_1_TFonly_annotation.tsv",
                     sep="\t",header=T)

# create node annotation for GTRD
GTRDTF0 <- GTRDTF0[!is.na(GTRDTF0$TF_SYMBOL),]
GTRDTF <- GTRDTF0[,c("TF_SYMBOL","TF_ENSEMBL")]
GTRDTarget <- unique(GTRD[,c("target_SYM","target_ensembl")])
targetnode <- setdiff(GTRDTarget$target_ensembl,GTRDTF$TF_ENSEMBL)
GTRDTarget <- unique(GTRDTarget[which(GTRDTarget$target_ensembl %in% targetnode),])
GTRDTF$type <- "tf"
GTRDTarget$type <- "target"
names(GTRDTF) <- c("SYMBOL","ENSEMBL","type")
names(GTRDTarget) <- c("SYMBOL","ENSEMBL","type")
annotGTRD <- rbind(GTRDTF,GTRDTarget)
write.table(annotGTRD,"annot_tables/2021_09_1_GTRDnoa_ensembl.tsv",
            sep = "\t", row.names = F, col.names = T, quote=F)
rm(GTRDTF)
rm(GTRDTarget)
rm(targetnode)
rm(GTRDTF0)

# create function for annotation and subsetting ----

subset_and_annotate <- function(celltype,edgefile,contrastname) {
  
  edgeTable <- read.delim(edgefile,sep="\t",header=T)
  FDRedgeTable <- edgeTable[which(edgeTable$adj.P.Val <= 0.1),]
  if(nrow(FDRedgeTable)==0){print("no significant edges\n")}
  if(nrow(FDRedgeTable)>0) {
    genesInFDRedge <- unique(unlist(strsplit(FDRedgeTable$edges, 
                                             split="_")))
    
    # write to file
    write.table(FDRedgeTable,
                file.path("GTRD_annotated_graphs",
                          paste0(celltype,"_",contrastname,"_glycolysisFDR01edgeweightLIMMA_table.tsv")),
                sep="\t",row.names = F,quote=F)
    
    # find these in GTRD
    FDR_GTRD_0 <- GTRD[which(GTRD$TF_ensembl %in% genesInFDRedge |
                               GTRD$target_ensembl %in% genesInFDRedge),]
    if(nrow(FDR_GTRD_0) == 0){print("no significant edges in GTRD\n")}
    
    if(nrow(FDR_GTRD_0)){
      # add logFC, adj pval and av exp for these edges
      # flip edges 
      FDRedgeTable$g1 <- lapply(strsplit(FDRedgeTable$edges,split="_"),"[[",1)
      FDRedgeTable$g2 <- lapply(strsplit(FDRedgeTable$edges,split="_"),"[[",2)
      FDRedgeTable$newedge <- paste(FDRedgeTable$g2,FDRedgeTable$g1, sep="_")
      
      # add a columns with pasted edges in GTRD
      FDR_GTRD_0$TF_target <- paste(FDR_GTRD_0$TF_ensembl,FDR_GTRD_0$target_ensembl,
                                    sep="_")
      
      # merge edges with TF_target
      oldedge_FDR <- merge(FDR_GTRD_0,FDRedgeTable,
                           by.x="TF_target",by.y="edges",
                           all = F)
      
      # merge newedge with TF_target
      newedge_FDR <- merge(FDR_GTRD_0,FDRedgeTable,
                           by.x="TF_target",by.y="newedge",
                           all = F)
      names(oldedge_FDR)[19] <- "edges"
      
      #newedge_FDR <- newedge_FDR[,names(oldedge_FDR)]
      
      tab = rbind(oldedge_FDR,newedge_FDR)
      #names(tab)
      tab <- unique(tab[,c("TF_SYM","target_SYM","TF_ensembl","target_ensembl",
                    "SiteCount","logFC","AveExpr","adj.P.Val","TF_target")])
      
      # create noa for these nodes
      noa_nodes <- union(tab$TF_ensembl,tab$target_ensembl)
      
      noa <- data.frame(ENSEMBL = noa_nodes,
                        TF = ifelse(noa_nodes %in% annotGTRD$ENSEMBL[which(annotGTRD$type=="tf")],"tf","target"))
      
      noa$nodes <- annotGTRD$SYMBOL[match(noa$ENSEMBL,annotGTRD$ENSEMBL)]
      
      # write edge and node tables
      write.table(tab,
                  file.path("GTRD_annotated_graphs",
                            paste0(celltype,"_",contrastname,"_glycolysisFDR_GTRD_edges.tsv")),
                  sep="\t",row.names = F,quote=F)
      write.table(noa,
                  file.path("GTRD_annotated_graphs",
                            paste0(celltype,"_",contrastname,"_glycolysisFDR_GTRD_nodes.tsv")),
                  sep="\t",row.names = F,quote=F)
      
    }
  }
  
  
}
  
# upload edgetables RAvhealthy ----
files_edgetables <- list.files(pattern="_edgeweightLIMMA_table_RAvhealthy.tsv",
                               recursive = T,full.names = T  )

celltypes <- sapply(strsplit(files_edgetables,split="/"),"[[",2)
contrastname <- "RAvhealthy"
# call function ----
test <- subset_and_annotate(celltype = celltypes[1],edgefile = files_edgetables[1], contrastname = contrastname)
for(i in 1:length(celltypes)) {
  subset_and_annotate(celltype = celltypes[i],edgefile = files_edgetables[i], contrastname = contrastname)
}

# upload edgetables IFXvUntreated ----
rm(files_edgetables)
rm(celltypes)
rm(contrastname)
files_edgetables <- list.files(pattern="_edgeweightLIMMA_table_IFXvUntreated.tsv",
                               recursive = T,full.names = T  )

celltypes <- sapply(strsplit(files_edgetables,split="/"),"[[",2)
contrastname <- "IFXvUntreated"
for(i in 1:length(celltypes)) {
  subset_and_annotate(celltype = celltypes[i],edgefile = files_edgetables[i], contrastname = contrastname)
}


# upload edgetables MTXvUntreated ----
rm(files_edgetables)
rm(celltypes)
rm(contrastname)
files_edgetables <- list.files(pattern="_edgeweightLIMMA_table_MTXvUntreated.tsv",
                               recursive = T,full.names = T  )

celltypes <- sapply(strsplit(files_edgetables,split="/"),"[[",2)
contrastname <- "MTXvUntreated"
for(i in 1:length(celltypes)) {
  subset_and_annotate(celltype = celltypes[i],edgefile = files_edgetables[i], contrastname = contrastname)
}




# upload edgetables TCZvUntreated ----
rm(files_edgetables)
rm(celltypes)
rm(contrastname)
files_edgetables <- list.files(pattern="_edgeweightLIMMA_table_TCZvUntreated.tsv",
                               recursive = T,full.names = T  )

celltypes <- sapply(strsplit(files_edgetables,split="/"),"[[",2)
contrastname <- "TCZvUntreated"
for(i in 1:length(celltypes)) {
  subset_and_annotate(celltype = celltypes[i],edgefile = files_edgetables[i], contrastname = contrastname)
}




