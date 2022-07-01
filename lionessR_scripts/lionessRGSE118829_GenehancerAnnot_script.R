## take unannotated graphs
## annotate with Genehancer
## create graphs for each celltype and condition
## save in new directory

# identify edges with FDR <= 0.1
# find TFs and targets from Genehancer in them
# subset TF-target edges
# assign direction to the edges
# identify DEGs in these.

# subsetting from FDR tables
# get list of genes in FDR table
# find edges in Genehancer with one or both nodes in this list
# save table
# create noa for these nodes

setwd(".")
dir.create("Genehancer_annotated_graphs")

# upload Genehancer network
Genehancer <- read.delim("./annot_tables/Genehancer/glygenes_TF_edgeAtt.csv",
                   header=T,sep=",")
Genehancer_nodes <- read.delim("./annot_tables/Genehancer/glygenes_TF_noa.csv",
                      sep=",",header=T)

# node annotation for Genehancer


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
                file.path("Genehancer_annotated_graphs",
                          paste0(celltype,"_",contrastname,"_glycolysisFDR01edgeweightLIMMA_table.tsv")),
                sep="\t",row.names = F,quote=F)
    
    # find these in Genehancer
    # One of the nodes in any gene hancer edge must be from the FDRedgetable
    FDR_Genehancer_0 <- Genehancer[which(Genehancer$TF_ensembl %in% genesInFDRedge |
                               Genehancer$glygenes_ensembl %in% genesInFDRedge),]
    if(nrow(FDR_Genehancer_0) == 0){print("no significant edges in Genehancer\n")}
    
    if(nrow(FDR_Genehancer_0)){
      # add logFC, adj pval and av exp for these edges
      # flip edges 
      FDRedgeTable$g1 <- lapply(strsplit(FDRedgeTable$edges,split="_"),"[[",1)
      FDRedgeTable$g2 <- lapply(strsplit(FDRedgeTable$edges,split="_"),"[[",2)
      FDRedgeTable$newedge <- paste(FDRedgeTable$g2,FDRedgeTable$g1, sep="_")
      
      # add a columns with pasted edges in Genehancer
      FDR_Genehancer_0$TF_target <- paste(FDR_Genehancer_0$TF_ensembl,FDR_Genehancer_0$glygenes_ensembl,
                                    sep="_")
      
      # merge edges of FDR_Genehancer with FDRedgetable for both types of edges
      oldedge_FDR <- merge(FDR_Genehancer_0,FDRedgeTable,
                           by.x="TF_target",by.y="edges",
                           all = F)
      
      # merge newedge with TF_target
      newedge_FDR <- merge(FDR_Genehancer_0,FDRedgeTable,
                           by.x="TF_target",by.y="newedge",
                           all = F)
      names(oldedge_FDR)[16] <- "edges"
      
      #newedge_FDR <- newedge_FDR[,names(oldedge_FDR)]
      
      tab = rbind(oldedge_FDR,newedge_FDR)
      #names(tab)
      tab <- unique(tab[,c("TF_ensembl","glygenes_ensembl","TF_symbol","glygenes_symbol",
                           "EnhancerSiteCount","logFC","AveExpr","adj.P.Val","TF_target","gene_pair")])
      
      # create noa for these nodes
      noa_nodes <- union(tab$TF_ensembl,tab$glygenes_ensembl)
      
      noa <- data.frame(ENSEMBL = noa_nodes,
                        TF = Genehancer_nodes$type[match(noa_nodes,Genehancer_nodes$ENSEMBL)])
      
      noa$nodes <- Genehancer_nodes$SYMBOL[match(noa$ENSEMBL,Genehancer_nodes$ENSEMBL)]
      
      # write edge and node tables
      write.table(tab,
                  file.path("Genehancer_annotated_graphs",
                            paste0(celltype,"_",contrastname,"_glycolysisFDR_Genehancer_edges.tsv")),
                  sep="\t",row.names = F,quote=F)
      write.table(noa,
                  file.path("Genehancer_annotated_graphs",
                            paste0(celltype,"_",contrastname,"_glycolysisFDR_Genehancer_nodes.tsv")),
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




