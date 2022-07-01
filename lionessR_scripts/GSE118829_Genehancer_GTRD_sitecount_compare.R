# compare genehancer and GTRD graphs
# compare binding site counts:
# bar plot of binding site per TF in GTRD and Genehancer

# upload GTRD and Genehancer tables
# for each edge get sitecount and enhancer sitecount
# plot 

library(ggplot2)
library(reshape2)

setwd("withoutNetDiffgt05")
GTRD_dir <- "./GTRD_annotated_graphs"
Genehancer_dir <- "./Genehancer_annotated_graphs"

plotsDir <- "./GTRD_Genehancer_comparison_plots"
outDir <- "./GTRD_Genehancer_comparison"
if(!dir.exists(plotsDir)){dir.create(plotsDir)}
if(!dir.exists(outDir)){dir.create(outDir)}

## get edge and node tables ----
# GTRD edge
GTRD_edges_list = list.files(GTRD_dir,
                             pattern = "_glycolysisFDR_GTRD_edges.tsv",
                             full.names = T)
names(GTRD_edges_list) <- sapply(strsplit(unlist(strsplit(GTRD_edges_list,
                                                          split="_glycolysisFDR_GTRD_edges.tsv")),
                                          split = "/"),"[[",9)


GTRD_edges_tab <- lapply(GTRD_edges_list,read.csv,sep="\t")

# GTRD noa

GTRD_nodes_list = list.files(GTRD_dir,
                             pattern = "_glycolysis_DE_FDRGTRD_nodes.tsv",
                             full.names = T)
names(GTRD_nodes_list) <- sapply(strsplit(unlist(strsplit(GTRD_nodes_list,
                                                          split="_glycolysis_DE_FDRGTRD_nodes.tsv")),
                                          split = "/"),"[[",9)


GTRD_nodes_tab <- lapply(GTRD_nodes_list,read.csv,sep="\t")

# Genehancer edge
Genehancer_edges_list = list.files(Genehancer_dir,
                                   pattern = "_glycolysisFDR_Genehancer_edges.tsv",
                                   full.names = T)
names(Genehancer_edges_list) <- sapply(strsplit(unlist(strsplit(Genehancer_edges_list,
                                                                split="_glycolysisFDR_Genehancer_edges.tsv")),
                                                split = "/"),"[[",9)


Genehancer_edges_tab <- lapply(Genehancer_edges_list,read.csv,sep="\t")

# Genehancer noa

Genehancer_nodes_list = list.files(Genehancer_dir,
                                   pattern = "_glycolysis_DE_FDRGenehancer_nodes.tsv",
                                   full.names = T)
names(Genehancer_nodes_list) <- sapply(strsplit(unlist(strsplit(Genehancer_nodes_list,
                                                                split="_glycolysis_DE_FDRGenehancer_nodes.tsv")),
                                                split = "/"),"[[",9)


Genehancer_nodes_tab <- lapply(Genehancer_nodes_list,read.csv,sep="\t")

celltype_list = sapply(strsplit(names(GTRD_edges_tab),
                                "_"),
                       function(x){paste(x[c(1:3)],collapse="_")})
condition_list = sapply(strsplit(names(GTRD_edges_tab),"_"),"[[",4)


## function to compare and plot ----
compare_and_plot <- function(GTRD_tab,GTRD_noa,Genehancer_tab,Genehancer_noa,celltype,condition){
  if(!dir.exists(file.path(outDir,celltype))) {dir.create(file.path(outDir,celltype))}
  if(!dir.exists(file.path(plotsDir,celltype))) {dir.create(file.path(plotsDir,celltype))}
  
  print(celltype)
  print(condition)
  ## for each table, get the TF-glygene edges alone ----
  GTRD_tfglygene <- GTRD_tab[which(GTRD_tab$target_ensembl %in% 
                                     GTRD_noa$ENSEMBL[which(GTRD_noa$TF == "target")]),]
  ## genehancer table only has TF-glygene edges
  
  ## get edges that are only in GTRD ----
  e_onlyGTRD <- setdiff(GTRD_tfglygene$TF_target,Genehancer_tab$TF_target)
  print("GTRD only edges: ")
  print(length(e_onlyGTRD))
  GTRD_only <- GTRD_tfglygene[which(GTRD_tfglygene$TF_target %in% 
                                      e_onlyGTRD),]
  write.csv(GTRD_only,
            file.path(outDir,celltype,
                      paste0(condition,"_GTRDonly_tfglygene_edges.csv")),
            row.names = F)
  ## get edges that are only in Genehancer ----
  e_onlyGenehancer <- setdiff(Genehancer_tab$TF_target,GTRD_tfglygene$TF_target)
  print("Genehancer only edges: ")
  print(length(e_onlyGenehancer))
  Genehancer_only <- Genehancer_tab[which(Genehancer_tab$TF_target %in% 
                                      e_onlyGenehancer),]
  write.csv(Genehancer_only,
            file.path(outDir,celltype,
                      paste0(condition,"_Genehanceronly_tfglygene_edges.csv")),
            row.names = F)
  ## get edges that are common to GTRD and Genehancer ----
  e_common <- intersect(Genehancer_tab$TF_target,GTRD_tfglygene$TF_target)
  print("common edges: ")
  print(length(e_common))
  Genehancer_GTRD_common <- merge(Genehancer_tab,GTRD_tfglygene,
                                  by = "TF_target",
                                  all = F)
  write.csv(Genehancer_GTRD_common,
            file.path(outDir,celltype,
                      paste0(condition,"_Genehancer_GTRD_common_tfglygene_edges.csv")),
            row.names = F)
  
  ## plot sitecounts for GTRD edges ----
  current_plotdir = file.path(plotsDir,celltype,paste(celltype,condition,sep="_"))
  if(!dir.exists(current_plotdir)){dir.create(current_plotdir)}
  
  ## split GTRD edgetable according to glygene ----
  GTRD_tfglygene_list <- split(GTRD_tfglygene,GTRD_tfglygene$target_ensembl)
  
  ## plot barplots of all TFs in each gene ----
  lapply(GTRD_tfglygene_list, function(df){
    n_tf = nrow(df)
    glygene = unique(df$target_SYM)
    print(glygene)
    #print(length(glygene))
    p<-ggplot(data=df, aes(x=TF_SYM, y=SiteCount)) +
      geom_bar(stat="identity") + 
      labs(title = paste(celltype, condition, glygene, "GTRD Sitecounts")) +
      xlab(label = "Transcription factor") +
      
      theme_bw() +
      theme(title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    ggsave(p,filename = file.path(current_plotdir,
                                  paste0(glygene,"_GTRD_SiteCounts.png")),
           width = 3+n_tf/2, height = 5, units = "in", limitsize = F,
           device = "png",dpi = 300)
  })

  ## split Genehancer edgetable according to glygene ----
  Genehancer_tfglygene_list <- split(Genehancer_tab,Genehancer_tab$glygenes_ensembl)
    
  ## plot barplots of all TFs in each gene ----
  lapply(Genehancer_tfglygene_list, function(df){
    n_tf = nrow(df)
    glygene = unique(df$glygenes_symbol)
    print(glygene)
    #print(length(glygene))
    p <- ggplot(data=df, aes(x=TF_symbol, y=EnhancerSiteCount)) + 
      geom_bar(stat="identity") + 
      labs(title = paste(celltype, condition, glygene, "Genehancer EnhancerSiteCount")) +
      xlab(label = "Transcription factor") +
      
      theme_bw() +
      theme(title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    ggsave(p,filename = file.path(current_plotdir,
                                  paste0(glygene,"_Genehancer_EnhancerSiteCounts.png")),
           width = 3+n_tf/2, height = 5, units = "in", limitsize = F,
           device = "png",dpi = 300)
    })
  
  ## Genehancer GTRD side-by-side ----
  ## split Genehancer_GTRD_common according to glygenes
  Genehancer_GTRD_common_list <- split(Genehancer_GTRD_common,Genehancer_GTRD_common$glygenes_ensembl) 
  
  ## melt each df using SiteCount and Enhancer stiecount
  ## plot both side-by-side -----
  lapply(Genehancer_GTRD_common_list, function(df){
    n_tf = nrow(df) 
    glygene = unique(df$glygenes_symbol)
    print(glygene)
    df_2 <- df[,c("TF_symbol","glygenes_symbol","SiteCount","EnhancerSiteCount")]
    names(df_2) <- c("TF_symbol","glygenes_symbol","GTRDSiteCount","GenehancerEnhancerSiteCount")
    df_long <- melt(df_2,value.name = "count")
    #print(head(df_long))
    
    p <- ggplot(data=df_long, aes(x=TF_symbol, y=count, fill=variable)) +
      geom_bar(stat="identity", width = 0.8, position=position_dodge())+
      geom_text(aes(label=count), vjust = 1,color="black",
                position = position_dodge(1), size=4)+
      scale_fill_manual(values=c('#999999','#E69F00'))+
      labs(title = paste(celltype, condition, glygene,"\n", "Genehancer EnhancerSiteCount and GTRD SiteCount")) +
      xlab(label = "Transcription factor") +
      
      theme_bw() +
      theme(title = element_text(size = 10),
            axis.text = element_text(size = 8),
            legend.title = element_blank(),
            legend.position = "bottom")
    ggsave(p,filename = file.path(current_plotdir,
                                  paste0(glygene,"_Genehancer_GTRD_sidebyside_Counts.png")),
           width = 3+n_tf/2, height = 5, units = "in", limitsize = F,
           device = "png",dpi = 300)
    
    
    
  })
  

  
  
  
}


## test the function ----
# res <- compare_and_plot(GTRD_tab = GTRD_edges_tab[[14]], GTRD_noa = GTRD_nodes_tab[[14]],
#                         Genehancer_tab = Genehancer_edges_tab[[14]], Genehancer_noa = Genehancer_nodes_tab[[14]],
#                         celltype = celltype_list[14], condition = condition_list[14])

## run over all tables ----
for (i in 1:length(GTRD_edges_tab)) {
  res = compare_and_plot(GTRD_tab = GTRD_edges_tab[[i]], 
                         GTRD_noa = GTRD_nodes_tab[[i]],
                         Genehancer_tab = Genehancer_edges_tab[[i]], 
                         Genehancer_noa = Genehancer_nodes_tab[[i]],
                         celltype = celltype_list[i], 
                         condition = condition_list[i])
  
}