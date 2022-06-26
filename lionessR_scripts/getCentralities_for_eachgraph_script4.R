## create graphs from differentially co-expressed GTRD annotated edges 
## calculate centrality measures

# upload edge and node tables
# create igraph object
# calculate centralities
# save centrality tables

setwd("/home/group_srivatsan01/shilpa/lionessR/GSE118829/withoutNetDiffgt05/GTRD_annotated_graphs")

library(igraph)
# files required;
# edgetable
# nodetable

edge_files <- list.files(pattern = "_glycolysisFDR_GTRD_edges.tsv",
                         full.names = T)
names(edge_files) <- gsub("\\./","",unlist(strsplit(edge_files,
                                                    "_glycolysisFDR_GTRD_edges.tsv"))) 

noa_files <- list.files(pattern = "_glycolysis_DE_FDRGTRD_nodes.tsv",
                        full.names = T)
names(noa_files) <- gsub("\\./","",unlist(strsplit(noa_files,
                                                   "_glycolysis_DE_FDRGTRD_nodes.tsv"))) 

celltypes = names(edge_files)

dir.create("centrality_results")

find_centralities <- function(celltype) {
  
  nodes <- read.delim(noa_files[celltype],sep="\t",header=T)
  #nodes <- nodes[,c("nodes","ENSEMBL","TF","DE","enzreg")]
  edges <- read.delim(edge_files[celltype],sep = "\t",header=T)
  edges <- edges[,c("TF_ensembl","target_ensembl",
                    "TF_SYM","target_SYM",
                    "SiteCount","logFC","AveExpr",
                    "adj.P.Val","TF_target")]
  if (nrow(edges) == 0){print(paste("No edges for",celltype))}
  cell_graph = graph_from_data_frame(d=edges,
                                     directed=T,
                                     vertices=nodes)
  # find in out and total degree of nodes in this graph ----------
  inDegree <- degree(cell_graph,mode="in")
  outDegree <- degree(cell_graph,mode="out")
  totalDegree <- degree(cell_graph,mode="total")
  # find centrality measures in graph -----
  inCloseness <- closeness(cell_graph,mode="in",weights = NA)
  outCloseness <- closeness(cell_graph,mode="out",weights = NA)
  totalCloseness <- closeness(cell_graph,mode="total",weights = NA)
  directedbtw <- betweenness(cell_graph,
                             directed = T,weights = NA)
  hs <- hub_score(cell_graph, weights = NA)$vector
  as <- authority_score(cell_graph, weights = NA)$vector
  # combine to one dataframe --------
  centrality_measures <- data.frame(node = V(cell_graph)$nodes,
                                    ENSEMBL = V(cell_graph)$name,
                                    DE = V(cell_graph)$DE,
                                    TF = V(cell_graph)$TF,
                                    inDegree = inDegree,
                                    outDegree = outDegree,
                                    totalDegree = totalDegree,
                                    inCloseness = inCloseness,
                                    outCloseness = outCloseness,
                                    totalCloseness = totalCloseness,
                                    #eigen = directedEigen,
                                    betweenness = directedbtw,
                                    hubScore = hs,
                                    authorityScore = as)
  # save -------
  write.table(centrality_measures,
              file.path("centrality_results",
                        paste0(celltype,"_lionessR_GTRD_centralities.tsv")),
              sep="\t",row.names = F,quote=F)
  
}
for(cell in 3:length(celltypes)) {
  find_centralities(celltype = celltypes[cell])
  
}
warnings()
