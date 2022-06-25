### Annotate genes in the dataset GSE118829 with KEGG pathway ID using AnnotationHub ----

# annotate maxMean genes for GSE118829 with KEGG pathIDs
library(AnnotationHub)
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
columns(orgdb)

setwd(".")
annot_dir <- "annot_tables"

# KEGGmet ----
KEGGmet <- read.delim(file.path(annot_dir,"2021_09_03_metabolic_pathways.tsv"),
                      sep="\t",header=T)
KEGGmet$pathwayclass <- gsub("Metabolism; ","",KEGGmet$pathwayclass)
KEGGmet$name <- gsub(" - Homo sapiens \\(human\\)","",KEGGmet$name)

metpathIDs <- KEGGmet$pathways

# maxMean and rlog ----
maxMean <- read.csv("GSE118829_STAR_deseq2_genes_countsmat_normalized.csv",sep=",",header=T)
rlog <- read.csv("GSE118829_STAR_deseq2_genes_rlog_countsmat.csv", sep = ",", header = T)

# annotate genes with KEGG pathway ID ----
path_mapping <- select(orgdb, 
                       keys=maxMean$Gene, 
                       columns=c("PATH","ENTREZID","SYMBOL"), keytype="ENSEMBL")

path_mapping_rlog <- select(orgdb, 
                       keys=rlog$Gene, 
                       columns=c("PATH","ENTREZID","SYMBOL"), keytype="ENSEMBL")


# remove unannotated rows ----
path_mapping_1 <- path_mapping[!(is.na(path_mapping$PATH)),]
path_mapping_1$PATH <- paste0("hsa",path_mapping_1$PATH) 
rm(path_mapping)

path_mapping_rlog1 <- path_mapping_rlog[!(is.na(path_mapping_rlog$PATH)),]
path_mapping_rlog1$PATH <- paste0("hsa",path_mapping_rlog1$PATH) 
rm(path_mapping_rlog)


# subset this from pathmapping ----
metPath_mapping <- path_mapping_1[which(path_mapping_1$PATH %in% metpathIDs),]
metPath_mapping_rlog <- path_mapping_rlog1[which(path_mapping_rlog1$PATH %in% metpathIDs),]


# convert to one gene per line -----
if (nrow(metPath_mapping)>0) {
  # since many to many mapping is there between genes and pathways,
  # make mapping gene centric. ie, one row for one gene with pathways comma separated
  t = split(metPath_mapping,f=metPath_mapping$ENTREZID)
  t[[1]]
  
  # add pathway name and category
  s <- lapply(t,function(x){sapply(x$PATH,grep,KEGGmet$pathways)})
  s[[1]]
  
  r <- lapply(s,function(x){KEGGmet[c(x),c("pathways","name","pathwayclass")]})
  
  # collapse pathway id, name and category into one vector each
  pathway_vectors <- lapply(t,function(x){paste(x$PATH,collapse=",")})
  
  pathwayname_vectors <- lapply(r,function(x){paste(x$name,collapse=";")})
  pathwayClass_vectors <- lapply(r,function(x){paste(x$pathwayclass,collapse=";")})
  pathwaycheck_vectors <- lapply(r,function(x){paste(x$pathways,collapse=";")})
  
  # create dataframe for pathway annotation
  metPath_mapping_1 <- data.frame(genes = names(pathway_vectors), 
                                  metPaths = as.character(unlist(pathway_vectors)),
                                  pathname = unlist(pathwayname_vectors),
                                  pathwayclass = unlist(pathwayClass_vectors))
  
  metPath_mapping_1$SYMBOL <- metPath_mapping$SYMBOL[match(metPath_mapping_1$genes,
                                                           metPath_mapping$ENTREZID)]
  metPath_mapping_1$ENSEMBL <- metPath_mapping$ENSEMBL[match(metPath_mapping_1$genes,
                                                            metPath_mapping$ENTREZID)]
  # save
  write.table(metPath_mapping_1,
              file.path(annot_dir,"GSE118829_deseqcountsmat_KEGGpathway_mapping.tsv"),
              sep="\t", row.names = F,quote = F)
  
}  

if (nrow(metPath_mapping_rlog)>0) {
  # since many to many mapping is there between genes and pathways,
  # make mapping gene centric. ie, one row for one gene with pathways comma separated
  t = split(metPath_mapping_rlog,f=metPath_mapping_rlog$ENTREZID)
  t[[1]]
  
  # add pathway name and category
  s <- lapply(t,function(x){sapply(x$PATH,grep,KEGGmet$pathways)})
  s[[1]]
  
  r <- lapply(s,function(x){KEGGmet[c(x),c("pathways","name","pathwayclass")]})
  
  # collapse pathway id, name and category into one vector each
  pathway_vectors <- lapply(t,function(x){paste(x$PATH,collapse=",")})
  
  pathwayname_vectors <- lapply(r,function(x){paste(x$name,collapse=";")})
  pathwayClass_vectors <- lapply(r,function(x){paste(x$pathwayclass,collapse=";")})
  pathwaycheck_vectors <- lapply(r,function(x){paste(x$pathways,collapse=";")})
  
  # create dataframe for pathway annotation
  metPath_mapping_1 <- data.frame(genes = names(pathway_vectors), 
                                  metPaths = as.character(unlist(pathway_vectors)),
                                  pathname = unlist(pathwayname_vectors),
                                  pathwayclass = unlist(pathwayClass_vectors))
  
  metPath_mapping_1$SYMBOL <- metPath_mapping_rlog$SYMBOL[match(metPath_mapping_1$genes,
                                                                metPath_mapping_rlog$ENTREZID)]
  metPath_mapping_1$ENSEMBL <- metPath_mapping_rlog$ENSEMBL[match(metPath_mapping_1$genes,
                                                                  metPath_mapping_rlog$ENTREZID)]
  # save
  write.table(metPath_mapping_1,
              file.path(annot_dir,"GSE118829_rlogcountsmat_KEGGpathway_mapping.tsv"),
              sep="\t", row.names = F,quote = F)
  
}  
