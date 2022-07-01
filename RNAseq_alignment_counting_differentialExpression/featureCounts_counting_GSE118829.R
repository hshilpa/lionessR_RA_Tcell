#-------------------------------------------------------------------------------
#---------------- FEATURECOUNTS SCRIPT -----------------------------------------
#-------------------------------------------------------------------------------

## Counting the aligned reads
## Counting using featureCounts and GRCh38 release 102 as reference
## alignment done with STAR 2.7.7a

### STEPS
# 0. set working directory as a folder in home directory, set log file
# 1. set paths for input files(sorted BAM) and GTF annotation 
# 2. set output directory for countsmatrix
# 3. create vectors of sorted BAM file names

# 4. 
#    a. count sorted BAM at transcript level
#    b. count transcript BAM at gene level
# 5. cleanup columns names and annotation
# 6. save count files to output directory

##---------------------------- LOAD LIBRARIES ----------------------------------


library(Rsubread)

##---------------------------- SET INPUT DIRS ----------------------------------

setwd(".")

# 1. set paths for input files(sorted BAM) and GTF annotation 
hg38_GTF = "/home/references/STAR_index/Homo_sapiens.GRCh38.102.gtf"
inDIR = "/home/RA_RNAseq/GSE118829/aligned_STAR"


##--------------------------- SET OUTPUT DIR AND SINK --------------------------

# 2. set output directory for countsmatrix
outDir = "./RA_RNAseq/GSE118829/counts"
sink(file.path(outDir,
               paste(gsub(" |-|:","_",Sys.time()),
                     "_featureCounts.log",sep="")))

##------------------------------------------------------------------------------

# 3. create vectors of sorted BAM file names
BAMfiles = list.files(inDIR)[grepl("Aligned.sortedByCoord.out",
                                   list.files(inDIR))]


##------------------------------------------------------------------------------
print(Sys.time())
cat("\n")
cat("------------------------------------\n")
cat("BAM FILES COUNTED: \n")
cat(BAMfiles, sep = "\n")
cat("------------------------------------\n\n")

##------------------------------------------------------------------------------

# 4.
#    a. count sorted BAM at transcript level
#    b. count sorted BAM at gene level


  infile = BAMfiles
  # path to BAM
  outBamfile = file.path(inDIR,BAMfiles)
  # extract samplename from filename
  sample = unlist(strsplit(infile,split = "Aligned.sortedByCoord.out.bam"))
  
  # sink
  cat("\n------------------------------------\n")
  cat("Sample counted: ",sample)
  cat("\n------------------------------------\n")
  
  
  #    a. count sorted BAM at transcript level
  tmp_transcript  <-  featureCounts(outBamfile, annot.ext = hg38_GTF, ignoreDup=FALSE, 
                                   isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,
                                   GTF.featureType="transcript", GTF.attrType="transcript_id",
                                   GTF.attrType.extra=c("gene_id","gene_name"),
                                   useMetaFeatures=FALSE)
  

  #    b. count sorted BAM at gene level
  tmp_gene  <-  featureCounts(outBamfile, annot.ext = hg38_GTF, ignoreDup=FALSE, 
                                    isGTFAnnotationFile=TRUE, isPairedEnd=FALSE,
                                    GTF.featureType="exon", GTF.attrType="gene_id",
                                    GTF.attrType.extra="gene_name",
                                    useMetaFeatures=TRUE)
 

geneCounts_df=as.data.frame(tmp_gene$counts)
transcriptCounts_df=as.data.frame(tmp_transcript$counts)
gene_annot=tmp_gene$annotation
transcript_annot=tmp_transcript$annotation
##------------------------------------------------------------------------------


# 5. cleanup columns names and annotation
colnames(transcriptCounts_df) = unlist(strsplit(colnames(transcriptCounts_df),
                                                split="Aligned.sortedByCoord.out.bam"))
colnames(geneCounts_df) = unlist(strsplit(colnames(geneCounts_df),
                                          split="Aligned.sortedByCoord.out.bam"))

transcriptCounts_df$TranscriptID = rownames(transcriptCounts_df)
transcriptCounts_df = transcriptCounts_df[,c((ncol(transcriptCounts_df)),1:(ncol(transcriptCounts_df)-1))]

geneCounts_df$GeneID = rownames(geneCounts_df)
geneCounts_df = geneCounts_df[,c((ncol(geneCounts_df)),1:(ncol(geneCounts_df)-1))]

transcriptCounts_df = merge(transcriptCounts_df,
                            transcript_annot, 
                            by.x = "TranscriptID", 
                            by.y = "GeneID")

geneCounts_df = merge(geneCounts_df,
                      gene_annot,
                      by.x = "GeneID",
                      by.y = "GeneID")

# sink
cat("\n------------------------------------\n")
cat("Trancripts counted: ", nrow(transcriptCounts_df),"\n")
cat("Genes counted: ", nrow(geneCounts_df))
cat("\n------------------------------------\n")

##------------------------------------------------------------------------------
# 6. save count files to output directory
# add Sys.time() to filename to avoid overwriting older files
# change filename if necessary
write.table(transcriptCounts_df,
            paste(outDir,"/",
                  gsub(" |-|:","_",Sys.time()),
                  "_transcriptCounts_STAR.tsv",sep=""),
            sep="\t", 
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
write.table(geneCounts_df,
            paste(outDir,"/",
                  gsub(" |-|:","_",Sys.time()),
                  "_geneCounts_STAR.tsv",sep=""),
            sep="\t", 
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

print(Sys.time())
cat("\n")
cat("------------------------------------\n")

sink()
##------------------------------------------------------------------------------
