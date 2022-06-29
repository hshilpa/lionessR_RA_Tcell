#-------------------------------------------------------------------------------
#---------------- TEMPLATE FOR STAR ALIGNEMNT SCRIPT ---------------------------
#-------------------------------------------------------------------------------
## Final Script for Mapping with STAR and hg38 version 2_0 2021/03/22

####!!!!!!!!!!!!! DO NOT MODIFY !!!!!!!!!!!!!!####

##------------------------------------------------------------------------------
## CHANGES
## 2021/03/22
## edited for single end reads
## previous version v_1_3 2021/03/11
## 2021/03/11
## edited path for tmp directory
## reduced threads to 3
## reduce threads further or adjust limitBAMsortRAM and outBAMsortingBinsN
## these settings are hitting the ulimit -n for users.
## discuss with Khan sir if adjustments also don't work.
## also consider sorting with Rsamtools.
## previous version v_1_2 2021/03/11
## 2021/03/10
## removed space from limitBAMsortRAM and outBAMsortingBinsN
## added path for tmp directory
## in STAR command creation, changed to paste0
## previous version v_1_1 2021/03/10
## 2021/03/10
## Earlier version did not allocate enough memory to sort BAM
## changed from default to --outBAMsortingBinsN 200 and 
## --limitBAMsortRAM 48000000000
## previous version 1_0 2021/03/04
##------------------------------------------------------------------------------

## Alignment of reads and counting the aligned reads
## alignment using STAR and GRCh38 release 102 as reference
## index file built with STAR 2.7.7a

### STEPS
# 0. input options for STAR in a vector
# 1. set paths for input files 
# 2. set output directory for countsmatrix and log file
# 3. create vectors of fwd and reverse file names
# 4. create a loop for the following:
#    a. create STAR command string using options and input fasta files
#    b. run STAR command in linux using system(<command string>)
#    c. get BAM file stats using samtools from linux
#    d. delete all objects created in the loop


##----------------------------------------------------------------------------##


# STAR --help
# Usage: 
# STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
# Spliced Transcripts Alignment to a Reference 
# (c) Alexander Dobin, 2009-2020

# STAR version=2.7.7a
# STAR compilation time,server,dir
# =Mon Dec 28 13:38:40 EST 2020 
# vega:/home/dobin/data/STAR/STARcode/STAR.master/source
# For more details see:
#  <https://github.com/alexdobin/STAR>
#  <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
### versions
#  versionGenome           2.7.4a
# string: earliest genome index version compatible with this STAR release. 
# Please do not change this value!

#-------------------------------------------------------------------------------
### Run Parameters
#-------------------------------------------------------------------------------

options = c(runMode = "--runMode alignReads ",
            # string: type of the run.
            runThreadN = "--runThreadN 3 ",
            # int: number of threads to run STAR
            ### Genome Parameters
            genomeDir = "--genomeDir /home/references/STAR_index/hg38_STAR_2.7.7a_ENSEMBL_primary_release102 ",
            #  string: path to the directory where genome files are stored 
            outBAMsortingBinsN = "--outBAMsortingBinsN 200 ",
            limitBAMsortRAM = "--limitBAMsortRAM 48000000000 ",
            ### Output: general
            outFileNamePrefix = "",
            # string: output files name prefix (including full or relative path).
            outTmpDir = "",
            ### Output: SAM and BAM
            outSAMtype = " --outSAMtype BAM SortedByCoordinate ",
            # strings: type of SAM/BAM output
            ### Quantification of Annotations
            quantMode = " --quantMode TranscriptomeSAM GeneCounts ") 
#    string(s): types of quantification requested
#------------------------------------------------------------------------------#

library(Rsamtools)


# set directory to save counts files


# setwd to a directory in your folder
setwd("./STAR_alignment")
# set temp file prefix
options["outTmpDir"] = paste(" --outTmpDir ",getwd(),"/tmp",sep="")


## GTF file for annotation
hg38_GTF = "/home/references/STAR_index/Homo_sapiens.GRCh38.102.gtf"

fastq_dir = "/home/RA_RNAseq/GSE118829/fastq"
outdir = "/home/RA_RNAseq/GSE118829/aligned_STAR"

dir.create(outdir)
#filtered and trimmed files
fwd = list.files(fastq_dir,pattern = ".fastq")



for (i in 1:length(fwd))
  
{
  read1 = fwd[i]
  
  
  # path to fastq infiles
  infile1 = file.path(fastq_dir,read1)
  
  
  # extract samplename from filename
  sample = unlist(strsplit(read1,split = ".fastq"))
  
  # set output file prefix
  options["outFileNamePrefix"] = paste(" --outFileNamePrefix ",file.path(outdir,sample),sep="")
  
  
  # record time
  sink(paste(getwd(),"/STARoutput_",sample,".txt",sep=""),append = TRUE)
  cat("\n------------------------\n")
  print(Sys.time())
  cat("ALIGNMENT AND COUNTING OF ", sample,"\n")
  
  
  # paste together the command
  STAR <- paste0("STAR ",
                 paste(options,collapse=""),
                 " --readFilesIn ",infile1, 
                 " > ",sample,".log")
  
  startTime = Sys.time()
  
  ## Use the following command to run directly on linux
  system(STAR)
  
  endTime1 = Sys.time()
  
  cat("\n")
  cat("Time to run STAR: ")
  t = endTime1 - startTime
  print(t)
  cat("\n")
  rm(t)
  rm(STAR)
  
  # get name of bam file
  outBamfile = paste0(outdir,"/",sample,"Aligned.sortedByCoord.out.bam")
  
  # samtools flagstat to get stats of created bam file
  samtools = paste("samtools flagstat ",
                   outBamfile,
                   " >",
                   paste(getwd(),"/",sample, "_stats.out",sep=""))
  
  system(samtools)
  
  endTime2 = Sys.time()
  
  cat("\n")
  cat("Time to check stats with Samtools: ")
  t = endTime2 - endTime1
  print(t)
  cat("\n")
  rm(t)
  rm(samtools)
  
  
  cat("Time taken for ",sample,": ")
  t = endTime2 - startTime
  print(t)
  cat("\n------------------------\n")
  sink()  
  
  rm(t)
  rm(endTime1)
  rm(endTime2)
  rm(startTime)
  
  rm(outBamfile)
  rm(read1)
  rm(infile1)
  rm(sample)
}

