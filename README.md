# lionessR_RA_Tcell
Code used in the paper titled "Altered transcriptional regulation of glycolysis in circulating CD8+ T  cells of rheumatoid arthritis patients"

The input for this analysis is the rlog transformed counts matrix from RNA seq data.
Transcription factors and glycolysis related genes are extracted from the matrix for the analysis.
Single sample edgeweights of all gene-pairs are calculted using lionessR.
The edgeweights are compared between untreated RA and healthy, and treated RA and untreated RA for each cell type.
Edges with significantly different edgeweights are selected and annotated with GTRD transcription factor - target information.
The selected edges are used to create networks with igraph and centralities are calculated.

## RNA sequencing differential expression
The fastq files were aligned using STAR (v.2.7.7a)
A wrapper script in R was used : STAR_alignment_singleEnd_GSE118829.R

Reads were counted using Rsubread (v2.4.3). The Rscript featureCounts_counting_GSE118829.R was used.

Differential expression was calculated using the DESeq2 package (v.1.36.0). The script used was deseq2_GSE118829.R
The counts were rlog transformed in this script for the lionessR analysis.

## lionessR analysis
Kegg pathway membership of each gene was annoataed using the script 
annotateKEGGmetabolicPaths_GSE118829_script0.R 

Glycolysis and transcription factor genes were subsetted, lionessR analysis and subsequent limma of single sample edgeweights were performed in the script lionessR_GSE118829_script1.R (lionessR v.1.10.0, limma v.3.52.1)

Subsetting significant edges and annotation with GTRD was done in the script annotateWith_GTRD_edges_GSE118829_script2.R

Differentially expressed nodes were annotated in annotateNodes_with_DEGs_script3.R

Centrality measures were calculated in getCentralities_for_eachgraph_script4.R (igraph 1.3.1)

## annot_tables
The table GSE118829_sampleDetails.csv gives the details of the samples used in the analysis

Key enzymes in glycolysis are given in the table glycolysis_keyEnzymes.tsv

The GTRD transcription factors are given in the table 2021_09_1_TFonly_annotation.tsv

Kegg metabolic pathways are given in the table 2021_09_03_metabolic_pathways.tsv

The GTRD TF target annotation is from v19 of GTRD.
The GeneHancer network edges and nodes are given in the files glygenes_TF_edgeAtt.csv and glygenes_TF_noa.csv respectively.
