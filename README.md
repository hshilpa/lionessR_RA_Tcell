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

