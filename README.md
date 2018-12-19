# Epigenetic Reprogramming in Human HSC Aging
Characterization of the aging human HSC epigenome

## Overview
We performed a comprehensive characterization of epigenomic changes during normal human hematopoietic stem cell (HSC) aging. The hematopoietic stem cell enriched population (Lineage-, CD34+, CD38-; HSCe) was isolated from young (18-30 yo) and aged (65-75) donors.  Using multiple (4-10) biological replicates per age group, we profiled H3K4me1, H3K4me3, H3K27me3 and H3K27ac using a low-input ChIP-seq protocol, mC using ERRBS, and the HSCe transcriptome using bulk as well as sc-RNA-seq. All data used for the below analyses have been deposited in GEO (accession number GSE104408) and sc-RNAseq data is available using accession number GSE104379. In addition, we also utilized previously published datasets (see 'Published data sets used for analysis' below).


## Analysis of histone modifications (ChIP-seq) in aging
ChIP-seq for H3K4me1, H3K4me3, H3K27me3, and H3K27ac in human HSCe was performed using 4-7 biological replicates per age group, per mark. For processing and alignment of ChIP-seq data, see 'Preprocessing_ChIPseq.txt'. For peak calling see 'ChIPseq_PeakCalling.txt'. For generating bigwig files that contain the fold-enrichment (IP/Input) (used to make tracks for UCSC and other downstream analyses), see 'FoldEnrichment_bw.py'. Functional annotation and annotation to the nearest TSS was performed using chipenrich (see 'ChIPenrich_Annotation.txt'). Briefly:
  - Illumina adapters were removed using cutadapt (v1.12)
  - fastQC was performed
  - Samples were aligned to hg19 using bowtie2 (v2.0.5)
  - Unique mapped reads were extracted using samtools (v1.2)
  - Reads were sorted and indexed using samtools (v1.3)
  - Replicates for each age group were pooled using samtools merge (v1.3)
  - Peak calling and differential peak calling were performed using macs2 (v2.0.10.20131216)
  - Read-normalized fold-enrichment(IP/Input) bigwig files were created for each replicate and pooled replicates using deepTools2 (see 'FoldEnrichment_bw.py').
  - Fold-change(Young/Aged) was calcualted for differential peaks using deepTools2 multiBigwigSummary and the R-package gtools
  - ChIP-seq peaks, enhancers, and bivalent promoters were annotated using the R-package chipenrich (v2.2.0)
  - The differential peak calling method was verified using permutation analysis to establish an FDR for MACS2 bdgdiff peak calling (see 'PermutationAnalysis.txt')

### k-means clustering of age-associated histone differences
k-means clustering was performed using all peaks that had significant (LLR > 3, absolute fold-change >1.5) H3K4me1, H3K4me3, H3K27ac, or H3K27me3 changes with age (n=37,058 peaks) and R version 3.5.1 (see 'kmClusteringOfHistoneChanges').
  - For the input matrix, scores for each histone modification, for each age group, were calculated using deeptools2 multiBigwigSummary with bigwig files that were read normalized and contained the fold-enrichment of the Pooled IP/Pooled Input (see 'FoldEnrichment_bw.py').
  - The fold-change of aged/young was calculated for each modification for each peak using gtools. Peaks with low-counts for a given histone modification (score aged <3 and score young <3) were filtered
  - K-means clustering of the fold-change values was then performed using base R
  - Genomic annotation of each peak was performed using the R-package Genomation (v1.2.2), with promoter regions being defined as +/- 3000 bp of the TSS. Genomation was also used to annotate peaks to bivalent promoters and active enhancers identified in young HSCe.
  - A heatmap of the clusters along with their annotation package was plotted using the R-package ComplexHeatmap
  
## Analysis of mC (ERRBS) in aging 
DNA methylation was profiled using ERRBS in young (n=7) and aged HSCe (n=5) (see 'ERRBS_Analysis.txt'). In order to plot a heatmap of the differentially methylated regions (DMR), the percent methylation of each DMR was calculated (see 'PercMethofDMR.py') for young and aged HSCe. 

- Data was aligned using the amp-errbs pipeline (https://code.google.com/archive/p/amp-errbs/) with:
  - Bismark (v0.4.1)
  - Cutadapt (v1.12)
  - Bowtie (v0.12.8)
 - Differential methylation analysis was performed using the R-packages methylKit (methylKit_0.9.4) and edmr (edmr_0.6.3.1)

 
## Epigenetic changes predisposing to AML
In order to determine if age-associated epigenetic changes predispose for myeloid malignancies, we utilized published data sets to identify clusters of differentially methylated regions (DMR) or H3K27ac peaks that are altered with normal HSCe aging and in acute myeloid leukemia (AML).

### DNA methylation
The percent methylation for each age-associated DMR (n=529) was calculated for each young and aged HSCe sample, and AML patient (n=119, see 'PercMethofDMR.py'). k-means clustering was then performed and boxplots generated for each cluster.
  - methCall data for the AML patients is from Glass et al., 2017 (GSE98350)
  - The optimal number of clusters was determined using the gap statistic method with the R-package factoExtra
   - k-means clustering of hypomethylated and hypermethylated DMR was performed using R
  
### H3K27ac
The log2(H3K27ac/Input) enrichment was calculated for AML blasts (n=52) and young and aged HSCe at regions with reduced H3K27ac with age. k-means clustering was then performed and boxplots generated for each cluster.
  - H3K27ac ChIP-seq data is from McKeown et al., 2017 (SRP103200)
  - Bigwig files containing the log2(H3K27ac/Input) were generated using deepTools2 bamCompare (version 2.5.6) (see 'bamCompare.py')
  - Enrichment for each differential peak was calculated using deepTools2 multiBigwigSummary (version 2.5.6)
  - The optimal number of clusters was determined using the R-package factoExtra and the gap statistic method. 
  
  
## Transcriptional changes with aging (RNA-seq)
RNA-seq was performed on young and aged HSCe (n=10 per age group). These libraries were 2nd-stranded paired-end 50bp and contained ERCC spike-ins. Additionally, RNA-seq was performed on CD34+ with shLMNA knockdown (n=8 biological replicates). These libraries were paired-end 1st-stranded, 75 bp and contain ERCC spike-ins. RNA-seq was also performed on CD34+ with sgKL6 knockout (n=6 replicates); these libraries were paired-end 1st-stranded, 75 bp and contained ERCC spike-ins.


### Data alignment and processing
For details on processing the HSCe libraries, see 'RNAseqPipeline_HSCe.py' for details on processisng the CD34+ libraries with shLMNA KD or sgKLF6 knockout, see 'RNAseqPipeline_LMNA_and_KLF6'. For all libraries the below steps were followed:
- Reads were trimmed and adapters removed using cutadapt (v1.12)
- Fastqc was performed 
- Reads were aligned to hg19 using the STAR aligner (v2.5.2b)
- Aligned files were sorted using Samtools (v1.3.1)
- Counts were generated using QoRTs (v1.0.7) and a .gtf file without rRNA, mitoRNA and tRNAs 


### Differential gene expression analysis
See 'DifferentialGeneExpressionAnalysis.txt'
  - Differential gene expression analysis was performed using DESeq2 (v1.10.1)
  - For the comparison of Aged vs. Young HSCe, a design that controlled for donor sex and batch while testing for differences between age group was used
  - For analysis of CD34+ cells with shLMNA knockdown, a design that controlled for donor while testing for differences in empty vector versus shLMNA was used.
  - For analysis of CD34+ cells with sgKLF6 KO, a design that controlled for donor while testing for differences in sgNT (non-targeting) vector versus sgKLF6 was used.
  - The Wald method was used for all analyses
  - The regularized log counts generated by DESeq2 were also used for downstream analysis
  
  
## Published data sets used for analysis:
- AML H3K27ac ChIPseq data: McKeown et al., 2017; SRA accession number SRP103200
- AML ERRBS mC data: Glass et al., 2017; GEO accession number GSE98350
- Gene annotations for Gencode v19 (Ensembl 74): Gencode, http://www.gencodegenes.org/releases/19.html
- Genome sequence for GRCh37.p13: Gencode, http://www.gencodegenes.org/releases/19.html
- ERCC spike in .fa and .gtf: ThermoFisher, https://www.thermofisher.com/order/catalog/product/4456739


 

