# ENCODE Splicing Analysis
[![DOI](https://zenodo.org/badge/611118810.svg)](https://zenodo.org/badge/latestdoi/611118810)

This repository contains the code necessary to execute the splicing noise analysis on the ENCODE shRNA knockdown samples. Using the metadata extracted from the [Encode Metadata Extraction](https://github.com/guillermo1996/ENCODE_Metadata_Extraction) repository and the script provided in [Splicing_Analysis_Script.R](https://github.com/guillermo1996/ENCODE_Splicing_Analysis/blob/main/Splicing_Analysis_Script.R), we measured the mis-splicing ratio for each annotated intron found across the different clusters.

Other than the scripts, two reports detailing the results are also provided:

* [Splicing Analysis results](https://guillermo1996.github.io/ENCODE_Splicing_Analysis/RMarkdown/Splicing_Analysis_Results.html): different metrics compared across all shRNA knockdown target genes are reported. Included are the comparison of the mis-splicing ratios between cases and control samples in order to measure a significant difference in the medians, as well as other studies regarding the percentage of mis-spliced annotated introns, the number of unique novel junctions and the percentage of reads that they represent.

* [Splicing Noise specific results](https://guillermo1996.github.io/ENCODE_Splicing_Analysis/RMarkdown/Splicing_Noise_Specific_Results.html): focusing on some of the RBP/NMD target genes, we studied the distributions of the distances between the novel splice site and the reference splice site. Also studied are the distribution of the [MaxEntScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) score results and the distribution of the mis-splicing ratios. In all situations, we focused on comparing different splice sites (acceptor site vs. donor site) and different clusters (cases vs. controls).

This repository is based on the [splicing-accuracy-manuscript](https://github.com/SoniaRuiz/splicing-accuracy-manuscript) repository from Sonia García Ruiz, and was developed in order to provide additional information regarding the effects of an RNA binding protein knockdown on the splicing accuracy.
