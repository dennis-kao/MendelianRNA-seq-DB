MendelianRNA-seq

Modification of Beryl Cummings scripts for damaging splicing events discovery in RNA-seq

1. Run bcbio rna-seq pipeline to get bam files.

2. Create a list of genes of interest (muscular or kidney), in the format
GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS
use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.