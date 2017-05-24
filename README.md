#MendelianRNA-seq

Modification of Beryl Cummings scripts for damaging splicing events discovery in RNA-seq

1. Run bcbio RNA-seq pipeline to get .bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS
use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

Some ready list are in data folder.

```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Run junction discovery script
```qsub ~/tools/MendelianRNA-seq/Analysis/rnaseq.splice_junction_discovery.pbs -v gene_list=kidney.glomerular.genes.list,bam_list=bam.list```

4. Run novel junction discovery script with the input file produced from step 3. and the name of the bamfile (without the extension) you want to discover novel junctions in
```input=All.kidney.glomerular.genes.list.splicing.txt sample=sampleName /MendelianRNA-seq/Analysis/NormalizeAndDiscoverNovelSpliceJunctions.sh```

Additional parameters can be specified:
a) minread, the minimum number of reads a site needs to have (default=10)
b) threshold, the minimum normalized read count a site needs to have (default=0.5)

Script will produce the following files in the same directory as the input file:
a) norm_All.kidney.glomerular.genes.list.splicing.txt - the input file with a column for normalized read counts
b) novel_sampleName_All.kidney.glomerular.genes.list.splicing.txt - the file from 1) with splice sites only seen in "sampleName" 
c) thresholdX.XX_novel_sampleName_All.kidney.glomerular.genes.list.splicing.txt - the file from 2) with splice sites with a normalized read count greater than X.XX 
	
	


