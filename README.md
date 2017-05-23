#MendelianRNA-seq

Modification of Beryl Cummings scripts for damaging splicing events discovery in RNA-seq

1. Run bcbio rna-seq pipeline to get bam files.

2. Create a list of genes of interest (muscular or kidney), in the format
GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS
use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

Some ready list are in data folder.

```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Run junction discovery script
```qsub ~/tools/MendelianRNA-seq/Analysis/rnaseq.splice_junction_discovery.pbs -v gene_list=kidney.glomerular.genes.list,bam_list=bam.list```

4. Run novel junction discovery script
```input=/path/to/All.kidney.glomerrlar.genes.list.splicing.txt sample=sampleName /MendelianRNA-seq/Analysis/NormalizeAndDiscoverNovelSpliceJunctions.sh```

Additional parameters can be specified. Read the .sh file.
