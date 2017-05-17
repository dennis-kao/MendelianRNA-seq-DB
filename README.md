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

4. Filter junction covered with >10 reads and found in 1 sample.
```cat All.kidney.glomerular.genes.list.splicing.txt | awk '{if (($6>10) && ($7==1)) print $0}' > novel_junctions.txt```

5. Add normalized read count values to the file produced from step 3.
```qsub ~/tools/MendelianRNA-seq/Analysis/rnaseq.normalize_splice_junction_values.pbs -v splice_file=All.kidney.glomerular.genes.list.splicing.txt,transcript_model=gencode.comprehensive.splice.junctions.txt,action=--normalize,$outputFileName=output.txt```

