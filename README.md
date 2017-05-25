# MendelianRNA-seq

#### Modification of Beryl Cummings scripts for damaging splicing events discovery in RNA-seq

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:

	GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS
	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Run the novel splice junction discovery script on the file produced in step 2 

	```qsub ~/tools/MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v gene_list=kidney.glomerular.genes.list,bam_list=bam.list,sample=sampleName```

Mandatory parameters:
a) gene_list - path to file produced in step 2
b) bam_list - path to a text file containing the names of all bam files used in the analysis, each on a seperate line

e.g.
control1.bam
control2.bam
control3.bam
findNovel.bam

c) sample - the name of the bam file you want to find novel junctions in, without the ".bam" extension. For example, if your file name is "findNovel.bam", then write "sample=findNovel"

Additional parameters can be specified:
a) minread, the minimum number of reads a site needs to have (default=10)
b) threshold, the minimum normalized read count a site needs to have (default=0.5)
	


