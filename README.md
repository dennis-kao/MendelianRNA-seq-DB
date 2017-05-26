# MendelianRNA-seq

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

## Steps

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
	
	```GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS```

	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Run the [novel splice junction discovery script](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/Analysis/rnaseq.novel_splice_junction_discovery.pbs)

	```qsub ~/tools/MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v gene_list=kidney.glomerular.genes.list,bam_list=bam.list,sample=sampleName```

	Mandatory parameters:
	1. gene_list, path to file produced in step 2
	2. bam_list, a text file containing the names of all bam files used in the analysis, each on a seperate line. For example:

		```
		control1.bam
		control2.bam
		control3.bam
		findNovel.bam
		```
		NOTE: bam files and bai files should be in the current working directory
	
	3. sample, the name of the bam file you want to find novel junctions in, without the ".bam" extension. For example, if your file name is "findNovel.bam", then write "sample=findNovel"

	Optional parameters:
	1. minread, the minimum number of reads a junction needs to have (default=10)
	2. threshold, the minimum normalized read count a junction needs to have (default=0.5)
	3. [transcript_model](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt), the absolute path to a text file containing a list of known canonical splice junctions (default=/home/dennis.kao/tools/MendelianRNA-seq/gencode.comprehensive.splice.junctions.txt)
	
## Output

The script outputs a single file with the name:

threshold**0.XX**\_novel\_**sampleName**\_norm\_**All.kidney.glomerular.genes.list**.splicing.txt

where 0.XX is the threshold value, sampName is the sample you want to discover novel junctions in, and All.kidney.glomerular.genes.list is the name of the input file.

The file contains text information in the format:

```GENE	CHROM:START-STOP	READ_COUNT SITES_ANNOTATED	NORM_READ_COUNT```

Here is a sample output:

```
PLEKHG5	1:6556900-6557380	11	One annotated	1.222
CNTN1	12:41303907-41312441	14	Both annotated	1.0
CNTN1	12:41330707-41331372	50	Both annotated	1.0
CNTN1	12:41327680-41330583	39	Both annotated	1.0
CNTN1	12:41302295-41303875	46	Both annotated	1.0
TRPV4	12:110222242-110224515	11	Both annotated	1.0
CFL2	14:35181674-35182098	10	Neither annotated	-
PMM2	16:8909579-8930060	10	Neither annotated	-
PMM2	16:8909579-8917904	17	Neither annotated	-
MYH2	17:10441098-10546261	10	Neither annotated	-
MYH2	17:10426966-10533715	37	Neither annotated	-
MYH2	17:10426947-10533478	42	Neither annotated	-
MYH2	17:10426959-10533502	57	Neither annotated	-
MYH2	17:10426961-10533718	10	Neither annotated	-
MYH2	17:10426961-10533480	13	Neither annotated	-
MYH2	17:10426949-10533504	98	Neither annotated	-
CAV3	3:8819222-8883028	27	Both annotated	1.0
TPM2	9:35683226-35684252	20	Neither annotated	-
```
## Footnotes

The transcript_model file _gencode.comprehensive.splice.junctions.txt_ is based off of gencode v19.
