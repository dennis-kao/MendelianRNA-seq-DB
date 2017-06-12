# MendelianRNA-seq

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq helps to discover novel splice sites in a sample given a list of bam files. SpliceJunctionDiscovery.py calls upon samtools to report the presence of introns given a list of regions of interest and summarizes these results for read count. NormalizeSpliceJunctionValues.py normalizes the read count of each site based on read support from nearby junctions. FilterSpliceJunctions.py then filters out any site that is of low quality and/or is present in samples other than the one being studied.

SpliceJunctionDiscovery.py usually takes the longest to execute because it calls upon samtools based on the number of samples * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the standalone script. This number should be equal to or less than the number of cores on your system.

## Steps

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
	
	```GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS```

	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Make and/or navigate to a directory containing all your .bam and corresponding .bai files. Run the [novel splice junction discovery script](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/Analysis/rnaseq.novel_splice_junction_discovery.pbs). NOTE: there should not be any .txt files present beforehand in order for SpliceJunctionDiscovery.py to run correctly.

	```qsub MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=kidney.glomerular.genes.list,bamList=bamlist.list,sample=sampleName```

	Mandatory parameters:
	1. transcriptFile, path to file produced in step 2
	2. sample, the name of the bam file you want to find novel junctions in, without the ".bam" extension. For example, if your file name is "findNovel.bam", then write "sample=findNovel"
	
	Optional parameters:
	1. bamList, a text file containing the names of all bam files used in the analysis, each on a seperate line. For example:

		```
		control1.bam
		control2.bam
		control3.bam
		findNovel.bam
		```
		
	2. minread, the minimum number of reads a junction needs to have (default=10)
	3. threshold, the minimum normalized read count a junction needs to have (default=0.5)
	4. [transcript_model](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt), the absolute path to a text file containing a list of known canonical splice junctions (default=/home/dennis.kao/tools/MendelianRNA-seq/gencode.comprehensive.splice.junctions.txt). This file is used in NormalizeSpliceJunctionValues.py.
	5. processes, the number of worker processes running in the background calling samtools. This the slowest step in the program. This number should be equal to or less than the number of cores on your machine. 
	
		For torque users: This number should also be equal to or less than the number specified for ppn in rnaseq.novel_splice_junction_discovery.pbs:

		
		#PBS -l walltime=10:00:00,nodes=1:ppn=10
	
## Output

The scripts output 2 files:

All.**kidney.glomerular.genes**.list, (where kidney.glomerular.genes is the name of your transcriptFile) 

which contains all splice site information pertaining to all samples,

and

threshold**X.XX**\_novel\_**sampleName**\_norm\_**All.kidney.glomerular.genes.list**.splicing.txt, (where X.XX is the threshold value, sampleName is the sample you want to discover novel junctions in, and All.kidney.glomerular.genes.list is the name of the input file) 

which contains splice sites only seen in sampleName that have a read count > minRead and a normalized read count > threshold.

The "threshold" file contains text information in the format:

```GENE	GENE-TYPE	CHROM:START-STOP	READ-COUNT	SAMPLES-SEEN	READ-COUNT:SAMPLE	SITES_ANNOTATED	NORM-READ-COUNT:SAMPLE```

Here is a sample output:

```
UNC13C	NEXON	15:54586263-54590009	22	1	22:SAMPLE       Neither annotated	-
UNC13C	NEXON	15:54685380-54707180	14	1	14:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54614294-54624241	19	1	19:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54914618-54915993	16	1	16:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54841890-54847630	19	1	19:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54527307-54528628	10	1	10:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54799393-54803951	14	1	14:SAMPLE	Neither annotated	-
UNC13C	NEXON	15:54624310-54625965	16	1	16:SAMPLE	Neither annotated	-
UROS	NEXON	10:127506885-127511598	14	1	14:SAMPLE	Neither annotated	-
UROS	NEXON	10:127505095-127506791	15	1	15:SAMPLE	Neither annotated	-
VSNL1	NEXON	2:17830893-17836464	12	1	12:SAMPLE	Neither annotated	-
VSNL1	NEXON	2:17773504-17830677	10	1	10:SAMPLE	Neither annotated	-
VSNL1	NEXON	2:17722186-17773337	10	1	10:SAMPLE	Neither annotated	-
VWC2L	NEXON	2:215338494-215342820	19	1	19:SAMPLE	Neither annotated	-
WDPCP	NEXON	2:64040842-64054756	14	1	14:SAMPLE	Neither annotated	-
YIPF7	NEXON	4:44631565-44637938	12	1	12:SAMPLE	Neither annotated	-
ZNF331	NEXON	19:54041746-54042470	40	1	40:SAMPLE	One annotated	40.0:SAMPLE
ZNF404	NEXON	19:44384289-44388108	12	1	12:SAMPLE	Neither annotated	-
```
## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.
