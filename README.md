# MendelianRNA-seq-DB

![alt text](https://github.com/dennis-kao/MendelianRNA-seq-DB/blob/master/SpliceJunctionSchema.png)

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq-DB is a tool to discover novel junctions in a list of bam files. SpliceJunctionDiscovery.py calls  samtools to report the presence of introns in a list of regions of interest and outputs their read counts. SpliceJunctionSummary.py reads the output and stores this information into a database. The database can then be querried in a number of ways to isolate or filter junctions in hopes of finding an aberrant splicing event which may then be used to diagnose a patient with a rare disease.

SpliceJunctionDiscovery.py usually takes the longest to execute because it calls upon samtools based on the number of samples * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the script. This number should be equal to or less than the number of threads on your system.

## Steps

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
	
	```GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS```

	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Make and/or navigate to a directory containing .bam and corresponding .bai files. Run the [novel splice junction discovery script](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/Analysis/rnaseq.novel_splice_junction_discovery.pbs).

	```qsub /home/MendelianRNA-seq-DB/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=kidney.glomerular.genes.list,bamList=bamlist.list```

	Mandatory parameters:
	1. transcript_file, path to file produced in step 2
	2. bam_list, a text file containing the names of all bam files used in the analysis, each on a seperate line. For example:

		```
		G65693.GTEX.8TY6-5R5T.2.bam
		G55612.GTEX.7G54-3GS8.1.bam
		G09321.GTEX.0EYJ-9E12.3.bam
		PATIENT.bam
		```
		
		All control bams should have the phrase "GTEX" in their file name. All other bams are considered to be patients.

	3. processes, the number of worker processes running in the background calling samtools. This the slowest step in the program. This number should be equal to or less than the number of cores on your machine. 
	
		For torque users: This number should also be equal to or less than the number specified for ppn in rnaseq.novel_splice_junction_discovery.pbs:

		
			#PBS -l walltime=10:00:00,nodes=1:ppn=10

4. Run SpliceJunctionDiscovery.py with --addGencode or --addGencodeWithFlanks to initally populate the database with gencode junctions. 

	```python3 SpliceJunctionDiscovery.py --addGencodeWithFlanks -transcript_model=gencode.comprehensive.splice.junctions.txt```
	
5. Run SpliceJunctionDiscovery.py with --addBAM to populate the database with junctions and read counts from your samples.

	```python3 SpliceJunctionDiscovery.py --addBAM -gene_list=kidney.glomerular.genes.list -processes=10 -bamlist=bamlist.list```

## Output

SpliceJunctionDiscovery.py generates a folder for each bam. Within this folder are text files containing summarized read counts for junctions pertaining to a specific gene.

By default the database is named SpliceJunction.db. There are 4 tables:

	1. SAMPLE_REF, a list of samples and their type (0 = GTEX or control, 1 = patient)
	2. JUNCTION_REF, a list of junctions and their frequency of appearances in samples
	3. JUNCTION_COUNTS, read counts of junctions specific to a sample
	4. GENE_REF, an annotation of junctions with genes, a single junction can map to multiple genes
	
A Python script to query the database will be added later on.

## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.

A gene can encompass partial or whole regions of other genes. This edge case has been accounted for in SpliceJunctionSummary.py and the mapping of a single junction to multiple genes has been done with the table GENE_REF.
