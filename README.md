# MendelianRNA-seq-DB

![alt text](./SpliceJunctionSchema.png)

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq-DB is a tool to discover splice sites in a collection of BAM files. It is a rewrite of [MendelianRNA-seq](https://github.com/dennis-kao/MendelianRNA-seq).

## Pipeline

SpliceJunctionDiscovery.py calls samtools to report the presence of introns from a list of regions of interest and outputs their read counts to text files. AddJunctionsToDatabase.py reads this output, performs gencode annotations and normalization, and stores the information into a database. FilterSpliceJunctions.py contains some pre-defined queries which can be used to filter junctions in hopes of finding an aberrant splicing event causative for disease.

SpliceJunctionDiscovery.py usually takes the longest to execute because it calls upon samtools based on the number of samples * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the script. AddJunctionsToDatabase.py is much faster and likely takes minutes to an hour for sample sizes less than 100. Querrying the database using FilterSpliceJunctions is probably the fastest step.

## Steps

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
	
	```GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS```

	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Make and/or navigate to a directory containing .bam and corresponding .bai files. Run the [novel splice junction discovery script](Analysis/rnaseq.novel_splice_junction_discovery.pbs).

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
		
		All control BAMs should have the phrase "GTEX" in their file name. BAMs without this condition are considered to be patients.

	3. processes, the number of worker processes running in the background calling samtools. This the slowest step in the program. This number should be equal to or less than the number of cores on your machine. 
	
		For torque users: This number should also be equal to or less than the number specified for ppn in [rnaseq.novel_splice_junction_discovery.pbs](Analysis/rnaseq.splice_junction_summary.pbs):

		
			#PBS -l walltime=10:00:00,nodes=1:ppn=10

4. Run SpliceJunctionDiscovery.py with --addGencode or --addGencodeWithFlanks to initally populate the database with gencode junctions. 

	```python3 SpliceJunctionDiscovery.py --addGencodeWithFlanks -transcript_model=gencode.comprehensive.splice.junctions.txt```
	
5. Run SpliceJunctionDiscovery.py with --addBAM to populate the database with junctions and read counts from your samples.

	```python3 SpliceJunctionDiscovery.py --addBAM -gene_list=kidney.glomerular.genes.list -processes=4 -bamlist=bamlist.list```

## Output

SpliceJunctionDiscovery.py generates a folder for each bam. Within this folder are text files containing summarized read counts for junctions pertaining to a specific gene.

By default the database is named SpliceJunction.db. There are 4 tables:

	1. SAMPLE_REF, a list of samples and their type (0 = GTEX or control, 1 = patient)
	2. JUNCTION_REF, a list of junctions and their frequency of appearances in samples
	3. JUNCTION_COUNTS, read counts of junctions in a sample
	4. GENE_REF, an annotation of junctions with genes, a single junction can map to multiple genes
	
Documentation on how to use FilterSpliceJunctions.py will be added later.

## Differences between MendelianRNA-seq-DB and Beryl Cumming's original MendelianRNA-seq

- SpliceJunctionDiscovery has been rewritten in Python and parallelized - decreasing processing time by a factor proprotional to the number of worker processes
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools
- All information produced by SpliceJunctionDiscovery is stored in a database instead of text files. This allows the user to utilize previously computed results instead of having to run the entire pipeline again when a new sample needs to be analyzed.
- The database has some new fields that can be used to filter junctions: 
	```
	n_patients_seen
	n_gtex_seen
	total_read_count
	total_patient_read_count
	total_gtex_read_count
	```
- Junction annotation now discriminates between START and STOP instead of 'ONE'. In addition, there is a new annotation, called 'EXON_SKIP' which denotes the event of exon skipping. This is done by checking the reference transcript_model to see if the start and stop positions belong to different junctions.
- Normalization of annotated junctions now considers read counts from all junctions which have at least one annotated junction as the denominator whereas before only "BOTH" annotated junctions were used
- From the gencode file, multipe junctions are generated to increase the definition for what is considered to be "annotated". The start and stop position of each junction both have a +/- 1 tolerance. The different combinations of these values (i.e. junction = start + 1, stop + 1) can be been in the function storeTranscriptModelJunctions() of [AddJunctionsToDatabase.py](Analysis/AddJunctionsToDatabase.py). This likely introduces a larger number of false positives for reported EXON_SKIP events. There is an option to not use flanking at all (--addGencode). 

## Citations

[Improving genetic diagnosis in Mendelian disease with transcriptome sequencing](http://stm.sciencemag.org/content/9/386/eaal5209)

Beryl Cumming's original scripts: [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq)

## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.

A gene can encompass partial or whole regions of other genes. This edge case has been accounted for in AddJunctionsToDatabase.py in two ways: 

	1. The mapping of a single junction to multiple genes has been done with the table GENE_REF
	2. If the script encounters the same junction in a sample more than once, it will utilize the result with the highest read count for read count and normalized read count and will discard the other.
