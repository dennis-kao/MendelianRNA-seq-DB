# MendelianRNA-seq-DB

![alt text](./SpliceJunctionSchema.png)

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq-DB is a tool to discover junctions in a collection of BAM files. It is a rewrite of the scripts found in the /Analysis folder of [MendelianRNA-seq](https://github.com/dennis-kao/MendelianRNA-seq) to support storing junction information in a database in addition to the parallel processing step.

The main benefit of using a database as opposed to text files is that sample bam files only have to be processed once. Analyzing a new sample does not require the user to re-run the pipeline on controls - instead, the user only needs to query the database. In addition, ram use stays relatively low because the results are stored on the disk and not in a Python dictionary (hash map) like the original scripts.

## Diagnosis methodology

MendelianRNA-seq-DB was initially developed to help researchers discover splice junctions causitive for rare Mendelian diseases. The methodology is as follows:

1. Generate 2 sets of splice junctions from a collection of .bam files. One set is considered to be "healthy" and the other is considered to be "disease"
2. Remove any shared splice junctions from the "disease" set since variants causitive for disease are likely not present in a "healthy" population (keep in mind we are dealing with rare diseases)
3. Remove splice sites from the "disease" set which have a low number of read counts and/or normalized read counts and thus can considered as noise
4. Priortize and analyze variants which pertain to regions in genes related to this disease or are consistent with the patient's phenotype
5. Priortize and analyze variants whose intronic regions share only one splice site with that of a known healthy transcript model<sup>*</sup>

<sup>*</sup>The most probable mutation event in an exon-intron-exon region is that which only alters one exon. This event is denoted as 'START' or 'STOP' in the database. 

## Pipeline details

SpliceJunctionDiscovery.py calls upon samtools to report the presence of introns in a list of regions of interest, summarizes their read counts, and writes this to a text file. AddJunctionsToDatabase.py reads this output, performs gencode annotations and normalization, and stores the information into a database. FilterSpliceJunctions.py contains some pre-defined queries which can be used to filter junctions in hopes of finding an aberrant splicing event causative for disease.

SpliceJunctionDiscovery.py usually takes the longest to execute because it calls upon samtools based on the number of samples * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the script. 

AddJunctionsToDatabase.py is much faster and likely takes minutes to an hour for sample sizes less than 100. Querrying the database using FilterSpliceJunctions is probably the fastest step taking seconds to execute.

## Required files

1. .bam (and .bai) files produced from an RNA-seq pipeline - All control or "healthy" .bams need to have the phrase 'GTEX' in their file name for read count logic to work properly. You need a sufficient number of high quality control BAMs so that you can filter out more splice junctions and discover those that are specific to a diseased sample. The [GTEx project](https://www.gtexportal.org/home/) is a good resource for control BAMs. These BAM files should all be from the same tissue due to tissue specific expression. A way to test for contaminated tissue samples has been described in the study cited below. Note that you can generate .bai files from .bam files using this line: ```parallel  samtools index ::: *.bam```

2. transcript_file - A text file containing a list of genes and their spanning chromosome positions that you want to discover junctions in:
	```
	GENE	ENSG	STRAND	CHROM	START	STOP	GENE_TYPE
	```
	You can use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that, or convert an existing .bed file using this bash line:
	```
	cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> gene.list
	```
	There is an included file which contains [all protein coding regions](all-protein-coding-genes-no-patches.list).
	
3. bamlist.list - A file containing the names of all the bams you want to discover junctions in. The file should quite simply be:
	
	
		G65693.GTEX.8TY6-5R5T.2.bam
		G55612.GTEX.7G54-3GS8.1.bam
		G09321.GTEX.0EYJ-9E12.3.bam
		PATIENT.bam
	
	
	An easy way to generate this file would be to navigate to a directory containing the .bam files you want to use and running this line: ```ls *.bam | grep '' > bamlist.list```

4. transcript_model - A text file containing a list of known canonical splice junctions. These will be used to evaluate a junction's annotation (none, one, both) and it's annotated normalization calculation. You can use your own, or use the [included file](gencode.comprehensive.splice.junctions.txt). This file contains junctions from gencode v19.

## Steps

1. Put bamlist.list, .bam files, .bai files in a new directory. Navigate to it. 
	NOTE: there should not be any .txt files present beforehand in order for SpliceJunctionDiscovery.py to run correctly.

2. For [Torque](http://www.adaptivecomputing.com/products/open-source/torque/) users there is a [PBS file](Analysis/rnaseq.novel_splice_junction_discovery.pbs) containing all the commands you need to run. Just change the "home" directory in the file to match where you placed the MendelianRNA-seq folder and run: 

	```qsub MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=transcript_file,bamList=bamlist.list,processes=10```
	
	For non-Torque users, SpliceJunctionDiscovery can be run from terminal:
	
	```python3 MendelianRNA-seq/Analysis/SpliceJunctionDiscovery.py -transcriptFile=$transcriptFile -bamList=$bamList -processes=$processes```
	
	Parameters:
	1. transcriptFile, path to file #2
	2. bamList, path to file #3
	3. processes, the number of worker processes running in the background calling samtools.This number should be equal to or less than the number of cores on your machine.
	
		For torque users: This number should also be equal to or less than the number specified for ppn in rnaseq.novel_splice_junction_discovery.pbs:
		
			#PBS -l walltime=10:00:00,nodes=1:ppn=10

3. Run AddJunctionsToDatabase.py with --addGencode to initally populate the database with gencode junctions. 

	```python3 AddJunctionsToDatabase.py --addGencode -transcript_model=gencode.comprehensive.splice.junctions.txt```
	
4. Run AddJunctionsToDatabase.py with --addBAM to populate the database with junctions and read counts from your samples.

	```python3 AddJunctionsToDatabase.py --addBAM -transcript_file=all-protein-coding-genes-no-patches.txt -processes=4 -bamlist=bamlist.list -flank=1```

	-flank is a parameter which enables flanking for gencode annotation. If a gencode junction was 1:100-300 and a junction in a sample was 1:99-299, the sample junction would be considered BOTH annotated. This is because both the start and stop positions fall within a +/- 1 range of the gencode junction.

5. At this point your database (SpliceJunction.db) has been populated with junction information from your samples. Now you can use FilterSpliceJunction.py to output junction information.

	To print out splice sites only seen in a "disease" sample and not in any GTEx sample use:

	```python3 FilterSpliceJunctions.py --sample SAMPLE_NAME MIN_READ_COUNT```
	
	It should be noted that 

	If you prefer to use awk and grep tools to filter splice sites and/or avoid writing your own database querries to perform more complex filters then use this to print out all junction information:

	```python3 FilterSpliceJunctions.py --all```

## Output

SpliceJunctionDiscovery.py generates a folder for each bam. Within this folder are text files containing summarized read counts for junctions pertaining to a specific gene. You can delete these folders once all the information has been added to the database.

By default the database is named SpliceJunction.db. There are 4 tables:

	1. SAMPLE_REF, a list of samples and their type (0 = GTEX or control, 1 = patient)
	2. JUNCTION_REF, a list of junctions and their frequency of appearances in samples
	3. JUNCTION_COUNTS, read counts of junctions in a sample
	4. GENE_REF, an annotation of junctions with genes, a single junction can map to multiple genes

Using one of the options of FilterSpliceJunctions.py will produce a text file containing junction information in the following format:

	gene	chromosome:start-stop	annotation	n_gtex_seen	n_patients_seen	total_patient_read_count	total_gtex_read_count	total_read_count	sample:read_count	sample:norm_read_count
	MT-ND1	MT:3540-3610	NONE	0	1	11	0	11	PATIENT.bam:11	PATIENT.bam:NULL
	AC002321.2	GL000201.1:4130-9415	NONE	1	1	32	4	36	PATIENT.bam:32	PATIENT.bam:NULL
	MT-CO1	MT:7276-13822	NONE	1	1	5	1	6	PATIENT.bam:5	PATIENT.bam:NULL
	MT-ATP6	MT:9234-9511	NONE	0	1	6	0	6	PATIENT.bam:6	PATIENT.bam:NULL
	AC002321.2	GL000201.1:9511-14322	START	1	1	70	2	72	PATIENT.bam:70	PATIENT.bam:NULL

## Differences between MendelianRNA-seq-DB and Beryl Cumming's original MendelianRNA-seq

- SpliceJunctionDiscovery has been rewritten in Python and parallelized - decreasing processing time by a factor proprotional to the number of worker processes
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools. As a result of improper parsing using bash tools, junction start and/or stop positions were not reported properly (e.x. 1:100-200*1D30 represents an alignment that should really be 1:100-230 or 1:100-231)
- Junction flanking in NormalizeSpliceJunctionValues.py has been fixed and now works. When flanking junctions were added to the set in the original make_annotated_junction_set(), individual characters in the string were added as opposed to the entire string itself (e.x. 1:100-200 gets added as '1', ':', '0', '2', '-')
- All information produced by SpliceJunctionDiscovery is stored in a database instead of text files. This allows the user to utilize previously computed results instead of having to run the entire pipeline again when a new sample needs to be analyzed.
- The database has some new fields that can be used to filter junctions: 
	```
	n_patients_seen
	n_gtex_seen
	total_read_count
	total_patient_read_count
	total_gtex_read_count
	```
- Junction annotation now discriminates between START and STOP instead of 'ONE'. In addition, there is a new annotation, called 'EXON_SKIP' which denotes the event of exon skipping. This is done by checking to see if the reported 3' and 5' positions from a sample's junction belong to different transcript_model junctions.
- Normalization of annotated junctions now considers read counts from all junctions which have at least one annotated splice site as the denominator whereas before only "BOTH" annotated junctions were used

## Citations

[Improving genetic diagnosis in Mendelian disease with transcriptome sequencing](http://stm.sciencemag.org/content/9/386/eaal5209)

Beryl Cumming's original scripts: [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq)

## Footnotes

A gene can encompass partial or whole regions of other genes. This edge case has been accounted for in AddJunctionsToDatabase.py in two ways: 

	1. The mapping of a single junction to multiple genes has been done with the table GENE_REF
	2. If the script encounters the same junction in a sample more than once, it will utilize the result with the highest read count for read count and normalized read count and will discard the other.

A +/- flanking region is considered when annotating the 5' and 3' positions of sample junctions to increase the number of annotated junctions. This value is specified by the -flank parameter (default 1). There is an option to not use flanking at all (-flank 0).
