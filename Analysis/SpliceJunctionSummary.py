#!/usr/bin/python3
#!/usr/bin/env bash

import os
import sys
import argparse
import multiprocessing
import subprocess
import sqlite3
import re
import gc
import psutil
from datetime import datetime
from sys import getsizeof

databasePath = ""

def connectToDB():
	conn = sqlite3.connect('SpliceJunction.db')
	cur = conn.cursor()

	return conn, cur

def commitAndClose(conn):
	conn.commit()
	conn.close()

def initializeDB():

	conn, cur = connectToDB()

	cur.execute('''PRAGMA journal_mode=WAL;''')

	# if not ('wal' in cur.fetchone()):
	# 	print("Could not set SQLite database to WAL mode. Exiting.")
	# 	exit(1)

	cur.execute('''PRAGMA foreign_keys = ON;''')

	cur.execute('''create table if not exists SAMPLE_REF (
		sample_name varchar(50) primary key, 
		type tinyint not null);''') # type = {0, 1} 
									# GTEX, patient

	cur.execute('''create table if not exists JUNCTION_REF (
		chromosome tinyint not null,
		start unsigned big int not null,
		stop unsigned big int not null,
		gencode_annotation tinyint not null,
		n_patients_seen unsigned big int default 0,
		n_gtex_seen unsigned big int default 0,
		total_read_count big int default 0,
		primary key (chromosome, start, stop));''') # gencode_annotation = {0, 1, 2, 3, 4}
													# none, only start, only stop, both, exon skipping

	cur.execute('''create table if not exists JUNCTION_COUNTS (
		bam_id integer not null,
		junction_id integer not null,
		read_count unsigned big int not null,
		norm_read_count float,
		foreign key(bam_id) references SAMPLE_REF(ROWID),
		foreign key(junction_id) references JUNCTION_REF(ROWID),
		primary key (bam_id, junction_id));''')

	cur.execute('''create table if not exists GENE_REF (
		gene varchar(30) not null,
		junction_id integer not null,
		foreign key(junction_id) references JUNCTION_REF(ROWID)
		primary key (gene, junction_id));''')

	commitAndClose(conn)

def getJunctionID(cur, chrom, start, stop, reads, bam_type):

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping
	# thus, gencode junctions will always have a gencode_annotation value of 3

	# check if start and stop are apart of an existing gencode annotation
	cur.execute('''select ROWID, gencode_annotation from JUNCTION_REF where 
		chromosome is ? and 
		start is ? and 
		stop is ?;''', (chrom, start, stop))
	res = cur.fetchone()

	# print(res)

	if res:
		ROWID, annotation = res
	else: # if no such junction determine annotation of new junction: novel junction, only one annotated or a case of exon skipping?
		
		cur.execute('''select * from JUNCTION_REF where 
			gencode_annotation is 3 and 
			chromosome is ? and 
			start is ?;''', (chrom, start))
		isStartAnnotated = cur.fetchone()

		cur.execute('''select * from JUNCTION_REF where 
			gencode_annotation is 3 and 
			chromosome is ? and 
			stop is ?;''', (chrom, stop))
		isStopAnnotated = cur.fetchone()

		if isStopAnnotated and isStartAnnotated:
			annotation = 4 # exon skipping
		elif isStopAnnotated:
			annotation = 2 # only stop
		elif isStartAnnotated:
			annotation = 1 # only start
		else:
			annotation = 0 # novel junction

		cur.execute('''insert into JUNCTION_REF (
			chromosome, 
			start, 
			stop, 
			gencode_annotation) 
			values (?, ?, ?, ?);''', (chrom, start, stop, annotation))
		
		ROWID = cur.lastrowid

	return ROWID, annotation

def makeSpliceDict(bam, gene_file):

	spliceDict = {}

	with open(gene_file, "r") as gf:
		for line in gf:

			chrom, start, stop, count = line.strip().split()
			uniqueSplice = (chrom, start, stop) 

			spliceDict[uniqueSplice] = int(count)

	return spliceDict

def normalizeReadCount(spliceDict, junction, annotation, annotated_counts):

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping

	chrom, start, stop = junction
	annotation = int(annotation)

	if (annotation == 0) or (annotation == 4): 
		return 'NULL'
	elif annotation == 1:
		key = makeStartString(chrom, start)
	elif annotation == 2:
		key = makeStopString(chrom, stop)
	elif annotation == 3:
		startString = makeStartString(chrom, start)
		stopString = makeStopString(chrom, stop)

		if annotated_counts[startString] > annotated_counts[stopString]:
			key = startString
		else:
			key = stopString

	res = round((float(spliceDict[junction]) / float(annotated_counts[key])), 3)

	return str(res)

def makeStartString(chrom, start):
	return ''.join([chrom,':','START',':',start])

def makeStopString(chrom, stop):
	return ''.join([chrom,':','STOP',':',stop])

def get_annotated_counts(spliceDict):

	count_dict = {}

	for junction in spliceDict:

		chrom, start, stop = junction

		startString = makeStartString(chrom, start)
		stopString = makeStopString(chrom, stop)

		if startString in count_dict:
			if count_dict[startString] < spliceDict[junction]:
				count_dict[startString] = spliceDict[junction]
		else:
			count_dict[startString] = spliceDict[junction]

		if stopString in count_dict:
			if count_dict[stopString] < spliceDict[junction]:
				count_dict[stopString] = spliceDict[junction]
		else:
			count_dict[stopString] = spliceDict[junction]

	return count_dict

# def summarizeGeneFile(poolArguement):

# 	bamList, gene = poolArguement
# 	conn, cur = connectToDB()

# 	for bam in bamList:

# 		bam_id = getBAMID(cur, bam)

# 		sample = bam[:-4]
# 		gene_file = ''.join([os.getcwd(), "/", sample, "/", gene, ".txt"])

# 		if not os.path.isfile(gene_file):
# 			continue

# 		spliceDict = makeSpliceDict(sample, gene_file)
# 		annotated_counts = get_annotated_counts(spliceDict)

# 		for junction in spliceDict:

# 			chrom, start, stop = junction

# 			junction_id, annotation = getJunctionID(cur, chrom, start, stop)

# 			# normalize read counts
# 			try:
# 				norm_read_count = normalizeReadCount(spliceDict, junction, annotation, annotated_counts)
# 			except ZeroDivisionError:
# 				print("Zero division error when normalizing %s:%s-%s in sample %s with annotation %d"%(chrom, start, stop, sample, annotation))
# 				norm_read_count = 'NULL'

# 			# complex insert select join query
# 			try:
# 				# print ("sample: " + sample)
# 				# print ("gene " + gene)
# 				# print ("bam_id, junction_id, read_count, norm_read_count " + ' '.join([str(bam_id), str(junction_id), str(spliceDict[junction]), str(norm_read_count)]))
				
# 				cur.execute('''insert into JUNCTION_COUNTS (
# 					bam_id,
# 					junction_id,
# 					read_count,
# 					norm_read_count) 
# 					values (?, ?, ?, ?);''', (bam_id, junction_id, spliceDict[junction], norm_read_count))
# 			except:
# 				# print ("============CRASH!!!!!!================")
# 				# print ("sample: " + sample)
# 				# print ("gene " + gene)
# 				# print ("bam_id, junction_id, read_count, norm_read_count " + ' '.join([str(bam_id), str(junction_id), str(spliceDict[junction]), str(norm_read_count)]))
# 				# print ("============CRASH!!!!!!================")
# 				# exit(1)
# 				pass
	
# 		del spliceDict, annotated_counts

# 	commitAndClose(conn)

def updateJunctionInformation(junction_id, gene, sample, new_read_count, bam_id, cur):

	# check if gene is related to junction, add it if not
	try:
		cur.execute('''insert into GENE_REF (gene, junction_id) values (?, ?);''', (gene, junction_id))
	except sqlite3.IntegrityError:
		pass

	# check if sample already has the junction in the database
	cur.execute('''select ROWID, read_count, norm_read_count from JUNCTION_COUNTS where junction_id is ? and bam_id is ?;''', (junction_id, bam_id))
	res = cur.fetchone()

	# if it is, check if new_reads > old_reads, update JUNCTION_REF and JUNCTION_COUNTS for the appropriate sample
	if res:
		sample_junction_id, old_read_count, norm_read_count = res

		if int(new_read_count) > int(old_read_count):

	# if not, add it with read counts and normalized read counts, increment n_samples_seen, increment n_times_seen, 
	else:
		

def getBAMID(cur, bam):
	cur.execute('''select ROWID, type from SAMPLE_REF where sample_name is ?;''', (bam, ))
	bam_id, bam_type = cur.fetchone()

	return bam_id, bam_type

def parallel_process_gene_files(num_processes, bam_files, gencode_file, gene_list):

	conn, cur = connectToDB()

	bamList = []
	poolArguements = []

	with open(bam_files, "r") as bf:
		for line in bf:

			bam = line.strip()

			try: # insert sample names into SAMPLE_REF
				if 'GTEX' in bam:
					cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 0);''', (bam, ))
				else: # sample is a patient
					cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 1);''', (bam, ))
			except sqlite3.IntegrityError as e:
				continue # if sample already in DB, don't process it

			bamList.append(bam) # only process samples if insert was successful

	commitAndClose(conn)

	# annotated_junction_set = make_annotated_junction_set(gencode_file)

	# print ("Size of AJS: " + str(getsizeof(annotated_junction_set)))

	with open(gene_list, "r") as gf:
		for line in gf:
			gene = line.strip().split()[0]

			poolArguements.append((bamList, gene))

	# print ("Creating a pool with " + str(num_processes) + " processes")
	# pool = multiprocessing.Pool(processes=int(num_processes),maxtasksperchild=100)
	# print ('pool: ' + str(pool))

	# pool.map(summarizeGeneFile, poolArguements)
	# pool.close()
	# pool.join

	for arg in poolArguements:

		bamList, gene = arg
		conn, cur = connectToDB()

		for bam in bamList:

			bam_id, bam_type = getBAMID(cur, bam)

			sample = bam[:-4]
			gene_file = ''.join([os.getcwd(), "/", sample, "/", gene, ".txt"])

			if not os.path.isfile(gene_file):
				continue

			spliceDict = makeSpliceDict(sample, gene_file)
			annotated_counts = get_annotated_counts(spliceDict)

			for junction in spliceDict:

				chrom, start, stop = junction
				reads = spliceDict[junction]

				junction_id, annotation = getJunctionID(cur, chrom, start, stop, reads, bam_type)

				# normalize read counts
				try:
					norm_read_count = normalizeReadCount(spliceDict, junction, annotation, annotated_counts)
				except ZeroDivisionError:
					print("Zero division error when normalizing %s:%s-%s in sample %s with annotation %d"%(chrom, start, stop, sample, annotation))
					norm_read_count = 'NULL'

				# complex insert select join query
				try:
					# print ("sample: " + sample)
					# print ("gene " + gene)
					# print ("bam_id, junction_id, read_count, norm_read_count " + ' '.join([str(bam_id), str(junction_id), str(spliceDict[junction]), str(norm_read_count)]))
					
					cur.execute('''insert into JUNCTION_COUNTS (
						bam_id,
						junction_id,
						read_count,
						norm_read_count) 
						values (?, ?, ?, ?);''', (bam_id, junction_id, spliceDict[junction], norm_read_count))
				except:
					# print ("============CRASH!!!!!!================")
					# print ("sample: " + sample)
					# print ("gene " + gene)
					# print ("bam_id, junction_id, read_count, norm_read_count " + ' '.join([str(bam_id), str(junction_id), str(spliceDict[junction]), str(norm_read_count)]))
					# print ("============CRASH!!!!!!================")
					# exit(1)
					pass
		
			del spliceDict, annotated_counts
			
		commitAndClose(conn)

def addTranscriptModelJunction(chrom, start, stop, gene, cur):

	try:
		cur.execute('''insert into JUNCTION_REF (chromosome, start, stop, gencode_annotation) values (?, ?, ?, ?);''', (chrom, start, stop, '3')) # 3 stands for both annotated
		junction_id = cur.lastrowid
	except sqlite3.IntegrityError:
		cur.execute('''select ROWID from JUNCTION_REF where chromosome is ? and start is ? and stop is ?;''', (chrom, start, stop))
		junction_id = cur.fetchone()[0]

	# assign the junction to the gene specified
	try:
		cur.execute('''insert into GENE_REF (gene, junction_id) values (?, ?);''', (gene, junction_id))
	except sqlite3.IntegrityError:
		pass # gencode file can have the same junctions coordinates, gene but differing transcript

def storeTranscriptModelJunctions(gencode_file, enableFlanking):

	conn, cur = connectToDB()

	print ('Started adding transcript_model junctions @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	with open(gencode_file, "r") as gf:
		for commitFreq, line in enumerate(gf):

			gene, gene2, plus, chrom, start, stop, gene_type = line.strip().split()

			start = int(start)
			stop = int(stop)

			if enableFlanking:

				# shifts the junction while maintaining the same distance between start and stop
				for offset in range(-1,2):
					startFlank = start + offset
					stopFlank = stop + offset
					addTranscriptModelJunction(chrom, startFlank, stopFlank, gene, cur)

				# generate junctions with the most extreme flanking regions of start and stop
				addTranscriptModelJunction(chrom, (start + 1), (stop - 1), gene, cur)
				addTranscriptModelJunction(chrom, (start - 1), (stop + 1), gene, cur)
				
			else:
				addTranscriptModelJunction(chrom, start, stop, gene, cur)

			if commitFreq % 500 == 0: #	yes this works
				conn.commit()

	commitAndClose(conn)

	print ('Finished adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

if __name__=="__main__":

	print ('SpliceJunctionSummary.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Summarize the read counts of the junctions reported by SpliceJunctionDiscovery.py')
	parser.add_argument('-transcript_model',help="Transcript model of canonical splicing, e.g. gencode v19. Default is set to /home/dennis.kao/tools/MendelianRNA-seq-DB/gencode.comprehensive.splice.junctions.txt",action='store',default = "/home/dennis.kao/largeWork/gene-lists/all-protein-coding-genes-no-patches.list")
	parser.add_argument('-processes',help='Number of worker processes to parse gene files, default=10.',default=10)
	parser.add_argument('-bamlist',help='A text file containing the names of bam files you want to discover splice junctions in each on a seperate line, default=bamlist.list',default='bamlist.list')
	parser.add_argument('-gene_list',help='A text file containing the names of all the genes you want to investigate, default=gene_list.txt',default='gene_list.txt')
	parser.add_argument('-db',help='The name of the database you are storing junction information in, default=SpliceJunction.db',default='SpliceJunction.db')
	mode_arguments = parser.add_mutually_exclusive_group(required=True)
	mode_arguments.add_argument('--addGencode',action='store_true',help='Populate the database with gencode junctions, this step needs to be done once before anything else')
	mode_arguments.add_argument('--addBAM',action='store_true',help='Add junction information from bamfiles found in the file bamlist.list')
	flanking = parser.add_mutually_exclusive_group(required=False)
	flanking.add_argument('--enableFlanking',action='store_true',help='Use with --addBAM. When transcript_model junctions are cross referenced for annotation or normalized, any start or stop position that falls within a +/- 1 nucleotide range is considered to be valid. Adds 5x as many transcript_model junctions and may introduce some false positives for exon skipping.')
	args=parser.parse_args()

	print ('Working in directory ' + str(os.getcwd()))

	#databasePath = os.getcwd() + '/' + args.db
	databasePath = args.db

	print (databasePath)

	initializeDB()

	if args.addGencode: 
		print ('Storing junctions from the transcript model file ' + args.transcript_model)
		storeTranscriptModelJunctions(args.transcript_model, args.enableFlanking)
	elif args.addBAM:
		print ('Storing junctions from bam files found in the file ' + args.bamlist)
		parallel_process_gene_files(args.processes, args.bamlist, args.transcript_model, args.gene_list)

	print ('SpliceJunctionSummary.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
