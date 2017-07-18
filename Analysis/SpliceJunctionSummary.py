#!/usr/bin/python3
#!/usr/bin/env bash

import os
import sys
import argparse
import multiprocessing
import subprocess
import sqlite3
import re
from datetime import datetime

databaseName = ""

def connectToDB():
	conn = sqlite3.connect(databaseName)
	cur = conn.cursor()

	return conn, cur

def commitAndClose(conn):
	conn.commit()
	conn.close()

def initializeDB():

	conn, cur = connectToDB()

	cur.execute('''create table if not exists SAMPLE_REF (
		sample_name varchar(50) primary key, 
		type tinyint not null);''') # type = {0, 1} 
									# GTEX, patient

	cur.execute('''create table if not exists JUNCTION_REF (
		chromosome tinyint not null,
		start unsigned big int not null,
		stop unsigned big int not null,
		gene varchar(20) not null,
		gencode_annotation tinyint not null,
		n_junctions_seen unsigned big int,
		n_times_seen unsigned big int,
		primary key (chromosome, start, stop));''') # gencode_annotation = {0, 1, 2, 3, 4}
													# none, only start, only stop, both, exon skipping

	cur.execute('''create table if not exists JUNCTION_COUNTS (
		sample_id integer not null,
		junction_id integer not null,
		read_count unsigned big int not null,
		norm_read_count float,
		foreign key(sample_id) references SAMPLE_REF(ROWID),
		foreign key(junction_id) references JUNCTION_REF(ROWID));''')

	cur.execute('''PRAGMA journal_mode=WAL;''')

	if not ('wal' in cur.fetchone()):
		print("Could not set SQLite database to WAL mode. Exiting.")
		exit(1)

	cur.execute('''PRAGMA foreign_keys = ON;''')

	commitAndClose(conn)

def getJunctionID(cur, gene, chrom, start, stop):

	# gencode junctions will always have a gencode_annotation value of 3

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping

	# check if start and stop are apart of an existing gencode annotation
	cur.execute('''select ROWID from JUNCTION_REF where 
		gencode_annotation is 3 and 
		chromosome is ? and 
		start is ? and 
		stop is ?;''', (chrom, start, stop))
	res = cur.fetchone()

	if res:
		annotation = 3
		return res, annotation
	else: # if no such junction, add it. Determine the annotation: novel junction, only one annotated or a case of exon skipping?

		# determine gencode_annotation
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

		# add new junction with gencode_annotation
		cur.execute('''insert into JUNCTION_REF (
			chromosome, 
			start, 
			stop, 
			gene, 
			gencode_annotation) 
			values (?, ?, ?, ?);''', (chrom, start, stop, gene, annotation))
		
		return cur.lastrowid, annotation

def makeSpliceDict(bam, gene_file):

	spliceDict = {}

	with open(gene_file, "r") as gf:
		for line in gf:

			chrom, start, stop = line.strip().split()
			uniqueSplice = (chrom, start, stop) 

			if uniqueSplice not in spliceDict:
				spliceDict[uniqueSplice] = 1
			else:
				spliceDict[uniqueSplice] += 1

	return spliceDict

def normalizeReadCount(spliceDict, start, stop, annotation):

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping

	if (annotation == 0) or (annotation == 4): 
		return 'NULL'
	elif annotation == 1:
		
	elif annotation == 2:

def makeAnnotatedCountDict(spliceDict):

	annotatedCountDict = {}

	for junction in spliceDict:

		chrom, start, stop = junction



def summarizeGeneFile(poolArguement):

	bamList, gene, conn = poolArguement
	cur = conn.cursor()

	for sample in bamList:

		gene_file = ''.join([os.getcwd(), "/", sample, "/", gene, ".txt"])

		if not os.path.isfile(gene_file):
			continue

		spliceDict = makeSpliceDict(sample, gene_file)
		annotatedPositions = makeAnnotatedPositions()
		annotatedCountDict = makeAnnotatedCountDict(spliceDict)

		for junction in spliceDict:

			chrom, start, stop = junction

			junction_id, annotation = getJunctionID(cur, gene, chrom, start, stop)

			cur.execute('''select ROWID from SAMPLE_REF where sample_name is ?;''', (sample))
			sample_id = cur.fetchone()

			# normalize read counts
			norm_read_count = normalizeReadCount(spliceDict, start, stop, annotation)

			# complex insert select join query
			cur.execute('''insert into JUNCTION_COUNTS (
				sample_id,
				junction_id,
				read_count,
				norm_read_count) 
				values (?, ?, ?, ?);''', (sample_id, junction_id, str(spliceDict[junction]), norm_read_count))
	
		del spliceDict, annotatedPositions, annotatedCountDict
		
	conn.commit()

def parallel_process_gene_files(num_processes, bam_files, transcriptFile):

	conn = sqlite3.connect(databaseName)

	with open(bam_files, "r") as bf:
		for bam in bf:
			bamList.append(bam)

			# insert sample names into SAMPLE_REF
			if 'GTEX' in bam:
				cur.execute('''insert into SAMPLE_REF (sample_name, type), values (?, 0);''', bam)
			else: # sample is a patient
				cur.execute('''insert into SAMPLE_REF (sample_name, type), values (?, 1);''', bam)

	with open(transcriptFile, "r") as tf:
		for line in tf:
			gene = line.strip().split()[0]

			poolArguements.append((bamList, gene, conn))

	print ("Creating a pool with " + numProcesses + " processes")
	pool = multiprocessing.Pool(int(numProcesses))
	print ('pool: ' + str(pool))

	conn.commit()

	pool.map(summarizeGeneFile, poolArguements)
	pool.close()
	pool.join()

	conn.close()

def storeGencodeJunctions(gencode_file):

	conn = sqlite3.connect(databaseName)
	cur = conn.cursor()

	print ('Started adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	with open(gencode_file, "r") as gf:
		for commitFreq, line in gf:

			chrom, start, stop, gene, transcript, gene_type = line.strip().split()

			cur.execute('''insert into JUNCTION_REF (chromosome, start, stop, gene, gencode_annotation) values (?, ?, ?, ?, ?);''', (chrom, start, stop, gene, '2'))

			if commitFreq % 500 == 0: #	yes this works
				conn.commit()

	commitAndClose(conn)

	print ('Finished adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

if __name__=="__main__":

	print ('SpliceJunctionSummary.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Summarize the read counts of the junctions reported by SpliceJunctionDiscovery.py')
	parser.add_argument('-processes',help='Number of worker processes to parse gene files, default=10.',default=10)
	parser.add_argument('-bamList',help='A text file containing the names of bam files you want to discover splice junctions in each on a seperate line',default='bamlist.list')
	parser.add_argument('-gencode', help='Specify a path to a file containing gencode junctions to be added to the database. This only needs to be run once.',default='')
	parser.add_argument('-db',help='Name of the database file including extension')
	args=parser.parse_args()

	print ('Working in directory' + str(os.getcwd()))

	databaseName = args.db

	initializeDB()

	if gencode: 
		storeGencodeJunctions(args.gencode)

	# parallel_process_gene_files(args.processes, args.bamList, args.transcriptFile)

	print ('SpliceJunctionSummary.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))