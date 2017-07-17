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

def connectToDB():
	conn = sqlite3.connect('SpliceJunctions.db')
	cur = conn.cursor()

	return conn, cur

def initializeDB():

	conn, cur = connectToDB()

	cur.execute('''create table if not exists SAMPLE_REF (
		sample_name varchar(50) primary key, 
		type tinyint not null);''') # type = {0, 1, 2} 
									# GTEX, patient, gencode_reference

	cur.execute('''create table if not exists JUNCTION_REF (
		chromosome tinyint not null,
		start unsigned big int not null,
		stop unsigned big int not null,
		gene varchar(20) not null,
		gencode_annotation tinyint not null,
		n_junctions_seen unsigned big int,
		n_times_seen unsigned big int,
		primary key (chromosome, start, stop));''') # gencode_annotation = {0, 1, 2, -1}
													# none, one, two, exon skipping

	cur.execute('''create table if not exists JUNCTION_COUNTS (
		read_count unsigned big int not null,
		norm_read_count float,
		sample_id integer not null,
		junction_id integer not null,
		foreign key(sample_id) references SAMPLE_REF(ROWID),
		foreign key(junction_id) references JUNCTION_REF(ROWID));''')

	cur.execute('''PRAGMA journal_mode=WAL;''')

	if not ('wal' in cur.fetchone()):
		print("Could not set SQLite database to WAL mode. Exiting.")
		exit(1)

	conn.commit()
	conn.close()

def makeUniqSpliceDict(gene_file):

	d = {}

	with open(gene_file, "r") as SpliceFile:
		for line in SpliceFile:

			elems = line.strip().split()

			gene, sample, chrom, junctionStart, junctionEnd, matchedExon, intronLength = elems
			uniqSplice = "%s %s %s %s"%(gene, chrom, junctionStart, junctionEnd)
			
			if uniqSplice not in d:
				d[uniqSplice] = {}

			if sample not in d[uniqSplice]:	
				d[uniqSplice][sample] = 1
			else:
				d[uniqSplice][sample] += 1

	return d	

def getJunctionID(cur, gene, chrom, start, stop):

	# NOTE: gencode junctions will always have a gencode_annotation value of 2

	# check if start and stop are apart of an existing gencode annotation
	cur.execute('''select ROWID from JUNCTION_REF where gencode_annotation is 2 and chromosome is ? and start is ? and stop is ?;''', (chrom, start, stop))
	res = cur.fetchone()

	if res:
		return res
	else: # if no such junction, add it. Determine the annotation: novel junction, only one annotated or a case of exon skipping?

		# determine gencode_annotation
		cur.execute('''select * from JUNCTION_REF where gencode_annotation is 2 and chromosome is ? and start is ?''', (chrom, start))
		isStartAnnotated = cur.fetchone()

		cur.execute('''select * from JUNCTION_REF where gencode_annotation is 2 and chromosome is ? and stop is ?''', (chrom, stop))
		isStopAnnotated = cur.fetchone()

		if isStopAnnotated and isStartAnnotated:
			annotation = -1 # exon skipping
		elif isStopAnnotated:
			annotation = 1
		elif isStartAnnotated:
			annotation = 1
		else:
			annotation = 0 # novel junction

		# add new junction with gencode_annotation
		cur.execute('''insert into JUNCTION_REF (chromosome, start, stop, gene, gencode_annotation) values (?, ?, ?, ?);''', (chrom, start, stop, gene, annotation))
		return cur.lastrowid

def summarizeGeneFile(gene_file):

	spliceDictionary = makeUniqSpliceDict(gene_file)

	conn, cur = connectToDB()

	for junction in spliceDictionary:

		gene, chrom, start, stop = junction.rstrip().split()
		
		junction_id = getJunctionID(cur, gene, chrom, start, stop)
		
		for sample in spliceDictionary[junction]:

			cur.execute('''insert into ''')

			# write results to DB
			item = "%s:%s"%(sample,spliceDictionary[junction][sample])

	conn.close()
	# os.system('rm ' + gene_file)

def parallel_process_gene_files(num_processes):

	gene_files = os.system('ls *.txt').splitlines()

	print ("Creating a pool with " + numProcesses + " processes")
	pool = multiprocessing.Pool(int(numProcesses))
	print ('pool: ' + str(pool))

	pool.map(summarizeGeneFile, gene_files)
	pool.close()
	pool.join()

def storeGencodeJunctions(gencode_file):

	conn, cur = connectToDB()

	print ('Started adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	with open(gencode_file, "r") as gf:
		for commitFreq, line in gf:

			chrom, start, stop, gene, transcript, gene_type = line.strip().split()

			cur.execute('''insert into JUNCTION_REF (chromosome, start, stop, gene, gencode_annotation) values (?, ?, ?, ?, ?);''', (chrom, start, stop, gene, '2'))

			if commitFreq % 500 == 0: #	yes this works
				conn.commit()

	conn.commit()
	conn.close()

	print ('Finished adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

if __name__=="__main__":

	print ('SpliceJunctionSummary.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Summarize the read counts of the junctions reported by SpliceJunctionDiscovery.py')
	parser.add_argument('-processes',help='Number of worker processes to parse gene files, default=10.',default=10)
	parser.add_argument('-gencode', help='Specify a path to a file containing gencode junctions to be added to the database. This only needs to be run once.',default='')
	#parser.add_argument('-db',help='Name of the database file including extension, default=SpliceJunctions.db',default='SpliceJunctions.db')
	args=parser.parse_args()

	print ('Working in directory' + str(os.getcwd()))

	initializeDB()

	if gencode: 
		storeGencodeJunctions(args.gencode)

	# parallel_process_gene_files(args.processes)

	print ('SpliceJunctionSummary.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))