#!/usr/bin/python3

import os
import sys
import sqlite3
from AddJunctionsToDatabase import connectToDB, commitAndClose

def applyOutputFormatSettings(conn, cur):
	cur.execute('.separator \t') #tabs for awk
	cur.execute('.header on')

	conn.commit()

def restoreDefaultSettings(cur):
	cur.execute('.seperator |')
	cur.execute('.header off')

def countGTEX(cur):

	counter = 0

	cur.execute('select SAMPLE_NAME from SAMPLE_REF;')

	res = cur.fetchall()

	for name in res:
		if 'GTEX' in name:
			counter += 1

	return counter

def writeToFile(res, file):

	with open(file, "a") as out:
		for line in cur.fetchall():
			out.write(line)

def sampleSpecificJunctions(cur, sample, min_read):

	count = str(countGTEX(cur))

	output = '_'.join([sample, 'rc' + str(min_read), 'specific', 'n_gtex_' + count])

	cur.execute('''select 
	sample_ref.sample_name,
	gene_ref.gene,
	junction_ref.chromosome,
	junction_ref.start,
	junction_ref.stop,
	junction_ref.gencode_annotation,
	junction_counts.read_count,
	junction_counts.norm_read_count,
	junction_ref.n_patients_seen,
	junction_ref.n_gtex_seen,
	junction_ref.total_read_count,
	junction_ref.total_patient_read_count
	from 
	sample_ref inner join junction_counts on junction_counts.bam_id = sample_ref.rowid
	inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
	left join gene_ref on junction_ref.rowid = gene_ref.junction_id
	where
	sample_ref.sample_name = ? and
	junction_counts.read_count >= ? and
	junction_ref.n_gtex_seen <= ?;''',
	(sample, min_read, 0))

def sampleSpecificJunctionsWithControlTolerance(cur, sample, min_read, max_n_gtex_seen, max_total_gtex_reads):

	output = '_'.join([str(sample), ('rc' + str(min_read)), ('maxGTEX' + str(max_n_gtex_seen)), ('maxGTEXrc' + str(max_total_gtex_reads))])

	if not max_n_gtex_seen:
		max_n_gtex_seen = 0

	if not max_total_gtex_reads:
		max_total_gtex_reads = 0

	if not min_read:
		min_read = 0

	cur.execute('''select 
		sample_ref.sample_name,
		gene_ref.gene,
		junction_ref.chromosome,
		junction_ref.start,
		junction_ref.stop,
		junction_ref.gencode_annotation,
		junction_counts.read_count,
		junction_counts.norm_read_count,
		junction_ref.n_patients_seen,
		junction_ref.n_gtex_seen,
		junction_ref.total_read_count,
		junction_ref.total_patient_read_count
		from 
		sample_ref inner join junction_counts on junction_counts.bam_id = sample_ref.rowid
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		left join gene_ref on junction_ref.rowid = gene_ref.junction_id
		where
		sample_ref.sample_name = ? and
		junction_counts.read_count >= ? and
		junction_ref.n_gtex_seen <= ? and
		junction_ref.total_gtex_reads <= ?;''',
		(sample, min_read, max_n_gtex_seen, max_total_gtex_reads))

	writeToFile(cur.fetchall(), output)

def printSamplesInDB():

	cur.execute('''select * from SAMPLE_REF;''')
	for line in cur.fetchall():
		print (line)

if __name__=="__main__":
	
	conn, cur = connectToDB()
	applyOutputFormatSettings(conn, cur)

	# sample = argv[2]
	# min_read = int(argv[3])
	# max_n_gtex_seen = int(argv[4])
	# max_total_gtex_reads = int(argv[5])

	if argv[1] == '--printsamples':
		printSamplesInDB()
	elif argv[1] == '--samplespecificjunctions':
		sampleSpecificJunctions(cur, argv[2], int(argv[3]))
	elif argv[1] == '--custom':
		sampleSpecificJunctionsWithControlTolerance(cur, argv[2], argv[3], argv[4], argv[5])

	restoreDefaultSettings(cur)
	commitAndClose(conn)