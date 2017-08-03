#!/usr/bin/python3

import os
import sys
import sqlite3
from AddJunctionsToDatabase import connectToDB, commitAndClose

def translateAnnotation(annotation):

	if annotation == 0:
		return 'NONE'
	elif annotation == 1:
		return 'START'
	elif annotation == 2:
		return 'STOP'
	elif annotation == 3:
		return 'BOTH'
	elif annotation == 4:
		return 'EXON_SKIP'

def tableHeader():
	header = ['gene', 'chromosome', 'start', 'stop', 'annotation', 'read_count', 'norm_read_count', 'n_patients_seen', 'n_gtex_seen', 'total_read_count', 'total_patient_read_count']

	return '\t'.join(header)

def countGTEX(cur):
	cur.execute('select count(*) from SAMPLE_REF where type = 0;') # 0 = GTEX, 1 = PATIENT

	num_gtex = cur.fetchone()[0]

	return num_gtex

def writeToFile(res, file):
	with open(file, "w") as out:

		out.write(tableHeader())

		for row in res:

			line = []

			for i, element in enumerate(row):
				if i == 4:
					line.append(translateAnnotation(element))
				else:
					line.append(str(element))

			out.write('\t'.join(line) + '\n')

def sampleSpecificJunctions(cur, sample, min_read):

	count = str(countGTEX(cur))

	output = '_'.join([sample, 'rc' + str(min_read), 'specific', 'n_gtex_' + count])

	cur.execute('''select 
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

	writeToFile(cur.fetchall(), output)

def customQuery(cur, sample, min_read, max_n_gtex_seen, max_total_gtex_reads):

	if not max_n_gtex_seen:
		max_n_gtex_seen = 0

	if not max_total_gtex_reads:
		max_total_gtex_reads = 0

	if not min_read:
		min_read = 0

	# if not min_n_patients_seen:
	# 	min_n_patients_seen = 0

	# if not min_total_patient_reads:
	# 	min_total_patient_reads = 0

	output = '_'.join([str(sample), ('rc' + str(min_read)), ('maxGTEX' + str(max_n_gtex_seen)), ('maxGTEXrc' + str(max_total_gtex_reads))])

	cur.execute('''select 
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

def deleteSample(cur, sample):

	cur.execute('select ROWID, type from sample_ref where sample_name = ?;', (sample, ))
	res = cur.fetchone()

	if not res:
		print ("Sample %s does not exist in the database!" % sample)
		exit(1)
	else:
		bam_id, bam_type = res

	cur.execute('select junction_id, read_count from JUNCTION_COUNTS where bam_id = ?;', (bam_id, ))

	for junction_id, read_count in cur.fetchall():

		if bam_type == 0:
			cur.execute('''update JUNCTION_REF set 
				n_gtex_seen = n_gtex_seen - 1,
				total_read_count = total_read_count - ?,
				total_gtex_read_count = total_gtex_read_count - ?
				where ROWID = ?;''', (read_count, read_count, junction_id))
		elif bam_type == 1:
			cur.execute('''update JUNCTION_REF set 
				n_patients_seen = n_patients_seen - 1,
				total_read_count = total_read_count - ?,
				total_patient_read_count = total_patient_read_count - ?
				where ROWID = ?;''', (read_count, read_count, junction_id))
		else:
			raise Exception ('FATAL ERROR - bam_id is not 0 or 1')

	cur.execute('''delete from JUNCTION_COUNTS where junction_id = ?;''', (junction_id, ))
	cur.execute('''delete from SAMPLE_REF where sample_name = ?;''', (sample, ))

	print ("Successfully deleted %s from database!" % sample)

def printSamplesInDB(cur):
	cur.execute('''select * from SAMPLE_REF;''')

	for line in cur.fetchall():
		print (line)

if __name__=="__main__":
	
	conn, cur = connectToDB()

	# sample = sys.argv[2]
	# min_read = int(sys.argv[3])
	# max_n_gtex_seen = int(sys.argv[4])
	# max_total_gtex_reads = int(sys.argv[5])

	if sys.argv[1] == '--printsamples':
		printSamplesInDB(cur)
	elif sys.argv[1] == '--sample':
		sampleSpecificJunctions(cur, sys.argv[2], int(sys.argv[3]))
	elif sys.argv[1] == '--custom':
		customQuery(cur, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
	elif sys.argv[1] == '--delete':
		deleteSample(cur, sys.argv[2])
	else:
		print('Invalid option. Use one of the following:')
		print('--printsamples')
		print('--sample	[SAMPLE]	[MIN_READ]')
		print('--custom [SAMPLE] [MIN_READ] [MAX_N_GTEX_SEEN] [MAX_TOTAL_GTEX_READS]')
		print('--delete [SAMPLE]')

	commitAndClose(conn)