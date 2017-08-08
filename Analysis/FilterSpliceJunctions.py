#!/usr/bin/python3

import os
import sys
import sqlite3
from AddJunctionsToDatabase import connectToDB, commitAndClose

def tableHeader():
	header = ['gene', 'chromosome:start-stop', 'annotation', 'n_gtex_seen', 'n_patients_seen', 'total_patient_read_count', 'total_gtex_read_count', 'total_read_count', 'sample:read_count', 'sample:norm_read_count']

	return '\t'.join(header)

def countGTEX(cur):
	cur.execute('select count(*) from SAMPLE_REF where type = 0;') # 0 = GTEX, 1 = PATIENT

	return cur.fetchone()[0]

def countPatients(cur):
	cur.execute('select count(*) from SAMPLE_REF where type = 1;') # 0 = GTEX, 1 = PATIENT

	return cur.fetchone()[0]

def writeToFile(res, file):
	with open(file, "w") as out:

		out.write(tableHeader())

		for row in res:
			out.write('\t'.join(str(element) for element in row) + '\n')

def sampleSpecificJunctions(cur, sample, min_read):

	count = str(countGTEX(cur))

	output = '_'.join([sample, 'specific', 'rc' + str(min_read), 'n_gtex_' + count])

	cur.execute('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.read_count),
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.norm_read_count)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id 
		where
		sample_ref.sample_name = ? and
		junction_counts.read_count >= ? and
		junction_ref.n_gtex_seen <= ?
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''',
		(sample, min_read, 0))

	writeToFile(cur.fetchall(), output)

def customSampleSpecificJunctions(cur, sample, min_read, max_n_gtex_seen, max_total_gtex_reads):

	if not max_n_gtex_seen:
		max_n_gtex_seen = 0

	if not max_total_gtex_reads:
		max_total_gtex_reads = 0

	if not min_read:
		min_read = 0

	output = '_'.join([str(sample), ('rc' + str(min_read)), ('maxGTEX' + str(max_n_gtex_seen)), ('maxGTEXrc' + str(max_total_gtex_reads))])

	cur.execute('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.read_count),
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.norm_read_count)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id 
		where
		sample_ref.sample_name = ? and
		junction_counts.read_count >= ? and
		junction_ref.n_gtex_seen <= ? and
		junction_ref.total_gtex_read_count <= ?
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''',
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
	cur.execute('''select sample_name,
	case
		when type = 0 then 'CONTROL'
		when type = 1 then 'PATIENT'
	END
	from SAMPLE_REF;''')

	for line in cur.fetchall():
		print('\t'.join(str(i) for i in line))

def makeSampleDict(cur):

	sampleDict = {}

	cur.execute('''select ROWID, sample from SAMPLE_REF;''')

	for bam_id, sample in cur.fetchall():
		sampleDict[bam_id] = sample

	return sampleDict

def printAllJunctions(cur):

	output = 'all_junctions_n_gtex_' + str(countGTEX(cur)) + '_n_paitents_' + str(countPatients(cur))

	cur.execute('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.read_count),
		group_concat(distinct sample_ref.sample_name||':'||junction_counts.norm_read_count)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id
		where junction_ref.total_read_count > 0
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''')

	writeToFile(cur.fetchall(), output)

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
		customSampleSpecificJunctions(cur, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
	elif sys.argv[1] == '--delete':
		deleteSample(cur, sys.argv[2])
	elif sys.argv[1] == '--all':
		printAllJunctions(cur)
	else:
		print('Invalid option. Use one of the following:')
		print('--printsamples')
		print('--sample	[SAMPLE]	[MIN_READ]')
		print('--custom [SAMPLE] [MIN_READ] [MAX_N_GTEX_SEEN] [MAX_TOTAL_GTEX_READS]')
		print('--delete [SAMPLE]')
		print('--all')

	commitAndClose(conn)