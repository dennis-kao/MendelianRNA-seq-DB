#!/usr/bin/python3

import os
import sys
import sqlite3

def connectToDB():

	conn = sqlite3.connect('SpliceJunction.db')
	cur = conn.cursor()

	return conn, cur

def applyOutputFormatSettings(cur):

	cur.execute('.separator \t') #tabs for awk
	cur.execute('.header on')

def restoreDefaultSettings(cur):

	cur.execute('.seperator |')
	cur.execute('.header off')

def printAll(res):

	for line in res:
		print (line)

def sampleSpecificJunctions(sample, min_read):

def affectedSpecificJunctions(cur, sample, min_read, max_n_gtex_seen, max_total_gtex_reads):

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

def printSamplesInDB():

	cur.execute('''select * from SAMPLE_REF;''')
	printAll(cur.fetchall())

if __name__=="__main__":
	