import sys
import os
import errno
import argparse
import multiprocessing
import subprocess
from subprocess import Popen, PIPE
from cigar import Cigar
from datetime import datetime

# Beryl Cummings' code, modified slightly
def printSplices(path, spliceDict):
	for key in spliceDict:
		chrom, junctionStart, junctionEnd = key.split(':')
		timesSeenInSample = str(spliceDict[key])

		with open(path, "a") as out:
			out.write("\t".join([str(chrom),str(junctionStart),str(junctionEnd),timesSeenInSample]+"\n"))

# e.x. cigar='3M1D40M20N'
def parseCIGARForIntrons(cigar):

	if 'N' in cigar:
		cigar = cigar.split('N')[0] + 'N' #remove all information after intron
	else:
		raise Exception('No intron detected')

	offset = 0
	matchedExon = 0
	intronLength = 0
	
	for c in list(Cigar(cigar).items()): # returns list of tuples : [(20, 'N')]
		if c[1] == 'N':
			intronLength += int(c[0])
		elif c[1] == 'D':
			offset += int(c[0])
		elif c[1] == 'I':
			offset -= int(c[0])
		elif c[1] == 'M':
			matchedExon += int(c[0])
		## soft clipping is ignored
		## hard clipping is ignored too

	return offset, matchedExon, intronLength

def intronDiscovery(poolArguement):

	bamFiles, gene, gene_type, chrom, start, stop, cwd = poolArguement

	print ('processing ' + gene)

	pos = ''.join([chrom, ':', start, '-', stop])

	for bam in bamFiles:

		spliceDict = {}
		geneFilePath = (cwd + "/" + bam[:-4] + "/" + gene + ".txt")

		try:
			exitcode, stdout, stderr = run(' '.join(['samtools view', bam, pos]))
		except Exception as e:
			print ('Exception message: ' + str(e))
			print ("Exception occured while running \"samtools view\" on " + bam + " for position " + pos + " Skipping.")
			continue

		if not stdout:
			#print ('No introns found for ' + gene + ' at ' + pos + ' in ' + bam)
			continue

		for line in stdout.splitlines():

			elems = line.decode().split()

			alignmentStart = int(elems[3])
			cigar = str(elems[5])
			alignmentScore = int(elems[1])
 
			if 'N' not in cigar:  	#only get introns
				continue

			if (alignmentScore >= 256):  	#only primary alignments
				continue

			if not ((alignmentStart > int(start)) and (alignmentStart < int(stop))):  	#check if alignment start is after known junction start but before known junction end 
				continue

			try:
				offset, matchedExon, intronLength = parseCIGARForIntrons(cigar)
			except Exception as e:
				print ('Error message: ' + str(e))
				print ('Error trying to parse CIGAR string: ' + cigar +  ' with the bam file ' + bam +  ' and the position: ' + pos + ' Skipping.')
				continue

			junctionStart = alignmentStart + matchedExon + offset
			junctionEnd = junctionStart + intronLength

			# Beryl Cummings' Code, taken from makeUniqSpliceDict()
			uniqueSplice = ' '.join([chrom, str(junctionStart), str(junctionEnd)])
			
			if uniqueSplice not in spliceDict:
				spliceDict[uniqueSplice] = {}

			if sample not in spliceDict[uniqueSplice]:	
				spliceDict[uniqueSplice][sample] = 1
			else:
				spliceDict[uniqueSplice][sample] += 1


		del stdout # saves ram in between samtool calls

		if spliceDict:
			printSpliceDict(geneFilePath, spliceDict)
			del spliceDict

	print ('finished ' + gene)

def processGenesInParallel(transcriptFile, bamList, numProcesses):

	cwd = os.getcwd()
	numProcesses = str(numProcesses) #handles ambigious Python behaviour: when set to default processes is an int, when specified as parameter processes is a string

	bamFiles = []
	poolArguements = []

	print ("Creating a pool with " + numProcesses + " processes")
	pool = multiprocessing.Pool(int(numProcesses))
	print ('pool: ' + str(pool))

	with open(bamList) as bl:
		for i in bl:

			i = i.rstrip("\n")

			bamLocation = os.getcwd() + '/' + i
			if not os.path.isfile(bamLocation):
				print ('bam file: ' + i + ' does not exist in CWD! Skipping.')
				continue

			os.system("mkdir " + bamLocation[:-4])
			bamFiles.append(i)

	with open(transcriptFile) as tf:
		for line in tf:

			elems = line.strip().split()
			try:
				gene, gene2, plus, chrom, start, stop, gene_type = elems #edit the transcript file so that you only deal with junction coordinates
			except Exception as e:
				print ('Error while parsing transcript file named: ' + str(transcriptFile) + "\n" + 'Error message: ' + str(e) + "\nExiting.")
				exit (3)

			poolArguements.append((bamFiles, gene, gene_type, chrom, start, stop, cwd))

	pool.map(intronDiscovery, poolArguements) # run the worker processes
	pool.close()
	pool.join()
	
if __name__=="__main__":

	print ('SpliceJunctionDiscover.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Discover splice junctions from a list of bam files')
	parser.add_argument('-transcriptFile',help="A list of positions that you want to discover junctions in",action='store',default = "/home/dennis.kao/largeWork/protein-coding-genes.list")
	parser.add_argument('-bamList',help='A text file containing the names of bam files you want to discover splice junctions in each on a seperate line',default='bamlist.list')
	parser.add_argument('-processes',help='number of processes to run multiple instances of: "samtools view", default=10',default=10)
	args=parser.parse_args()

	print ('Working in directory' + str(os.getcwd()))
	print ('Transcript file is ' + str(args.transcriptFile))
	print ('Identifying splice junction is ' + str(args.bamList))

	processGenesInParallel(args.transcriptFile, args.bamList, args.processes)
	
	# transcriptFile = str(args.transcriptFile).rsplit('/')[-1] #remove paths

	print ('SpliceJunctionDiscover.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
