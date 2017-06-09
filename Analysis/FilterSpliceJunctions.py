import sys
import os
import argparse

def filter(file, sampleName, minReadCount, minNormReadCount):

	minReadCount = int(minReadCount)
	minNormReadCount = int(minNormReadCount)

	with open(file) as f:
		for line in f:

			gene, gene_type, pos, ntimes, nsamp, samptimes_sorted, tag, annotation = line.strip().split("\t")

			if (sampleName not in line) or (int(nsamp) != 1): # only get junctions specific to sample
				continue

			if int(ntimes) < minReadCount: # minRead filtration
				continue

			if (tag == 'Neither annotated') or ('*' in line): # if neither sites are annotated or normalized read count could not be determined, print the line anyways
				print line
				continue

			normReadCount = int(tag.split(':')[0]) # tag = normReadCount:sample

			if normReadCount < minNormReadCount: # normalization filtration
				continue

			print line

if __name__=="__main__":

	# print 'FilterSpliceJunction.py started working on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")

	parser = argparse.ArgumentParser(description = 'Filter out splice junctions based on read count, normalized read count and specificity to a sample')
	parser.add_argument('-file',help="The path to the file produced by NormalizeSpliceJunctionValues.py")
	parser.add_argument('-sample',help="The name of the bam file you want to discover novel junctions in without the /".bam/" extension")
	parser.add_argument('-minRead',help="The minimum read count a junction needs to have to be considered, default=10", default=10)
	parser.add_argument('-minNormRead',help="The minimum normalized read count a junction needs to have to be considered, default=0.5", default=0.5)
	args=parser.parse_args()

	# print 'Working in directory' + str(os.getcwd())
	# print 'Splice file is ' + str(args.transcriptFile)

	filter(args.file, args.sample, args.minRead, args.minNormRead)

	# input = str(args.file).rsplit('/')[-1]
	# output = 'threshold' + str(args.normRead) + '_novel_' + sample + '_norm_' + input

	# print 'Output file is: ' + output
	# print 'FilterSpliceJunction.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")