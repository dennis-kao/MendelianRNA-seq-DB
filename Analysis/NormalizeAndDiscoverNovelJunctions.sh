#
#	NormalizeAndDiscoverNovelJunctions.sh
#	May 18th, 2017
#	dennis.kao@sickkids.ca
#
#	
#	Usage: input=input.txt sample=sample ./NormalizeAndDiscoverNovelJunctions.sh
#	OR		input=input.txt sample=sample mindread=10 minnormread=0.5 transcript_model=gencode.txt ./NormalizeAndDiscoverNovelJunctions.sh
#
#
#	Input - file generated by SpliceJunctionDiscovery.py or rnaseq.splice_junction_discovery.pbs
#	Output - Generates 2 files in the same directory as the input file: 
# 		1) normalized_input.txt: input file with an additional column for normalized read counts
# 		2) novel_normalized_input.txt: file from 1) that has been filtered for a minimum read count, minimum normalized read count and only has junctions specific to sample

beryl_home=~/tools/MendelianRNA-seq
baseDir=`pwd`

#	User parameters
if [ -z "$input" ];
	then
		echo "ERROR - Specify an input file"
		exit 1
fi

if [ -z "$sample" ];
	then
		echo "ERROR - Specify the sample name used in the input file"
		exit 2
fi

if [ -z "$minread" ];
	then
		minread=10
fi

if [ -z "$minnormread" ];
	then
		minnormread=0.5
fi

if [ -z "$transcript_model" ];
	then
		transcript_model="$beryl_home/gencode.comprehensive.splice.junctions.txt"
fi

inputFileName=`basename $input`
output="normalized_$inputFileName"
outputFilePath=`dirname $input`

#	Actual computation
echo "1. Normalizing read counts in $input"
~/tools/MendelianRNA-seq/Analysis/NormalizeSpliceJunctionValues.py -transcript_model=$transcript_model -splice_file=$input --normalize > $outputFilePath/$output 
echo "Output file is called $output"

echo "2. Filtering for minimum read count, minimum normalized read count and novel junctions"
echo | cat $output | grep $sample | awk '{ if ($5 == 1 && $4 >= $(minread)) print $0}' > $outputFilePath/novel_$output
echo "Output file is called novel_$output"

echo "DONE - NormalizeAndDiscoverNovelJunctions.sh ran successfully"
echo "Output files can be found in: $outputFilePath"