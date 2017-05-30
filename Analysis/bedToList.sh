#	bedToList
#
#	May 29th, 2017
#	dennis.kao@sickkids.ca
#
#	Usage: input=/path/to/gene_list.bed ./bedToList.sh
#
#	Input: a .bed file containing information about the location of genes in this format:
#	"CHROMOSOME	START	STOP	NAME"
#	Output: a .list file to be used as the gene_list parameter in SpliceJunctionDiscovery.sh
#

if [ -z "$input" ];
	then
		echo "ERROR - Please specify an input file" 
		exit 1
fi

output=$input".list"

echo "Input file: "$input
cat $input | awk '{print $4"\t"$4"\t""+""\t"$1"\t"$2"\t"$3"\t""FILLER"}' > $output
echo "Output file: "$output
