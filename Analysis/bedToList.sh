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

if [ -z "$1" ];
	then
		echo "ERROR - Please specify an input file" 
		exit 1
fi

output=$1".list"

echo "Input file: "$1
cat $1 | awk '{print $4"\t"$4"\t""+""\t"$1"\t"$2"\t"$3"\t""NEXON"}' > $output
echo "Output file: "$output
