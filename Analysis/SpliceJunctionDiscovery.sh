#This version (v6) has the sed command to account for softclipped reads

# what is transriptFile - processed from Gencode ?

# some change

beryl_home=/home/naumenko/tools/MendelianRNA-seq

#transcriptFile=${beryl_home}/NMD_genes.list

transcriptFile=$1
bamList=$2

baseDir=`pwd`

IFS=$'\n'

groupname=`basename $transcriptFile | sed 's/gencode.v19.genes.filtered.group//'`

echo "Started on" $(date)
echo "Working in directory $baseDir"
echo "Transcript file is $transcriptFile"
echo "Identifying splice junction is $bamList"

echo -e "Gene\tType\tChrom\tStart\tEnd\tNTimesSeen\tNSamplesSeen\tSamples:NSeen" > All.${groupname}.splicing.txt
for line in `cat $transcriptFile`
do
    start=`echo $line | cut -f5`
    stop=`echo $line | cut -f6`
    chrom=`echo $line | cut -f4`
    gene=`echo $line | cut -f1`
    gene_type=`echo $line | cut -f7`
    base=$baseDir/$gene
    pos=$chrom":"$start"-"$stop
    echo 'processing' $gene
    for i in `cat $bamList`
    do
    sample=`echo $i | awk -F"/" '{print $NF}' | cut -d "." -f 1`


    #$6 ~ /N/ - splitted read - CIGAR string
    #$5=60 mapping quality threshold is very high, something is wrong with our alignment, we have 255 mostly 
    #int($2) < 256 - not a secondary alignment
    #$4 - mapping start, $6 - CIGAR
    #parsing CIGAR MATCH GAP
    #soft clipping
    #output: gene gene_type(number of exons) sample chr junction_start junction_end exonic_part intron_length
    

    # only get split reads, no secondary reads, print mapping start, space, then CIGAR string, print everything before and not including 'N', 
    # replace all instances of 'M' with a space, outright remove soft clipping information if its there, print resulting string - $4 can be empty or $3 is CIGAR remnant and $4 is 

    # output up until final awk line is one of the two: 
    # 1) start, matching nucleotides, intron length EX) 19670857 56 752
    # 2) start, matching nucleotides, #D## - deletion and matching nucleotides, intron length EX) 19670857 56 1D15 752

    samtools view ${i} ${pos} | awk '($6 ~ /N/)' | awk 'int($2)< 256' | awk -v sta=$start -v sto=$stop '$4>sta&&$4<sto {print $4,$6}' | \
    cut -d "N" -f 1 | tr 'M' ' ' |  sed -r 's/[0-9]+S//' | awk -v s=$sample -v ge=$gene -v t=$gene_type -v chr=$chrom '{print ge,t,s,chr,$1+$2-1,$1+$2+$3,$2,$3,$4}' >> ${base}.splicing.txt
    done
    echo $gene 'is complete'
    ~/work/tools/bcbio/anaconda/bin/python ${beryl_home}/Analysis/SpliceJunctionSummary.py <  ${base}.splicing.txt >> All.${groupname}.splicing.txt
    #rm ${base}.splicing.txt
done

echo "Finished on" $(date)

