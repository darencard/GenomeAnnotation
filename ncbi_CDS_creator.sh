#!/usr/bin/env zsh

usage()
{
cat << EOF
ncbi_CDS_creator.sh

Version 1.0 (16 February, 2016)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (daren.card@gmail.com).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script creates a CDS file from a NCBI genome using the feature table. The only
dependencies are rsync (should be installed on all Unix systems) and the NCBI Entrez
Direct E-utilities (www.ncbi.nlm.nih.gov/books/NBK179288). The script can automatically
download the feature table for the target genome and will use it to gather CDS sequences
from NCBI databases using Entrez. The user only has to specify the type of taxa (-t) using
the list of options and the name of the organism (-n) in the format "Genus_species". Note that
the name must adhere to the naming convention adopted by NCBI. Alternatively, the user can
pass in the feature table file if they have already obtained it from NCBI by specifying
the -f flag with the file (and leaving off the -n flag). The user can specify an
output directory (-o) where all intermediate, temporary and final files will be placed.
Output is a multi-FASTA sequence of the full genome CDSs.

Note: this script will apparently not work for all NCBI genomes and failed on those tested
outside of the 'refseq' directory at ftp://ftp.ncbi.nlm.nih.gov/genomes/ (i.e., those in 'all').
Use with caution and if an error arises, it is probably not an easy fix!

zsh ncbi_CDS_creator.sh -n <Genus_species> [-t <taxon> | -f <feature_table>] [-o <output> -h]

OPTIONS:
        -h      usage information and help (this message)
        -n      the name of the organism (following NCBI convention) as "Genus_species" (e.g., Homo_sapiens)
	-t	taxa type - either 'archaea', 'bacteria', 'fungi', 'invertebrate', 'plant', 'protozoa', 
		'vertebrate_mammalian', 'vertebrate_other', or 'viral'
	-f	feature table from NCBI for the genome
        -o      directory for output files (files are named automatically) [.]
EOF
}

NAME=
TAXON=
FEAT=
OUTPUT="."

while getopts "hn:t:f:o:" OPTION
do
        case $OPTION in
                h)
                        usage
                        exit 1
                        ;;
		n)
			NAME=$OPTARG
			;;
		t)
			TAXON=$OPTARG
			;;
		f)
			FEAT=$OPTARG
			;;
		o)
			OUTPUT=$OPTARG
			;;
	esac
done

if [[ -z $NAME ]]
then
	usage
	exit 1
fi

if [[ -n $TAXON && -n $FEAT ]]
then
	usage
	exit 1
fi

if [ ! -d "$OUTPUT" ]; then
	mkdir ${OUTPUT}
fi

SAFE_NAME_1=$(echo $NAME | cut -d "_" -f 1)
SAFE_NAME_2=$(echo $NAME | cut -d "_" -f 2)

if [[ -z $FEAT && -n $TAXON ]]; then

# use rsync to download the feature table for the target genome
# should work for all refseq genomes
echo "\n**Downloading the genome feature table**"
rsync -av rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/${TAXON}/${NAME}/representative/\*/\*_feature_table.txt.gz ${OUTPUT}

# from downloaded feature table, parse CDS feature lines, eliminate redundancy using geneID field, and export accessions
echo "\n**Extracting CDS metadata**"
# if file was downloaded
zcat ${OUTPUT}/*${SAFE_NAME_1}*${SAFE_NAME_2}*_feature_table.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}$1=="CDS"{print $16,$11,$13,$18}' | \
sort -k1,1 -k4,4nr -k3,3 | awk 'BEGIN{FS="\t";OFS="\t";gene=0}{if(gene!=$1){print $3};gene=$1}END{if(gene!=$1){print $3","}}' \
> ${OUTPUT}/${NAME}_id_list.tmp

fi

if [[ -n $FEAT && -z $TAXON ]]; then
# from supplied feature table, parse CDS feature lines, eliminate redundancy using geneID field, and export accessions
echo "\n**Extracting CDS metadata**"
# if file was supplied
zcat ${FEAT} | awk 'BEGIN{FS="\t";OFS="\t"}$1=="CDS"{print $16,$11,$13,$18}' | \
sort -k1,1 -k4,4nr -k3,3 | awk 'BEGIN{FS="\t";OFS="\t";gene=0}{if(gene!=$1){print $3};gene=$1}END{if(gene!=$1){print $3","}}' \
> ${OUTPUT}/${NAME}_id_list.tmp

fi

#split list of accessions into smaller batches (not sure why... inherited from hacked script)
split -l 1000 ${OUTPUT}/${NAME}_id_list.tmp ${OUTPUT}/${NAME}_tmp.

# for each batch of accessions, remove new lines and separate by only commas
foreach file (`ls $OUTPUT/${NAME}_tmp*`)
  cat ${file} | sed 's/$/,/' | perl -pe 'chop' | sed 's/$/\n/' >> ${OUTPUT}/${NAME}_chunks.tmp
end

echo "\n**Creating CDS sequences**"

# for each accession, run efetch to extract CDS sequence and append to the output file
foreach chunk (`cat ${OUTPUT}/${NAME}_chunks.tmp`)
  efetch -db nuccore -id "${chunk}" -format fasta_cds_na >> ${OUTPUT}/${NAME}.transcript_cds.fa
end

# remove all temporary intermediate files
rm -f ${OUTPUT}/*tmp*

echo "\n**Done! See ${OUTPUT}/${NAME}.transcript_cds.fa for full genome CDS**"
