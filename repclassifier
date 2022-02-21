#!/usr/bin/env bash

usage()
{
cat << EOF
repclassifier

Version 1.1 (2022-02-04)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (daren.card@gmail.com).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script attempts to classify unknown repeat elements based on databases or
libraries of known repeat elements. It is designed to annotate unknown repeats 
resulting from RepeatModeler, but could be used in other circumstances. User 
supplies the unknown elements, a clade for classification using Repbase, and a 
library of known elements. All libraries must be in standard FASTA format. 
The user should also supply an output directory that will store the results of
the run. When naming the output directory, note that this script is meant to be 
run iteratively. Finally, the user should also specify an existing library file 
of known repeat elements to append the newly classified elements to - normally 
this will be the same file supplied as a library of known elements.

Dependencies include RepeatMasker, which is used to perform the repeat searches 
(it is a handy way of automating BLAST), bioawk, and seqkit.

This script perform classification two ways:
1. If a clade name within the Repbase database is supplied using the -d option, 
the script will first use the supplied Repclass clade to attempt to classify unknowns. 
Remaining unknowns from that classification are then used to run a subsequent round 
of classification using the library of known elements provided.
2. If no Repbase database clade is provided using the -d option, the script will only 
perform a single round of classification using the library of known elements provided. 

Therefore, normally, this script should be run iteratively/sequentially for several 
rounds. In each round, known elements from previous rounds appended to existing known 
elements can be used to classify repeats. Repeating this leads to an increasing number 
of classified elements over rounds. Thus, one may wish to name output directories 
with sequential names like "round1 - roundN" to keep track of all sets of results. 
Note: If running this script iteratively, be sure to turn off the Repbase search after
the first run to avoid unnecessary computation.

This script uses RepeatMasker with default settings to perform the classification. 
When an unknown element is matched in part or whole against a single family/subfamily 
from the known repeat database or library, that family/subfamily is used as the 
classifier. If an unknown element is matched in part against more than one family/
subfamily from the known repeat database or library, one of three things happens:
1. If all matches are to the same exact subfamily, that subfamily is used for
classification.
2. If matches are to elements in different subfamilies of the same family, the
shared family is used for classification and no subfamily classification is given.
3. If matches are to elements in different families, the elements are regarded as 
chimeric and are not classified - they are simply left as unknowns.

There are several output files that will be contained in a directory supplied with -d:
1. Standard RepeatMasker outputs are included in subdirectories "known_match" for searches 
against the known repeat library and "database_match" for searches against the repeat 
database. Round 1 will have both subdirectories and subsequent rounds will only have the 
"known_match" subdirectory. See RepeatMasker documentation for details of these outputs.
2. The "full_match" subdirectory contains concatenated versions of the .cat and 
.out files of any known or database searches. It also contains repeats that are classified 
unambiguously based on subfamilies (#1 above; subfamily_unambiguous_classified_elements.txt), 
classified unambiguously based on families (#2 above; family_unambiguous_classified_elements.txt), 
and apparantly chimeric elements (#3 above; chimeric_elements.txt). The unambiguous elements 
are combined together into combined_classified_elements.txt, which represents all successfully 
classified repeats for a given round.
3. New .known and .unknown files provide libraries of all classified elements (newly classified 
elements appended to file supplied with -a option) and all elements that remain unclassified, 
respectively. These can be passed into subsequent rounds of classification with this program.

repclassifier -u <unknown_repeats> -k <known_repeats> -a <known_repeats> -o <output_dir>
	 [ -d <repbase_clade> -t <threads> -h ]

OPTIONS:
  -h		usage information and help (this message)
  -u		unknown repeat sequences in FASTA format
  -k		known repeat sequences in FASTA format
  -d		name of species or clade from Repbase to use for annotation
  -a		existing library of known repeat sequences in FASTA format
		to append newly classified elements to; file must exist - 
		even if empty; usually same file as set for '-k'
  -o		name for output directory for results
  -t		number of computer threads to use in BLAST [1]
EOF
}

THREAD=1

while getopts "hu:k:d:a:o:t:" OPTION
do
    case $OPTION in
	help)
		usage
		exit 1
		;;
	u)
		UNKNOWN=$OPTARG
		;;
	k)
		KNOWN=$OPTARG
		;;
	d)
		DATABASE=$OPTARG
		;;
	a)
		APPEND=$OPTARG
		;;
	o)
		OUTPUT=$OPTARG
		;;
	t)
		THREAD=$OPTARG
		;;
	?)
		usage
		exit
		;;
    esac
done

if [[ -z ${UNKNOWN} ]] || [[ -z ${OUTPUT} ]] || [[ -z ${KNOWN} ]]
then
	usage
	exit 1
fi


if [[ -z ${DATABASE} ]]
then

  # create output directories
  mkdir -p ${OUTPUT}/known_match

  # IGNORE: determine input round (${ROUND} - 1)
  # INROUND=`echo -e "${ROUND}" | awk '{ print $1 - 1 }'`

  # IGNORE: create "round_N" symlink
  # ln -s ${UNKNOWN} ${INROUND}.input

  # create "input.fasta" symlink
  ln -s ${UNKNOWN} input.fasta

  # perform known repeat masking
  cmd1="RepeatMasker -pa ${THREAD} -e ncbi -nolow -dir ${OUTPUT}/known_match -lib ${KNOWN} input.fasta"
  #echo ${cmd1}
  eval ${cmd1}

  # create combined output directory
  mkdir -p ${OUTPUT}/full_match

  # combine all .out and .cat.gz files together
  # note that normal header from .out files is not retained!
  cat ${OUTPUT}/known_match/*.out | awk '{ if ($1 != "SW" && $1 != "score") print $0 }' | sed '/^$/d' | sed '/^[[:space:]]*$/d' > ${OUTPUT}/full_match/${OUTPUT}.out

  cat ${OUTPUT}/known_match/*.cat > ${OUTPUT}/full_match/${OUTPUT}.cat

  # summarize repeat classifications
  # subfamily non-ambiguous elements - elements in the same subfamily, which means they are totally nonambiguous (attach full info)
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  cut -f 1 | sort | uniq |
  while read unknown;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${unknown}" | cut -f 4 | sort | uniq |
  awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR == 1) print unknown, $1 }';
  done > ${OUTPUT}/full_match/subfamily_unambiguous_classified_elements.txt

  # family non-ambiguous elements - elements in same family even if subfamilies are not the same (attach family info)
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  cut -f 1 | sort | uniq |
  while read unknown;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${unknown}" | cut -f 4 | sort | uniq |
  awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR > 1) print unknown }';
  done |
  while read element;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${element}" | cut -f 4 |
  awk -v OFS="\t"  -v element="${element}" -F "/|-" '{ print element, $1, $2 }';
  done |
  sort | uniq -c |
  awk '{ print $2 }' | uniq -u |
  while read clean;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${clean}" | cut -f 4 |
  awk -v OFS="\t"  -v clean="${clean}" -F "/|-" '{ print clean, $1"-"$2 }';
  done |
  uniq > ${OUTPUT}/full_match/family_unambiguous_classified_elements.txt

  # combine classified elements
  cat ${OUTPUT}/full_match/subfamily_unambiguous_classified_elements.txt ${OUTPUT}/full_match/family_unambiguous_classified_elements.txt > ${OUTPUT}/full_match/combined_classified_elements.txt

  # chimeric elements
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep -v -f <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  cut -f 1 | sort | uniq |
  while read element;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${element}" | cut -f 4 | paste -s -d "," |
  awk -v OFS="\t" -v element="${element}" '{ print element, $1 }';
  done > ${OUTPUT}/full_match/chimeric_elements.txt

  # create FASTA output of unknowns that remain unknown
  cat input.fasta | bioawk -c fastx '{ print $name }' |
  grep -wv -f <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  seqkit grep -f - input.fasta > ${OUTPUT}/${OUTPUT}.unknown

  # create FASTA output of unknowns that were classified into known category
  cat ${APPEND} \
  <(cat input.fasta | bioawk -c fastx '{ print $name }' |
  grep -wf <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  while read seq;
  do
  replace=`cat ${OUTPUT}/full_match/combined_classified_elements.txt | grep "${seq}" | cut -f 2`;
  seqkit grep -p ${seq} input.fasta | seqkit fx2tab |
  awk -F "\t|#" -v replace="${replace}" -v OFS="\t" '{ print $1"#"replace, $3 }' | seqkit tab2fx -w0;
  done) > ${OUTPUT}/${OUTPUT}.known

  # compress outputs to save space
  find ${OUTPUT} -name "*.cat" | while read file; do gzip -f ${file}; done
  find ${OUTPUT} -name "*.out" | while read file; do gzip -f ${file}; done
  find ${OUTPUT} -name "*.masked" | while read file; do gzip -f ${file}; done

  # clean up "round_N" symlink
  rm input.fasta

else

  # create output directories
  mkdir -p ${OUTPUT}/database_match ${OUTPUT}/known_match

  # create "input.fasta" symlink
  ln -s ${UNKNOWN} input.fasta

  # perform database repeat masking
  cmd1="RepeatMasker -pa ${THREAD} -e ncbi -nolow -dir ${OUTPUT}/database_match -species ${DATABASE} input.fasta"
  #echo ${cmd1}
  eval ${cmd1}

  # find FASTA for input of next search
  DBOUTPUT=`find ${OUTPUT}/database_match -name "*.masked"`

  # perform known repeat masking
  cmd2="RepeatMasker -pa ${THREAD} -e ncbi -nolow -dir ${OUTPUT}/known_match -lib ${KNOWN} ${DBOUTPUT}"
  #echo ${cmd2}
  eval ${cmd2}

  # create combined output directory
  mkdir -p ${OUTPUT}/full_match

  # combine all .out and .cat.gz files together
  # note that normal header from .out files is not retained!
  cat ${OUTPUT}/database_match/*.out ${OUTPUT}/known_match/*.out | awk '{ if ($1 != "SW" && $1 != "score") print $0 }' | sed '/^$/d' | sed '/^[[:space:]]*$/d' > ${OUTPUT}/full_match/${OUTPUT}.out

  cat ${OUTPUT}/database_match/*.cat ${OUTPUT}/known_match/*.cat > ${OUTPUT}/full_match/${OUTPUT}.cat

  # summarize repeat classifications
  # subfamily non-ambiguous elements - elements in the same subfamily, which means they are totally nonambiguous (attach full info)
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  cut -f 1 | sort | uniq |
  while read unknown;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${unknown}" | cut -f 4 | sort | uniq |
  awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR == 1) print unknown, $1 }';
  done > ${OUTPUT}/full_match/subfamily_unambiguous_classified_elements.txt

  # family non-ambiguous elements - elements in same family even if subfamilies are not the same (attach family info)
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  cut -f 1 | sort | uniq |
  while read unknown;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${unknown}" | cut -f 4 | sort | uniq |
  awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR > 1) print unknown }';
  done |
  while read element;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${element}" | cut -f 4 |
  awk -v OFS="\t"  -v element="${element}" -F "/|-" '{ print element, $1, $2 }';
  done |
  sort | uniq -c |
  awk '{ print $2 }' | uniq -u |
  while read clean;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${clean}" | cut -f 4 |
  awk -v OFS="\t"  -v clean="${clean}" -F "/|-" '{ print clean, $1"-"$2 }';
  done |
  uniq > ${OUTPUT}/full_match/family_unambiguous_classified_elements.txt

  # combine classified elements
  cat ${OUTPUT}/full_match/subfamily_unambiguous_classified_elements.txt ${OUTPUT}/full_match/family_unambiguous_classified_elements.txt > ${OUTPUT}/full_match/combined_classified_elements.txt

  # chimeric elements
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep -v -f <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  cut -f 1 | sort | uniq |
  while read element;
  do
  cat ${OUTPUT}/full_match/${OUTPUT}.out |
  awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' |
  awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' |
  sort -k1,1 -k5,5nr |
  grep "${element}" | cut -f 4 | paste -s -d "," |
  awk -v OFS="\t" -v element="${element}" '{ print element, $1 }';
  done > ${OUTPUT}/full_match/chimeric_elements.txt

  # create FASTA output of unknowns that remain unknown
  cat input.fasta | bioawk -c fastx '{ print $name }' |
  grep -wv -f <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  seqkit grep -f - input.fasta > ${OUTPUT}/${OUTPUT}.unknown

  # create FASTA output of unknowns that were classified into known category
  cat ${APPEND} \
  <(cat input.fasta | bioawk -c fastx '{ print $name }' |
  grep -wf <(cat ${OUTPUT}/full_match/combined_classified_elements.txt | cut -f 1) |
  while read seq;
  do
  replace=`cat ${OUTPUT}/full_match/combined_classified_elements.txt | grep "${seq}" | cut -f 2`;
  seqkit grep -p ${seq} input.fasta | seqkit fx2tab |
  awk -F "\t|#" -v replace="${replace}" -v OFS="\t" '{ print $1"#"replace, $3 }' | seqkit tab2fx -w0;
  done) > ${OUTPUT}/${OUTPUT}.known

  # compress outputs to save space
  find ${OUTPUT} -name "*.cat" | while read file; do gzip -f ${file}; done
  find ${OUTPUT} -name "*.out" | while read file; do gzip -f ${file}; done
  find ${OUTPUT} -name "*.masked" | while read file; do gzip -f ${file}; done

  # clean up "input.fasta" symlink
  rm input.fasta

fi
