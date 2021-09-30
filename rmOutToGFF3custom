#!/usr/bin/env bash
usage()
{
cat << EOF
rmOutToGFF3custom
Version 1.1 (2021-09-30)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script converts the .out file from RepeatMasker to a GFF3 file. Note that the output
is probably not perfect GFF3, so beware with downstream applications. This script
emulates the rmOutToGFF3.pl script supplied with RepeatMasker but provides a fuller ID
("target=") for each element in column 9 of the GFF. This ID includes the matching element,
like rmOutToGFF3.pl, but also includes the repeat family: in the format <Family>/<Element>.
This change is because many matching elements produced from RepeatModeler have IDs that
provide no information about repeat family classification. Output is written to standard
output (SDOUT).

This script requires requires awk, which should be available on any standard Unix system.

rmOutToGFF3custom -o <RM.out> [-h] > <name.gff3>

OPTIONS:
        -h		usage information and help (this message)
	-o		RepeatMasker .out file
EOF
}

while getopts "ho:" OPTION
do
        case $OPTION in
                help)
                        usage
                        exit 1
                        ;;
		o)
			RMOUT=$OPTARG
			;;
	esac
done

if [[ -z $RMOUT ]]
then
	usage
	exit 1
fi

cat <(echo "##gff-version 3") \
<(cat ${RMOUT} | tail -n +4 | \
awk -v OFS="\t" '{ if ($12 ~ /)/) print $5, "RepeatMasker", "dispersed_repeat", $6, $7, $1, $9, ".", "Target="$11"/"$10" "$14" "$13; \
else print $5, "RepeatMasker", "dispersed_repeat", $6, $7, $1, $9, ".", "Target="$11"/"$10" "$12" "$13 }' | \
awk -v OFS="\t" '{ if ($7 == "C") print $1, $2, $3, $4, $5, $6, "-", $8, $9; else print $0 }' | \
sort -k1,1 -k4,4n -k5,5n)
