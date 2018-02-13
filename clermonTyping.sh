#!/usr/bin/env bash

#################################################
########## Clermont Typing pipeline #############
#################################################
# From a set a contigs in fasta format:
# 1] Launch mash for getting phylogroup
# 2] Make a blast db
# 3] Launch blast on a primers fasta file
# 4] Launch in silicco PCR for getting phylogroup
# 5] Reportings tools and send e-mail to user
#
# Current version : 1.0.0 (Oct. 2017)
#
# Contact: johann.beghain@inserm.fr

MY_PATH="`dirname \"$0\"`"
#Default threshold = 0 (disabled)
THRESHOLD=0
#Default name = date
DATE=$( date "+%F_%H%M%S")
NAME=analysis_$DATE
#Global variables
THREADS=4
#BLAST settings
PRIMERS="${MY_PATH}/data/primers.fasta"
PERC_IDENTITY=90
BLAST_TASK='blastn'
#MASH settings
DB_MASH="${MY_PATH}/data/mash/mash_reference.msh"

function usage(){
	printf "Script usage :\n"
	printf "\t-h					: print this message and exit\n"
	printf "\t--fasta					: fasta contigs name(s). Can be separated by an arobase (@) value\n"
	printf "\t--name					: name for this analysis (optional)\n"
	printf "\t--threshold				: Option for ClermontTyping, do not use contigs under this size (optional)\n"
}

function mash_analysis(){
	echo "============== Running mash ================"
	${MY_PATH}/bin/mash screen -w $DB_MASH $FASTA >$WORKING_DIR/${FASTA_NAME}_mash_screen.tab
	#Old mash analysis
	#mash sketch -p 3 -k 32 -s 5000 -o mash_strain $FASTA
	#mash paste mash_matrix.msh $DB_MASH mash_strain.msh
	#mash dist -t mash_matrix.msh mash_matrix.msh >mash_matrix.txt
	#result=`./mash_get_phylogroup.py --genome $file --matrix mash_matrix.txt --annotation $ANNOTATION_FILE`
	#echo "$NAME	$result"
}

function blast_analysis(){
	echo "============== Making blast db ================"
	echo "makeblastdb -in $FASTA -input_type fasta -out $WORKING_DIR/db/$NAME -dbtype nucl"
	makeblastdb -in $FASTA -input_type fasta -out $WORKING_DIR/db/$FASTA_NAME -dbtype nucl
	echo "============== Running blast =================="
	blastn -query $PRIMERS -perc_identity $PERC_IDENTITY -task $BLAST_TASK -word_size 6 -outfmt 5 -db $WORKING_DIR/db/$FASTA_NAME -out $WORKING_DIR/$FASTA_NAME.xml
}

function report_calling(){
	# rscript = path to clermontReport.R
	# clermont_out = path to clermonTyping output 
	# namus = report name
	# out_dir = self explanatory!
	echo "============= Generating report ==============="
	rscript=$1
	shift
	clermont_out=$1
	shift
	namus=$1
	shift
	out_dir=$1

	# echo "$rscript ; $clermont_out ; $namus ; $out_dir"

	modif_script=${out_dir}/${namus}.R
	cp ${rscript} ${modif_script}

	sed -i "s:TARTAMPION:$clermont_out:g" "${modif_script}"

	Rscript --slave -e "library(markdown); sink('/dev/null');rmarkdown::render('${modif_script}')"
}

OPTS=$( getopt -o h -l fasta,threshold,name: -- "$@" )
if [ $? != 0 ]
then
    exit 1
fi

while [[ $# -gt 1 ]]
do
    case "$1" in
        -h)
        usage
        exit 0
        ;;
        --fasta) 
		FASTAS="$2";
        shift
        ;;
        --name) 
		NAME="$2";
        shift
        ;;
        --threshold) 
		THRESHOLD="$2";
        shift
        ;;
        --) shift; break;;
    esac
    shift
done

if [ -z $FASTAS ]
then
	echo "Missing the contigs file. Option --fasta"
	usage
	exit 1
fi

echo "You asked for a Clermont typing analysis named $NAME of phylogroups on $FASTAS with a threshold under $THRESHOLD."

if [ ! -d $NAME ]
then
	mkdir $NAME
fi
CURRENT_DIR=`pwd`
WORKING_DIR=$CURRENT_DIR/$NAME
#Analysis of each fasta file
IFS='@' read -ra ARRAY_FASTA <<< "$FASTAS"
for FASTA in "${ARRAY_FASTA[@]}"; do
	if [ -f $FASTA ] && [ ! -z $FASTA ]
	then
		#Rename file
		BASE_NAME_FASTA=`basename $FASTA`
		FASTA_NAME=${BASE_NAME_FASTA%£*}
		echo "Analysis of ${FASTA_NAME}"
		cp $FASTA $WORKING_DIR/${FASTA_NAME}
		FASTA=$WORKING_DIR/$FASTA_NAME
		##### Step 1: MASH analysis #####
		# Generate ${FASTA_NAME}_mash_screen.tab
		mash_analysis
		##### Step 2: Blast #############
		# Generate ${FASTA_NAME}.xml
		blast_analysis
		##### Step 3: ClermonTyping #####
		echo "============== ClermonTyping =================="
		results=`${MY_PATH}/bin/clermont.py -x ${WORKING_DIR}/${FASTA_NAME}.xml -s $THRESHOLD`
		echo "$FASTA_NAME	$results	${FASTA_NAME}_mash_screen.tab" >> $WORKING_DIR/${NAME}_phylogroups.txt
	else
		echo "$FASTA doesn't exists"
	fi
done
##### Step 4: Reporting #########
report_calling "${MY_PATH}/bin/clermontReport.R" $WORKING_DIR/${NAME}_phylogroups.txt $NAME $WORKING_DIR


exit 0
