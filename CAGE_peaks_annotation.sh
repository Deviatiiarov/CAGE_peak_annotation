#!/bin/bash

# cut -f2,13 ./example/ncbiRefSeq.txt > ./example/ncbiRefSeqToGeneName.txt
# awk 'BEGIN{OFS="\t"};{print $4, $4}'  ./example/Entrez_gene_galGal6.bed > ./example/Entrez_gene_galGal6_names.txt
# echo -e ./example/augustusGene.txt'\n'./example/ensGene.txt'\t'./example/ensemblToGeneName.txt'\n'./example/ncbiRefSeq.txt'\t'./example/ncbiRefSeqToGeneName.txt'\n'./example/genscan.txt > ./example/ucsc_path.txt
# echo -e ./example/Entrez_gene_galGal6.bed'\t'./example/Entrez_gene_galGal6_names.txt > ./example/bed_path.txt

# chmod a+x CAGE_peaks_annotation.sh
# ./CAGE_peaks_annotation.sh -c ./example/rDPI_merged_chicken6.bed -b ./example/bed_path.txt -u ./example/ucsc_path.txt -o ./example/output -t ./example/temp -l 1000 -a Y

# sed '/^#/ d' /data5/ruslan/chicken/gg6_mapped/level2_merge/gg6_merged_level2.osc | tail -n +2 | awk 'BEGIN {OFS="\t"};{print $2, $3, $4, $1, $7,$5, $6}' > /data5/ruslan/chicken/gg6_mapped/level2_merge/gg6_merged_level2.bed
# ./CAGE_peaks_annotation.sh -c /data5/ruslan/chicken/gg6_mapped/level2_merge/gg6_merged_level2.bed -b ./example/bed_path.txt -u ./example/ucsc_path.txt -o /data5/ruslan/chicken/gg6_mapped/ann_auto/out2 -t ./example/temp -l 1000 -a Y


echoerr() { echo "$@" 1>&2; }

usage() {
	echoerr
	echoerr ${0##/*}
	echoerr "  Ruslan D (2019), ruselusalbus@gmail.com"
	echoerr "  CAGE peaks annotation"
	echoerr "Usage:"
	echoerr "  ${0##/} OPTIONS -c FILE -u FILE -b FILE -o PATH"
	echoerr
	echoerr "Required:"
    echoerr "  [-c PATH]     A single BED6 file for CAGE clusters"
	echoerr "  [-u PATH]     A single file with paths to ucsc files"
	echoerr "  [-o PATH]     Path where to put the output files"
	echoerr "  [-b BED]      A single file with paths to BED6 files"
	echoerr
	echoerr "  -u and -b files could be a tab separated files with two"
    echoerr "  columns, where second column provides paths to transcripts"
	echoerr "  annotation files, for example:"
    echoerr
    echoerr "  XM_025152649.1  LOC112532827"
	echoerr "  XM_025152753.1  LOC107049475"
	echoerr "  XM_025152751.1  LOC107049475"
    echoerr 
	echoerr "OPTIONS:"
	echoerr 
	echoerr "  [-t PATH]     Temp directory"
	echoerr "  [-l NUM]      Length of up and downstream extension area"
	echoerr "                For example 200,1000,1200"
	echoerr "                Default is 50 and 500 bp"
	echoerr "  [-a y/n]      Write input files for TssClassifier"
	echoerr "                Default is N"
	echoerr
	exit
}

#BINDIR="`dirname "$0"`"/bin
#if [ ! -d $BINDIR ]; then
#	echoerr "cannot locate bin directory (../bin relative the path to this script)"
#	exit
#fi

#echoerr $BINDIR

## show usage if '-h' or  '--help' is the first argument or no argument is given
case $1 in
	""|"-h"|"--help") usage ;;
esac

TOME="N"
TOME_input=("y" "Y" "n" "N")
while getopts c:u:o:b:t:l:a: opt; do
	case ${opt} in
		c) CAGE_BED=${OPTARG};;
        u) FILES=${OPTARG};;
		o) OPATH=${OPTARG};;
		b) BED=${OPTARG};;
		t) TEMP=${OPTARG};;
		l) IFS=',' read -a LEN <<< "$OPTARG" ;;
		a) TOME=${OPTARG};;
		*) usage;;
	esac
done

RMTMP=1

## check the parameters
if [ "$CAGE_BED" == "" ]; then echoerr "c parameter is required"; usage; fi
if [ "$CAGE_BED" != "" ] && [ ! -e $CAGE_BED ] ; then echoerr "No such file ${CAGE_BED}"; usage; fi

if [ "$FILES" == "" ] && [ "$BED" == "" ]; then echoerr "f or b parameter is required"; usage; fi
if [ "$FILES" != "" ] && [ ! -e $FILES ] ; then echoerr "No such file ${FILES}"; usage; fi
if [ "$BED" != "" ] && [ ! -e $BED ]; then echoerr "No such file ${BED}"; usage; fi

if [ "$OPATH" == "" ]; then echoerr "o parameter is required"; usage; fi
if [ "$OPATH" != "" ] && [ ! -d $OPATH ]; then mkdir $OPATH; fi
if [ "$TEMP" == "" ]; then TEMP=`mktemp -d`; fi
if [ "$TEMP" != "" ]; then RMTMP=0; fi
if [ "$TEMP" != "" ] && [ ! -d $TEMP ]; then mkdir $TEMP; fi

if  [[ ! $OPATH == */ ]]; then
	OPATH=${OPATH}/
fi
if  [[ ! $TEMP == */ ]]; then
	TEMP=${TEMP}/
fi

if [ "$LEN" != "" ] && [[ ! $LEN =~ ^[0-9]+$ ]]; then echoerr "length must be an integer"; usage; fi
if [[ ! " ${TOME_input[@]} " =~ " ${TOME} " ]]; then echoerr "a must be y/n or Y/N"; usage; fi

echoerr ${TEMP}
echoerr "check #1 passed"
n_files=0
n_bed=0
number_of_ann=0
array=(50 500)


if [ "$LEN" != "" ] && [[ $LEN =~ ^[0-9]+$ ]]; then array2=("${array[@]}" "${LEN[@]}"); fi
if [ "$LEN" == "" ]; then array2=${array[@]}; fi

final_array=$(echo ${array2[*]}| tr " " "\n" | sort -n | tr "\n" " ")
echoerr ${final_array[@]}
MAX=$(echo ${final_array[*]}| tr " " "\n" | sort -n | tail -n 1)
echoerr $MAX

awk 'BEGIN{OFS="\t"};{$2=$7;$3=$7+1; print}' $CAGE_BED > ${TEMP}CAGE_1bp_BED.txt # change DPI tab to 1 bp size!
awk 'BEGIN{OFS="\t"};{$2=$7-50;$3=$7+50; print}' $CAGE_BED | sort -k1,1 -k2,2n | awk '$2>0'  > ${TEMP}CAGE_100bp_BED.txt
CAGE_BED_DIM=$(awk '{print NF}' $CAGE_BED | sort -nu | tail -n 1)
echoerr $CAGE_BED_DIM

CAGE_BED_OUT=$(basename -- ${CAGE_BED})

awk 'BEGIN {OFS="\t"};{print $4, $1, $2, $3, $6, $7}' $CAGE_BED | sort -k 1b,1 > ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i '1i CLUSTER_ID\tCHROMOSOME\tSTART\tEND\tSTRAND\tTSS' ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt

### UCSC
number_of_file=0
if [ "$FILES" != "" ] && [ -e $FILES ] ; then
n_files=$(wc -l $FILES)
for FILE in  $(cut -f1 $FILES ); do
((number_of_file+=1))
# loop N2 for numbers
for n_len in ${final_array[@]}; do 
# prepare input file (UCSC)
RES_FILE=${TEMP}$(basename -- ${FILE%.txt*}_${n_len}bp.txt)
awk -v len=$n_len 'BEGIN{OFS="\t"};{if($4=="+"){print $3, ($5-len), ($5 + len),$2,$1,$4,$5,$6}}; {if($4=="-"){print $3, ($6-len), ($6+len), $2,$1,$4,$5,$6}}'  $FILE | awk 'BEGIN{OFS="\t"};{if($2 < 0) $2=0;print}' | sort -k1,1 -k2,2n > $RES_FILE; 
# intersection
echoerr $RES_FILE
bedtools intersect -s -wa -wb -a ${TEMP}CAGE_1bp_BED.txt -b $RES_FILE > ${RES_FILE%.txt*}_CAGEpeaks.txt
RES_FILE_DIM=$(awk '{print NF}' $RES_FILE | sort -nu | tail -n 1)
all_col=$((CAGE_BED_DIM + RES_FILE_DIM))
column_of_gene_start=$((CAGE_BED_DIM+RES_FILE_DIM-1))
column_of_gene_end=$((CAGE_BED_DIM+RES_FILE_DIM))
awk '$6=="+"' ${RES_FILE%.txt*}_CAGEpeaks.txt | awk -v n_col=$column_of_gene_start 'BEGIN{OFS="\t"};{print $0, $n_col - $7}' | awk -v col=$((all_col+1)) 'BEGIN{OFS="\t"};{gsub("-", "", $col)}1' | sort --key $((all_col+1)) -n | awk '!seen[$4]++' | awk -v end=$column_of_gene_end '$7<$end' >  ${RES_FILE%.txt*}_CAGEpeaks_plus.txt
awk '$6=="-"' ${RES_FILE%.txt*}_CAGEpeaks.txt | awk -v n_col=$column_of_gene_end 'BEGIN{OFS="\t"};{print $0, $7 - $n_col}' | awk -v col=$((all_col+1)) 'BEGIN{OFS="\t"};{gsub("-", "", $col)}1' | sort --key $((all_col+1)) -n | awk '!seen[$4]++' | awk -v start=$column_of_gene_start '$7>$start' >  ${RES_FILE%.txt*}_CAGEpeaks_minus.txt
cat ${RES_FILE%.txt*}_CAGEpeaks_plus.txt ${RES_FILE%.txt*}_CAGEpeaks_minus.txt | cut -f 4,$((CAGE_BED_DIM+4)),$((all_col+1)) | sort -k 1b,1 > ${RES_FILE%.txt*}_CAGEpeaks.txt 
# add gene id
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' ${RES_FILE%.txt*}_CAGEpeaks.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp.txt;
header=$(awk -v name=$(basename -- ${RES_FILE%.txt*}) 'NR==1 {print ($0, name)}' ${TEMP}Annotation_tmp.txt);
mv ${TEMP}Annotation_tmp.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i "1s/.*/$header/" ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
# add dist
if [ "$n_len" == "$MAX"  ]; then 
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$3;next}{ print $0, a[$1]}' ${RES_FILE%.txt*}_CAGEpeaks.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp.txt;
header=$(awk -v name=$(basename -- ${RES_FILE%_*}_dist) 'NR==1 {print ($0, name)}' ${TEMP}Annotation_tmp.txt);
mv ${TEMP}Annotation_tmp.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i "1s/.*/$header/" ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
# check annotation
if [ "$(awk -v n_row=$number_of_file 'NR == n_row {print $2}' $FILES)" != "" ]; then
ANN_FILE=$(awk -v n_row=$number_of_file 'NR == n_row {print $2}' $FILES)
Annotated_DIM=$(( $(awk '{print NF}' ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt | sort -nu | tail -n 1) - 1 ))
cut -f 1,$Annotated_DIM ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt |  awk 'NF==2' | awk 'BEGIN{OFS="\t"};{print $2, $1}' > ${RES_FILE%.txt*}_0_TempName.txt
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' $ANN_FILE  ${RES_FILE%.txt*}_0_TempName.txt  | awk 'BEGIN{OFS="\t"};{print $2,$3,$1}' | awk 'NF==3{print $0}'  > ${TEMP}${number_of_file}_$(basename -- ${FILE%.txt*}_${n_len}bp)_GeneName.txt
fi;
fi;
done;
done; 
fi

#echoerr $( wc -l ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt)

### BED
unset  number_of_file
if [ "$BED" != "" ] && [ -e $BED ] ; then
n_bed=$(wc -l $BED)
for FILE in  $(cut -f1 $BED ); do 
((number_of_file+=1))
# loop N2 for numbers
for n_len in ${final_array[@]}; do 
# prepare input file (BED6)
RES_FILE=${TEMP}$(basename -- ${FILE%.bed*}_${n_len}bp.txt)
awk -v len=$n_len 'BEGIN{OFS="\t"};{if($6=="+"){print $1, ($2-len), ($2 + len),$4,$5,$6,$2,$3}}; {if($6=="-"){print $1, ($3-len), ($3+len), $4,$5,$6,$2,$3}}'  $FILE | awk 'BEGIN{OFS="\t"};{if($2 < 0) $2=0;print}' | sort -k1,1 -k2,2n > $RES_FILE; 
echoerr $RES_FILE
# intersection 
bedtools intersect -s -wa -wb -a ${TEMP}CAGE_1bp_BED.txt -b $RES_FILE > ${RES_FILE%.txt*}_CAGEpeaks.txt
RES_FILE_DIM=$(awk '{print NF}' $RES_FILE | sort -nu | tail -n 1)
all_col=$((CAGE_BED_DIM + RES_FILE_DIM))
column_of_gene_start=$((CAGE_BED_DIM+RES_FILE_DIM-1))
column_of_gene_end=$((CAGE_BED_DIM+RES_FILE_DIM))
awk '$6=="+"' ${RES_FILE%.txt*}_CAGEpeaks.txt | awk -v n_col=$column_of_gene_start 'BEGIN{OFS="\t"};{print $0, $n_col - $7}' | awk -v col=$((all_col+1)) 'BEGIN{OFS="\t"};{gsub("-", "", $col)}1' | sort --key $((all_col+1)) -n | awk '!seen[$4]++' | awk -v end=$column_of_gene_end '$7<$end' >  ${RES_FILE%.txt*}_CAGEpeaks_plus.txt
awk '$6=="-"' ${RES_FILE%.txt*}_CAGEpeaks.txt | awk -v n_col=$column_of_gene_end 'BEGIN{OFS="\t"};{print $0, $7 - $n_col}' | awk -v col=$((all_col+1)) 'BEGIN{OFS="\t"};{gsub("-", "", $col)}1' | sort --key $((all_col+1)) -n | awk '!seen[$4]++' | awk -v start=$column_of_gene_start '$7>$start' >  ${RES_FILE%.txt*}_CAGEpeaks_minus.txt
cat ${RES_FILE%.txt*}_CAGEpeaks_plus.txt ${RES_FILE%.txt*}_CAGEpeaks_minus.txt | cut -f 4,$((CAGE_BED_DIM+4)),$((all_col+1)) | sort -k 1b,1 > ${RES_FILE%.txt*}_CAGEpeaks.txt 
# add gene id
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' ${RES_FILE%.txt*}_CAGEpeaks.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp.txt;
header=$(awk -v name=$(basename -- ${RES_FILE%.txt*}) 'NR==1 {print ($0, name)}' ${TEMP}Annotation_tmp.txt);
mv ${TEMP}Annotation_tmp.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i "1s/.*/$header/" ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
# add dist
if [ "$n_len" == "$MAX"  ]; then 
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$3;next}{ print $0, a[$1]}' ${RES_FILE%.txt*}_CAGEpeaks.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp.txt;
header=$(awk -v name=$(basename -- ${RES_FILE%_*}_dist) 'NR==1 {print ($0, name)}' ${TEMP}Annotation_tmp.txt);
mv ${TEMP}Annotation_tmp.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i "1s/.*/$header/" ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
# check annotation
if [ "$(awk -v n_row=$number_of_file 'NR == n_row {print $2}' $BED)" != "" ]; then
ANN_FILE=$(awk -v n_row=$number_of_file 'NR == n_row {print $2}' $BED)
Annotated_DIM=$(( $(awk '{print NF}' ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt | sort -nu | tail -n 1) - 1 ))
cut -f 1,$Annotated_DIM ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt |  awk 'NF==2' | awk 'BEGIN{OFS="\t"};{print $2, $1}' > ${RES_FILE%.txt*}_0_TempName.txt
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' $ANN_FILE ${RES_FILE%.txt*}_0_TempName.txt | awk 'BEGIN{OFS="\t"};{print $2,$3,$1}' | awk 'NF==3{print $0}'  > ${TEMP}${number_of_file}_$(basename -- ${FILE%.bed*}_${n_len}bp)_GeneName.txt
fi;
fi;
done;
done;
fi

#echoerr $( wc -l ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt)

if [ "$(awk 'NF==2' $FILES | wc -l )" > 0 ] || [ "$(awk 'NF==2' $BED | wc -l )" > 0  != "" ]; then
Annotated_DIM=$(awk '{print NF}' ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt | sort -nu | tail -n 1)

for ann in ${TEMP}*_GeneName.txt; do
((number_of_ann+=1))
awk 'BEGIN {OFS="\t"};FNR==NR{a[$1]=$2;next}{ print $0, a[$1]}' $ann ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp.txt;
header=$(awk -v name="GeneName" 'NR==1 {print ($0, name)}' ${TEMP}Annotation_tmp.txt);
mv ${TEMP}Annotation_tmp.txt ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
sed  -i "1s/.*/$header/" ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt
done;

Annotated_DIM2=$((Annotated_DIM+number_of_ann))      #######################################################
cut -f $((Annotated_DIM+1))-$Annotated_DIM2 ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt | sed 's/\t/;/g' | sed -e 's/\(;\)*$//g' |   awk 'BEGIN{OFS="\t"};{gsub(".*;","")};{print $0}' > ${TEMP}Annotation_tmp.txt
cut -f 1-$Annotated_DIM ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt > ${TEMP}Annotation_tmp2.txt # 
paste -d '\t' ${TEMP}Annotation_tmp2.txt ${TEMP}Annotation_tmp.txt > ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt

#echoerr $( wc -l ${OPATH}${CAGE_BED_OUT%.bed*}_annotated.txt)
fi

##################################


for i in ${TEMP}*_50bp_CAGEpeaks.txt; do
 cut -f1 $i >> ${TEMP}Pos_TSS_id.txt;
done

case ${TOME} in
  "n"|"N")
    (
	echoerr "skipping TssClass files writing"
       );;
  "y"|"Y")
    (
	sort -u ${TEMP}Pos_TSS_id.txt | grep -wFf - ${TEMP}CAGE_100bp_BED.txt > ${OPATH}${CAGE_BED_OUT%.bed*}_100bp_pos.bed
	sort -u ${TEMP}Pos_TSS_id.txt | grep -vwFf - ${TEMP}CAGE_100bp_BED.txt > ${OPATH}${CAGE_BED_OUT%.bed*}_100bp_neg.bed
        );;
  *)
    usage;;
esac


if [ $RMTMP -eq 1 ]; then rm -rf $TEMP; fi
rm -rf 0 # wth is this zero file?

# END