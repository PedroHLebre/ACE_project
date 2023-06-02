#!/bin/bash
cd ~/normalization
#2019 sequences
for NAME in ~/ACE_2019/*R1*fastq
do
	CURRENT=$(echo $NAME | sed s/"\/home\/aceproject\/ACE_2019\/"/""/ )
	CURRENT=$(echo $CURRENT | grep -Eo ^S[0-9]+ )
	/home/databases/jamstec/bbmap/reformat.sh in1=~/ACE_2019/${CURRENT}_R1_trim.fastq in2=~/ACE_2019/${CURRENT}_R2_trim.fastq out=${CURRENT}_filtered_R1_R2.fa minlength=140
done
#2020 sequences
for NAME in ~/acetrimmedfastqc/*R1*fastq
do
	FILEA=$(echo $NAME)
	FILEB=$(echo $NAME | sed s/"R1.trim"/"R2.trim"/ )
	CURRENT=$(echo $NAME | grep -Eo S[0-9]+)
	/home/databases/jamstec/bbmap/reformat.sh in1=$FILEA in2=$FILEB out=${CURRENT}_filtered_R1_R2.fa minlength=140
done
