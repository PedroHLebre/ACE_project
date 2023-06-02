#!/bin/bash

cd ~/normalization/results
for NAME in *_hmmsearch_result
do
	CURRENT=$(echo $NAME | grep -Eo ^S[0-9]+)
	Rscript ~/normalization/scripts/Normalisation.R $NAME ../${CURRENT}_filtered_R1_R2.fa ACE_key_file 1e-10 ${CURRENT}_ACE_normalisation_result
done
