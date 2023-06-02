#!/bin/bash

cd ~/normalization
for NAME in *_filtered_R1_R2.fa
do
	CURRENT=$(echo $NAME | grep -Eo ^S[0-9]+)
	transeq -sequence $NAME -outseq ${CURRENT}_subsampled_aa.fa -frame 6 -clean >> transeq.log
	hmmsearch --tblout ${CURRENT}_hmmsearch_result -E 1e-06 --cpu 8 ACE_hmm_profiles.hmm ${CURRENT}_subsampled_aa.fa
	rm ${CURRENT}_subsampled_aa.fa
	mv ${CURRENT}_hmmsearch_result ./results/
done
