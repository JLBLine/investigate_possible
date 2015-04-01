#!/bin/sh
./investigate_possibles.py --matched_cats=comp_v5,vlssr,mrc,sumss,nvss \
	--pref_cats=nvss,sumss,comp_v1,mrc,vlssr \
	--cat_freqs=182.435,74,408,843,1400 \
	--input_bayes=KGS_compv5-v-m-s-n_split-eyeball.txt \
	--output=extras_test