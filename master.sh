#!/bin/bash

#STEP ONE
bash prepare_data.sh

#STEP TWO
python get_transcripts.py

#STEP THREE
Rscript simulate_sample.R

#STEP FOUR
bash do_quantification.sh

#STEP FIVE
RScript compare_quantification.R
