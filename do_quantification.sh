#!/bin/bash

#index the transcriptome for both algorithms
salmon index -p 20 -i salmon_index -t data/select_transcripts.fa
kallisto index -i kallisto_index data/select_transcripts.fa

for i in {1..10}
    do
        #run salmon
        salmon quant -p 20 -i salmon_index --biasCorrect -l 'IU' \
            -1 data/sample_01_1.fasta -2 data/sample_01_2.fasta \
            -o sim_test_salmon_${i}
    done

#run kallisto
kallisto quant -i kallisto_index -b 10 -o sim_test_kallisto \
    data/sample_01_1.fasta data/sample_01_2.fasta
#extract the kallisto bootstraps
kallisto h5dump -o sim_test_kallisto_bs sim_test_kallisto/abundance.h5
