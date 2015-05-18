#!/bin/bash
grep -Po 'ENSG\d{11}' extdata/Homo_sapiens.GRCh38.rel79.cdna.all.fa | sort | uniq > data/gene_ids.txt
#RANGE is how many gene ids in transcriptome fasta
RANGE=$(wc -l data/gene_ids.txt | cut -f1 -d ' ')

for i in {1..250}
    do
        #Select a random gene
        number=$(perl -e "print int(rand(${RANGE}))")
        head -n ${number} data/gene_ids.txt | tail -1 >> data/select_genes.txt
    done
