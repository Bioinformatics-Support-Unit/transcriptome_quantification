from Bio import SeqIO
import re

fh = open('data/select_genes.txt', 'r')
genes = fh.readlines()
genes = map(str.rstrip, genes)
records = SeqIO.parse('Homo_sapiens.GRCh38.rel79.cdna.all.fa', 'fasta')
keep_records = []
for record in records:
    desc = record.description
    pattern = '(ENSG\d{11})'
    m = re.search(pattern, desc)
    gene_id = m.groups()[0]
    if gene_id in genes:
        keep_records.append(record)

fh.close()
out = open('data/select_transcripts.fa', 'w')
SeqIO.write(keep_records, out, 'fasta')
out.close()
