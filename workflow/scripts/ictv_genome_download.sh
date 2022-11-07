#!/usr/bin/env bash

cd ${snakemake_params[ictv_dir]}

# Download the protein → taxid association and filter for viruses
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.EXTRA.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"

cat nucl_gb.accession2taxid.gz nucl_wgs.accession2taxid.EXTRA.gz nucl_wgs.accession2taxid.gz > nucl.accession2taxid.FULL.gz
gunzip nucl.accession2taxid.FULL.gz

# filter to keep only viral proteins
awk '{print $2}' nucl.accession2taxid.FULL \
    | sort -u \
    | taxonkit --data-dir ncbi-taxdump lineage \
    | rg "\tViruses;" \
    | awk '{print $1}' \
    > virus_taxid.list

csvtk grep -t -f 2 -P virus_taxid.list nucl.accession2taxid.FULL > virus.accession2taxid

rm nucl.accession2taxid.FULL

# Find the ICTV-compliant genomes and write a new table with the ICTV taxids
./get_ictv_taxids.py

# Download and filter NR proteins
# aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz"

cat viral.*.1.genomic.fna.gz > nt.gz

# Create a list containing the accessions of the proteins of ICTV viruses
cut -f 1 virus.accession2taxid.ictv > virus.accession.txt

# Filter the NR proteins to keep the proteins encoded by ICTV viruses
seqkit grep -j 4 -f virus.accession.txt nt.gz | seqkit seq -i -w 0 -o nt.virus.fna.gz

rm nr.gz

# Filter the NR virus taxid table
seqkit fx2tab -n -i nt.virus.fna.gz > nr.virus.list.txt
csvtk grep -t -H -f 1 -P nr.virus.list.txt virus.accession2taxid.ictv > nr.virus.accession2taxid.ictv
