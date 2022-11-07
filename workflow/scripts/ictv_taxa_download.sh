#!/usr/bin/env bash

cd ${snakemake_params[ictv_dir]}

aria2c -x 4 -o ictv.xlsx "https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/13426/download"

# convert xlsx to tsv
csvtk xlsx2csv ictv.xlsx \
    | csvtk csv2tab \
    | sed 's/\xc2\xa0/ /g' \
    | csvtk replace -t -F -f "*" -p "^\s+|\s+$" \
    > ictv.tsv

# choose columns, and remove duplicates
csvtk cut -t -f "Realm,Subrealm,Kingdom,Subkingdom,Phylum,Subphylum,Class,Subclass,Order,Suborder,Family,Subfamily,Genus,Subgenus,Species" ictv.tsv \
    | csvtk uniq -t -f "Realm,Subrealm,Kingdom,Subkingdom,Phylum,Subphylum,Class,Subclass,Order,Suborder,Family,Subfamily,Genus,Subgenus,Species" \
    | csvtk del-header -t \
    > ictv.taxonomy.tsv

# create a file that will store all the ICTV taxa names
csvtk cut -t -H -f 1,3,5,7,9,11,13,15 ictv.taxonomy.tsv \
    | sed 's/\t/\n/g' \
    | awk '!/^[[:blank:]]*$/' \
    | sort -u \
    > ictv.names.txt

# Download the NCBI taxdump
aria2c -x 4 "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
mkdir ncbi-taxdump
tar zxfv taxdump.tar.gz -C ncbi-taxdump
rm taxdump.tar.gz

# use taxonkit create-taxdump to create a custome taxdump for ICTV
taxonkit create-taxdump -K 1 -P 3 -C 5 -O 7 -F 9 -G 11 -S 13 -T 15 \
--rank-names "realm","kingdom","phylum","class","order","family","genus","species" ictv.taxonomy.tsv \
--out-dir ictv-taxdump --data-dir ncbi-taxdump

# execute fix_taxdump.py to make taxids sequential for compatibility with MMseqs2
python ./fix_taxdump.py

# Download the protein → taxid association and filter for viruses
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz"

gunzip prot.accession2taxid.FULL.gz

# filter to keep only viral proteins
awk '{print $2}' prot.accession2taxid.FULL \
    | sort -u \
    | taxonkit --data-dir ncbi-taxdump lineage \
    | rg "\tViruses;" \
    | awk '{print $1}' \
    > virus_taxid.list

csvtk grep -t -f 2 -P virus_taxid.list prot.accession2taxid.FULL > virus.accession2taxid

rm prot.accession2taxid.FULL

# Find the ICTV-compliant proteins and write a new table with the ICTV taxids
./get_ictv_taxids.py

# Download and filter NR proteins
# aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz"
aria2c -x 4 "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.faa.gz"

cat viral.*.protein.faa.gz > nr.gz

# Create a list containing the accessions of the proteins of ICTV viruses
cut -f 1 virus.accession2taxid.ictv > virus.accession.txt

# Filter the NR proteins to keep the proteins encoded by ICTV viruses
seqkit grep -j 4 -f virus.accession.txt nr.gz | seqkit seq -i -w 0 -o nr.virus.faa.gz

rm nr.gz

# Filter the NR virus taxid table
seqkit fx2tab -n -i nr.virus.faa.gz > nr.virus.list.txt
csvtk grep -t -H -f 1 -P nr.virus.list.txt virus.accession2taxid.ictv > nr.virus.accession2taxid.ictv

