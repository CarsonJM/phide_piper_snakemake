# Download informative vogdb HMM info from MiuVIG paper and save as a .tsv file
# Saved as vogdb_taxonomic_hmms.tsv
# Download vogdb HMM profiles
wget -O /home/carsonjm/resources/taxonomic_hmms/vogdb_lca.tsv.gz http://fileshare.csb.univie.ac.at/vog/latest/vog.lca.tsv.gz
cd /home/carsonjm/resources/taxonomic_hmms
gunzip vogdb_lca.tsv.gz


# Download viphog database
f http://ftp.ebi.ac.uk/pub/databases/metagenomics/vpHMM_database.tar.gz
cd /home/carsonjm/resources/taxonomic_hmms
tar -xvf vpHMM_database.tar.gz

# Download viphog metadata
wget -O /home/carsonjm/resources/taxonomic_hmms/additional_data_vpHMMs_v2.tsv http://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/additional_data_vpHMMs_v2.tsv
