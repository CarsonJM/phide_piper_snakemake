# make directory for spacer files
mkdir /home/carsonjm/resources/crispr_spacers

# download uhgg crispr spacers from ftp site
wget -O /home/carsonjm/resources/crispr_spacers/uhgg_spacers.fna https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/uhgg_spacers.fna

# download uhgg metadata
wget -O /home/carsonjm/resources/crispr_spacers/uhgg_spacers_metadata.tsv http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv

# download crispropendb spacers
wget -O /home/carsonjm/resources/crispr_spacers/PhageHostIdentifier_DBfiles.zip http://crispr.genome.ulaval.ca/dash/PhageHostIdentifier_DBfiles.zip
unzip /home/carsonjm/resources/crispr_spacers/PhageHostIdentifier_DBfiles.zip

# Downloaded crispropendb metadata from http://crispr.genome.ulaval.ca/ by
# 1. Opening "Download" tab
# 2. Selecting "Complete version with 11 767 783 spacers" from "Download the database" dropdown
# 3. Selecting "CSV" from "Download format" dropdown
# 4. Pressing "Prepare file" button
# 5. Pressing "Download" button
# 6. Moved file to /home/carsonjm/resources/crispr_spacers/SpacersDB.csv.zip
unzip /home/carsonjm/resources/crispr_spacers/SpacersDB.csv.zip -d /home/carsonjm/resources/crispr_spacers/

