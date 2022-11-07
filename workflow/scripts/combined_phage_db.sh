### Benler 2020 (https://doi.org/10.1186/s40168-021-01017-w) ###
# download benler phage genomes (from https://ftp.ncbi.nih.gov/pub/yutinn/benler_2020/gut_phages/)
cd /gscratch/scrubbed/carsonjm/resources/benler_phage

cat phage_genomes | xargs -n 1 -P 10 wget -P phages -q -w 5 -T 5
cat phages/* > benler_2020_phages.fna

# all phages are complete so this database is ready


### Nayfach 2020 (https://doi.org/10.1186/s40168-021-01017-w)
# checkv database downloaded by CheckV
# all phages are complete so this database is ready


### Tisza 2021 (https://doi.org/10.1073/pnas.2023202118)
# download chvd (from https://zenodo.org/record/4498884#.Y07qNtfMJhk)
cd /gscratch/scrubbed/carsonjm/resources/

wget https://zenodo.org/record/4498884/files/CHVD_clustered_mash99_v1.tar.gz?download=1
mv CHVD_clustered_mash99_v1.tar.gz?download=1 CHVD_clustered_mash99_v1.tar.gz
tar -xvf CHVD_clustered_mash99_v1.tar.gz

wget https://zenodo.org/record/4498884/files/HV3_table2_master_table.xlsx?download=1
mv HV3_table2_master_table.xlsx?download=1 HV3_table2_master_table.xlsx

# extract only high quality genomes
# 1) retain only genomes that are > 90% complete
python /gscratch/stf/carsonjm/CarsonJM/phide_piper/workflow/scripts/filter_chvd.py


### Camarillo-Guerrero 2021 (https://doi.org/10.1016/j.cell.2021.01.029)
# download gpd (from http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/)
cd /gscratch/scrubbed/carsonjm/resources/gpd

wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_metadata.tsv

wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_sequences.fa.gz
gunzip GPD_sequences.fa.gz

# extract only high quality genomes
# 1) retain only genomes that are > 90% complete
python /gscratch/stf/carsonjm/CarsonJM/phide_piper/workflow/scripts/filter_gpd.py


### Gregory 2020 (https://doi.org/10.1016/j.chom.2020.08.003)
# download gvd (from https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Gregory_and_Zablocki_GVD_Jul2020)
cd /gscratch/scrubbed/carsonjm/resources/gvd

wget https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Gregory_and_Zablocki_GVD_Jul2020/GVD_Viral_Populations/GVDv1_viralpopulations.fna.tar.gz
tar -xvf GVDv1_viralpopulations.fna.tar.gz

# extract only high quality genomes
# 1) rerun genomes through Checkv
# 2) retain only genomes that are > 90% complete


### Nayfach 2021 (https://doi.org/10.1038/s41564-021-00928-6)
# download mgv (from https://portal.nersc.gov/MGV/)
cd /gscratch/scrubbed/carsonjm/resources/mgv_db

wget https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08.tar.gz
tar -xvf MGV_v1.0_2021_07_08.tar.gz

# extract only high quality genomes
# 1) retain only genomes that are > 90% complete
python /gscratch/stf/carsonjm/CarsonJM/phide_piper/workflow/scripts/filter_mgv.py


# Li 2022 (https://doi.org/10.1016/j.isci.2022.104418)
# download ovd (from https://github.com/grc145/Temp)
cd /gscratch/scrubbed/carsonjm/resources/oral_virome_database

wget https://github.com/grc145/Temp/raw/master/OVD-genomes.fa.bz2

# extract only high quality genomes
# 1) retain only genomes that are > 90% complete
python /gscratch/stf/carsonjm/CarsonJM/phide_piper/workflow/scripts/filter_ovd.py


# download img vr4 (from https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html)
# extract genomes > 90% complete
python /gscratch/stf/carsonjm/CarsonJM/phide_piper/workflow/scripts/filter_imgvr.py

cat /gscratch/scrubbed/carsonjm/resources/benler_phage/benler_2020_phages.fna \
/gscratch/scrubbed/carsonjm/resources/checkv/checkv-db-v1.4/genome_db/checkv_reps.fna \
/gscratch/scrubbed/carsonjm/resources/chvd/CHVD_hq.fasta \
/gscratch/scrubbed/carsonjm/resources/gpd/GPD_hq_sequences.fa \
/gscratch/scrubbed/carsonjm/resources/imgvr_6/IMGVR_hq_nucleotides.fna \
/gscratch/scrubbed/carsonjm/resources/mgv_db/MGV_v1.0_2021_07_08/mgv_hq_contigs.fna \
/gscratch/scrubbed/carsonjm/resources/oral_virome_db/OVD_hq_genomes.fa > /gscratch/scrubbed/carsonjm/resources/hq_virus_db/combined_hq_genomes.fa