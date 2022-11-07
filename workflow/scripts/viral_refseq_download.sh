# make directory to store viral refseq files in
mkdir /home/carsonjm/resources/viralrefseq212
mkdir /home/carsonjm/resources/viralrefseq212/protein
mkdir /home/carsonjm/resources/viralrefseq212/genomic

# download viral refseq proteins from NCBI
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz -P /home/carsonjm/resources/viralrefseq212/protein -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz -P /home/carsonjm/resources/viralrefseq212/protein -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz -P /home/carsonjm/resources/viralrefseq212/protein -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.faa.gz -P /home/carsonjm/resources/viralrefseq212/protein -w 5 -T 5

# concatenate viral refseq proteins
cat /home/carsonjm/resources/viralrefseq212/protein/viral.*.protein.faa.gz > /home/carsonjm/resources/viralrefseq212/protein/viral.combined.protein.faa.gz

# gunzip viral refseq proteins
gunzip /home/carsonjm/resources/viralrefseq212/protein/viral.combined.protein.faa.gz

# download viral refseq genomes from NCBI
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz -P /home/carsonjm/resources/viralrefseq212/genomic -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz -P /home/carsonjm/resources/viralrefseq212/genomic -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz -P /home/carsonjm/resources/viralrefseq212/genomic -w 5 -T 5
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz -P /home/carsonjm/resources/viralrefseq212/genomic -w 5 -T 5

# concatenate viral refseq genomes
cat /home/carsonjm/resources/viralrefseq212/genomic/viral.*.1.genomic.fna.gz > /home/carsonjm/resources/viralrefseq212/genomic/viral.combined.1.genomic.fna.gz

# download refseq catalog
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release212.catalog.gz -P /home/carsonjm/resources/viralrefseq212 -w 5 -T 5

# gunzip refseq catalog
gunzip /home/carsonjm/resources/viralrefseq212/RefSeq-release212.catalog.gz