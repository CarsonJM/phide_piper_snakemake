import pandas as pd
from Bio import SeqIO

# load checkv quality report
checkv = pd.read_csv(str(snakemake.input.checkv), sep='\t')

# retain HQ genomes
checkv_complete = checkv[checkv['completeness'] >= snakemake.params.min_completeness]

# read clustering output
clusters = open(str(snakemake.input.clusters), 'r')
votu_reps= []
for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    votu_reps.append(centroid)

votu_reps_set = set(votu_reps)

# retain only hq genomes that are cluster representatives
hq_genomes = checkv_complete[checkv_complete['contig_id'].isin(votu_reps_set)]
hq_genomes_set = set(hq_genomes['contig_id'])

# filter sequences to retain hq genome reps
hq_viruses = []

# read in the fasta
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    # retain contigs in hq set
    if record.id in hq_genomes_set:
        # add sequences passing filter to list
        hq_viruses.append(record)

# save saved sequences to specified file
SeqIO.write(hq_viruses, str(snakemake.output), "fasta")