import pandas as pd
from Bio import SeqIO

# read viral identification report
report = pd.read_csv(str(snakemake.input.viral_report))

# filter based on config paramters
# extract mgv virus sequences
if snakemake.params.run_mgv: 
    mgv_report = set(report[report['MGV_viral'] == 'Viral']['contig_id'])
else:
    mgv_report = []

if snakemake.params.run_vf:
    vf_report = set(report[(report['VirFinder_score'] >= snakemake.params.vf_score)]['contig_id'])
else:
    vf_report = []

if snakemake.params.run_vs:
    report['VirSorter_cat'] = report.apply(lambda x: x.Category_number if x.Category_text == 'complete_phage' else x.Category_number + 3, axis=1)
    vs_report = set(report[report['VirSorter_cat'].isin(snakemake.params.vs_cat)]['contig_id'])
else:
    vs_report = []

if snakemake.params.run_vs2:
    vs2_report = set(report[report['VirSorter2_max_score'] >= snakemake.params.vs2_score]['contig_id'])
else:
    vs2_report = []

if snakemake.params.run_dvf:
    dvf_report = set(report[(report['DeepVirFinder_score'] >= snakemake.params.dvf_score)]['contig_id'])
else:
    dvf_report = []

if snakemake.params.run_vb:
    vb_report = set(report[report['VIBRANT_viruses'].notnull()]['contig_id'])
else:
    vb_report = []

if snakemake.params.run_genomad:
    genomad_report = set(report[(report['virus_score'] >= snakemake.params.genomad_score) & (report['fdr'] <= snakemake.params.genomad_fdr)]['contig_id'])
else:
    genomad_report = []


if snakemake.params.run_external:
    external_report = set(report[(report['identity'] >= snakemake.params.min_mash_score) & (report['shared-hashes'] >= snakemake.params.min_mash_hashes) & (report['median-multiplicity'] >= snakemake.params.min_mash_multiplicity)]['query-id'])
else:
    external_report = []

# list to store all combined sequences
combined_sequences = []

# extract viral contigs
for record in SeqIO.parse(str(snakemake.input.contigs), "fasta"):
    record.id = snakemake.params.assembly + '_' + record.id
    if record.id in mgv_report:
        combined_sequences.append(record)
    elif record.id in vf_report:
        combined_sequences.append(record)
    elif record.id in vs_report:
        combined_sequences.append(record)
    elif record.id in vs2_report:
        combined_sequences.append(record)
    elif record.id in dvf_report:
        combined_sequences.append(record)
    elif record.id in vb_report:
        combined_sequences.append(record)
    elif record.id in genomad_report:
        combined_sequences.append(record)
    elif record.id in external_report:
        combined_sequences.append(record)

# prophage = []

# for record in SeqIO.parse(str(snakemake.input.gn_prophage), "fasta"):
#     record.id = snakemake.params.assembly + '_' + record.id
#     prophage.append(record.id.partition('|')[0])
#     combined_sequences.append(record)

# for record in SeqIO.parse(str(snakemake.input.vb_prophage), "fasta"):
#     record.id = snakemake.params.assembly + '_' + record.id
#     name = (record.id.partition('_fragment')[0])
#     if name in prophage:
#         continue
#     else:
#         combined_sequences.append(record)

SeqIO.write(combined_sequences, str(snakemake.output), "fasta")