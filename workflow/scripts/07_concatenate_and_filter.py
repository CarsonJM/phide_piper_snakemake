import Bio.SeqIO
import os

genome_ids = set([_[1:].split()[0].rsplit('_', 1)[0] for _ in open(str(snakemake.input.viruses)) if _[0]=='>'])

msa_catted = dict([(id, '') for id in genome_ids])
msa_count = dict([(id, 0) for id in genome_ids])
for index, file in enumerate(os.listdir(str(snakemake.params.in_dir))):
	marker_id = file.split('.')[0]
	genome_to_gene = {}
	inpath = '%s/%s' % (str(snakemake.params.in_dir), file)
	for r in Bio.SeqIO.parse(inpath, format='clustal'):
		genome_id = r.id.rsplit('_', 1)[0]
		genome_to_gene[genome_id] = r
	marker_length = len(genome_to_gene[genome_id].seq)
	for genome_id in genome_ids:
		if genome_id not in genome_to_gene:
			msa_catted[genome_id] += "-" * marker_length
		else:
			gene = genome_to_gene[genome_id]
			msa_catted[genome_id] += str(gene.seq).replace('.', '-')
			msa_count[genome_id] += 1

genome_to_gaps = {}
for genome_id in genome_ids:
	seq = msa_catted[genome_id]
	perc_gaps = round(100*seq.count('-')/float(len(seq)), 2)
	genome_to_gaps[genome_id] = perc_gaps

hits_dir = str(snakemake.params.out_dir)
if not os.path.exists(hits_dir):
	os.makedirs(hits_dir)

with open('%sconcat.faa' % str(snakemake.params.out_dir), 'w') as f:
	for genome_id, seq in msa_catted.items():
		perc_gaps = genome_to_gaps[genome_id]
		num_hits = msa_count[genome_id]
		if perc_gaps < snakemake.params.max_percent_gaps and msa_count[genome_id] >= snakemake.params.min_hits:
			f.write('>%s percent_gaps=%s num_hits=%s\n' % (genome_id, perc_gaps, num_hits))
			f.write('%s\n' % seq)
