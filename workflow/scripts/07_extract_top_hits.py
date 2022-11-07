import Bio.SeqIO
import os

faa = dict([(r.id, str(r.seq).rstrip('*')) for r in Bio.SeqIO.parse(str(snakemake.input.viruses), 'fasta')])

hits_dir = str(snakemake.params.out_dir)
if not os.path.exists(hits_dir):
	os.makedirs(hits_dir)

hits = {}
for line in open(str(snakemake.input.hmmsearch), 'r'):
	if line[0]=='#': continue
	r = line.split()
	if float(r[4]) > snakemake.params.max_evalue: continue
	elif r[0] not in hits:
		hits[r[0]] = r
	elif float(r[5]) > float(hits[r[0]][5]):
		hits[r[0]] = r

marker_to_hits = {}
for r in hits.values():
	marker_id = r[2]
	if marker_id not in marker_to_hits:
		marker_to_hits[marker_id] = []
	gene = faa[r[0]]
	marker_to_hits[marker_id].append([r[0], gene])

for marker_id in marker_to_hits:
	outpath = (str(snakemake.params.out_dir) + marker_id + ".faa")
	out = open(outpath, 'w')
	for gene_id, seq in marker_to_hits[marker_id]:
		out.write('>'+gene_id+'\n'+seq+'\n')
	out.close()

fp = open(str(snakemake.output), 'x')
fp.close()