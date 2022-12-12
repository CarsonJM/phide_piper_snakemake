import os
import pandas as pd
from Bio import SeqIO

# create report to add results to
virus_report = pd.DataFrame(columns = ['contig_id'])

# extract mgv virus sequences
if snakemake.params.run_mgv or snakemake.params.run_virfinder:
    if os.stat(str(snakemake.input.mgv_results)).st_size != 0:
        mgv_results = pd.read_csv(str(snakemake.input.mgv_results), sep='\t')
        mgv_viruses = pd.read_csv(str(snakemake.input.mgv_viruses), sep='\t')
        mgv_results.rename(columns = {'vfr_score':'VirFinder_score', 'vfr_pvalue':'VirFinder_pvalue'}, inplace=True)
        mgv_results['MGV_viral'] = mgv_results.apply(lambda x: 'Viral' if x.contig_id in set(mgv_viruses['contig_id']) else 'Non-viral', axis=1)
        virus_report = virus_report.merge(mgv_results[['contig_id', 'length', 'gene_len', 'cds_density', 'switch_rate', 'VirFinder_score', 'vpfs', 'pfams', 'MGV_viral']], on='contig_id', how='outer')
        virus_report.contig_id = virus_report.contig_id.astype(str)

# extract virsorter virus sequences
if snakemake.params.run_virsorter:
    if os.stat(str(snakemake.input.vs_results)).st_size != 0:
        vs_results = pd.read_csv(str(snakemake.input.vs_results), sep="\t", header=None,
        names = ['VirSorter_contig_id', 'Contig_genes', 'Fragment', 'Fragment_genes', 'Category_text', 'Category_number',
        'Phage hallmark genes', 'Phage gene enrichment sig', 'Non-Caudovirales phage gene enrichment sig',
        'Pfam depletion sig', 'Uncharacterized enrichment sig', 'Strand switch depletion sig', 'Short genes enrichment sig'])
        vs_key = pd.read_csv(str(snakemake.input.vs_translation), sep="\t", header=None, names=['contig_id', 'VirSorter_contig_id'])
        vs_results['VirSorter_contig_id'] = vs_results['VirSorter_contig_id'].str.partition('-circular')[0]
        vs_key_results = vs_key.merge(vs_results, on='VirSorter_contig_id', how='outer')
        vs_key_results.contig_id = vs_key_results.contig_id.astype(str)
        virus_report = virus_report.merge(vs_key_results, on='contig_id', how='outer')
        virus_report.contig_id = virus_report.contig_id.astype(str)

# extract virsorter2 virus sequences
if snakemake.params.run_virsorter2:
    if os.stat(str(snakemake.input.vs2_results)).st_size != 0:
        vs2_results = pd.read_csv(str(snakemake.input.vs2_results), sep="\t")
        vs2_results['contig_id'] = vs2_results['seqname'].str.partition('||')[0]
        vs2_results.rename(columns = {'max_score':'VirSorter2_max_score'}, inplace=True)
        virus_report = virus_report.merge(vs2_results[['contig_id', 'VirSorter2_max_score', 'max_score_group', 'hallmark', 'viral', 'cellular']], on='contig_id', how='outer')
        virus_report.contig_id = virus_report.contig_id.astype(str)

# extract vibrant virus sequences
if snakemake.params.run_vibrant:
    if os.stat(str(snakemake.input.vb_results)).st_size != 0:
        vb_results = pd.read_csv(str(snakemake.input.vb_results), sep="\t", header=None)
        vb_results = vb_results.astype(str)
        if vb_results[0].str.contains('_fragment').any():
            vb_results[['contig_id', 'fragment']] = vb_results[0].str.split('_fragment', n=1, expand=True)
        else:
            vb_results['contig_id'] = vb_results[0]
        vb_results.rename(columns = {0:'VIBRANT_viruses'}, inplace=True)
        virus_report = virus_report.merge(vb_results[['contig_id', 'VIBRANT_viruses']], on='contig_id', how='outer')

# extract deepvirfinder virus sequences
if snakemake.params.run_deepvirfinder:
    if os.stat(str(snakemake.input.dvf_results)).st_size != 0:
        dvf_results = pd.read_csv(str(snakemake.input.dvf_results), sep="\t")
        dvf_results.rename(columns = {'name': 'contig_id', 'score':'DeepVirFinder_score', 'pvalue':'DeepVirFinder_pvalue'}, inplace=True)
        dvf_results.contig_id = dvf_results.contig_id.astype(str)
        virus_report = virus_report.merge(dvf_results[['contig_id', 'DeepVirFinder_score']], on='contig_id', how='outer')

# extract diamond virus sequences
if snakemake.params.run_genomad:
    if os.stat(str(snakemake.input.genomad_results)).st_size != 0:
        genomad_results = pd.read_csv(str(snakemake.input.genomad_results), sep='\t')
        if len(genomad_results['seq_name'].str.split('|', expand=True)) > 1:
            genomad_results['contig_id'] = genomad_results['seq_name'].str.split('|', expand=True)[0]
        virus_report = virus_report.merge(genomad_results[['contig_id', 'topology', 'n_genes', 'genetic_code', 'virus_score', 'fdr', 'n_hallmarks', 'marker_enrichment']], on='contig_id', how='outer')

# extract diamond virus sequences
if snakemake.params.run_external:
    if os.stat(str(snakemake.input.external_results)).st_size != 0:
        external_results = pd.read_csv(str(snakemake.input.external_results), sep='\t', names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-id', 'query-comment'])
        external_results['shared-hashes'] = external_results['shared-hashes'].str.partition('/')[0].astype(int)
        external_results = external_results[external_results['identity'] >= 0.99]
        virus_report = pd.concat([virus_report, external_results], axis=0, ignore_index=True)

virus_report['assembly'] = snakemake.params.assembly
virus_report['contig_id'] = virus_report['assembly'] + '_' + virus_report['contig_id']

if snakemake.params.run_external:
    virus_report['query-id'] = virus_report['assembly'] + '_' + virus_report['query-id']


virus_report.to_csv(str(snakemake.output), index=False)

