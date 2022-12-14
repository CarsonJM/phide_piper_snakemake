import os
import pandas as pd
from Bio import SeqIO

# extract diamond virus sequences
if snakemake.params.run_genomad:
    if os.stat(str(snakemake.input.genomad_results)).st_size != 0:
        genomad_results = pd.read_csv(str(snakemake.input.genomad_results), sep='\t')
        genomad_results.rename(columns = {'seq_name':'vls_id'}, inplace=True)

# extract diamond virus sequences
if snakemake.params.run_external:
    if os.stat(str(snakemake.input.external_results)).st_size != 0:
        external_results = pd.read_csv(str(snakemake.input.external_results), sep='\t', names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'vls_id', 'query-comment'])
        external_results['shared-hashes'] = external_results['shared-hashes'].str.partition('/')[0].astype(int)
        external_results = external_results[(external_results['identity'] >= snakemake.params.min_mash_score) & (external_results['shared-hashes'] >= snakemake.params.min_mash_hashes) & (external_results['median-multiplicity'] >= snakemake.params.min_mash_multiplicity)]
        virus_report = pd.concat([genomad_results, external_results], axis=0, ignore_index=True)

virus_report['assembly'] = snakemake.params.assembly
virus_report['vls_id'] = virus_report['assembly'] + '_' + virus_report['vls_id']

virus_report.to_csv(str(snakemake.output), index=False)

