import pandas as pd
import os
import shutil

# read ccp77
ccp77 = pd.read_csv(str(snakemake.input.ccp77), sep='\t')
ccp77_set = set(ccp77['VOGid'])

# list all files in vogdb
vogdb_hmms = os.listdir(str(snakemake.params.vogdb_dir))

# copy files in ccp77
for hmm in vogdb_hmms:
    vogdb = hmm.partition('.')[0]
    if vogdb in ccp77_set:
        shutil.copy(str(snakemake.params.vogdb_dir) + hmm, str(snakemake.params.ccp77_dir))
    else:
        continue
