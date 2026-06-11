#####################
# import statements #
#####################
import os
import pandas as pd
import sys
import tqdm
import RNA
from collections import defaultdict
from tabulate import tabulate
from collections import defaultdict
import re
import gzip
from Bio import SeqIO
from Bio import Seq
import pysam
import pickle
from pandarallel import pandarallel
from functools import reduce
import random
import gc

#############
# functions #
#############

def build_seed_patterns(mirbase_df):
    patterns = {}
    
    for mirna_id, mirna_row in mirbase_df.iterrows():
        patterns[mirna_id] = {}
        sixmer_pattern = mirna_row['6-mer'].replace('N', '[ACGT]')
        
        patterns[mirna_id]['mirna_name'] = mirna_row['name']
        for seed_type in ['8-mer', '7-mer-A1', '7-mer-m8', '6-mer']:
            patterns[mirna_id][seed_type] = mirna_row[seed_type]
        patterns[mirna_id]['6-mer_pattern'] = re.compile(sixmer_pattern)
    
    return patterns

###############
# main script #
###############
if __name__ == '__main__':        
    #--- Step 1: define variables ---
    # from args
    ref_db_dir = sys.argv[1]

    # other directory paths
    miRBase_dir = os.path.join(ref_db_dir, 'miRBase')

    #---Step 2: Load miRBase Data ----
    print(f'Loading miRBase Sequence Data...')
    mirbase_out_path = os.path.join(miRBase_dir, 'mirbase.parquet')

    # check if mirbase df was built previously, if not build it
    if os.path.isfile(mirbase_out_path):
        print('mirbase_df built previously...\n')
        mirbase_df = pd.read_parquet(mirbase_out_path)
        mirbase_df['6mer_pattern'] = mirbase_df['6-mer'].apply(lambda x: x.replace('N', '[ACTG]'))
        mirbase_df['6mer_pattern'] = mirbase_df['6mer_pattern'].apply(re.compile)
    else:
        # Parse the miRBase file and filter for human miRNAs
        print('mirbase_df not Built Previously, Parsing:')
        miRBase_raw_path = os.path.join(miRBase_dir, 'miRNA.dat.gz')
        mirbase_df_data = {
            'accession': [],
            'name': [],
            'seq': [],
            'evidence': [],
            'experiment': []
        }

        with gzip.open(miRBase_raw_path, 'rt') as f:
            for record in tqdm.tqdm(SeqIO.parse(f, 'embl'), desc='Parsing miRBase Raw Data:'):
                # Human miRNAs start with 'hsa-' in their ID
                if record.name.startswith('hsa-'):
                    for rec_feat in record.features:
                        if rec_feat.type.lower() == 'mirna':
                            if len(rec_feat.qualifiers['accession']) == 1:
                                mirna_acc = rec_feat.qualifiers['accession'][0]
                            else:
                                mirna_acc = ','.join(rec_feat.qualifiers['accession'])
                                print(f'WARNING: {record.id} has miRNA feature with multiple acessions')

                            if len(rec_feat.qualifiers['product']) == 1:
                                mirna_name = rec_feat.qualifiers['product'][0]
                            else:
                                mirna_name = ','.join(rec_feat.qualifiers['product'])
                                print(f'WARNING: {record.id} has miRNA feature with multiple products')
                            
                            mirbase_df_data['accession'].append(    mirna_acc)
                            mirbase_df_data['name'].append(         mirna_name)
                            mirbase_df_data['seq'].append(          str(rec_feat.location.extract(record.seq).back_transcribe()).upper())
                            mirbase_df_data['evidence'].append(     ','.join(rec_feat.qualifiers.get('evidence', [])))
                            mirbase_df_data['experiment'].append(   ','.join(rec_feat.qualifiers.get('experiment', [])))

        mirbase_df = pd.DataFrame(mirbase_df_data)
        mirbase_df = mirbase_df.groupby(['accession', 'name', 'seq'], as_index=False).agg({
            'evidence': lambda x: ','.join(x.unique()),
            'experiment': lambda x: ','.join(filter(None, x.unique()))  # filters out empty strings
        })
        mirbase_df = mirbase_df.set_index('accession')
        mirbase_df['seed']          = mirbase_df['seq'].apply(lambda x: x[1:8])
        mirbase_df['8-mer']         = mirbase_df['seq'].apply(lambda x: str(Seq.Seq(x[1:8]).reverse_complement()) + 'A')
        mirbase_df['7-mer-A1']      = mirbase_df['seq'].apply(lambda x: str('N' + Seq.Seq(x[1:7]).reverse_complement()) + 'A')
        mirbase_df['7-mer-m8']      = mirbase_df['seq'].apply(lambda x: str(Seq.Seq(x[1:8]).reverse_complement()) + 'N')
        mirbase_df['6-mer']         = mirbase_df['seq'].apply(lambda x: str('N' + Seq.Seq(x[1:7]).reverse_complement()) + 'N')

        column_order = [
            'name',
            'seq',
            'seed',
            '8-mer',
            '7-mer-A1',
            '7-mer-m8',
            '6-mer',
            'evidence',
            'experiment'
        ]
        
        mirbase_df = mirbase_df[column_order]
        mirbase_df.to_parquet(mirbase_out_path, compression='gzip')

        # garbage collection
        del mirbase_df_data
        gc.collect()


    print(f'miRbase Total Entries: {len(mirbase_df)}')
    print(f'Loaded miRBase Sequence Data.\n')

    #--- Step 3: Build miRBase Seed Pattern Dict ---
    print(f'Build miRBase Seed Pattern Dict...')
    seed_patterns_path = os.path.join(miRBase_dir, 'miRBase.seed-patterns.pkl.gz')
    if os.path.isfile(seed_patterns_path):
        print(f'\nPattern Dictionary Built Previously, Loading...')
        with gzip.open(seed_patterns_path, 'rb') as f:
                seed_patterns = pickle.load(f)
        print(f'Loaded Pattern Dictionary for {len(seed_patterns)} Target Sites.\n')
    else:
        print(f'\nBuilding Pattern Dictionary for {len(mirbase_df)} Target Sites...')
        seed_patterns = build_seed_patterns(mirbase_df=mirbase_df)
        with gzip.open(seed_patterns_path, 'wb') as f:
            pickle.dump(seed_patterns, f)
        print(f'Built Pattern Dictionary for {len(seed_patterns)} Target Sites.\n')
    print(f'Built miRBase Seed Pattern Dict.')