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

def import_raw_miRTarBase(csv_path):
    rename_dict = {
        'miRTarBase ID':'mirtar_id', 
        'miRNA':'mirna_name',
        'Species (miRNA)':'mirna_species',
        'Target Gene':'gene_name',
        'Target Gene (Entrez ID)':'gene_id_entrez',
        'Species (Target Gene)':'gene_species',
        'Experiments':'experiments',
        'Support Type':'support',
        'References (PMID)':'ref_pmid',
        'Target Site':'target_site'
    }

    dtype_dict = {
        'Target Gene (Entrez ID)':'int'
    }

    df = pd.read_csv(
        csv_path,
        dtype= dtype_dict,
        compression='infer'
    )

    df = df.rename(columns=rename_dict)

    return df

def convert_to_nt_pattern(seq, pattern_dict):
    pattern = seq
    for wc_char in pattern_dict.keys():
        if wc_char in pattern:
            pattern = pattern.replace(wc_char, pattern_dict[wc_char])
    return pattern

def clean_mirtar_tgt_sites(mirtar_tgt_sites_raw_path):
    mirtar_tgt_sites_df = import_raw_miRTarBase(mirtar_tgt_sites_raw_path)
    
    non_hsa_idx = mirtar_tgt_sites_df[mirtar_tgt_sites_df['mirna_species'] != 'hsa'].index
    mirtar_tgt_sites_df = mirtar_tgt_sites_df.drop(non_hsa_idx)
    
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x.strip().upper())
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: ''.join(x.split(' ')))
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x.replace('U', 'T'))
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x.replace('-', 'N'))
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x.replace("5'", ''))
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x.replace("3'", ''))
    mirtar_tgt_sites_df['target_site'] = mirtar_tgt_sites_df['target_site'].apply(lambda x: x[1:] if x[0] == '.' else x)

    group_by = mirtar_tgt_sites_df.groupby(by=['target_site'])
    mirtar_tgt_sites_df = group_by['mirtar_id'].apply(lambda x: ','.join(sorted(list(x.unique())))).reset_index()
    mirtar_tgt_sites_df['mirna_name'] = group_by['mirna_name'].apply(lambda x: ','.join(sorted(list(x.unique())))).reset_index(drop=True)
    mirtar_tgt_sites_df['mirna_species'] = group_by['mirna_species'].apply(lambda x: ','.join(sorted(list(x.unique())))).reset_index(drop=True)
    mirtar_tgt_sites_df['gene_name'] = group_by['gene_name'].apply(lambda x: ','.join(sorted(list(x.unique())))).reset_index(drop=True)
    mirtar_tgt_sites_df['gene_id_entrez'] = group_by['gene_id_entrez'].apply(lambda x: ','.join(sorted(list(set([str(v) for v in x]))))).reset_index(drop=True)
    mirtar_tgt_sites_df['support'] = group_by['support'].apply(lambda x: ','.join(sorted(list(x.unique())))).reset_index(drop=True)
    
    nt_list = ['A', 'C', 'G', 'T']
    expanded_nt_list = ['A', 'C', 'G', 'T', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']
    non_nt_mask = ~(mirtar_tgt_sites_df['target_site'].apply(lambda x: pd.Series(list(x)).isin(expanded_nt_list).all()))
    non_nt_idx = mirtar_tgt_sites_df[non_nt_mask].index
    mirtar_tgt_sites_df = mirtar_tgt_sites_df.drop(non_nt_idx)

    pattern_dict ={
        'W':'[AT]',
        'S':'[CG]',
        'M':'[AC]',
        'K':'[GT]',
        'R':'[AG]',
        'Y':'[CT]',
        'B':'[CGT]',
        'D':'[AGT]',
        'H':'[ACT]',
        'V':'[ACG]',
        'N':'[ACGT]'
    }
    mirtar_tgt_sites_df['target_site_pattern'] = mirtar_tgt_sites_df['target_site'].apply(convert_to_nt_pattern, pattern_dict=pattern_dict)

    return mirtar_tgt_sites_df

def build_target_site_pattens(mirtar_tgt_sites_df):
    patterns = {}

    for mirtar_idx, mirtar_row in tqdm.tqdm(mirtar_tgt_sites_df.iterrows(), total=len(mirtar_tgt_sites_df), desc='Building Pattern Dictionary for miRTarBase Target Sites'):
        patterns[mirtar_idx] = {}
        patterns[mirtar_idx]['mirtar_id'] = mirtar_row['mirtar_id']
        patterns[mirtar_idx]['mirna_name'] = mirtar_row['mirna_name']
        patterns[mirtar_idx]['tgt_site_pattern'] = mirtar_row['target_site_pattern']
    
    return patterns

###############
# main script #
###############
if __name__ == '__main__':        
    #--- Step 1: define variables ---
    # from args
    ref_db_dir = sys.argv[1]

    # other directory paths
    miRTarBase_dir = os.path.join(ref_db_dir, 'miRtarBase')

    #---Step 2: Load miRTarBase Target Sites Data ----
    print(f'Loaded miRTarBase Target Sites Data...')
    mirtar_tgt_sites_pq_path = os.path.join(miRTarBase_dir, 'miRTarBase.tgt-sites.hsa.parquet')
    if os.path.isfile(mirtar_tgt_sites_pq_path):
        mirtar_tgt_sites_df = pd.read_parquet(mirtar_tgt_sites_pq_path)
        mirtar_tgt_sites_df['target_site_pattern'] = mirtar_tgt_sites_df['target_site_pattern'].apply(re.compile)
    else:
        mirtar_tgt_sites_raw_path = os.path.join(miRTarBase_dir, 'MicroRNA_Target_Sites.csv.gz')
        mirtar_tgt_sites_df = clean_mirtar_tgt_sites(mirtar_tgt_sites_raw_path=mirtar_tgt_sites_raw_path)

        mirtar_tgt_sites_df.to_parquet(mirtar_tgt_sites_pq_path, compression='gzip')
        mirtar_tgt_sites_df['target_site_pattern'] = mirtar_tgt_sites_df['target_site_pattern'].apply(re.compile)
    print(f'Loaded miRTarBase Target Sites Data, {len(mirtar_tgt_sites_df)} Entries.\n')

    #---Step 3: Build miRTarBase Target Site Patterns Dict ----
    print(f'Building Pattern Dict for miRTarBase tgt sites...')
    tgt_site_patterns_path = os.path.join(miRTarBase_dir, 'miRTarBase.tgt-sites.patterns.pkl.gz')
    if os.path.isfile(tgt_site_patterns_path):
        print(f'\nPattern Dictionary Built Previously, Loading...')
        with gzip.open(tgt_site_patterns_path, 'rb') as f:
                mirtar_patterns = pickle.load(f)
        print(f'Loaded Pattern Dictionary for {len(mirtar_patterns)} Target Sites.\n')
    else:
        print(f'\nBuilding Pattern Dictionary for {len(mirtar_tgt_sites_df)} Target Sites...')
        mirtar_patterns = build_target_site_pattens(mirtar_tgt_sites_df=mirtar_tgt_sites_df)
        with gzip.open(tgt_site_patterns_path, 'wb') as f:
            pickle.dump(mirtar_patterns, f)
        print(f'Built Pattern Dictionary for {len(mirtar_patterns)} Target Sites.\n')
    print(f'Built Pattern Dict for miRTarBase tgt sites.')