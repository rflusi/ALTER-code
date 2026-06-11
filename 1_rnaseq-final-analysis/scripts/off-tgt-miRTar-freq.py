#####################
# import statements #
#####################

import os
import pandas as pd
import sys
import tqdm
import numpy as np
from collections import defaultdict
from tabulate import tabulate
from collections import defaultdict
from functools import reduce

#############
# functions #
#############
def load_sample_map(proj_dir):
    sample_map_df = pd.read_csv(os.path.join(proj_dir,'sample-map.csv'))
    
    non_wt_mask = ~sample_map_df['condition'].isin(sample_map_df['wt_condition'])
    non_wt_conditions = list(sample_map_df.loc[non_wt_mask, 'condition'].unique())
    non_wt_conditions.sort()

    return sample_map_df, non_wt_conditions

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

def mirtarbase_by_gene(mirtar_path):
    # -------------------------------------------------
    # Input:    1. mirtar_path - path to the downloaded csv of human mirtarbase
    # Output:   a df of humn mirtar base grouped by gene
    # -------------------------------------------------
    df = import_raw_miRTarBase(mirtar_path)

    df = df[df['mirna_species'] == 'hsa'].reset_index(drop=True)
    df = df.drop('mirna_species', axis=1)
    df = df.drop('gene_species', axis=1)
    df['ref_pmid'] = df['ref_pmid'].apply(lambda x: str(int(x)) if not pd.isna(x) else '')

    group_by = df.groupby(['gene_name', 'gene_id_entrez'])
    grouped_df = group_by['mirtar_id'].apply(lambda x: ','.join(x)).reset_index()
    grouped_df['mirna_name'] = list(group_by['mirna_name'].apply(lambda x: ','.join(x)))
    grouped_df['experiments'] = list(group_by['experiments'].apply(lambda x: ','.join(x)))
    grouped_df['support'] = list(group_by['support'].apply(lambda x: ','.join(x)))
    grouped_df['ref_pmid'] = list(group_by['ref_pmid'].apply(lambda x: ','.join(x)))

    return grouped_df

###############
# main script #
###############
if __name__ == '__main__':
    tqdm.tqdm.pandas()
    
    #--- Step 1: define variables ---
    # from args
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    var_dir = sys.argv[3]
    ref_db_dir = sys.argv[4]
    # target SNP as a tuple (ref, alt)
    target_snp = (sys.argv[5][0], sys.argv[5][1])

    # other directory paths
    miRTarBase_dir = os.path.join(ref_db_dir, 'miRtarBase')
    biomart_dir = os.path.join(ref_db_dir, 'biomart')
    align_dir = os.path.join(var_dir, 'transcript-align')
    align_snp_dir = os.path.join(align_dir, f'{target_snp[0]}{target_snp[1]}')
    expression_dir = os.path.join(proj_dir, '5_expression')
    salmon_dir = os.path.join(expression_dir, 'salmon')
    salmon_analysis_dir = os.path.join(salmon_dir, 'analysis')
    hit_dir = os.path.join(var_dir, 'hit-id')
    out_dir = os.path.join(hit_dir, 'miRTarBase-freq')

    os.makedirs(out_dir, exist_ok=True)

    sample_map_df, non_wt_conditions = load_sample_map(proj_dir=proj_dir)

    overlap_path = os.path.join(hit_dir, f'{experiment}.overlap.{target_snp[0]}{target_snp[1]}.parquet')
    overlap_df = pd.read_parquet(overlap_path)
    non_out_overlap_cols = pd.Series(overlap_df.columns)

    #---Step 1: Import miRTarBase Data----
    print(f'\nLoading miRTarBase Gene Targets Data...')
    miRTarBase_gene_path = os.path.join(miRTarBase_dir, 'miRTarBase-by-gene.parquet')
    if os.path.isfile(miRTarBase_gene_path):
        print(f'Gene Level db Created Previously, Loading...')
        miRTarBase_gene_df = pd.read_parquet(miRTarBase_gene_path)
    else:
        print(f'Building Gene Level db...')
        miRTarBase_csv_path = os.path.join(miRTarBase_dir, 'hsa_MTI.csv.gz')
        miRTarBase_gene_df = mirtarbase_by_gene(miRTarBase_csv_path)
        miRTarBase_gene_df = miRTarBase_gene_df.set_index('gene_id_entrez')
        miRTarBase_gene_df.to_parquet(miRTarBase_gene_path)
    print (f'Loaded miRTarBase Gene Targets Data.\n')

    #---Step 2: Import Gene Level Raw Expression Data----
    print (f'Loading Gene Expression Data...')
    gene_raw_path = os.path.join(salmon_analysis_dir, 'salmon-gene-raw-counts.tsv.gz')
    gene_raw_df = pd.read_csv(gene_raw_path, sep='\t', compression='infer', index_col='Unnamed: 0')
    print (f'Loaded Gene Expression Data.\n')

    print (f'Determining Expressed Genes (DP >= 20)...\n')
    expressed_mask_dict = {}
    for condition in non_wt_conditions:
        sample_mask = sample_map_df['condition'] == condition
        samples = sample_map_df.loc[sample_mask, 'sample']
        sample_cols = [f'{sample}_reads' for sample in samples]
        expressed_mask = (gene_raw_df[sample_cols] >= 20).all(axis=1)
        expressed_mask_dict[condition] = expressed_mask

        print(f'{len(gene_raw_df[expressed_mask])}\tExpressed Genes for {condition}')

    biomart_gene_path = os.path.join(biomart_dir, 'biomart_gene.parquet')
    entrez_lookup = pd.read_parquet(biomart_gene_path)['gene_id_entrez']
    entrez_lookup = entrez_lookup.apply(lambda x: [int(val) for val in x.strip().split(',') if val != ''])
    print (f'\nDetermined Expressed Genes (DP >= 20).\n')

    #---Step 3: Calculate Frequency of Gene Subsets in miRTarBase ----
    print (f'Determining Frequency of miRNA Target Genes...\n')
    
    gene_id_lists = overlap_df['gene_id'].apply(lambda x: [id for id in x.split(',')])                          # split multiple genes at site into list
    gene_id_lists = gene_id_lists.apply(lambda x: [id.split('.')[0] for id in x])                               # de-version gene ids for entrez lookup
    gene_id_lists = gene_id_lists.apply(lambda x: [entrez_lookup[id] for id in x if id in entrez_lookup.index]) # lookup matching entrez id for miRTarBase search, these are in lists as some ensembl ids match multipl entrez ids
    gene_id_lists = gene_id_lists.apply(lambda x: reduce(lambda y,z: y+z, x) if len(x) > 0 else x)              # combine entrez id lists if there are multiple 
    gene_id_lists = gene_id_lists.apply(lambda x: sorted(list(set(x))))                                         # sort entrez id list
    overlap_df['gene_id_entrez'] = gene_id_lists.apply(lambda x: ','.join([str(id) for id in x]))               # convert to csv string 
        
    entrez_path = os.path.join(hit_dir, f'{experiment}.overlap.{target_snp[0]}{target_snp[1]}.entrez.parquet')
    overlap_cols = pd.Series(overlap_df.columns)
    out_overlap_col_mask = ~(overlap_cols.isin(non_out_overlap_cols))
    out_overlap_cols = overlap_cols[out_overlap_col_mask]
    overlap_df.to_parquet(entrez_path, compression='gzip')                                                     # update overlap_df file

    mirna_tgt_idx = [
        'expressed_genes',
        'off_tgt_hits',
        'unique_off_tgt_hits',
    ]

    freq_cols = defaultdict(lambda: [])
    count_cols = defaultdict(lambda: [])
    entrez_out_dict = {}
    unique_entrez_out_dict = {}
    for condition in non_wt_conditions:
        expressed_genes_df = gene_raw_df[expressed_mask_dict[condition]].copy() # expressed gene set all genes that passed depth filter for the condition
        
        off_tgt_df = overlap_df[overlap_df[f'{condition}_hit']].copy()              # get entries that were hits for the condition
        off_tgt_df = off_tgt_df['gene_id'].apply(lambda x: x.strip().split(','))    # split csv into list
        off_tgt_df = list(set(sum(off_tgt_df.tolist(), [])))                         # combine all lists and keep only unique gene_ids (consistent with expressed dataset)
        off_tgt_df = pd.DataFrame({                                                 # for analysis a dataframe with a single column where each row is one gene id
            'gene_id': off_tgt_df
        })
        
        # same as above but additional filtering for off targets that didnt appear in any other condition before gene_id collection
        unique_off_tgt_df = overlap_df[overlap_df[f'{condition}_hit']].copy()   
        other_conditions = pd.Series(non_wt_conditions)[pd.Series(non_wt_conditions) != condition].to_list()
        unique_mask = ~(unique_off_tgt_df[[f'{other_condition}_hit' for other_condition in other_conditions]].any(axis=1))
        unique_off_tgt_df = unique_off_tgt_df[unique_mask]
        unique_off_tgt_df = unique_off_tgt_df['gene_id'].apply(lambda x: x.strip().split(','))
        unique_off_tgt_df = list(set(sum(unique_off_tgt_df.tolist(), [])))
        unique_off_tgt_df = pd.DataFrame({
            'gene_id': unique_off_tgt_df
        })

        for df_list_idx, df_to_count in enumerate([expressed_genes_df, off_tgt_df, unique_off_tgt_df]):
            df_to_count['gene_id_entrez'] = df_to_count['gene_id'].apply(lambda x: x.strip().split('.')[0])
            df_to_count['gene_id_entrez'] = df_to_count['gene_id_entrez'].apply(lambda x: entrez_lookup[x] if x in entrez_lookup.index else [])
            entrez_series = pd.Series(list(set(sum(df_to_count['gene_id_entrez'].tolist(), []))))

            total_count = len(entrez_series)
            tgt_count = len(entrez_series[entrez_series.isin(miRTarBase_gene_df.index)])
            
            if total_count != 0:
                freq_cols[condition].append(tgt_count/total_count)
            else:
                freq_cols[condition].append(np.nan)

            count_cols[f'{condition}_total'].append(total_count)
            count_cols[f'{condition}_mirna_tgt'].append(tgt_count)

            if (df_list_idx == 1) or (df_list_idx == 2):
                na_idx = df_to_count[df_to_count['gene_id_entrez'].apply(lambda x: len(x) == 0)].index
                df_to_count = df_to_count.drop(na_idx)
                
                if not df_to_count.empty:
                    out_df_rows = df_to_count.apply(lambda x:
                        [[x['gene_id'], entrez_id]
                        for entrez_id in x['gene_id_entrez']],
                        axis=1
                    ).sum()
                    out_df = pd.DataFrame(
                        out_df_rows,
                        columns=['gene_id', 'gene_id_entrez']
                    )
                else:
                    out_df = df_to_count

                if (df_list_idx == 1):
                    entrez_out_dict[condition] = out_df
                elif (df_list_idx == 2):
                    unique_entrez_out_dict[condition] = out_df
                
    freq_cols = dict(freq_cols)
    total_cols = dict(count_cols)

    miRTarBase_freq_df = pd.DataFrame(freq_cols, index=mirna_tgt_idx)
    miRTarBase_freq_df = miRTarBase_freq_df.T
    miRTarBase_count_df = pd.DataFrame(count_cols, index=mirna_tgt_idx)
    miRTarBase_count_df = miRTarBase_count_df.T

    print(tabulate(miRTarBase_freq_df, headers='keys', tablefmt='fancy_grid'))
    print('')
    print(tabulate(miRTarBase_count_df, headers='keys', tablefmt='fancy_grid'))
    print (f'\nDetermined Frequency of miRNA Target Genes.\n')

    #---Step 4: Write Outputs ----
    xlsx_out_path = f'{hit_dir}/{experiment}.miRTar-freq.{target_snp[0]}{target_snp[1]}.xlsx'
    with pd.ExcelWriter(xlsx_out_path, engine='openpyxl') as writer:
        miRTarBase_freq_df.to_excel(writer, sheet_name=f'{target_snp[0]}{target_snp[1]}-freq', index=True)
        miRTarBase_count_df.to_excel(writer, sheet_name=f'{target_snp[0]}{target_snp[1]}-count', index=True)

    for condition, out_df in entrez_out_dict.items():
        pq_path = os.path.join(out_dir, f'{experiment}.{condition}.off-tgt-entrez-all.{target_snp[0]}{target_snp[1]}.parquet')
        out_df.to_parquet(pq_path, compression='gzip')
    
    for condition, out_df in unique_entrez_out_dict.items():
        pq_path = os.path.join(out_dir, f'{experiment}.{condition}.off-tgt-entrez-unique.{target_snp[0]}{target_snp[1]}.parquet')
        out_df.to_parquet(pq_path, compression='gzip')