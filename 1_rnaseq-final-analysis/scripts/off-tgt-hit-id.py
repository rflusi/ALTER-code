#####################
# import statements #
#####################

import os
import pandas as pd
import tqdm as tqdm
from functools import reduce
import sys
import numpy as np
from tabulate import tabulate
from pandarallel import pandarallel
from itertools import combinations

#############
# functions #
#############
def generate_filter_list(var_df, sample_map, target_snp, mutect=False):
    ###################################################################################################
    # Purpose: filter variant table for likely off targets, filtering strategy similar to CURE paper  #
    #          and all filtering is done within matched replicates. Filtering summary:                #
    #          1. var    - sample was called as a variant by HapplotypeCaller                         #
    #          2. VOI    - sample variant matches the variant of interest                             #
    #          3. DP     - read depth >= 20                                                           #
    #          4. GQ     - depth by quality >=20 (corresponds to 99% confidence)                      #
    #          5. non-wt - pct_ref in the matched wt sample is >99%, meaning the SNP was introduced   #
    #                      by editing                                                                 #
    # Inputs: 1. var_df - the annotated varaints table                                                #
    #         2. sample_map - the dataframe matching sample identifiers to biological conditions      #
    #         3. replicates - the number of replicates                                                #
    #         4. wt_condition - the reference biological condition                                    #
    #         5. target_snp - the snp being analyzed in this notebook                                 #
    # Output: a list of filter names and a dictionary with the following structure, every sample id   #
    #         in the map is a key:                                                                    #
    #         filter_dict                                                                          #
    #         |_ key: <sample_id>                                                                     #
    #            |_value: dictionary of filter masks for <sample id>                                  #
    #              |_key: <filter name>                                                               #
    #                |_value: a filter mask for var_df matching <filter name>                         # 
    ###################################################################################################

    filter_dict = {}
    combined_filter_dict = {}
    filter_names = [
        'var',
        'VOI',
        'DP',
        'GQ',
        'non-wt',
    ]

    for map_idx, map_row in sample_map.iterrows():
        rep = map_row['rep']
        sample = map_row['sample']
        wt_condition = map_row['wt_condition']

        filter_dict[sample] = {}
        
        # var filter components
        filter_list = [
            var_df[f'{sample}.binGT'] != '0/0',
            var_df[f'{sample}.binGT'] != '0|0',
            var_df[f'{sample}.binGT'] != './.',
            var_df[f'{sample}.binGT'] != '.|.',
            ~(var_df[f'{sample}.binGT'].apply(pd.isna)),
        ]
        filter_dict[sample]['var'] = reduce(lambda x,y: x&y, filter_list)

        # VOI filter components
        filter_list = [
            var_df[f'{sample}_ref'] == target_snp[0],
            (var_df[f'{sample}_all_1'] == target_snp[1]) | (var_df[f'{sample}_all_2'] == target_snp[1])
        ]
        filter_dict[sample]['VOI'] = reduce(lambda x,y: x&y, filter_list)

        # DP and GP filters
        filter_dict[sample]['DP'] = var_df[f'{sample}.DP'] >= 20
        if mutect:
            filter_dict[sample]['GQ'] = pd.Series([True]*len(var_df), index=var_df.index)
        else:
            var_df[f'{sample}.GQ'] >= 20

        filter_dict[sample]['GQ'] = var_df[f'{sample}.GQ'] >= 20

        # non-wt filter
        search_mask = sample_map['condition'] == wt_condition
        search_mask = search_mask & (sample_map['rep'] == rep)
        wt_sample = sample_map.loc[search_mask, 'sample'].iloc[0]
        filter_dict[sample]['non-wt'] = var_df[f'{wt_sample}_pct_ref'] > 99

        # combined filters
        combined_filter_names = []
        combined_filter_dict[sample] = {}
        for filter_idx in range(1,len(filter_names)+1):
            filters_names_to_combine = filter_names[:filter_idx]
            filters_to_combine = [filter_dict[sample][f] for f in filters_names_to_combine]

            combined_filter_name = '_'.join(filters_names_to_combine)
            combined_filter = reduce(lambda x,y: x&y, filters_to_combine)

            combined_filter_names.append(combined_filter_name)
            combined_filter_dict[sample][combined_filter_name] = combined_filter
        


    return filter_dict, filter_names, combined_filter_dict, combined_filter_names

def calc_sample_mean(df_row, col_type, condition, sample_map_df):
    # Input:
    #   - df_row: row of a var_df or similar
    #   - col_type: target column type that is formulaically used with sample i.e '_pct_snp' for '{sample}_pct_snp
    #   - condition: target condition from the sample map csv
    #   - sample_map_df: sample map dataframe 
    # Output:
    #   - Average value for target column type
    
    condition_mask = sample_map_df['condition'] == condition
    samples = sample_map_df.loc[condition_mask, 'sample']
    
    if not samples.empty:
        rep_values = df_row[[f'{sample}{col_type}' for sample in samples]]
        return rep_values.mean()
    else:
        return np.nan
    
def determine_gene_info(gtf_hit_list, info_type):
    if (len(gtf_hit_list) == 0):
        gene_info =''
    else:
        gene_info = set()
        for gtf_dict in gtf_hit_list:
            if info_type == 'gene_name':
                gene_info.add(gtf_dict['gene'][1])
            elif info_type == 'gene_id':
                gene_info.add(gtf_dict['gene'][0])
        
        gene_info = list(gene_info)
        gene_info.sort()
        gene_info = ','.join(gene_info)
    return gene_info


###############
# main script #
###############
if __name__ == '__main__':
    #--- Step 1: define variables ---
    tqdm.tqdm.pandas()
    pandarallel.initialize(progress_bar=True, verbose=2, nb_workers=8)

    # from args
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    var_dir = sys.argv[3]
    # target SNPs as a list of stings for example we are interested in C-to-T (C-to-U) so ['CT']
    target_snp = (sys.argv[4][0], sys.argv[4][1])
    # mutect mode or no
    mutect = True if sys.argv[5] == 'True' else False

    # directories
    hit_dir = os.path.join(var_dir, 'hit-id')
    out_dir = os.path.join(hit_dir,f'{target_snp[0]}{target_snp[1]}')

    os.makedirs(hit_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # sample map
    sample_map_df = pd.read_csv(os.path.join(proj_dir,'sample-map.csv'))
    # max number of biological replicates
    replicates = len(sample_map_df['rep'].unique())

    print('Loading var_df....')
    var_pq_path = os.path.join(var_dir, f'{experiment}.variants.annotated.parquet')
    var_df = pd.read_parquet(var_pq_path)
    print('var_df Loaded.\n')

    # --- Step 2: Filter for hits ---
    print('Filtering Variants for Off-Target Hits....\n')
    # build summary table showing number of entries passing each sucessive filter
    filter_dict, filter_names, combined_filter_dict, combined_filter_names = generate_filter_list(var_df=var_df, sample_map=sample_map_df, target_snp=target_snp)

    filter_counts_df = pd.DataFrame()

    # add 'subset' column to filter counts_df that name every filtered subset then make it the index
    subset_list = ['total_entries']
    for bio_condition in sample_map_df['condition'].unique():
        for filt_name in combined_filter_names:
            subset_list.append(f'{bio_condition}_{filt_name}')
    filter_counts_df['subset'] = subset_list
    filter_counts_df = filter_counts_df.set_index('subset', drop=True)

    # for each biological replicate tabulate the hits after each filter then add the column of counts to filter_counts_df
    for rep in range(1, replicates + 1):
        count_col_list = []
        count_col_list.append(len(var_df))

        for bio_condition in sample_map_df['condition'].unique():
            sample_mask = sample_map_df['condition'] == bio_condition
            # calculate hits for each replicate, if sample groups have different rep numbers columns will have null values after final rep
            if rep in sample_map_df.loc[sample_mask, 'rep'].values:
                sample_mask = sample_mask & (sample_map_df['rep'] == rep)
                sample = sample_map_df.loc[sample_mask, 'sample'].iloc[0]
                sample_filters = combined_filter_dict[sample]
                
                for filt_name in combined_filter_names:
                    count_col_list.append(len(var_df[sample_filters[filt_name]]))
            else:
                count_col_list = count_col_list + ([np.nan]*len(sample_filters))

        filter_counts_df[f'r{rep}'] = count_col_list

    # the last column in filter_counts_df will hold counts of var_df entries where all replicates passed a given filter
    count_col_list = []
    count_col_list.append(len(var_df))

    for bio_condition in sample_map_df['condition'].unique():
        condition_mask = sample_map_df['condition'] == bio_condition
        
        for filt_name in combined_filter_names:
            rep_filters = []
            
            for sample in sample_map_df.loc[condition_mask, 'sample']:
                # combined filter for the sample is added to rep_filters
                rep_filters.append(combined_filter_dict[sample][filt_name])

            # add a count of how many entries passed all replicate filters
            pass_count = len(var_df[reduce(lambda x,y: x&y, rep_filters)])
            count_col_list.append(pass_count)
            
    filter_counts_df['intersection'] = count_col_list
    print(tabulate(filter_counts_df, headers='keys', tablefmt='fancy_grid'))
    print('')
    print('Off-Target Hits Identified.\n')

    # --- Step 3: Write Hit Outputs ---
    print('Writing Outputs...\n')

    # get mask for all final hits
    final_filt_name = combined_filter_names[-1]
    non_wt_mask = ~(sample_map_df['condition'].isin(sample_map_df['wt_condition']))
    non_wt_conditions = sample_map_df.loc[non_wt_mask, 'condition'].unique()

    all_hits_mask = pd.Series([False]*len(var_df), index=var_df.index)
    for condition in non_wt_conditions:
        sample_mask = sample_map_df['condition'] == condition
        condition_hits_mask = [combined_filter_dict[sample][final_filt_name] for sample in sample_map_df.loc[sample_mask, 'sample']]
        condition_hits_mask = reduce(lambda x,y: x&y, condition_hits_mask)
        all_hits_mask = all_hits_mask | condition_hits_mask 

    # calculate mean pct_snp and pct_ref for all conditions
    for condition in sample_map_df['condition'].unique():
        for col_type in ['_pct_ref', '_pct_snp']:
            print(f'Calculating mean{col_type} for {condition}\n')
            mean_col = f'{condition}_mean{col_type}'
            var_df.loc[all_hits_mask , mean_col] = var_df[all_hits_mask].parallel_apply(
                calc_sample_mean, 
                col_type=col_type, condition=condition, sample_map_df=sample_map_df,
                axis=1
            )

    # add hit mask columns
    for condition in non_wt_conditions:
        condition_mask = sample_map_df['condition'] == condition
        bio_condition_samples = sample_map_df.loc[condition_mask, 'sample'].to_list()
        
        hit_filter = [combined_filter_dict[sample][final_filt_name] for sample in bio_condition_samples]
        hit_filter = reduce(lambda x,y: x&y, hit_filter)
        hit_col = f'{condition}_hit'

        var_df[hit_col] = False
        var_df.loc[hit_filter, hit_col] = True

    # add columns with legible gene identifiers
    var_df.loc[all_hits_mask, 'gene_id'] = var_df.loc[all_hits_mask, 'gtf_hits'].apply(determine_gene_info, info_type='gene_id')
    var_df.loc[all_hits_mask, 'gene_name'] = var_df.loc[all_hits_mask, 'gtf_hits'].apply(determine_gene_info, info_type='gene_name')

    # for each biological condition and each filter output a parquet of all entries passing the filter these will be grouped into subdirectories for each condition
    # the count matrix will be output as a tsv in the top level output directory
    for bio_condition in sample_map_df['condition'].unique():    
        print(f'{bio_condition}:\n')

        bio_condition_out_dir = os.path.join(out_dir, bio_condition)
        os.makedirs(bio_condition_out_dir, exist_ok=True)

        condition_mask = sample_map_df['condition'] == bio_condition
        bio_condition_samples = sample_map_df.loc[condition_mask, 'sample'].to_list()

        for filt_idx, combined_filter_name in enumerate(combined_filter_names):
            if combined_filter_name != 'var':
                
                rep_filters = [combined_filter_dict[sample][combined_filter_name] for sample in bio_condition_samples]
                combined_filter = reduce(lambda x,y: x&y, rep_filters)

                if len(var_df[combined_filter]) > 0:
                    out_pq_path = os.path.join(bio_condition_out_dir, f'{bio_condition}.{combined_filter_name}.parquet')
                    var_df[combined_filter].to_parquet(out_pq_path, compression='gzip')

                    if filt_idx == len(combined_filter_names) - 1:
                        wt_condition = sample_map_df.loc[condition_mask, 'wt_condition'].iloc[0]
                        wt_condition_mask = sample_map_df['condition'] == wt_condition
                        wt_samples = sample_map_df.loc[wt_condition_mask,'sample'].to_list()
                        
                        out_pq_path = os.path.join(var_dir, f'{experiment}.{bio_condition}.variants.{target_snp[0]}{target_snp[1]}-hits.parquet')
                        var_df[combined_filter].to_parquet(out_pq_path, compression='gzip')

                        out_col_list = ['chrom', 'pos', 'ref', 'alt']
                        for col_type in ['.DP', '_pct_ref', '_pct_snp']:
                            for sample in (wt_samples + bio_condition_samples):
                                out_col_list.append(f'{sample}{col_type}')

                        for col_type in ['_pct_ref', '_pct_snp']:
                            for condition in [wt_condition, bio_condition]:
                                mean_col = f'{condition}_mean{col_type}'
                                out_col_list.append(mean_col)

                        out_tsv_path = os.path.join(out_dir, f'{experiment}.{bio_condition}.variants.{target_snp[0]}{target_snp[1]}-hits.tsv')
                        var_df.loc[combined_filter, out_col_list].to_csv(out_tsv_path, sep='\t', index=False)

                print(f'\tFilter:     \t{combined_filter_name}')
                print(f'\ttsv entries:\t{len(var_df[combined_filter])}\n')
            
            out_tsv_path = os.path.join(var_dir, f'filter-counts-{target_snp[0]}{target_snp[1]}.tsv')
            filter_counts_df.to_csv(out_tsv_path, sep='\t', index=False, float_format='%.2f')

    # write overlap df that contains all hits with hit mask cols and retains var_df idx and an intersection .xlsx with counts
    overlap_df = var_df.loc[all_hits_mask].copy()

    intersection_df_dict = {}
    for condition in non_wt_conditions:
        other_conditions = pd.Series(non_wt_conditions)[~(pd.Series(non_wt_conditions) == condition)].to_list()

        intersection_df_idx = []
        for set_size in range(1,len(other_conditions)+1):
            size_combinations = list(combinations(other_conditions, set_size))
            combinations_csv = [','.join(sorted(list(c))) for c in size_combinations]
            intersection_df_idx += combinations_csv

        intersection_df_idx = [condition] + intersection_df_idx
        col_vals = []
        for combination in intersection_df_idx:
            intersection_conditions = set(combination.split(',') + [condition])
            intersection_conditions = list(intersection_conditions)
            intersection_conditions.sort()
            
            non_intersection_mask = ~(pd.Series(non_wt_conditions).isin(intersection_conditions))
            non_intersection_conditions = list(pd.Series(non_wt_conditions)[non_intersection_mask])
            
            intersection_mask = overlap_df[[f'{c}_hit' for c in intersection_conditions]].all(axis=1)
            intersection_mask = intersection_mask & ~(overlap_df[[f'{c}_hit' for c in non_intersection_conditions]]).any(axis=1)

            col_vals.append(len(overlap_df[intersection_mask]))
        
        col_vals.append(sum(col_vals))
        intersection_df_idx.append('Total')
        
        intersection_df = pd.DataFrame({condition:col_vals}, index=intersection_df_idx)
        intersection_df = intersection_df.sort_values(by=[condition], ascending=False)
        empty_intersection_mask = intersection_df[condition] == 0
        intersection_df = intersection_df.drop(intersection_df[empty_intersection_mask].index)

        intersection_df_dict[condition] = intersection_df

    overlap_df_out_path = f'{hit_dir}/{experiment}.overlap.{target_snp[0]}{target_snp[1]}.parquet'
    overlap_df.to_parquet(overlap_df_out_path)

    xlsx_out_path = f'{hit_dir}/{experiment}.overlap.{target_snp[0]}{target_snp[1]}.xlsx'
    with pd.ExcelWriter(xlsx_out_path, engine='openpyxl') as writer:
        for condition, intersection_df in intersection_df_dict.items():
            intersection_df.to_excel(writer, sheet_name=f'{condition}-{target_snp[0]}{target_snp[1]}', index=True)

    print('Outputs Written.\n')