#####################
# import statements #
#####################
import os
import pandas as pd
import sys
import tqdm
from tabulate import tabulate
import re
import gzip
import pysam
import pickle
from pandarallel import pandarallel
from functools import reduce
import random
import gc

#############
# functions #
#############

def load_sample_map(proj_dir):
    sample_map_df = pd.read_csv(os.path.join(proj_dir,'sample-map.csv'))
    
    non_wt_mask = ~sample_map_df['condition'].isin(sample_map_df['wt_condition'])
    non_wt_conditions = list(sample_map_df.loc[non_wt_mask, 'condition'].unique())
    non_wt_conditions.sort()

    return sample_map_df, non_wt_conditions

def get_norm_seq_context(transcript_row, pos_col, offset, transcript_fa):
    transcript_contig = transcript_row['transcript_contig']
    transcript_pos = transcript_row[pos_col]

    if transcript_pos == -1:
        seq_context = ''
    else:
        transcript_seq = transcript_fa.fetch(reference=transcript_contig)

        seq_context_start = max(1, transcript_pos-offset)
        seq_context_end = min(len(transcript_seq), transcript_pos+offset)

        start_missing = max(0, -1*(transcript_pos-offset-1))
        end_missing = max(0, transcript_pos+offset-len(transcript_seq))

        seq_context = transcript_seq[seq_context_start-1:seq_context_end]
        seq_context = f'{''.join(['-']*start_missing)}{seq_context}{''.join(['-']*end_missing)}'
    
    return seq_context

def find_random_seq_position(transcript_contig, transcript_fa, seq_motif, central_base_idx):
    
    transcript_seq = transcript_fa.fetch(reference=transcript_contig).upper()
    pos = random.randint(0 + central_base_idx, len(transcript_seq) - (len(seq_motif) - central_base_idx))
    if seq_motif != 'N':
        loop_break = 0
        while transcript_seq[pos-central_base_idx:pos+(len(seq_motif)-central_base_idx)] != seq_motif:
            pos = random.randint(0, len(transcript_seq))
            loop_break += 1
            if loop_break > 10000:
                break
        if loop_break > 10000:
            if 'N' in seq_motif:
                seq_pattern = seq_motif.replace('N', '[ACGT]')
            else:
                seq_pattern = seq_motif
            
            seq_pattern = re.compile(seq_pattern)
            matches = list(seq_pattern.finditer(transcript_seq))
            if len(matches) > 0:
                pos = matches[random.randint(0, len(matches)-1)].start() + central_base_idx
            else:
                pos = -2
    return pos + 1 # dna positions are index 1

###############
# main script #
###############
if __name__ == '__main__':        
    tqdm.tqdm.pandas()
        
    #--- Step 1: define variables ---
    # from args
    experiment_list = sys.argv[1].split(',')
    proj_dir = sys.argv[2]
    exp_dir_list = sys.argv[3].split(',')
    ref_dir = sys.argv[4]
    resource_dir = sys.argv[5]

    # other directory paths
    exp_dir_map ={experiment_list[i]:exp_dir_list[i] for i in range(len(experiment_list))}
    out_dir = os.path.join(resource_dir, 'expr-transcripts')

    os.makedirs(out_dir, exist_ok=True)

    # ref file paths
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')
    transcript_contig_path = os.path.join(
        exp_dir_list[0], '4_variants', 'transcript-align', 'transcript-id-to-contig.pkl.gz'
    )

    # sample map
    non_wt_set = set()
    sample_map_df_list = []
    for experiment, exp_dir in exp_dir_map.items():
        sample_map_df, non_wt_conditions = load_sample_map(proj_dir=exp_dir)
        sample_map_df_list.append(sample_map_df)
        non_wt_set = non_wt_set | set(non_wt_conditions)
    sample_map_df = pd.concat(sample_map_df_list)
    non_wt_conditions = sorted(non_wt_set)
    print(sample_map_df)
    print(non_wt_conditions)

    # output paths
    transcript_expression_pq_path = os.path.join(out_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.parquet') 
    random_pos_pq_path = os.path.join(out_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.random.parquet')
    print(random_pos_pq_path)

    #---Step 2: Import Transcript Level Raw Expression Data----
    
    print (f'\nLoading and Merging Transcript Expression Data for {', '.join(experiment_list)}...')

    merge_keys = set()
    expr_df_list = []
    for experiment, exp_dir in exp_dir_map.items():
        transcript_raw_path = os.path.join(exp_dir, '5_expression', 'salmon', 'analysis', 'salmon-transcript-raw-counts.tsv.gz')
        transcript_expr_df = pd.read_csv(transcript_raw_path, sep='\t', compression='infer', index_col='Unnamed: 0')
        # transcript_expr_df = transcript_expr_df.set_index('transcript_id', drop=True)
        if len(merge_keys) == 0:
            merge_keys = set(transcript_expr_df.columns)
        else:
            merge_keys = merge_keys.intersection(set(transcript_expr_df.columns))
        expr_df_list.append(transcript_expr_df)
    
    transcript_expr_df = reduce(
        lambda x,y: pd.merge(left=x, right=y, on=sorted(merge_keys), how='outer'),
        expr_df_list
    )
    transcript_expr_df = transcript_expr_df.set_index('transcript_id')
    for col in transcript_expr_df.columns:
        if ('_reads' in col) or ('_tpm' in col):
            transcript_expr_df[col] = transcript_expr_df[col].fillna(0)

    keep_cols = ['gene_id', 'gene_name']

    print (f'\nLoaded and Merged Transcript Expression Data for {', '.join(experiment_list)}.')

    print (f'Determining Expressed Transcripts (DP >= 20)...\n')
    expressed_cols = list(f'{c}_expressed' for c in non_wt_conditions)
    if not reduce(lambda x,y: x and y, [(c in transcript_expr_df.columns) for c in expressed_cols]):
        for condition in non_wt_conditions:
            expressed_col = f'{condition}_expressed' 

            sample_mask = sample_map_df['condition'] == condition
            samples = sample_map_df.loc[sample_mask, 'sample']
            sample_cols = [f'{sample}_reads' for sample in samples]
            expressed_mask = (transcript_expr_df[sample_cols] >= 20).all(axis=1)

            transcript_expr_df[expressed_col] = expressed_mask
        
        not_expressed_mask = ~(transcript_expr_df[expressed_cols].any(axis=1))
        not_expressed_idx = transcript_expr_df[not_expressed_mask].index

        keep_cols += expressed_cols
        drop_cols = pd.Series(transcript_expr_df.columns)[~(pd.Series(transcript_expr_df.columns).isin(keep_cols))].to_list()

        print (f'{len(transcript_expr_df)} Initial Transcripts.')
        print (f'{len(not_expressed_idx)} Transcripts not Expressed in Any Condition, Dropped.\n')
        transcript_expr_df = transcript_expr_df.drop(not_expressed_idx)
        transcript_expr_df = transcript_expr_df.drop(columns=drop_cols)

    # add contig column for expressed transcripts
    if ('transcript_contig' not in transcript_expr_df.columns) or ('transcript_len' not in transcript_expr_df.columns):
        transcript_fa = pysam.FastaFile(transcript_fa_path)
        with gzip.open(transcript_contig_path, 'rb') as f:
            transcript_contig_dict = pickle.load(f)

        contig_series = pd.Series((['']*len(transcript_expr_df)), index=transcript_expr_df.index)
        len_series = pd.Series(([int(-1)]*len(transcript_expr_df)), index=transcript_expr_df.index)
        
        for transcript_contig in transcript_fa.references:
            transcript_id = transcript_contig.strip().split('|')[0]
            if transcript_id in transcript_expr_df.index:
                contig_series[transcript_id] = transcript_contig
                len_series[transcript_id] = int(len(transcript_fa.fetch(reference=transcript_contig)))
        
        transcript_expr_df['transcript_contig'] = contig_series
        transcript_expr_df['transcript_len'] = len_series
        
        no_contig_mask = transcript_expr_df['transcript_contig'].apply(lambda x: x == '')
        no_len_mask = transcript_expr_df['transcript_len'].apply(lambda x: x == -1)
        no_contig_idx = transcript_expr_df[no_contig_mask | no_len_mask].index

        print (f'{len(no_contig_idx)} Transcripts with No Contig Found, Dropped.\n')

        transcript_expr_df = transcript_expr_df.drop(no_contig_idx)

        # garbage collect
        del contig_series
        del len_series
        del transcript_contig_dict
        del transcript_fa
        gc.collect()

    for condition in non_wt_conditions:
        print(f'{len(transcript_expr_df[transcript_expr_df[f'{condition}_expressed']])}\tExpressed Transcripts for {condition}')
    print (f'\nDetermined Expressed Transcripts (DP >= 20).\n')

    #---Step 3: Get Random Positions in Expressed Transcripts for Background Controls----
    print(f'Getting Random Positions in Expressed Transcripts...\n')
    transcript_fa = pysam.Fastafile(transcript_fa_path)
    random_seq_motifs = {
        'N':0,
        'C':0,
        'TC':1,
        'TCG':1
    }

    if os.path.isfile(random_pos_pq_path):
        print('Existing Random pos Columns Loaded.')
        random_pos_df = pd.read_parquet(random_pos_pq_path)
    else:
        print('Building New Random pos Columns...')
        random_pos_df = pd.DataFrame()

    new_entry_mask = ~(pd.Series(transcript_expr_df.index).isin(random_pos_df.index)) # apply to transcript_expr_df to find rows not already in random_pos_df
    new_entry_idx = list(pd.Series(transcript_expr_df.index)[new_entry_mask]) # list of indexes to be added to random_pos_df
    if len(new_entry_idx) > 0:
            new_entry_row_df = pd.DataFrame(index=new_entry_idx)
            random_pos_df = pd.concat([random_pos_df, new_entry_row_df])
    print(f'{len(new_entry_idx)}/{len(transcript_expr_df)} transcript_ids In transcript_expr_df not Present in random_pos_df')
    for seq_motif, central_base_idx in random_seq_motifs.items():
        random_col = f'random_{seq_motif}_pos'
        
        print(f'Getting Random {seq_motif} Positions in Expressed Transcripts...')

        # if column is a motif that hasn't had random positions pulled already calculate an entirely new column
        if random_col not in random_pos_df.columns:
            random_pos_df.loc[transcript_expr_df.index, random_col] = transcript_expr_df['transcript_contig'].progress_apply(
                find_random_seq_position,
                transcript_fa=transcript_fa,
                seq_motif=seq_motif, central_base_idx=central_base_idx
            ) # restricted to rows in transcript_expr_df so the lengths match and thats all we care about for this experiment
        else:
            na_mask = random_pos_df[random_col].apply(pd.isna) # check where the missing values are
            na_idx = pd.Series(random_pos_df[na_mask].index) # get indicies with na values
            in_transcr_expr_mask = na_idx.isin(transcript_expr_df.index)
            na_idx = na_idx[in_transcr_expr_mask] # make sure we dont try and calculate these values from an entry that doesn't exist in transcript_expr_df
            random_pos_df.loc[na_idx, random_col] = transcript_expr_df.loc[na_idx, 'transcript_contig'].progress_apply(
                find_random_seq_position,
                transcript_fa=transcript_fa,
                seq_motif=seq_motif, central_base_idx=central_base_idx
            ) # fill the na values
        print(f'Got Random {seq_motif} Positions in Expressed Transcripts.')

    for seq_motif, central_base_idx in random_seq_motifs.items():
        random_col = f'random_{seq_motif}_pos'
        print(f'{random_col} Sequence Check:')
        
        check_df = transcript_expr_df[['transcript_contig']].join(random_pos_df[[random_col]])
        seq_check = check_df.apply(
            lambda x: transcript_fa.fetch(
                reference=x['transcript_contig'],
                start=int(x[random_col]) - 1 - central_base_idx,
                end=int(x[random_col]) + (len(seq_motif) - central_base_idx - 1)
            ) if x[random_col] != -1 else '', axis=1
        ).value_counts()

        print(tabulate(pd.DataFrame(seq_check), headers='keys', tablefmt='fancy_grid'))

        # garbage collection
        del seq_check
        del check_df
        gc.collect()

    print(f'Got Random Positions in Expressed Transcripts.\n')

    # #---Step 4: Write Output ----
    transcript_expr_df.to_parquet(transcript_expression_pq_path, compression='gzip')
    random_pos_df.to_parquet(random_pos_pq_path, compression='gzip')