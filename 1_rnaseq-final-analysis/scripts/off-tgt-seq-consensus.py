#####################
# import statements #
#####################

import os
import pandas as pd
import sys
import pysam
import tqdm
from Bio import Seq
from functools import reduce
import numpy as np
import RNA
from collections import defaultdict
from tabulate import tabulate

#############
# functions #
#############
def load_sample_map(proj_dir):
    sample_map_df = pd.read_csv(os.path.join(proj_dir,'sample-map.csv'))
    
    non_wt_mask = ~sample_map_df['condition'].isin(sample_map_df['wt_condition'])
    non_wt_conditions = list(sample_map_df.loc[non_wt_mask, 'condition'].unique())
    non_wt_conditions.sort()

    return sample_map_df, non_wt_conditions

def get_norm_seq_context(transcript_row, offset, transcript_fa):
    transcript_contig = transcript_row['transcript_contig']
    transcript_pos = transcript_row['transcript_pos']

    transcript_seq = transcript_fa.fetch(reference=transcript_contig)

    seq_context_start = max(1, transcript_pos-offset)
    seq_context_end = min(len(transcript_seq), transcript_pos+offset)

    start_missing = max(0, -1*(transcript_pos-offset-1))
    end_missing = max(0, transcript_pos+offset-len(transcript_seq))

    seq_context = transcript_seq[seq_context_start-1:seq_context_end]
    seq_context = f'{''.join(['-']*start_missing)}{seq_context}{''.join(['-']*end_missing)}'
    
    return seq_context

def check_seq_series(seq_series, idx_tuple=(None, None), RNA=False):
    if len(seq_series.apply(len).value_counts()) != 1:
        print('[ERROR] requires all sequences to be the same length')
        return pd.Series(), -1, []

    seq_series = seq_series.copy()
    na_mask = ~seq_series.apply(pd.isna)
    seq_series[na_mask] = seq_series[na_mask].apply(str.upper)
    seq_len = len(seq_series[na_mask].iloc[0])

    if idx_tuple == (None, None):
        idx_tuple = (0,seq_len)

    seq_series = seq_series.apply(lambda x: x[idx_tuple[0]:idx_tuple[1]])

    if RNA:
        seq_series = seq_series.apply(lambda x: str(Seq.Seq(x).transcribe()))
        base_list = ['A', 'C', 'G', 'U', '-']
    else:
        base_list = ['A', 'C', 'G', 'T', '-']
    
    return seq_series, seq_len, base_list

# Given a sequence return a dictionary
# keys = bases
# values = a list of 1s and 0s corresponding to where each base occurs in the sequence
def count_bases(seq):
    count_dict = {}
    for seq_idx, base in enumerate(list(seq)):
        if base not in count_dict.keys():
            count_dict[base] = [0]*len(seq)
        count_dict[base][seq_idx] += 1

    return count_dict

# given a series of sequences return a dataframe with total counts for each base at each position
# columns = bases
# rows = column base count at index
def calc_base_counts(seq_series, idx_tuple=(None,None), RNA=False):
    
    seq_series, seq_len, base_list = check_seq_series(seq_series=seq_series, idx_tuple=idx_tuple, RNA=RNA)

    count_dict_series = seq_series.apply(count_bases)
    count_df = pd.DataFrame()

    for base in base_list:
        total_base_counts = count_dict_series.apply(lambda x: np.array(x.get(base, [0]*seq_len)))
        total_base_counts = reduce(lambda x,y: np.add(x,y), total_base_counts.values)
        count_df[base] = total_base_counts
    
    return count_df


# given a series of sequences return a dataframe with a base frequency matrix 
# columns = bases
# rows = column base count at index
def calc_base_frequency(seq_series, idx_tuple=(None,None), RNA=False, logo=False):

    seq_series, seq_len, base_list = check_seq_series(seq_series=seq_series, idx_tuple=idx_tuple, RNA=RNA)

    count_df = calc_base_counts(seq_series=seq_series, RNA=RNA)
    
    if logo:
        freq_df = count_df.drop(columns=['-'])
        base_list.remove('-')
        
    total_series = freq_df[base_list].sum(axis=1)
    for base_col in freq_df.columns:
        freq_df[base_col] = freq_df[base_col]/total_series
    
    return count_df, freq_df

def calc_struct_frequency(struct_series):
    struct_len = struct_series.apply(len).iloc[0]
    
    struct_count_series = struct_series.apply(count_bases)
    char_set = struct_count_series.apply(lambda x: list(x.keys()))
    char_set = list(set(char_set.sum()))
    char_set.sort()

    count_df = pd.DataFrame()
    for struct_char in char_set:
        total_struct_counts = struct_count_series.apply(lambda x: np.array(x.get(struct_char, [0]*struct_len)))
        total_struct_counts = reduce(lambda x,y: np.add(x,y), total_struct_counts.values)
        count_df[struct_char] = total_struct_counts
    
    freq_df = count_df/len(struct_series)

    return count_df, freq_df

###############
# main script #
###############
if __name__ == '__main__':
    tqdm.tqdm.pandas()
    
    #--- Step 1: define variables ---
    # from args
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    ref_dir = sys.argv[3]
    # target SNP as a tuple (ref, alt)
    target_snp = (sys.argv[4][0],sys.argv[4][1])
    # offset
    offset = int(sys.argv[5])

    # based on args
    var_dir = os.path.join(proj_dir, '4_variants')
    align_dir = os.path.join(var_dir, 'transcript-align')
    align_snp_dir = os.path.join(align_dir, f'{target_snp[0]}{target_snp[1]}')
    expression_dir = os.path.join(proj_dir, '5_expression')
    salmon_dir = os.path.join(expression_dir, 'salmon')
    salmon_analysis_dir = os.path.join(salmon_dir, 'analysis')
    seq_consensus_dir = os.path.join(var_dir, 'seq-context')
    seq_snp_dir = os.path.join(seq_consensus_dir, f'{target_snp[0]}{target_snp[1]}-o{offset}')
    
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')
    genome_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.genome.fa.gz')
    genome_gtf_path = os.path.join(ref_dir, f'{ref_genome_name}.sorted.gtf.gz')
    
    sample_map_df, non_wt_conditions = load_sample_map(proj_dir=proj_dir)
    
    #--- Step 1: import transcript dfs ---
    print(f'Importing Transcript dfs for {experiment}...')
    
    transcript_df_dict = defaultdict(
        lambda: defaultdict(dict)
    )

    for condition in non_wt_conditions:
        transcript_df_path = os.path.join(align_snp_dir, f'{experiment}.{condition}.transcripts.{target_snp[0]}{target_snp[1]}.parquet')
        transcript_df = pd.read_parquet(transcript_df_path)
        transcript_df_dict[condition]['path'] = transcript_df_path
        transcript_df_dict[condition]['df'] = transcript_df

    transcript_df_dict=dict(transcript_df_dict)
    
    print(f'Imported Transcript dfs for {experiment}.\n')

    #--- Step 2: get sequence context and analyze for consensus---
    output_category = ['seq', 'struct']
    output_types = ['enrich', 'freq', 'count']
    output_dict = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(dict)
        )
    )
    for condition, condition_dict in transcript_df_dict.items():
        transcript_df = condition_dict['df']
        non_out_transcript_cols = pd.Series(transcript_df.columns)
        
        # get sequence context based on a given offset
        print(f'Getting Sequence Context (Offset={offset}) for {len(transcript_df)} Rows in {condition} transcript_df...')
        
        offset_col = f'o{offset}_seq_context'
        transcript_df[offset_col] = transcript_df.progress_apply(get_norm_seq_context, offset=offset, transcript_fa=pysam.FastaFile(transcript_fa_path), axis=1, result_type='reduce')
            
        print(f'\nSequence Context Length Value Counts (Offset={offset}):')
        print(tabulate(pd.DataFrame(transcript_df[offset_col].apply(len).value_counts()), headers='keys', tablefmt='fancy_grid'))
        print(f'\nSequence Context Central Base Counts (Offset={offset}):')
        print(tabulate(pd.DataFrame(transcript_df[offset_col].apply(lambda x: x[offset]).value_counts()), headers='keys', tablefmt='fancy_grid'))
        
        print(f'\nGot Sequence Context (Offset={offset}) for {condition}.\n')

        # calculate secondary structures for seq contexts
        print(f'Calculating Sequence Context Secondary Structure (Offset={offset}) for {condition}...')
        offset_col = f'o{offset}_struct_context'
        transcript_df[offset_col] = transcript_df[f'o{offset}_seq_context'].progress_apply(lambda x: RNA.fold(str(Seq.Seq(x).transcribe()))[0])
        
        print(f'\nStructure Value Counts (Offset={offset}):')
        print(tabulate(pd.DataFrame(transcript_df[offset_col].value_counts()).head(), headers='keys', tablefmt='fancy_grid'))

        print(f'\nCalculated Sequence Context Secondary Structure (Offset={offset}) for {condition}.\n')

        # Calculate base frequencies for consensus logo
        print(f'Calculating Base Frequencies (Offset={offset}) for {condition}...')
        base_list = ['A', 'C', 'G', 'U']
        not_na_mask = ~(transcript_df[f'o{offset}_seq_context'].apply(pd.isna)) & ~(transcript_df[f'o{offset}_seq_context'] == '')
        
        if not transcript_df[not_na_mask].empty:
            base_count_df, base_freq_df = calc_base_frequency(transcript_df.loc[not_na_mask, f'o{offset}_seq_context'], RNA=True, logo=True)
            
            base_enrich_df = pd.DataFrame()
            for base in base_list:
                base_enrich_df[base] = (base_freq_df[base]-0.25)/0.25
                base_enrich_df[base] = base_enrich_df[base].apply(lambda x: (max(0,x)))
            
            base_count_df = base_count_df.set_index(base_count_df.index-offset)
            base_freq_df = base_freq_df.set_index(base_freq_df.index-offset)
            base_enrich_df = base_enrich_df.set_index(base_enrich_df.index-offset)
            
            output_dict['seq']['count'][condition] = base_count_df
            output_dict['seq']['freq'][condition] = base_freq_df
            output_dict['seq']['enrich'][condition] = base_enrich_df
        else:
            output_dict['seq']['count'][condition] = pd.DataFrame()
            output_dict['seq']['freq'][condition] = pd.DataFrame()
            output_dict['seq']['enrich'][condition] = pd.DataFrame()

        print(f'Calculated Base Frequencies (Offset={offset}) for {condition}.\n')

        # Calculate structural frequencies for consensus
        print(f'Calculating Structure Frequencies (Offset={offset}) for {condition}...')
        not_na_mask = ~(transcript_df[f'o{offset}_struct_context'].apply(pd.isna)) & ~(transcript_df[f'o{offset}_struct_context'] == '')

        if not transcript_df[not_na_mask].empty:
            struct_count_df, struct_freq_df = calc_struct_frequency(transcript_df[f'o{offset}_struct_context'])
            
            baseline_freq = (1/len(struct_freq_df.columns))
            struct_enrich_df = (struct_freq_df - baseline_freq)/baseline_freq
            struct_enrich_df[struct_enrich_df < 0] = 0

            struct_count_df = struct_count_df.set_index(struct_count_df.index-offset)
            struct_freq_df = struct_freq_df.set_index(struct_freq_df.index-offset)
            struct_enrich_df = struct_enrich_df.set_index(struct_enrich_df.index-offset)

            output_dict['struct']['count'][condition] = struct_count_df
            output_dict['struct']['freq'][condition] = struct_freq_df
            output_dict['struct']['enrich'][condition] = struct_enrich_df
        else:
            output_dict['struct']['count'][condition] = pd.DataFrame()
            output_dict['struct']['freq'][condition] = pd.DataFrame()
            output_dict['struct']['enrich'][condition] = pd.DataFrame()
        
        print(f'Calculated Structure Frequencies (Offset={offset}) for {condition}.\n')

        print(f'Writing Transcript Context File for {condition}...')
        out_transcript_path = os.path.join(align_snp_dir, f'{experiment}.{condition}.transcripts.{target_snp[0]}{target_snp[1]}.o{offset}.context.parquet')
        transcript_cols = pd.Series(transcript_df.columns)
        out_transcript_col_mask = ~(transcript_cols.isin(non_out_transcript_cols))
        out_transcript_cols = transcript_cols[out_transcript_col_mask]
        
        transcript_df.to_parquet(out_transcript_path, compression='gzip')
        transcript_df_dict[condition]['df'] = transcript_df
        print(f'Wrote Transcript Context File for {condition}.')

    output_dict = dict(output_dict)

    #--- Step 3: Write Outputs ---
    excel_path = os.path.join(var_dir, f'{experiment}.transcripts.{target_snp[0]}{target_snp[1]}.o{offset}-consensus.xlsx')
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        conditions = list(transcript_df_dict.keys())
        conditions.sort()

        for category in output_category:
            for out_type in output_types:
                for condition in conditions:
                    df = output_dict[category][out_type][condition]
                    df.to_excel(writer, sheet_name=f'{condition}_{category}_{out_type}'[:31], index=True)

                    out_path = os.path.join(seq_snp_dir, f'{experiment}.{condition}.seq.o{offset}-{category}-{out_type}.parquet')
                    df.to_parquet(out_path, compression='gzip')
                    
