#####################
# import statements #
#####################

import os
import pandas as pd
import sys
import tqdm
from collections import defaultdict
from tabulate import tabulate
import re
import gzip
from Bio import Seq
import pysam
import pickle
from pandarallel import pandarallel
import random
import gc
import numpy as np
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

def get_norm_seq_context(transcript_row, pos_col, offset, transcript_fa):
    transcript_contig = transcript_row['transcript_contig']
    transcript_pos = transcript_row[pos_col]

    if pd.isna(transcript_pos) or transcript_pos == -1:
        seq_context = ''
    else:
        transcript_pos = int(transcript_pos)  # coerce float64 from left-merge dtype promotion
        transcript_seq = transcript_fa.fetch(reference=transcript_contig)

        seq_context_start = max(1, transcript_pos-offset)
        seq_context_end = min(len(transcript_seq), transcript_pos+offset)

        start_missing = max(0, -1*(transcript_pos-offset-1))
        end_missing = max(0, transcript_pos+offset-len(transcript_seq))

        seq_context = transcript_seq[seq_context_start-1:seq_context_end]
        seq_context = f'{''.join(['-']*start_missing)}{seq_context}{''.join(['-']*end_missing)}'

    return seq_context

def find_tgt_site_matches(seq_context, mirtar_patterns):
    matches = []

    if seq_context != '':
        for mirtar_idx, tgt_site_info in mirtar_patterns.items():
            pattern = tgt_site_info['tgt_site_pattern']
            for match in pattern.finditer(seq_context):
                matches.append({
                    'mirtar_idx': mirtar_idx,
                    'mirtar_id': tgt_site_info['mirtar_id'],
                    'mirna_name': tgt_site_info['mirna_name'],
                    'position': (match.start(), match.end()),
                    'matched_seq': match.group(),
                })
    
    return matches

def find_seed_matches(seq_context, seed_patterns):
    matches = []
    if seq_context != '':
        for mirna_id, mirna_info in seed_patterns.items():
            
            for match_idx, match in enumerate(mirna_info['6-mer_pattern'].finditer(seq_context)):
                seed_type = ''
                if match.group()[7] != 'A':
                    if match.group()[0] != mirna_info['8-mer'][0]:
                        seed_type       = '6-mer'
                        complementary   = (match.start()+1, match.start()+7)
                    else:
                        seed_type       = '7-mer-m8'
                        complementary   = (match.start()+1, match.start()+8)
                elif match.group()[0] != mirna_info['8-mer'][0]:
                    seed_type       = '7-mer-A1'
                    complementary   = (match.start(), match.start()+7)
                else:
                    seed_type       = '8-mer'
                    complementary   = (match.start(), match.start()+8)

                matches.append({
                    'mirna_id': mirna_id,
                    'mirna_name': mirna_info['mirna_name'],
                    'seed_type': seed_type,
                    'position': match.start(),
                    'matched_seq': match.group(),
                    'complementary': complementary
                })
    
    return matches

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

def add_random_sample(sample_size, random_sample_dict, transcript_expr_df, motif):
    random_sample_sizes = list(pd.Series(random_sample_dict.keys()).sort_values(ascending=True))

    if sample_size not in random_sample_sizes:
        if len(random_sample_dict) == 0:
            expr_idx_set = set()
            pick_from = list(transcript_expr_df.index)
        else:
            next_largest_idx = -1
            next_smallest_idx = -1
            for size_idx, random_sample_size in enumerate(random_sample_sizes):
                if sample_size > random_sample_size:
                    next_smallest_idx = size_idx
                else:
                    next_largest_idx = size_idx
                    break
            
            if next_smallest_idx == -1:
                expr_idx_set = set()
            else:
                next_smallest_sample = random_sample_sizes[next_smallest_idx]
                expr_idx_set = set(random_sample_dict[next_smallest_sample])
            
            if next_largest_idx == -1:
                pick_from = list(transcript_expr_df.index)
            else:
                next_largets_sample = random_sample_sizes[next_largest_idx]
                pick_from = random_sample_dict[next_largets_sample]

        random_pos_col = f'random_{motif}_pos'
        while len(expr_idx_set) < sample_size:
            idx_iloc = random.randint(0,len(pick_from)-1)
            random_idx = pick_from[idx_iloc]
            if transcript_expr_df.loc[random_idx, random_pos_col] != -1:
                expr_idx_set.add(random_idx)
        
        sample_entry = list(expr_idx_set)
        sample_entry.sort()
        random_sample_dict[sample_size] = sample_entry
        
    return random_sample_dict

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

def count_seed_type_hits(matches, seed_types=None, mirna_ids=None):
    count = 0
    target_matches = matches
    
    if mirna_ids:
        target_matches = [m for m in target_matches if m['mirna_id'] in mirna_ids]
    if seed_types:
        target_matches = [m for m in target_matches if m['seed_type'] in seed_types]
    
    count = len(target_matches)
    return count

def count_positional_matches(matches, offset, mode='mirtarbase', seed_types=None, mirna_ids=None):
    if len(matches) ==0:
        return np.zeros((offset * 2) + 1, dtype=int)

    if mode == 'mirbase':
        if mirna_ids:
            matches = [m for m in matches if m['mirna_id'] in mirna_ids]
        if seed_types:
            matches = [m for m in matches if m['seed_type'] in seed_types]
        target_positions = [m['complementary'] for m in matches]
    else:
        target_positions = [m['position'] for m in matches]
    
    
    alignment_array = np.zeros((offset * 2) + 1, dtype=int)
    for start, end in target_positions:
        alignment_array[start:end] += 1
    
    return alignment_array

def sum_positional_matches(match_series, offset, mode='mirtarbase', seed_types=None, mirna_ids=None):
    count_array =  np.zeros((offset * 2) + 1, dtype=int)

    for matches in match_series:
        count_array = np.add(
            count_array,
            count_positional_matches(
                matches=matches,
                offset=offset,
                mode=mode,
                seed_types=seed_types,
                mirna_ids=mirna_ids
            )
        )
    
    return count_array    

###############
# main script #
###############
if __name__ == '__main__':
    tqdm.tqdm.pandas()
    pandarallel.initialize(progress_bar=True, verbose=2, nb_workers=8)
        
    #--- Step 1: define variables ---
    # from args
    experiment_list = sys.argv[1].split(',')
    proj_dir = sys.argv[2]
    exp_dir_list = sys.argv[3].split(',')
    ref_dir = sys.argv[4]
    ref_db_dir = sys.argv[5]
    # target SNP as a tuple (ref, alt)
    target_snp = (sys.argv[6][0], sys.argv[6][1])
    offset = int(sys.argv[7])
    resource_dir = sys.argv[8]
    agg_res_dir = sys.argv[9]

    # other directory paths
    exp_dir_map = {experiment_list[i]:exp_dir_list[i] for i in range(len(experiment_list))}
    miRTarBase_dir          = os.path.join(ref_db_dir, 'miRtarBase')
    miRBase_dir             = os.path.join(ref_db_dir, 'miRBase')
    expr_resource_dir       = os.path.join(resource_dir, 'expr-transcripts')
    expr_mirna_resource_dir = os.path.join(expr_resource_dir, 'mirna-align')
    overlap_resource_dir    = os.path.join(resource_dir, 'transcript-overlap')

    # os.makedirs(out_dir, exist_ok=True)
    os.makedirs(agg_res_dir, exist_ok=True)
    os.makedirs(overlap_resource_dir, exist_ok=True)
    os.makedirs(expr_mirna_resource_dir, exist_ok=True)

    # # ref file paths
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')
    tgt_site_patterns_path = os.path.join(miRTarBase_dir, 'miRTarBase.tgt-sites.patterns.pkl.gz')
    seed_patterns_path = os.path.join(miRBase_dir, 'miRBase.seed-patterns.pkl.gz')

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

    # analysis file paths and information
    overlap_df_list = []
    non_out_transcript_cols = set()
    merge_keys = set()
    for experiment, exp_dir in exp_dir_map.items():
        transcript_overlap_path = os.path.join(
            exp_dir, '4_variants', 'transcript-align', f'{experiment}.transcripts.{''.join(target_snp)}.parquet'
        )
        transcript_overlap_df = pd.read_parquet(transcript_overlap_path)
        transcript_overlap_df = transcript_overlap_df.drop(columns='var_df_idx')
        overlap_df_list.append(transcript_overlap_df)
        col_set = set(transcript_overlap_df.columns)
        if len(merge_keys) == 0:
            merge_keys = col_set
        else:
            merge_keys = merge_keys.intersection(col_set)
    
    transcript_overlap_df = reduce(
        lambda x,y: pd.merge(left=x, right=y, on=sorted(merge_keys), how='outer'),
        overlap_df_list
    )
    for col in transcript_overlap_df.columns:
        if ('_hit' in col) or ('_expressed' in col):
            na_mask = transcript_overlap_df[col].apply(pd.isna)
            transcript_overlap_df.loc[na_mask, col] = False

    transcript_expr_path = os.path.join(expr_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.parquet')
    random_pos_pq_path = os.path.join(expr_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.random.parquet')
    random_sample_path = os.path.join(expr_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.random.pkl.gz')
    random_seq_motifs = {
        'N':0,
        'C':0,
        'TC':1,
        'TCG':1
    }
    
    transcript_expr_df = pd.read_parquet(transcript_expr_path)
    random_pos_df = pd.read_parquet(random_pos_pq_path)
    transcript_expr_df = pd.merge(left=transcript_expr_df, right=random_pos_df, left_index=True, right_index=True, how='left')
    non_out_transcript_expr_cols = set(transcript_expr_df.columns)

    # out paths
    transcript_overlap_out_path = os.path.join(overlap_resource_dir, f'{'.'.join(experiment_list)}.transcripts.{''.join(target_snp)}.mirna.o{offset}.parquet')
    transcript_expr_out_path = os.path.join(expr_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.mirna.o{offset}.parquet')

    #---Step 2: Calculate Sequence Context if Not Done Previously ----
    print(f'\nGetting Sequence Context for {''.join(target_snp)} Off-Targets (Offset = {offset})...')
    offset_col = f'o{offset}_seq_context'

    transcript_fa = pysam.FastaFile(transcript_fa_path)
    transcript_overlap_df[offset_col] = transcript_overlap_df.progress_apply(
        get_norm_seq_context, pos_col='transcript_pos', offset=offset, transcript_fa=transcript_fa,
        axis=1
    )

    # garbage collection
    del transcript_fa
    gc.collect()

    print(f'\nSequence Context Length Value Counts (Offset={offset}):')
    print(tabulate(pd.DataFrame(transcript_overlap_df[offset_col].apply(len).value_counts()), headers='keys', tablefmt='fancy_grid'))
    print(f'\nSequence Context Central Base Counts (Offset={offset}):')
    print(tabulate(pd.DataFrame(transcript_overlap_df[offset_col].apply(lambda x: x[offset]).value_counts()), headers='keys', tablefmt='fancy_grid'))

    print(f'\nGot Sequence Context for {''.join(target_snp)} Off-Targets (Offset = {offset}).\n')

    #---Step 3: Calculate Frequency of miRTarBase Target Sites in Off Target Transcripts ----
    print(f'Finding miRTarBase Target Site Matches in {''.join(target_snp)} Off-Targets (Offset = {offset})...\n')
    with gzip.open(tgt_site_patterns_path, 'rb') as f:
            mirtar_patterns = pickle.load(f)

    match_col = f'o{offset}_mirtarbase_matches'
    transcript_overlap_df[match_col] = transcript_overlap_df[f'o{offset}_seq_context'].parallel_apply(
        find_tgt_site_matches, mirtar_patterns=mirtar_patterns
    )

    total_matches = transcript_overlap_df[match_col].apply(len).sum()
    print(f'{total_matches} Matches in Off-Targets.\n')

    print(f'Found miRTarBase Target Site Matches in {''.join(target_snp)} Off-Targets (Offset = {offset})...\n')

    # garbage collection
    del mirtar_patterns
    gc.collect()

    #---Step 4: Calculate Frequency of miRBase Target Sites in Off Target Transcripts ----
    # Align seed motifs to seq_contexts
    print(f'Finding miRBase Seed Matches in {''.join(target_snp)} Off-Targets (Offset = {offset})...\n')
    with gzip.open(seed_patterns_path, 'rb') as f:
        seed_patterns = pickle.load(f)

    match_col = f'o{offset}_mirbase_matches'
    transcript_overlap_df[match_col] = transcript_overlap_df[f'o{offset}_seq_context'].parallel_apply(
        find_seed_matches, seed_patterns=seed_patterns
    )

    total_matches = transcript_overlap_df[match_col].apply(len).sum()
    print(f'{total_matches} Matches in Off-Targets.\n')

    print(f'Found miRBase Seed Matches in {''.join(target_snp)} Off-Targets (Offset = {offset}).\n')

    #---Step 5: Calculate Frequencies in Expressed Transcripts As Background Control----
    print(f'Align miRNA Databases to Expressed Transcripts...\n')
    # load expressed transcripts info and random control positions and merge to one df
    with gzip.open(tgt_site_patterns_path, 'rb') as f:
        mirtar_patterns = pickle.load(f)

    for seq_motif in random_seq_motifs.keys():
        random_col = f'random_{seq_motif}_pos'
        mirtar_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.mirtarbase.{seq_motif}.o{offset}.parquet')
        mirbase_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.mirbase.{seq_motif}.o{offset}.parquet')

        if (not os.path.isfile(mirtar_path)) or (not os.path.isfile(mirbase_path)):
            print(f'\n Getting Seq Contexts for Random {seq_motif}...\n')
            transcript_fa = pysam.FastaFile(transcript_fa_path)
            seq_context = transcript_expr_df.progress_apply(get_norm_seq_context, pos_col=random_col, transcript_fa=transcript_fa, offset=offset, axis=1)
            
            if not os.path.isfile(mirtar_path):
                print(f'\nAligning miRTarBase Target Sites to Random {seq_motif}...')
                mirtar_col = f'random_{seq_motif}_o{offset}_mirtarbase_matches'
                hits_col = seq_context.parallel_apply(
                    find_tgt_site_matches, mirtar_patterns=mirtar_patterns
                )
                pd.DataFrame({mirtar_col:hits_col}).to_parquet(mirtar_path, compression='gzip')

                # garbage collection
                del hits_col
                gc.collect()
            
            if (not os.path.isfile(mirbase_path)):
                print(f'\nAligning miRBase Target Sites to Random {seq_motif}...')
                mirbase_col = f'random_{seq_motif}_o{offset}_mirbase_matches'
                hits_col = seq_context.parallel_apply(
                    find_seed_matches, seed_patterns=seed_patterns
                )
                pd.DataFrame({mirbase_col:hits_col}).to_parquet(mirbase_path, compression='gzip')
                
                # garbage collection
                del hits_col
                gc.collect()

            # garbage collection
            del transcript_fa
            del seq_context
            gc.collect()

    # garbage collection
    del random_pos_df
    del mirtar_patterns
    del seed_patterns
    gc.collect()

    print(f'Aligned miRNA Databases to Expressed Transcripts.\n')

    #---Step 6: Guide Alignment----
    print('Aligning Editor Guides...')
    # get guides for each condition and build an info df
    condition_guides = sample_map_df.set_index(zip(
        sample_map_df['condition'],
        sample_map_df['guide'],
        sample_map_df['guide_seed'],
    )).index.unique()
    condition_guides = [list(g) for g in list(condition_guides)]

    guide_info_df = pd.DataFrame(
        condition_guides,
        columns=['condition', 'guide', 'seed']
    )
    guide_info_df['name'] = ''
    no_guide_idx = guide_info_df[(guide_info_df['guide'].apply(pd.isna)) | (guide_info_df['guide'] == '')].index
    guide_info_df = guide_info_df.drop(no_guide_idx)
    guide_info_df = guide_info_df.set_index('condition', drop=True)

    # add info to match mirbase processing
    guide_info_df['8-mer']         = guide_info_df['seed'].apply(lambda x: str(Seq.Seq(x[0:7]).reverse_complement()) + 'A')
    guide_info_df['7-mer-A1']      = guide_info_df['seed'].apply(lambda x: str('N' + Seq.Seq(x[0:6]).reverse_complement()) + 'A')
    guide_info_df['7-mer-m8']      = guide_info_df['seed'].apply(lambda x: str(Seq.Seq(x[0:7]).reverse_complement()) + 'N')
    guide_info_df['6-mer']         = guide_info_df['seed'].apply(lambda x: str('N' + Seq.Seq(x[0:6]).reverse_complement()) + 'N')

    # build pattern dict like mirbase
    guide_patterns = build_seed_patterns(guide_info_df)

    # align to off targets
    match_col_name = f'o{offset}_guides_matches'
    transcript_overlap_df[match_col_name] = transcript_overlap_df[f'o{offset}_seq_context'].progress_apply(find_seed_matches, seed_patterns=guide_patterns)
        
    print('Aligned Editor Guides.\n')

    # align to expressed transcripts
    for seq_motif in random_seq_motifs.keys():
        guide_align_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.guides.{seq_motif}.o{offset}.parquet')
        if not os.path.isfile(guide_align_path):
            random_col = f'random_{seq_motif}_pos'
            transcript_fa = pysam.FastaFile(transcript_fa_path)
            seq_context = transcript_expr_df.progress_apply(get_norm_seq_context, pos_col=random_col, transcript_fa=transcript_fa, offset=offset, axis=1)
            
            match_col_name = f'random_{seq_motif}_o{offset}_guides_matches'
            match_col = seq_context.progress_apply(find_seed_matches, seed_patterns=guide_patterns)
            pd.DataFrame({match_col_name:match_col}).to_parquet(guide_align_path, compression='gzip')

            # garbage collection
            del transcript_fa
            del seq_context
            del match_col
            gc.collect()

    #---Step 7: Select Random Samples for Controls with matched sample sizes----
    print(f'Selecting Random Expressed Transcripts for Control Populations...')
    rs_dict_update = False

    if os.path.isfile(random_sample_path):
        with gzip.open(random_sample_path, 'rb') as f:
            random_sample_dict = pickle.load(f)

        # Prune IDs no longer expressed in the current experiment (handles cross-experiment dict reuse)
        expr_idx_set = set(transcript_expr_df.index)
        for motif in list(random_sample_dict.keys()):
            for old_size in list(random_sample_dict[motif].keys()):
                stored_ids = random_sample_dict[motif][old_size]
                valid_ids = sorted(set(stored_ids).intersection(expr_idx_set))
                if len(valid_ids) < len(stored_ids):
                    rs_dict_update = True
                    del random_sample_dict[motif][old_size]
                    new_size = len(valid_ids)
                    if new_size > 0:
                        random_sample_dict[motif][new_size] = valid_ids
    else:
        random_sample_dict = dict()

    hit_cols = [[f'{condition}_hit', f'{condition}_expressed'] for condition in non_wt_conditions]
    sample_sizes = [len(transcript_overlap_df[transcript_overlap_df[hit_col].all(axis=1)]) for hit_col in hit_cols]

    for motif in random_seq_motifs.keys():
        if motif not in list(random_sample_dict.keys()):
            random_sample_dict[motif] = {}
        
        for sample_size in sample_sizes:
            if sample_size not in list(random_sample_dict[motif].keys()):
                rs_dict_update = True
                random_sample_dict[motif] = add_random_sample(
                    sample_size=sample_size,
                    random_sample_dict=random_sample_dict[motif],
                    transcript_expr_df=transcript_expr_df,
                    motif=motif
                )

    if rs_dict_update:
        with gzip.open(random_sample_path, 'wb') as f:
            pickle.dump(random_sample_dict, f)

    print(f'Selected Random Expressed Transcripts for Control Populations.\n')

    # garbage collection
    del random_sample_dict
    gc.collect()

    #---Step 8: Final Calculations----
    # count total hits per transcript in a new column for all datasets
    print('Counting Total Hits for Each Transcript...\n')
    # off targets
    print(f'Counting Hits in {''.join(target_snp)} Off-Targets (Offset = {offset})...')
    for db in ['mirtarbase', 'mirbase', f'guides']:
        match_col = f'o{offset}_{db}_matches'
        total_col = f'o{offset}_{db}_total'

        transcript_overlap_df[total_col] = transcript_overlap_df[match_col].progress_apply(len)

        if (db == 'mirbase') or (db == f'guides'):
            seed_types = ['8-mer', '7-mer-m8']
            total_col = f'o{offset}_{db}_full-seed_total'
            
            transcript_overlap_df[total_col] = transcript_overlap_df[match_col].progress_apply(count_seed_type_hits, seed_types=seed_types)

    print(f'Counted Hits in {''.join(target_snp)} Off-Targets (Offset = {offset}).\n')

    # expressed transcripts
    print(f'Counting Hits in Expressed Transcript Background (Offset = {offset})...')
    for db in ['mirtarbase', 'mirbase', f'guides']:
        for seq_motif in random_seq_motifs.keys():
            print(f'\nRandom {seq_motif} {db}...')

            match_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.{db}.{seq_motif}.o{offset}.parquet')
            total_col = f'random_{seq_motif}_o{offset}_{db}_total'
            match_col = f'random_{seq_motif}_o{offset}_{db}_matches'
            
            match_data = pd.read_parquet(match_path)
            transcript_expr_df[total_col] = match_data[match_col].progress_apply(len)


            if (db == 'mirbase') or (db == f'guides'):
                seed_types = ['8-mer', '7-mer-m8']
                total_col = f'random_{seq_motif}_o{offset}_{db}_full-seed_total'

                transcript_expr_df[total_col] = match_data[match_col].progress_apply(count_seed_type_hits, seed_types=seed_types)

            # garbage collection
            del match_data
            gc.collect()

    print(f'Counted Hits in Expressed Transcript Background (Offset = {offset}).\n')
    print('Counted Total Hits for Each Transcript.\n')

    # build final tables of total hit counts and hits/transcript
    print('Tabulating Counts and Frequencies...')
    print(f'\nTotal Hits...')
    with gzip.open(random_sample_path, 'rb') as f:
        random_sample_dict = pickle.load(f)

    count_df_idx = ['total_matches'] + [f'random_{m}_matches' for m in random_seq_motifs.keys()]
    freq_df_idx = ['matches_per_transcript'] + [f'random_{m}_matches_per_transcript' for m in random_seq_motifs.keys()]

    count_df_cols = defaultdict(lambda: defaultdict(list))
    freq_df_cols = defaultdict(lambda: defaultdict(list))
    seed_types = ['8-mer', '7-mer-m8']
    for condition in non_wt_conditions:
        hit_filter = transcript_overlap_df[f'{condition}_hit'] & transcript_overlap_df[f'{condition}_expressed'] 
        sample_size = len(transcript_overlap_df[hit_filter])
        
        for db in ['mirtarbase', 'mirbase', 'mirbase_full-seed']:
            match_total_col = f'o{offset}_{db}_total'
            total_matches = transcript_overlap_df.loc[hit_filter, match_total_col].sum()
            if sample_size != 0:
                match_freq = total_matches/sample_size 
            else:
                match_freq = 0

            count_df_cols[db][condition].append(total_matches)
            freq_df_cols[db][condition].append(match_freq)
        
            for seq_motif in random_seq_motifs.keys():
                random_idx = random_sample_dict[seq_motif][sample_size]
                match_total_col = f'random_{seq_motif}_o{offset}_{db}_total'
                total_matches = transcript_expr_df.loc[random_idx, match_total_col].sum()
                if sample_size != 0:
                    match_freq = total_matches/sample_size 
                else:
                    match_freq = 0

                count_df_cols[db][condition].append(total_matches)
                freq_df_cols[db][condition].append(match_freq)

        match_col = f'o{offset}_guides_matches'
        total_matches = transcript_overlap_df.loc[hit_filter, match_col].apply(count_seed_type_hits, mirna_ids=[condition]).sum()
        if sample_size != 0:
            match_freq = total_matches/sample_size 
        else:
            match_freq = 0
        count_df_cols['guides'][condition].append(total_matches)
        freq_df_cols['guides'][condition].append(match_freq)

        total_matches = transcript_overlap_df.loc[hit_filter, match_col].apply(count_seed_type_hits, mirna_ids=[condition], seed_types=seed_types).sum()
        if sample_size != 0:
            match_freq = total_matches/sample_size 
        else:
            match_freq = 0
        count_df_cols['guides_full-seed'][condition].append(total_matches)
        freq_df_cols['guides_full-seed'][condition].append(match_freq)

    for seq_motif in random_seq_motifs.keys():
        match_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.guides.{seq_motif}.o{offset}.parquet')
        match_col = f'random_{seq_motif}_o{offset}_guides_matches'
        match_data = pd.read_parquet(match_path)[match_col]
        
        for condition in non_wt_conditions:
            hit_filter = transcript_overlap_df[f'{condition}_hit'] & transcript_overlap_df[f'{condition}_expressed'] 
            sample_size = len(transcript_overlap_df[hit_filter])
            random_idx = random_sample_dict[seq_motif][sample_size]

            total_matches = match_data[random_idx].apply(count_seed_type_hits, mirna_ids=[condition]).sum()
            if sample_size != 0:
                match_freq = total_matches/sample_size 
            else:
                match_freq = 0
            count_df_cols['guides'][condition].append(total_matches)
            freq_df_cols['guides'][condition].append(match_freq)

            total_matches = match_data[random_idx].apply(count_seed_type_hits, mirna_ids=[condition], seed_types=seed_types).sum()
            if sample_size != 0:
                match_freq = total_matches/sample_size 
            else:
                match_freq = 0
            count_df_cols['guides_full-seed'][condition].append(total_matches)
            freq_df_cols['guides_full-seed'][condition].append(match_freq)
        
        # garbage collection
        del match_data
        gc.collect()

    count_df_cols = dict(count_df_cols)
    freq_df_cols = dict(freq_df_cols)

    xlsx_out_path = os.path.join(agg_res_dir, f'{'.'.join(experiment_list)}.mirna-align.{target_snp[0]}{target_snp[1]}.o{offset}.xlsx')
    with pd.ExcelWriter(xlsx_out_path, engine='openpyxl') as writer:
        for db in ['mirtarbase', 'mirbase', 'mirbase_full-seed', f'guides', f'guides_full-seed']:
            count_df = pd.DataFrame(count_df_cols[db], index=count_df_idx).T
            freq_df = pd.DataFrame(freq_df_cols[db], index=freq_df_idx).T

            count_df.to_excel(writer, sheet_name=f'{db}-count'[:31], index=True)
            freq_df.to_excel(writer, sheet_name=f'{db}-freq'[:31], index=True)

            print(f'\n{db}:')
            print(tabulate(count_df, headers='keys', tablefmt='fancy_grid'))
            print(tabulate(freq_df, headers='keys', tablefmt='fancy_grid'))

    print(f'Done Total Hits.\n')

    print(f'Positional...')
    count_df_cols = defaultdict(lambda: defaultdict(dict))
    freq_df_cols = defaultdict(lambda: defaultdict(dict))
    for condition in tqdm.tqdm(non_wt_conditions, total=len(non_wt_conditions), desc=f'Counting Positional Matches for Off-Targets:'):
        hit_filter = transcript_overlap_df[f'{condition}_hit'] & transcript_overlap_df[f'{condition}_expressed'] 
        
        for db in ['mirtarbase', 'mirbase']:
            match_col = f'o{offset}_{db}_matches'
            pos_counts = sum_positional_matches(
                match_series=transcript_overlap_df.loc[hit_filter, match_col],
                offset=offset,
                mode=db,
            )
            count_df_cols[condition][db]    = pos_counts
            freq_df_cols[condition][db]     = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

            if db == 'mirbase':
                pos_counts = sum_positional_matches(
                    match_series=transcript_overlap_df.loc[hit_filter, match_col],
                    offset=offset,
                    mode=db,
                    seed_types=seed_types
                )
                count_df_cols[condition][f'{db}_full-seed'] = pos_counts
                freq_df_cols[condition][f'{db}_full-seed']     = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

        match_col = f'o{offset}_guides_matches'
        pos_counts = sum_positional_matches(
                match_series=transcript_overlap_df.loc[hit_filter, match_col],
                offset=offset,
                mode='mirbase',
                mirna_ids=[condition]
            )
        count_df_cols[condition][f'guides']    = pos_counts
        freq_df_cols[condition][f'guides']     = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

        pos_counts = sum_positional_matches(
                match_series=transcript_overlap_df.loc[hit_filter, match_col],
                offset=offset,
                mode='mirbase',
                seed_types=seed_types,
                mirna_ids=[condition]
            )
        count_df_cols[condition][f'guides_full-seed'] = pos_counts
        freq_df_cols[condition][f'guides_full-seed']  = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

    for seq_motif in tqdm.tqdm(random_seq_motifs.keys(), total=len(random_seq_motifs), desc=f'Counting Positional Alignment for Random Background:'):
        for db in ['mirtarbase', 'mirbase']:
            match_path = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.{db}.{seq_motif}.o{offset}.parquet')
            match_col = f'random_{seq_motif}_o{offset}_{db}_matches'
            match_data = pd.read_parquet(match_path)[match_col]
            
            for condition in non_wt_conditions:
                hit_filter = transcript_overlap_df[f'{condition}_hit'] & transcript_overlap_df[f'{condition}_expressed'] 
                sample_size = len(transcript_overlap_df[hit_filter])
                random_idx = random_sample_dict[seq_motif][sample_size]

                pos_counts = sum_positional_matches(
                    match_series=match_data[random_idx],
                    offset=offset,
                    mode=db
                )
                count_df_cols[condition][f'random_{seq_motif}_{db}'] = pos_counts
                freq_df_cols[condition][f'random_{seq_motif}_{db}']  = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

                if db == 'mirbase':
                    pos_counts = sum_positional_matches(
                        match_series=match_data[random_idx],
                        offset=offset,
                        mode=db,
                        seed_types=seed_types
                    )
                    count_df_cols[condition][f'random_{seq_motif}_{db}_full-seed'] = pos_counts
                    freq_df_cols[condition][f'random_{seq_motif}_{db}_full-seed']  = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts
            
            # garbage collection
            del match_data
            gc.collect()

        db = f'guides'
        match_path  = os.path.join(expr_mirna_resource_dir, f'{'.'.join(experiment_list)}.transcripts.expressed.{db}.{seq_motif}.o{offset}.parquet')
        match_col   = f'random_{seq_motif}_o{offset}_{db}_matches'
        match_data  = pd.read_parquet(match_path)[match_col]
        
        for condition in non_wt_conditions:
            hit_filter = transcript_overlap_df[f'{condition}_hit'] & transcript_overlap_df[f'{condition}_expressed'] 
            sample_size = len(transcript_overlap_df[hit_filter])
            random_idx = random_sample_dict[seq_motif][sample_size]

            pos_counts = sum_positional_matches(
                match_series=match_data[random_idx],
                offset=offset,
                mode='mirbase',
                mirna_ids=[condition]
            )
            count_df_cols[condition][f'random_{seq_motif}_{db}'] = pos_counts
            freq_df_cols[condition][f'random_{seq_motif}_{db}']  = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

            pos_counts = sum_positional_matches(
                match_series=match_data[random_idx],
                offset=offset,
                mode='mirbase',
                mirna_ids=[condition],
                seed_types=seed_types
            )
            count_df_cols[condition][f'random_{seq_motif}_{db}_full-seed'] = pos_counts
            freq_df_cols[condition][f'random_{seq_motif}_{db}_full-seed']  = pos_counts/pos_counts.sum() if pos_counts.sum() > 0 else pos_counts

        # garbage collection
        del match_data
        gc.collect()

    xlsx_out_path = os.path.join(agg_res_dir, f'{'.'.join(experiment_list)}.mirna-align.{target_snp[0]}{target_snp[1]}.o{offset}.pos.xlsx')
    with pd.ExcelWriter(xlsx_out_path, engine='openpyxl') as writer:
        for condition in non_wt_conditions:
            count_df = pd.DataFrame(count_df_cols[condition], index=list(range(-offset, offset+1)))
            freq_df = pd.DataFrame(freq_df_cols[condition], index=list(range(-offset, offset+1)))

            count_df.to_excel(writer, sheet_name=f'{condition}-count'[:31], index=True)
            freq_df.to_excel(writer, sheet_name=f'{condition}-freq'[:31], index=True)

            print(f'\n{condition}:')

    transcript_overlap_df.to_parquet(transcript_overlap_out_path)

    expr_out_cols = list(set(transcript_expr_df.columns).difference(non_out_transcript_expr_cols))
    transcript_expr_df[expr_out_cols].to_parquet(transcript_expr_out_path)

    print(f'Done Positional.')
    print('Tabulated Counts and Frequencies.\n')