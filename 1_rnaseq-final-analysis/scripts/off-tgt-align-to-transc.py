#####################
# import statements #
#####################

import os
import pandas as pd
import sys
import pysam
import pickle
import gzip
import tqdm
from tabulate import tabulate
from Bio import Seq
import re
from pandarallel import pandarallel
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

def var_row_to_trans_df_rows(var_row, transcript_contig_dict, sample_map_df):

    chrom = var_row['chrom']
    pos = int(var_row['pos'])
    strand = var_row['strand']

    conditions = sorted(list(sample_map_df['condition'].unique()))
    non_wt_mask = ~sample_map_df['condition'].isin(sample_map_df['wt_condition'])
    non_wt_conditions = list(sample_map_df.loc[non_wt_mask, 'condition'].unique())
    non_wt_conditions.sort()
    
    transcript_df_rows = []
    for hit_dict in var_row['gtf_hits']:
        transcript_id = hit_dict['transcript']
        exon_id = hit_dict['exon']
        
        try:
            transcript_contig = transcript_contig_dict[transcript_id]
            
            transcript_df_row = [var_row.name, chrom, pos, strand, var_row['ref'], var_row['alt']]
            transcript_df_row += [var_row[f'{condition}_mean_pct_ref'] for condition in conditions]
            transcript_df_row += [var_row[f'{condition}_mean_pct_snp'] for condition in conditions]
            transcript_df_row += [var_row[f'{condition}_hit'] for condition in non_wt_conditions]
            transcript_df_row += [hit_dict['gene'][1], hit_dict['gene'][0], transcript_id, exon_id, transcript_contig]

            transcript_df_rows.append(transcript_df_row)
        except KeyError:
            pass

    return transcript_df_rows

def gtf_rec_from_line(gtf_line):
    gtf_rec = {}
    gtf_list = gtf_line.strip().split('\t')
    
    gtf_rec['chrom'] = gtf_list[0]
    gtf_rec['source'] = gtf_list[1]
    gtf_rec['feat_type'] = gtf_list[2]
    gtf_rec['start'] = int(gtf_list[3])
    gtf_rec['end'] = int(gtf_list[4])
    gtf_rec['strand'] = gtf_list[6]
    gtf_rec['info'] = {}

    gtf_info =  gtf_list[8].strip().split('; ')
    for item in gtf_info:
        key = item.split(' ')[0].strip()
        val = item.split(' ')[1].strip()
        val = val.replace(';','')

        if '"' in val:
            val = val.split('"')[1].strip()
        else:
            val = int(val)

        gtf_rec['info'][key] = val

    return gtf_rec    

def get_exon_info(transcript_row, gtfgz):

    # gtfgz = pysam.TabixFile(gtf_path)

    exon_id = transcript_row['exon_id']
    chrom = transcript_row['chrom']
    pos = transcript_row['pos']

    for gtf_line in gtfgz.fetch(reference=chrom, start=pos-1, end=pos):
        if gtf_line.strip().split('\t')[2].strip() == 'exon':
            gtf_rec = gtf_rec_from_line(gtf_line=gtf_line)
            if gtf_rec['info']['exon_id'] == exon_id:
                return pd.Series([gtf_rec['info']['exon_number'], gtf_rec['start'], gtf_rec['end']])
            
    # # garbage collection
    # del gtfgz
    # gc.collect()
            
    return pd.Series([-1, -1, -1])

def get_transcript_pos(transcript_row, genome_fa, transcript_fa):
    # genome_fa = pysam.Fastafile(genome_fa_path)
    # transcript_fa = pysam.Fastafile(transcript_fa_path)
    
    chrom =transcript_row['chrom']
    pos = transcript_row['pos']
    exon_start = transcript_row['exon_start']
    exon_end = transcript_row['exon_end']
    transcript_contig = transcript_row['transcript_contig']
    strand = transcript_row['strand']

    exon_seq = genome_fa.fetch(reference=chrom, start=exon_start-1, end=exon_end)
    exon_seq = str(Seq.Seq(exon_seq).reverse_complement()) if strand == '-' else exon_seq 
    transcript_seq = transcript_fa.fetch(reference=transcript_contig)

    exon_matches = list(re.finditer(exon_seq, transcript_seq))
    if len(exon_matches) < 1:
        transcript_pos = -1
    elif len(exon_matches) > 1:
        transcript_pos = -2
    else:
        exon_pos = pos - (exon_start - 1)
        exon_pos = exon_end - (pos - 1) if strand == '-' else exon_pos

        for match in exon_matches:
            transcript_exon_start = match.start() + 1
        
        transcript_pos = exon_pos + (transcript_exon_start - 1)

    # # garbage collection
    # del genome_fa
    # del transcript_fa
    # gc.collect()

    return transcript_pos

###############
# main script #
###############
if __name__ == '__main__':
    tqdm.tqdm.pandas()
    pandarallel.initialize(progress_bar=True, verbose=2, nb_workers=8)
        
    #--- Step 1: define variables ---
    # from args
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    var_dir = sys.argv[3]
    ref_dir = sys.argv[4]
    # target SNP as a tuple (ref, alt)
    target_snp = (sys.argv[5][0],sys.argv[5][1])

    # directories
    hit_dir = os.path.join(var_dir, 'hit-id')
    align_dir = os.path.join(var_dir, 'transcript-align')
    out_dir = os.path.join(align_dir, f'{target_snp[0]}{target_snp[1]}')
    expression_dir = os.path.join(proj_dir, '5_expression')
    salmon_dir = os.path.join(expression_dir, 'salmon')
    salmon_analysis_dir = os.path.join(salmon_dir, 'analysis')

    # ref_fles
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')
    genome_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.genome.fa.gz')
    genome_gtf_path = os.path.join(ref_dir, f'{ref_genome_name}.sorted.gtf.gz')

    sample_map_df, non_wt_conditions = load_sample_map(proj_dir=proj_dir)

    overlap_path = os.path.join(hit_dir, f'{experiment}.overlap.{''.join(target_snp)}.parquet')
    overlap_df = pd.read_parquet(overlap_path)

    transcript_raw_counts_path= os.path.join(salmon_analysis_dir, 'salmon-transcript-raw-counts.tsv.gz')
    transcript_raw_counts_df = pd.read_csv(transcript_raw_counts_path, sep='\t', compression='gzip')
    transcript_raw_counts_df = transcript_raw_counts_df.set_index('transcript_id', drop=True)

    #---Step 2: Create Transcript df ---
    # build dict of transcript contigs
    transcript_contig_path = os.path.join(align_dir, 'transcript-id-to-contig.pkl.gz')
    with gzip.open(transcript_contig_path, 'rb') as f:
        transcript_contig_dict = pickle.load(f)
    print(f'\nLoaded Previously Built Transcript Contig Lookup for {len(transcript_contig_dict)} Transcripts.\n')

    print(f'Building Initial Transcript df from Overlap Data...')

    transcript_df_cols = ['var_df_idx', 'chrom', 'pos',  'strand', 'ref', 'alt']
    transcript_df_cols += [f'{c}_mean_pct_ref' for c in sample_map_df['condition'].unique()] 
    transcript_df_cols += [f'{c}_mean_pct_snp' for c in sample_map_df['condition'].unique()] 
    transcript_df_cols += [f'{c}_hit' for c in non_wt_conditions]
    transcript_df_cols += ['gene_name', 'gene_id', 'transcript_id', 'exon_id', 'transcript_contig']

    transcript_df = overlap_df.parallel_apply(var_row_to_trans_df_rows, transcript_contig_dict=transcript_contig_dict, sample_map_df=sample_map_df, axis=1)
    transcript_df = transcript_df.sum()
    transcript_df = pd.DataFrame(transcript_df, columns=transcript_df_cols)

    print(f'{len(transcript_df)} Rows in Initial Transcript df')
    print(f'Built Initial Transcript df.\n')

    print(f'Filter Transcript df for Transcripts with DP >= 20...')
    expressed_cols = {}
    for condition in non_wt_conditions:
        condition_mask = sample_map_df['condition'] == condition
        condition_samples = sample_map_df.loc[condition_mask, 'sample']
        
        print(f'Determining {condition} Expressed Transcripts:')
        transcript_expressed_mask = transcript_df['transcript_id'].parallel_apply(
            lambda x:
                all(transcript_raw_counts_df.loc[x, [f'{sample}_reads' for sample in condition_samples]] >=20) 
                if x in transcript_raw_counts_df.index else False
        )
        print('')

        expressed_cols[f'{condition}_expressed'] = transcript_expressed_mask

    transcript_df = pd.concat([transcript_df, pd.DataFrame(expressed_cols)], axis=1)

    none_expressed_mask = ~(transcript_df[[col for col in expressed_cols.keys()]].any(axis=1))
    print(f'{len(transcript_df[none_expressed_mask])} Non-Expressed Transcripts in Initial Transcript df.')

    transcript_df = transcript_df.drop(transcript_df[none_expressed_mask].index)
    print(f'Filtered Transcript df for Transcripts with DP >= 20.\n')

    print(f'Get Exon Info...')
    gtfgz = pysam.TabixFile(genome_gtf_path)
    exon_info_cols = ['exon_num', 'exon_start', 'exon_end']
    transcript_df[exon_info_cols] = transcript_df.progress_apply(
            get_exon_info, gtfgz=gtfgz,
            axis=1
        )
    no_exon_found_mask = transcript_df['exon_num'] == -1
    print(f'Dropping {len(transcript_df[no_exon_found_mask])} Rows With No Exon Info...')

    transcript_df = transcript_df.drop(transcript_df[no_exon_found_mask].index)
    print(f'Got Exon Info.\n')

    print(f'Get Transcript pos...')
    genome_fa = pysam.FastaFile(genome_fa_path)
    transcript_fa = pysam.FastaFile(transcript_fa_path)
    transcript_df['transcript_pos'] = transcript_df.progress_apply(
            get_transcript_pos, genome_fa=genome_fa, transcript_fa=transcript_fa,
            axis=1
        )

    no_match_mask = transcript_df['transcript_pos'] == -1
    multi_match_mask = transcript_df['transcript_pos'] == -2
    print(f'{len(transcript_df[no_match_mask])} Rows where Exon was not Found in Transcript.')
    print(f'{len(transcript_df[multi_match_mask])} Rows where Exon was Found in Transcript Multiple Times.')

    # I've never seen this print >10 so I'm not writing code to handle these cases
    transcript_df = transcript_df.drop(transcript_df[no_match_mask | multi_match_mask].index)
    print(f'Got Transcript pos.\n')

    print(f'Sanity Checks...')
    print(f'transcript_df[ref] Value Counts:')
    print(transcript_df['ref'].value_counts())
    print(f'\ntranscript_df[alt] Value Counts:')
    print(transcript_df['alt'].value_counts())

    print(f'\nBase Identities Using Genomic Position:')
    print(transcript_df.progress_apply(lambda x: genome_fa.fetch(reference=x['chrom'], start=x['pos']-1,end=x['pos']) if x['strand'] == '+' else
                                str(Seq.Seq(genome_fa.fetch(reference=x['chrom'], start=x['pos']-1,end=x['pos'])).reverse_complement()), axis=1).value_counts())

    print(f'\nBase Identities Using Transcript Position:')
    print(transcript_df.progress_apply(lambda x: transcript_fa.fetch(reference=x['transcript_contig'], start=x['transcript_pos']-1, end=x['transcript_pos']), axis=1).value_counts())

    print(f'Sanity Checked for {condition}.\n')

    print(tabulate(transcript_df.head(), headers='keys', tablefmt='fancy_grid'))
    print('')

    #---Step 3: Write Outputs ---
    out_pq_path = os.path.join(align_dir, f'{experiment}.transcripts.{''.join(target_snp)}.parquet')
    transcript_df.to_parquet(out_pq_path, compression='gzip')
    for condition in non_wt_conditions:
        hit_col = f'{condition}_hit'
        expressed_col = f'{condition}_expressed'

        condition_mask = transcript_df[hit_col] & transcript_df[expressed_col]

        out_path = os.path.join(out_dir, f'{experiment}.{condition}.transcripts.{target_snp[0]}{target_snp[1]}.parquet')
        transcript_df[condition_mask].to_parquet(out_path, compression='gzip')