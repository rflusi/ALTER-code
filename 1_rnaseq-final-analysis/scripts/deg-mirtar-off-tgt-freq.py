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

###############
# main script #
###############
if __name__ == '__main__':
    #--- Step 1: define variables ---
    # from args
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    var_dir = sys.argv[3]
    ref_dir = sys.argv[4]
    ref_db_dir = sys.argv[5]
    deseq_fc_cutoff = int(sys.argv[6])
    deseq_padj_cutoff = float(sys.argv[7])
    snps = sys.argv[8].split(',')

    print(snps)

    # other directory paths
    miRTarBase_dir = os.path.join(ref_db_dir, 'miRtarBase')
    miRBase_dir = os.path.join(ref_db_dir, 'miRBase')
    biomart_dir = os.path.join(ref_db_dir, 'biomart')

    align_dir = os.path.join(var_dir, 'transcript-align')
    # align_snp_dir = os.path.join(align_dir, f'{target_snp[0]}{target_snp[1]}')
    hit_dir = os.path.join(var_dir, 'hit-id')

    expression_dir = os.path.join(proj_dir, '5_expression')

    salmon_dir = os.path.join(expression_dir, 'salmon')
    salmon_analysis_dir = os.path.join(salmon_dir, 'analysis')
    transc_expr_mirna_align_dir = os.path.join(salmon_dir, 'mirna-alignment')

    deseq_dir = os.path.join(expression_dir, 'deseq2')
    deseq_results_dir = os.path.join(deseq_dir, '3_results')
    result_tables_dir = os.path.join(deseq_results_dir, '2_result-tables')
    deg_tables_dir = os.path.join(result_tables_dir, '2_deg')

    # ref file paths
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')
    genome_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.genome.fa.gz')
    genome_gtf_path = os.path.join(ref_dir, f'{ref_genome_name}.sorted.gtf.gz')
    gene_map_path = os.path.join(ref_dir, f'{ref_genome_name}.gene-map.parquet')

    transcript_contig_path = os.path.join(align_dir, 'transcript-id-to-contig.pkl.gz')
    mirtar_gene_path = os.path.join(miRTarBase_dir, 'miRTarBase-by-gene.parquet')
    mirtar_tgt_sites_pq_path = os.path.join(miRTarBase_dir, 'miRTarBase.tgt-sites.hsa.parquet')
    tgt_site_patterns_path = os.path.join(miRTarBase_dir, 'miRTarBase.tgt-sites.patterns.pkl.gz')
    mirbase_out_path = os.path.join(miRBase_dir, 'mirbase.parquet')
    seed_patterns_path = os.path.join(miRBase_dir, 'miRBase.seed-patterns.pkl.gz')

    # sample map
    sample_map_df, non_wt_conditions = load_sample_map(proj_dir=proj_dir)

    # analysis file paths and information
    raw_gene_path = os.path.join(salmon_analysis_dir, 'salmon-gene-raw-counts.tsv.gz')

    #---Step 2: Import Gene Level Raw Expression Data----
    print (f'Loading Gene Expression Data...')
    gene_raw_path = os.path.join(salmon_analysis_dir, 'salmon-gene-raw-counts.tsv.gz')
    gene_raw_df = pd.read_csv(gene_raw_path, sep='\t', compression='infer', index_col='Unnamed: 0')
    print (f'Loaded Gene Expression Data.\n')

    print (f'Determining Expressed Genes (DP >= 20)...\n')
    raw_count_mask_dict = {}
    expressed_mask_dict = {}
    for condition in non_wt_conditions:
        sample_mask = sample_map_df['condition'] == condition
        samples = sample_map_df.loc[sample_mask, 'sample']
        sample_cols = [f'{sample}_reads' for sample in samples]
        
        raw_mask = (gene_raw_df[sample_cols] >= 0).any(axis=1)
        raw_count_mask_dict[condition] = raw_mask
        expressed_mask = (gene_raw_df[sample_cols] >= 10).all(axis=1)
        expressed_mask_dict[condition] = expressed_mask

        print(f'{len(gene_raw_df[expressed_mask])}\tExpressed Genes for {condition}')

    biomart_gene_path = os.path.join(biomart_dir, 'biomart_gene.parquet')

    entrez_lookup = pd.read_parquet(biomart_gene_path)['gene_id_entrez']
    entrez_lookup = entrez_lookup.apply(lambda x: [int(val) for val in x.strip().split(',') if val != ''])

    ensembl_lookup = pd.DataFrame(entrez_lookup.explode(ignore_index=False)).reset_index(drop=False)
    ensembl_lookup = ensembl_lookup.groupby('gene_id_entrez')['gene_id'].apply(lambda x: [val for val in x.unique()])

    print (f'\nDetermined Expressed Genes (DP >= 10).\n')

    #--- Step 3: Get Gene Subsets ---
    print(f'Collecting Gene Subsets...')

    gene_map_df = pd.read_parquet(gene_map_path)
    mirtar_df = pd.read_parquet(mirtar_gene_path)
    subsets = [
            'All Genes',
            'Raw Count Genes',
            'Expressed Genes',
            'Assayed Genes',
            'DEGs',
            'Up-regulated',
            'Down-regulated'
        ]

    gene_subsets = {}
    for condition in non_wt_conditions: 
        gene_subsets[condition] = {}
        gene_subsets[condition]['ensembl'] = {}
        gene_subsets[condition]['entrez'] = {}


        drop_na_path    = os.path.join(deg_tables_dir, f'{condition}.2_de-results.drop-na.tsv.gz')
        drop_na_df      = pd.read_csv(drop_na_path, sep='\t', compression='infer')
        
        deg_path        = os.path.join(deg_tables_dir, f'{condition}.3_de-results.padj{deseq_padj_cutoff}-lfc{deseq_fc_cutoff}.tsv')
        deg_df          = pd.read_csv(deg_path, sep= '\t', compression='infer')
        upreg_mask      = deg_df['log2FoldChange'] > 0
        downreg_mask    = deg_df['log2FoldChange'] < 0

        gene_subsets[condition]['ensembl']['All Genes']        = set(gene_map_df.index)
        gene_subsets[condition]['ensembl']['Raw Count Genes']  = set(gene_raw_df.loc[raw_count_mask_dict[condition], 'gene_id'])
        gene_subsets[condition]['ensembl']['Expressed Genes']  = set(gene_raw_df.loc[expressed_mask_dict[condition], 'gene_id'])
        gene_subsets[condition]['ensembl']['Assayed Genes']    = set(drop_na_df['gene_id'])
        gene_subsets[condition]['ensembl']['DEGs']             = set(deg_df['gene_id'])
        gene_subsets[condition]['ensembl']['Up-regulated']     = set(deg_df.loc[upreg_mask, 'gene_id'])
        gene_subsets[condition]['ensembl']['Down-regulated']   = set(deg_df.loc[downreg_mask, 'gene_id'])

        for subset_name, gene_subset in gene_subsets[condition]['ensembl'].items():
            dever_subset = set([id.split('.')[0] for id in gene_subset])
            
            entrez_lookup_idx = list(dever_subset.intersection(set(entrez_lookup.index)))
            entrez_subset = set(entrez_lookup[entrez_lookup_idx].explode())

            gene_subsets[condition]['entrez'][subset_name] = entrez_subset

    print(f'Collected Gene Subsets.\n')

    #--- Step 4: calculate frequency of miRTarBase targets in DEGs ---
    print(f'Calculating Frequency of miRTarBase Targets in Gene Subsets...')

    intersection_dict = {result_type:{} for result_type in ['mirtar', 'off_tgts']}
    intersection_dict['mirtar'] = {}
    intersection_dict['mirtar']['ensembl'] = {}
    intersection_dict['mirtar']['entrez'] = {}

    mirtar_out_dict = {c:{} for c in non_wt_conditions}
    mirtar_summary_rows = []
    for condition in non_wt_conditions:
        intersection_dict['mirtar']['ensembl'][condition] = {}
        intersection_dict['mirtar']['entrez'][condition] = {}
        
        out_dict = {col:[] for col in ['subset', 'total_count', 'mirna_count', 'pct_mirna']}
        out_dict['subset'] = subsets
        
        mirtar_summary_row = [condition]
        for subset_name in subsets:
            entrez_set = gene_subsets[condition]['entrez'][subset_name] #convert to set of entrez ids
            mirtar_set = entrez_set.intersection(set(mirtar_df.index)) # set intersection with mirtarbase, entrez ids

            total_count = len(entrez_set)
            mirtar_count = len(mirtar_set)

            out_dict['total_count'].append(total_count)
            out_dict['mirna_count'].append(mirtar_count)
            out_dict['pct_mirna'].append((mirtar_count/total_count)*100)

            mirtar_summary_row.append((mirtar_count/total_count)*100)

            mirtar_set_ensembl = set(ensembl_lookup[list(mirtar_set)].explode()) # mirtar set converted to de-versioned ensembl
            intersection_dict['mirtar']['ensembl'][condition][subset_name] = mirtar_set_ensembl
            intersection_dict['mirtar']['entrez'][condition][subset_name]  = mirtar_set

        out_df = pd.DataFrame(out_dict)
        mirtar_out_dict[condition] = out_df

        mirtar_summary_rows.append(mirtar_summary_row)
        
        print(f'\n{condition}')
        print(tabulate(out_df, headers='keys', tablefmt='fancy_grid'))


    mirtar_summary_df = pd.DataFrame(mirtar_summary_rows, columns=(['condition'] + subsets))
    print(f'\nmiRTarBase Frequency Summary:')
    print(tabulate(mirtar_summary_df, headers='keys', tablefmt='fancy_grid'))
    print(f'\nCalculated Frequency of miRTarBase Targets in Gene Subsets.\n')

    #--- Step 5: calculate frequency of off-target edits in DEGs ---
    print(f'Calculating Frequency of Off-Target Edits in Gene Subsets...\n')

    # prep off tgt sets and store in gene_subsets
    for condition in non_wt_conditions:
        # add keys to gene subsets for off tgt sets
        gene_subsets[condition]['ensembl']['off_tgt'] = {s:{} for s in snps}
        gene_subsets[condition]['entrez']['off_tgt'] = {s:{} for s in snps}
    for snp in snps:
        # load overlap df
        overlap_path = os.path.join(hit_dir, f'{experiment}.overlap.{snp}.parquet')
        overlap_df = pd.read_parquet(overlap_path)

        for condition in non_wt_conditions:
            # pull hits for condition out of overlap_df
            hit_mask = overlap_df[f'{condition}_hit']
            all_off_tgt_ensembl = overlap_df.loc[hit_mask, 'gene_id'].apply(lambda x: x.split(',')) # check for entries with multiple ids
            all_off_tgt_ensembl_dever = all_off_tgt_ensembl.apply(lambda x: [id.split('.')[0] for id in x]) # de-version
            
            all_off_tgt_ensembl = set(all_off_tgt_ensembl.explode()) # convert to set of unique ids
            gene_subsets[condition]['ensembl']['off_tgt'][snp] = all_off_tgt_ensembl # add key for off tgt sets of snp
            
            all_off_tgt_ensembl_dever = set(all_off_tgt_ensembl_dever.explode()) # convert to set of unique de-versioned ensembl ids
            all_off_tgt_ensembl_dever = all_off_tgt_ensembl_dever.intersection(set(entrez_lookup.index)) # restrict to those with matching entrez id
            all_off_tgt_entrez = set(entrez_lookup[list(all_off_tgt_ensembl_dever)].explode()) # translate to set of entrez ids
            gene_subsets[condition]['entrez']['off_tgt'][snp] = all_off_tgt_entrez # add key for off tgt sets of snp

    intersection_dict['off-tgt'] = {} # add key to intersection_dict for off tgt sets

    off_tgt_out_dict = {s:{} for s in snps} # dict for future out_df cols
    off_tgt_summary_dict = {}
    for snp in snps:
        # add keys to intersection_dict for off tgt sets of snp
        intersection_dict['off-tgt'][snp] = {}
        intersection_dict['off-tgt'][snp]['ensembl'] = {}
        intersection_dict['off-tgt'][snp]['entrez'] = {}

        off_tgt_summary_rows = []
        for condition in non_wt_conditions:
            intersection_dict['off-tgt'][snp]['ensembl'][condition] = {}
            intersection_dict['off-tgt'][snp]['entrez'][condition] = {}

            # prep out_dict for condition
            out_dict = {col:[] for col in ['subset', 'total_count', 'off-tgt_count', 'pct_off-tgt']}
            out_dict['subset'] = subsets
            
            off_tgt_summary_row = [condition]
            for subset_name in subsets:
                ensembl_set = gene_subsets[condition]['ensembl'][subset_name] # convert to set of entrez ids
                entrez_set = gene_subsets[condition]['entrez'][subset_name] # convert to set of entrez ids

                off_tgt_set_ensembl = ensembl_set.intersection(set(gene_subsets[condition]['ensembl']['off_tgt'][snp])) # set intersection with off-tgts, ensembl ids
                off_tgt_set_entrez = ensembl_set.intersection(set(gene_subsets[condition]['entrez']['off_tgt'][snp])) # set intersection with off-tgts, entrez ids

                total_count = len(ensembl_set)
                off_tgt_count = len(off_tgt_set_ensembl)

                out_dict['total_count'].append(total_count)
                out_dict['off-tgt_count'].append(off_tgt_count)
                out_dict['pct_off-tgt'].append((off_tgt_count/total_count)*100)

                off_tgt_summary_row.append((off_tgt_count/total_count)*100)

                intersection_dict['off-tgt'][snp]['ensembl'][condition][subset_name] = off_tgt_set_ensembl
                intersection_dict['off-tgt'][snp]['entrez'][condition][subset_name] = off_tgt_set_entrez

            out_df = pd.DataFrame(out_dict)
            off_tgt_out_dict[snp][condition] = out_df

            off_tgt_summary_rows.append(off_tgt_summary_row)
            
            print(f'\n{condition} {snp} Off-Targets:')
            print(tabulate(out_df, headers='keys', tablefmt='fancy_grid'))

        off_tgt_summary_df = pd.DataFrame(off_tgt_summary_rows, columns=(['condition'] + subsets))
        off_tgt_summary_dict[snp] = off_tgt_summary_df
        print(f'\nSummary, {snp} Off-Targets Frequency:')
        print(tabulate(off_tgt_summary_df, headers='keys', tablefmt='fancy_grid'))
        print(f'\nCalculated Frequency of miRTarBase Targets in Gene Subsets.\n')

    #--- Step 5: write output ---
    print(f'Writing Outputs...')
    # out paths
    mirtar_out_path = os.path.join(deseq_results_dir, f'{experiment}.degs.mirna-freq.xlsx')
    with pd.ExcelWriter(mirtar_out_path, engine='openpyxl') as writer:
        mirtar_summary_df.to_excel(writer, sheet_name='summary', index=False)
        for condition, out_df in mirtar_out_dict.items():
            out_df.to_excel(writer, sheet_name=condition, index=False)
    print(f'\nmiRTarBase Data at:')
    print(f'{mirtar_out_path}\n')

    for snp in snps:
        off_tgt_out_path = os.path.join(deseq_results_dir, f'{experiment}.degs.off-tgt-freq.{snp}.xlsx')
        off_tgt_summary_df = off_tgt_summary_dict[snp]
        with pd.ExcelWriter(off_tgt_out_path, engine='openpyxl') as writer:
            off_tgt_summary_df.to_excel(writer, sheet_name='summary', index=False)
            for condition, out_df in off_tgt_out_dict[snp].items():
                out_df.to_excel(writer, sheet_name=condition, index=False)
        print(f'\n{snp} Off-Target Data at:')
        print(f'{off_tgt_out_path}\n')
    
    print(f'Wrote Outputs.\n')