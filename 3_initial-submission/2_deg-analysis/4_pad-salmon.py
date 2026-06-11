###########

# imports #
###########
import json
import pandas as pd
import os
import pysam
from functools import reduce

#############
# functions #
#############

def pad_raw_counts(samples, expr_dir, quant_file_name):
    # -------------------------------------------------
    # Input:    sample list, expression directory, name of the quantification file to be padded, without extension
    #
    # Output:   determines a full list of contigs present in any sample and pads all sample data with rows for 
    #           the missing contigs where all counts are 0.
    # -------------------------------------------------
    print(f'\n[RUN] Padding Quantification Files: {quant_file_name}.sf\n')

    if quant_file_name == 'quant':
        merge_col_list = ['Name', 'Length']
    elif quant_file_name == 'quant.genes':
        merge_col_list = ['Name']
    
    all_contig_col = pd.DataFrame()
    for sample in samples:
        sample_raw_counts_path = os.path.join(expr_dir, sample, f'{quant_file_name}.sf')
        sample_contigs_col = pd.read_csv(sample_raw_counts_path, sep='\t')
        sample_contigs_col = sample_contigs_col[merge_col_list]

        if all_contig_col.empty:
            all_contig_col = sample_contigs_col
        else:
            all_contig_col = pd.merge(left=all_contig_col, right=sample_contigs_col, on=merge_col_list, how='outer')
    
    # all_contig_col is now a single column with all contigs present in any sample

    for sample in samples:
        print(f'\n\t{sample}:\n')

        sample_raw_counts_path = os.path.join(expr_dir, sample, f'{quant_file_name}.sf')
        sample_raw_counts_df = pd.read_csv(sample_raw_counts_path, sep='\t')
        sample_padded_df = pd.merge(left=all_contig_col, right=sample_raw_counts_df, on=merge_col_list, how='outer')

        display_mask = sample_padded_df['NumReads'].apply(pd.isna)
        print(f'\t\tContigs added:\n')
        print(sample_padded_df.loc[display_mask, 'Name'])
        print('\n')

        sample_padded_df = sample_padded_df.fillna(0)
        for padded_col in sample_padded_df.columns:
            test_mask = sample_padded_df[padded_col].apply(pd.isna)
            if not sample_padded_df[test_mask].empty:
                print(f'\n\t\t[WARNING] nan values in {padded_col}')

        sample_out_dir = os.path.join(expr_dir, 'padded-results', sample)
        os.makedirs(sample_out_dir, exist_ok=True)
        padded_out_path = os.path.join(sample_out_dir, f'{quant_file_name}.sf')
        sample_padded_df.to_csv(padded_out_path, sep='\t', index=False)

    print(f'\n[DONE] Padding Quantification Files: {quant_file_name}.sf\n')
    return None

def pull_gene_name_with_gene_id (gene_id, transcript_gene_map):
    gene_id_mask = transcript_gene_map['gene_id'] == gene_id
    gene_name = transcript_gene_map.loc[transcript_gene_map[gene_id_mask].index[0], 'gene_name']
    return gene_name

def pull_transcript_ids_with_gene_id (gene_id, transcript_gene_map):
    gene_id_mask = transcript_gene_map['gene_id'] == gene_id
    transcript_ids = transcript_gene_map.loc[gene_id_mask, 'transcript_ids']
    transcript_ids = ','.join(list(transcript_ids.values))
    return transcript_ids
    
def build_raw_count_df(samples, transcript_gene_map, gene_name_map, expr_dir, quant_file_name):
    # -------------------------------------------------
    # Input:    - sample list
    #           - transcript to gene map, a df read from a tsv file with transcript_id, gene_id, and gene_name
    #           - gene_id to gene_name map, a df mapping this information
    #           - expression directory
    #           - name of the quantification file to be padded, without extension
    #
    # Output:   aggregates raw count data into a single df for reading and visualization
    # -------------------------------------------------
    print(f'\n[RUN] Building combined raw count df: {quant_file_name}.sf\n')

    padded_dir = os.path.join(expr_dir, 'padded-results')

    if quant_file_name == 'quant':
        merge_col_list = ['transcript_id']
    elif quant_file_name == 'quant.genes':
        merge_col_list = ['gene_id']

    raw_count_df = pd.DataFrame()
    for sample in samples:
        id_col = 'transcript_id' if quant_file_name == 'quant' else 'gene_id'

        rename_dict = {
            'Name':id_col,
            'Length':f'{sample}_length',
            'EffectiveLength':f'{sample}_eff_length',
            'TPM':f'{sample}_tpm',
            'NumReads':f'{sample}_reads',
        }

        sample_raw_counts_path = os.path.join(padded_dir, sample, f'{quant_file_name}.sf')
        sample_raw_count_df = pd.read_csv(sample_raw_counts_path, sep='\t')
        sample_raw_count_df = sample_raw_count_df.rename(columns=rename_dict)
        
        if raw_count_df.empty:
            raw_count_df = sample_raw_count_df[[id_col, f'{sample}_tpm', f'{sample}_reads']]
        else:
            raw_count_df = pd.merge(left=raw_count_df, right=sample_raw_count_df[[id_col, f'{sample}_tpm', f'{sample}_reads']], on=merge_col_list, how='outer')

    for sample in samples:
        nan_mask = raw_count_df[f'{sample}_reads'].apply(pd.isna)
        if not raw_count_df[nan_mask].empty:
            print(f'\n\t[WARNING] nan values found:\n')
            print(raw_count_df.loc[nan_mask, [id_col, f'{sample}_reads']])

    if quant_file_name == 'quant':
        raw_count_df['gene_id'] = raw_count_df['transcript_id'].apply(lambda x: transcript_gene_map.loc[x, 'gene_id'])
        raw_count_df['gene_name'] = raw_count_df['transcript_id'].apply(lambda x: transcript_gene_map.loc[x, 'gene_name'])

    elif quant_file_name == 'quant.genes':
        test_mask = raw_count_df['gene_id'].apply(lambda x: 'ENST' in x)

        if reduce(lambda x,y: x or y, test_mask.values):
            print(print(f'\n\t[WARNING] ENST contigs found in gene quant files\n'))

        raw_count_df['gene_name'] = raw_count_df['gene_id'].apply(lambda x: gene_name_map.loc[x, 'gene_name'])
        # raw_count_df['transcript_id'] = raw_count_df['gene_id'].apply(pull_transcript_ids_with_gene_id, transcript_gene_map=transcript_gene_map)
        
    print(raw_count_df)
    print('\n')

    print(f'\n[DONE] Building combined raw count df: {quant_file_name}.sf\n')
    
    return raw_count_df


#############
# main body #
#############
if __name__ == '__main__':
    
    ##########################
    # User-Defined Variables #
    ##########################
    # - you should only need to make changes in this section

    # full path to tutorial workflow directory
    proj_dir = ''
    # full path to transcript_id map made during ref file curation
    transc_map_path = ''

    # create output directories
    sample_map_path = os.path.join(proj_dir, 'sample-map.tsv')
    salmon_results_dir = os.path.join(proj_dir, 'salmon-results')
    padded_dir = os.path.join(salmon_results_dir, 'padded-results')
    os.makedirs(padded_dir, exist_ok=True)

    combined_dir = os.path.join(salmon_results_dir, 'combined-results')
    os.makedirs(combined_dir, exist_ok=True)

    deseq2_dir = os.path.join(salmon_results_dir, 'deseq2')
    os.makedirs(deseq2_dir, exist_ok=True)

    deseq2_input_dir = os.path.join(deseq2_dir, '1_inputs')
    os.makedirs(deseq2_input_dir, exist_ok=True)

    # read in sample and transcript maps
    sample_map_df = pd.read_csv(sample_map_path, sep='\t')
    samples = list(sample_map_df['sample'])
    
    transc_map_df = pd.read_csv(transc_map_path, sep='\t', compression='infer', index_col='transcript_id')

    # create a map that maps gene_id to gene_name
    gene_map_df = pd.Series(zip(transc_map_df['gene_id'], transc_map_df['gene_name']))
    gene_map_df = gene_map_df.unique()
    gene_id_list = []
    gene_name_list = []
    for gene_id, gene_name in gene_map_df:
        gene_id_list.append(gene_id)
        gene_name_list.append(gene_name)
    gene_map_df = pd.DataFrame()
    gene_map_df.index = gene_id_list
    gene_map_df['gene_name'] = gene_name_list
    
    # pad sample salmon counts with contigs not in each individaual sample reference genome
    # padded files written to padded_dir
    pad_raw_counts(samples=samples, expr_dir=salmon_results_dir, quant_file_name='quant')
    pad_raw_counts(samples=samples, expr_dir=salmon_results_dir, quant_file_name='quant.genes')

    # build combined raw count dfs from salmon transcript and gene level results
    transcript_raw_count_df = build_raw_count_df(samples=samples, transcript_gene_map=transc_map_df, gene_name_map=gene_map_df, expr_dir=salmon_results_dir, quant_file_name='quant')
    gene_raw_count_df = build_raw_count_df(samples=samples, transcript_gene_map=transc_map_df, gene_name_map=gene_map_df, expr_dir=salmon_results_dir, quant_file_name='quant.genes')
    
    # write combined dfs to combined_dir
    transcript_raw_tsv_path = os.path.join(combined_dir, 'transcript-raw.tsv')
    transcript_raw_count_df.to_csv(transcript_raw_tsv_path, index=False, sep='\t')

    gene_raw_tsv_path = os.path.join(combined_dir, 'gene-raw.tsv')
    gene_raw_count_df.to_csv(gene_raw_tsv_path, index=False, sep='\t')

    deseq_sample_map_path = os.path.join(deseq2_input_dir, 'sample-map.tsv')
    sample_map_df.to_csv(deseq_sample_map_path, sep='\t', index=False)