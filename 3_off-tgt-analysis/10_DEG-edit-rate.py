###########
# imports #
###########
import pandas as pd
import numpy as np
import os
import openpyxl
#############
# functions #
#############
def calc_off_tgt_edit_rates (biocondition, off_tgt_df, sig_deg_df, drop_na_df, sample_map_df):
    off_tgt_df['gene_id'] = off_tgt_df['gene_id'].fillna('')
    off_tgt_df['mean_pct_snp'] = 0

    search_mask = sample_map_df['condition'] == biocondition
    for sample in sample_map_df.loc[search_mask, 'sample']:
        off_tgt_df['mean_pct_snp'] = off_tgt_df['mean_pct_snp'] + off_tgt_df[f'{sample}_pct_snp']
    off_tgt_df['mean_pct_snp'] = off_tgt_df['mean_pct_snp']/len(sample_map_df.loc[search_mask, 'sample'])

    print('Off Targets:\n')
    print(off_tgt_df.head())
    print(f'\nSignificant DEGs:\n')
    print(sig_deg_df.head())

    deg_upreg_mask = sig_deg_df['log2FoldChange'] > 0
    deg_downreg_mask = sig_deg_df['log2FoldChange'] < 0

    off_tgt_gene_count_df = pd.DataFrame()
    off_tgt_gene_set = off_tgt_df['gene_id'].apply(lambda x: x.split(','))
    off_tgt_gene_set = off_tgt_gene_set.sum()
    off_tgt_gene_set = set(off_tgt_gene_set)
    off_tgt_gene_count_df['gene_id'] = list(off_tgt_gene_set)
    off_tgt_gene_count_df = off_tgt_gene_count_df.set_index(off_tgt_gene_count_df['gene_id'])
    off_tgt_gene_count_df['off_tgt_count'] = 0

    for off_tgt_gene_id in off_tgt_gene_set:
        search_mask = off_tgt_df['gene_id'].apply(lambda x: off_tgt_gene_id in x)
        off_tgt_gene_count_df.loc[off_tgt_gene_id, 'off_tgt_count'] = len(off_tgt_df[search_mask])

    deg_merge_df = pd.merge(left=sig_deg_df.reset_index(drop=True), right=off_tgt_gene_count_df.reset_index(drop=True),on=['gene_id'] , how='inner')
    deg_merge_df['off_tgt_count'] = deg_merge_df['off_tgt_count'].fillna(0)

    up_merge_df = pd.merge(left=sig_deg_df[deg_upreg_mask].reset_index(drop=True), right=off_tgt_gene_count_df.reset_index(drop=True),on=['gene_id'] , how='inner')
    up_merge_df['off_tgt_count'] = up_merge_df['off_tgt_count'].fillna(0)

    down_merge_df = pd.merge(left=sig_deg_df[deg_downreg_mask].reset_index(drop=True), right=off_tgt_gene_count_df.reset_index(drop=True),on=['gene_id'] , how='inner')
    down_merge_df['off_tgt_count'] = down_merge_df['off_tgt_count'].fillna(0)

    summary_df = pd.DataFrame()
    summary_df['subset'] = [
        'Assayed Genes',
        'DEGs',
        'Up-regulated DEGs',
        'Down-regulated DEGs',
    ]
    summary_df['gene_count'] = [
        len(drop_na_df),
        len(sig_deg_df),
        len(sig_deg_df[deg_upreg_mask]),
        len(sig_deg_df[deg_downreg_mask])
    ]
    summary_df['off_tgt_count'] = [
        len(off_tgt_df),
        sum(deg_merge_df['off_tgt_count']),
        sum(up_merge_df['off_tgt_count']),
        sum(down_merge_df['off_tgt_count'])
    ]
    summary_df['off_tgt_pct'] = round((summary_df['off_tgt_count']/summary_df['gene_count'])*100,2)

    off_tgt_deg_mask = off_tgt_df['gene_id'].apply(lambda x: len(set(x.split(',')).intersection(set(sig_deg_df['gene_id']))) > 0)
    off_tgt_upreg_mask = off_tgt_df['gene_id'].apply(lambda x: len(set(x.split(',')).intersection(set(sig_deg_df.loc[deg_upreg_mask, 'gene_id']))) > 0)
    off_tgt_downreg_mask = off_tgt_df['gene_id'].apply(lambda x: len(set(x.split(',')).intersection(set(sig_deg_df.loc[deg_downreg_mask, 'gene_id']))) > 0)

    pct_snp_df = off_tgt_df[['chrom', 'pos', 'ref', 'alt', 'gene_id', 'gene_name', 'mean_pct_snp']].copy()
    pct_snp_df['deg_pct_snp'] = off_tgt_df['mean_pct_snp']
    pct_snp_df['upreg_pct_snp'] = off_tgt_df['mean_pct_snp']
    pct_snp_df['downreg_pct_snp'] = off_tgt_df['mean_pct_snp']
    pct_snp_df.loc[~off_tgt_deg_mask, 'deg_pct_snp'] = 0
    pct_snp_df.loc[~off_tgt_upreg_mask, 'upreg_pct_snp'] = 0
    pct_snp_df.loc[~off_tgt_downreg_mask, 'downreg_pct_snp'] = 0
    for pct_snp_col in ['deg_pct_snp', 'upreg_pct_snp', 'downreg_pct_snp']:
        pct_snp_df[pct_snp_col] = pct_snp_df[pct_snp_col].apply(lambda x: np.nan if x == 0 else x)


    print(f'\nTotal Expressed Genes:              \t{len(drop_na_df)}')
    print(f'Total Off-Tgt Hits:                   \t{len(off_tgt_df)}')
    print(f'Unique Genes in Off-Tgt Hits:         \t{len(off_tgt_gene_set)}')
    print(f'Total Significant DEG Count:          \t{len(sig_deg_df)}')
    print(f'Significant DEGs Upregulation Count:  \t{len(sig_deg_df[deg_upreg_mask])}')
    print(f'Significant DEGs Downregulation Count:\t{len(sig_deg_df[deg_downreg_mask])}\n')
    print(summary_df)
    print('\n')

    return summary_df, pct_snp_df

#############
# main body #
#############
if __name__ == '__main__':
    #------Manually defined variables------
    # - you should only need to change this section
    
    # full path to ALTER-code/5_tutorial-workflows/1_degs-off-tgts
    proj_dir = '/Users/rflusi/Library/CloudStorage/OneDrive-Stanford/my-documents/work-and-school/2022_Stanford/2022-02-24_banik-lab/experimental/RNAseq/1_pipeline-templates/manuscript/ALTER-code/5_tutorial-workflows/1_degs-off-tgts'
    # number of replicates
    replicates = 3
    # ref condition
    ref_condition = '01_transfection.control'
    

    #------Auto-defined variables------
    salmon_dir = os.path.join(proj_dir, 'salmon-results')
    deseq2_dir = os.path.join(salmon_dir, 'deseq2')
    proc_dir = os.path.join(proj_dir, 'init-processing')
    off_tgt_dir = os.path.join(proc_dir, 'off-tgt-analysis', f'{replicates}-reps')
    
    hit_filt_name = 'var_VOI_DP_GQ_non-wt'
    deg_dir = os.path.join(deseq2_dir, '3_results', '2_result-tables')
    deg_edit_dir = os.path.join(off_tgt_dir, 'deg-analysis')
    sample_map_path = os.path.join(deseq2_dir, '1_inputs', 'sample-map.tsv')

    #------Make dirs------
    os.makedirs(deg_edit_dir, exist_ok=True)

    #------Script------
    summary_df_dict = {}
    sample_map_df = pd.read_csv(sample_map_path, sep='\t')
    for biocondition in sample_map_df['condition'].unique():
        if biocondition != ref_condition:
            off_tgt_path = os.path.join(off_tgt_dir, biocondition, f'{biocondition}-{hit_filt_name}.tsv.gz')
            off_tgt_df = pd.read_csv(off_tgt_path, sep='\t')

            sig_deg_path = os.path.join(deg_dir, f'{biocondition}_de-results-padj0.01-lfc1.tsv')
            sig_deg_df = pd.read_csv(sig_deg_path, sep='\t')
            sig_deg_df = sig_deg_df.set_index(sig_deg_df['gene_id'])

            drop_na_path = os.path.join(deg_dir, f'{biocondition}_de-results-drop-na.tsv')
            drop_na_df = pd.read_csv(drop_na_path, sep='\t')

            summary_df, pct_snp_df = calc_off_tgt_edit_rates(biocondition=biocondition,
                                                 off_tgt_df=off_tgt_df, sig_deg_df=sig_deg_df, drop_na_df=drop_na_df, sample_map_df=sample_map_df
                                                 )
            
            out_dir = os.path.join(deg_edit_dir, biocondition)
            os.makedirs(out_dir, exist_ok=True)

            summary_df.to_csv(os.path.join(out_dir, 'deg-off-tgts-summary.tsv'), sep='\t', index=False)
            pct_snp_df.to_csv(os.path.join(out_dir, 'deg-off-tgts-pct-snp.tsv'), sep='\t', index=False, float_format='%.2f')

            summary_df_dict[biocondition] = summary_df

    excel_path = os.path.join(deg_edit_dir, 'deg-edit-rates.xlsx')
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        all_pct_df = pd.DataFrame(columns=['condition', "Assayed Genes",'DEGs',"Up-regulated DEGs", "Down-regulated DEGs"])
        nex_idx = 0
        for biocondition, summary_df in summary_df_dict.items():
            sheet_name = biocondition.split('_')[1]
            summary_df.to_excel(writer, sheet_name=sheet_name, index=False)
            all_pct_df.loc[nex_idx] = [sheet_name] + list(summary_df['off_tgt_pct'])
            nex_idx += 1
            
        all_pct_df.to_excel(writer, sheet_name='pct_summary', index=False)