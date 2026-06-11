#####################
# import statements #
#####################
import os
import pandas as pd
import numpy as np
import tqdm as tqdm
from Bio.Seq import Seq
import pysam
import math
from functools import reduce
import sys
import tqdm
from tabulate import tabulate
import tempfile
import shutil

#############
# functions #
#############
def import_VariantsToTable_tsv(path, samples):
    ###################################################################################################
    # Purpose: import variant table output from GATK, standardize datatypes and adjust column headers #
    # Inputs: 1. path - string, path of GATK output tsv file                                          #
    #         2. samples - list of strings with sample names                                          #
    # Output: formatted pandas dataframe created from the tsv data                                    #
    ###################################################################################################

    dtype_dict = {'CHROM':'str', 'POS':'int', 'REF':'str', 'ALT':'str', 'QUAL':'float', 'FILTER':'str'}
    for sample in samples:
        dtype_dict[f'{sample}.GT'] = 'str'
        dtype_dict[f'{sample}.AD'] = 'str'
        dtype_dict[f'{sample}.DP'] = 'float'
        dtype_dict[f'{sample}.GQ'] = 'float'
    df = pd.read_csv(path, sep='\t', dtype=dtype_dict, compression='infer')
    rename_dict = {}
    for column in df.columns:
        if '.' not in column:
            rename_dict[column] = column.lower()
    df = df.rename(columns=rename_dict)

    return df

def strand_id(var_row, gtfgz):
    ###################################################################################################
    # Purpose: take variant information from a row of the variant dataframe and obtain relevant       #
    #          feature information from the reference gtf                                             #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. gtfgz - gtfgz file imported with pysam.TabixFile()                                   #
    #         3. feat_priority - a list of feature types in priority order for assigning strands (ie. #
    #                            if feat_priority=['exon','transcript'] and a variant location has    #
    #                            transcript features in both strands but exon features in only one,   #
    #                            the exon strand will be assigned)                                    #
    # Output: values for strand, gene_name, gene_id, transcript_id, exon_id, and warning notes        #
    ###################################################################################################
    
    # the warning notes will be kept as a list then combined into a string at the end
    warning_str = []
    hit_dict_list = []
    strand = ''

    try:
        chrom = var_row['chrom']
        query_start = var_row['pos']
        query_end = var_row['pos'] + 1 

        # build a dict of every gtf feature at pos
        pos_gtf_recs = {
            'exon':[],
            'transcript':[],
            'gene':[]
        }

        for gtfgz_row in gtfgz.fetch(chrom, query_start, query_end):
            gtfgz_rec = gtfgz_row.strip().split('\t')
            feat_type = gtfgz_rec[2].strip()
            if feat_type in ['exon', 'transcript', 'gene']:
                pos_gtf_recs[feat_type].append(gtfgz_rec)

        # assemble a list of hits hits are only related to the feature type of highest priority exon > transcript > gene
        if len(pos_gtf_recs['exon']) > 0:
            target_feature = 'exon'
            
        elif len(pos_gtf_recs['transcript']) > 0:
            warning_str.append('no exons')
            target_feature = 'transcript'

        elif len(pos_gtf_recs['gene']) > 0:
            warning_str.append('no exons')
            warning_str.append('no transcripts')
            target_feature = 'gene'
        else:
            warning_str.append('no features')
            target_feature = ''

        if target_feature != '':
            for gtf_rec in pos_gtf_recs[target_feature]:
                hit_dict = {}
                hit_dict['strand'] = gtf_rec[6].strip()
                hit_dict['exon'] = gtf_rec[8].split('exon_id')[1].split('"')[1].strip() if target_feature == 'exon' else ''
                hit_dict['transcript'] = gtf_rec[8].split('transcript_id')[1].split('"')[1].strip() if target_feature in ['exon', 'transcript'] else ''
                hit_dict['gene'] = (
                    gtf_rec[8].split('gene_id')[1].split('"')[1].strip(),
                    gtf_rec[8].split('gene_name')[1].split('"')[1].strip()
                )
                hit_dict_list.append(hit_dict)
        
            # determine if all hits are in the same strand
            strand_set = set()
            for hit_dict in hit_dict_list:
                strand_set.add(hit_dict['strand'])
            strand = ''.join(strand_set) if len(strand_set) == 1 else 'ambiguous'
    except ValueError:
        warning_str.append('ValueError in gtf')
    
    return strand, hit_dict_list, ','.join(warning_str)
    
def strand_sequences(var_row, samples):
    ###################################################################################################
    # Purpose: adjust variant and GT sequences if the coding strand is -                              #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. samples - a list of sample names                                                     #
    # Output: values for ref, alt and all GT columns with appropriate sequences                       #
    ###################################################################################################
    
    # if the coding strand is + simply return the existing sequence data unchanged
    if var_row['strand'] != '-':
        return var_row['ref'], var_row['alt'], *[var_row[f'{sample}.GT'] for sample in samples]
    elif var_row['strand'] == '-':
        init_ref = var_row['ref']
        init_alt = var_row['alt']

        # get the reverse complement sequences for ref and alt, check for multiple values in each and handle that if necessary
        rt_list = []
        for alt in [init_ref, init_alt]: # I know the variable naming here is confusing, just roll with it
            if ('|' not in alt) and (',' not in alt):
                alt_rc = str(Seq(alt).reverse_complement()).strip()
            elif ('|' in alt) and (',' in alt):
                alt_list = alt.split('|')
                alt_list = [single_alt.strip() for single_alt in alt_list]

                for i in range(len(alt_list)):
                    if ',' in alt_list[i]:
                        alt_sublist = alt_list[i].split(',')
                        alt_sublist = [single_alt.strip() for single_alt in alt_sublist]
                        new_sublist = []
                        for single_alt in alt_sublist:
                            new_sublist.append(str(Seq(single_alt).reverse_complement()).strip())
                        new_sublist = ','.join(new_sublist)
                        alt_list[i] = new_sublist
                    else:
                        alt_list[i] = str(Seq(alt_list[i]).reverse_complement().strip())

                alt_rc = '|'.join(alt_list)
            else:
                alt_delim = '|' if ('|' in alt) else ','

                alt_list = alt.split(alt_delim)
                alt_list = [single_alt.strip() for single_alt in alt_list]
                alt_list = [(str(Seq(single_alt).reverse_complement().strip())) for single_alt in alt_list]
                
                alt_rc = alt_delim.join(alt_list)
                
            rt_list.append(alt_rc)

        # get the reverse complement sequences for GT values
        gt_rcs = []
        for sample in samples:
            gt = var_row[f'{sample}.GT']
            if not pd.isna(gt):    
                if '/' in gt:
                    gt_delim = '/'
                    alls = gt.split(gt_delim)
                elif '|' in gt:
                    gt_delim = '|'
                    alls = gt.split(gt_delim)
                else:
                    print(gt)
                gt_rc = []
                for all in alls:
                    gt_rc.append(str(Seq(all).reverse_complement().strip()))
                gt_rc = gt_delim.join(gt_rc)
                gt_rcs.append(gt_rc)
            else:
                gt_rc = './.'
                gt_rcs.append(gt_rc)
            
        return rt_list[0], rt_list[1], *gt_rcs
    
def id_sample_nts(var_row, sample, all_samples):
    ###################################################################################################
    # Purpose: add sample specific ref and alt columns                                                #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. sample - the relevant sample                                                         #
    #         3. all_samples - a list of all samples                                                  #
    # Output: for each sample, values for the reference nucleotide and genotype nucleotides 1 and 2   #
    ###################################################################################################
    if var_row[f'{sample}.GT'].strip() == './.':
        return '','',''
    else:
        init_ref = var_row['ref'].strip()
        init_gt = var_row[f'{sample}.GT'].strip()
        
        gt_delim = '|' if ('|' in init_gt) else '/'
        sample_alls = init_gt.split(gt_delim)
        sample_alls = [sample_all.strip() for sample_all in sample_alls]

        if (',' not in init_ref) and ('|' not in init_ref):
            return init_ref, sample_alls[0], sample_alls[1]
        elif ('|' in init_ref) or (',' in init_ref):
            grouped_refs = init_ref.split('|')
            grouped_refs = [refs.strip() for refs in grouped_refs]
            grouped_refs = [refs.split(',') for refs in grouped_refs]

            if len(grouped_refs) == len(all_samples):
                sample_ref = grouped_refs[all_samples.index(sample)]
                sample_ref = ','.join(sample_ref)
                return sample_ref, sample_alls[0], sample_alls[1]
            
            refs_list = []
            for ref_list in grouped_refs:
                for ref in ref_list:
                    refs_list.append(ref.strip())
            
            init_alt = var_row['alt'].strip()
            alts_list = init_alt.split('|')
            alts_list = [alt.strip() for alt in alts_list]
            alts_list = ','.join(alts_list)
            alts_list = alts_list.split(',')
            alts_list = [alt.strip() for alt in alts_list]
            
            ref_alt_overlap = set(refs_list).intersection(set(alts_list))
            if len(ref_alt_overlap) > 0:
                var_samples = []
                for samp in all_samples:
                    if var_row[f'{samp}.GT'] != './.':
                        var_samples.append(samp)

                if len(grouped_refs) == len(var_samples):
                    sample_ref = grouped_refs[var_samples.index(sample)]
                    sample_ref = ','.join(sample_ref)
                elif (sample_alls[0] not in ref_alt_overlap) or (sample_alls[1] not in ref_alt_overlap):
                    sample_ref = ''
                    for sample_all in sample_alls:
                        if sample_all not in ref_alt_overlap:
                            if sample_all in refs_list:
                                if (sample_ref != '') and (sample_ref != sample_all):
                                    print(f'[WARNING] Conflicting refs:\tindex:{var_row.name}ref:{init_ref}\talt:{init_alt}\tGT:{init_gt}')
                                    return '', sample_alls[0], sample_alls[1]
                                sample_ref = sample_all
                else:
                    sample_ref = ''

            else:
                sample_ref = ''
                for sample_all in sample_alls:
                    if sample_all in refs_list:
                        if (sample_ref != '') and (sample_ref != sample_all):
                            print(f'[WARNING] Conflicting refs:\tindex:{var_row.name}ref:{init_ref}\talt:{init_alt}\tGT:{init_gt}')
                            return '', sample_alls[0], sample_alls[1]
                        sample_ref = sample_all
            return sample_ref, sample_alls[0], sample_alls[1]
        
def check_in_clinvar(var_row, clinvar_vcf, samples):
    ###################################################################################################
    # Purpose: add columns with ClinVar data to each variant entry                                    #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. clinvar_vcf - ClinVar vcf file read with pysam.VariantFile()                         #
    #         3. samples - a list of all samples                                                      #
    # Output: 1. in_clinvar - boolean of whether the variant appears in clinvar                       #
    #         2. clinvar_ids - all clinvar ids associated with the variant                            #
    #         3. clinvar_sig - the variant significance as annotated in clinvar                       #
    #         4. clinvar_disease - diseases associated with the variant in clinvar                    #
    ###################################################################################################

    chrom = var_row['chrom']
    pos = var_row['pos']
    strand = var_row['strand']

    in_clinvar = False
    clinvar_ids = set()
    clinvar_sig = set()
    clinvar_disease = []

    for sample in samples:
        ref = var_row[f'{sample}_ref']
        sample_alts = set()
        for sample_all in [var_row[f'{sample}_all_1'], var_row[f'{sample}_all_2']]:
            if (sample_all != ref):
                sample_alts.add(sample_all)

        # check whether the variant is in clinvar and pull the related clinvar information if so
        if len(sample_alts) > 0:
            try:
                for clinvar_rec in clinvar_vcf.fetch(chrom, pos-1, pos): # fetch in 0 index vcf is 1 index
                    if clinvar_rec.alts != None:
                        if (clinvar_rec.alleles[0] == ref) and (len(sample_alts.intersection(set(clinvar_rec.alts))) > 0):
                            in_clinvar = True
                            clinvar_ids.add(clinvar_rec.id)
                            if 'CLNSIG' in clinvar_rec.info.keys():
                                clinvar_sig.add(','.join(clinvar_rec.info['CLNSIG']))
                            if 'CLNDN' in clinvar_rec.info.keys():
                                clinvar_disease.append(','.join(clinvar_rec.info['CLNDN']))

            except ValueError:
                pass

    clinvar_ids = ','.join(clinvar_ids)
    clinvar_sig = ','.join(clinvar_sig)
    clinvar_disease = ','.join(clinvar_disease)

    return in_clinvar, clinvar_ids, clinvar_sig, clinvar_disease

def get_bin_gt(var_row, sample):
    ###################################################################################################
    # Purpose: convert all genotypes to binary genotyping (0 for ref 1 for alt)                       #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. sample - the sample for which the genotype should be converted                       #
    # Output: The binary genotype value for the specified sample                                      #
    ###################################################################################################
    
    sample_gt = var_row[f'{sample}.GT'].strip()
    if sample_gt == './.':
        return './.'
    else:
        sample_ref = var_row[f'{sample}_ref'].strip()

        gt_delim = ''
        for delim in ['/', '|']:
            gt_delim = delim if (delim in sample_gt) else gt_delim
        
        sample_alls = sample_gt.split(gt_delim)
        sample_alls = [sample_all.strip() for sample_all in sample_alls]
        bin_gt = ['.', '.']

        for i in range(len(sample_alls)):
            bin_gt[i] = str(0) if (sample_alls[i] == sample_ref) else str(1)
        
        return gt_delim.join(bin_gt)

def get_value_count_col(count_index, value_col_name, target_col, df_slice):
    ###################################################################################################
    # Purpose: return a 1 column dataframe that can be used to build a value count table across       #
    #           multiple columns with the same set of values                                          #
    # Inputs: 1. count_idex - list of values to be counted                                            #
    #         2. value_col_name - column header for the output column                                 #
    #         3. target_col - name of the column whose values will be counted                         #
    #         4. df_slice - dataframe with the target column                                          #
    # Output: A dataframe with one column that is titled value_col_name and contains counts, in order #
    #         for the values in count_index                                                           #
    ###################################################################################################

    df_counts = df_slice[target_col].value_counts()
    df_count_list = []

    for count_idx in count_index:
        if count_idx in df_counts.index:
            df_count_list.append(int(df_counts[count_idx]))
        else:
            df_count_list.append(0)
            
    return pd.DataFrame({value_col_name:df_count_list})

def calc_pct_ref(ad_str):
    ###################################################################################################
    # Purpose: calculate the % of reads matching the reference nucleotide for a given AD value        #
    # Inputs: 1. ad_str - an AD value                                                                 #
    # Output: a %_ref value which is nan if the ad value is also nan                                  #
    ###################################################################################################
    
    if not pd.isna(ad_str):
        ads_list = ad_str.strip().split(',')
        ads_list = [int(var_ad) for var_ad in ads_list]
        total_ad = sum(ads_list)
        if total_ad > 0:
            ref_ad = ads_list[0]
            pct_ref = round(float((ref_ad/total_ad)*100), 2)
        else:
            pct_ref = np.nan

        return pct_ref
    else:
        return np.nan

def calc_pct_snp(var_row, sample, target_snps):
    ###################################################################################################
    # Purpose: calculate the % of reads matching the target SNP if the sample shows the target SNP    #
    # Inputs: 1. var_row - a subset of the dataframe row passed with pd.apply(axis=1)                 #
    #         2. sample - the sample for which to target %_snp                                        #
    #         3. target_snps - a list of the target SNPs                                              #
    # Output: a %_snp value, this will be nan if the sample does not have a variant that matches a    #
    #         target SNP                                                                              #
    ###################################################################################################

    ad_str = var_row[f'{sample}.AD']
    sample_ref = var_row[f'{sample}_ref']
    if (not pd.isna(ad_str)):
        all_1 = var_row[f'{sample}_all_1']
        all_2 = var_row[f'{sample}_all_2']
        target_snp = target_snps[0]

        ads_list = ad_str.strip().split(',')
        ads_list = [int(var_ad) for var_ad in ads_list]
        total_ad = sum(ads_list)

        if total_ad > 0:
            if len(ads_list) == 2:
                pct_snp = round(float((ads_list[1]/total_ad)*100),2)
            # if there are more than 2 alt nucleotides for a given sample it is ignored, this generally occurs for a miniscule fraction of sample entries
            elif len(ads_list) > 3:
                pct_snp = np.nan
            # if there are 3 ADs check whether the alts match or not and calculate accordingly
            else:
                if all_1 == all_2:
                    pct_snp = round(float((sum(ads_list[1:])/total_ad)*100),2)
                else:
                    pct_snp = round(float((ads_list[([all_1, all_2].index(target_snp[1]) + 1)]/total_ad)*100),2)
        else:
            pct_snp = np.nan
    else:
        pct_snp = np.nan

    return pct_snp

def hist_table(var_df, target_cols, bucket_lims, mode):
    ###################################################################################################
    # Purpose: create a histogram table of values in multiple dataframe columns                       #
    # Inputs: 1. var_df - the dataframe holding variant data                                          #
    #         2. target_cols - the columns whose values should be counted for the table               #
    #         3. bucket_lims - a list of numerical limits for the histogram buckets                   #
    #         4. mode - if mode == 'log' the beckets will be made on log scale and values will be     #
    #                   incremented by 1 so that values of 0 are on scale                             #
    # Output: a dataframe containing the histogram table                                              #
    ###################################################################################################
    if ('log' in mode) and bucket_lims[0] == 0:
        bucket_lims[0] = 1
    
    buckets = []
    for i in range(len(bucket_lims[1:])):
        buckets.append(f'{bucket_lims[i]}-{bucket_lims[i+1]}') 

    hist_df = pd.DataFrame({
        f'bucket':buckets
    })

    for target_col in target_cols:
        bucket_counts = []
        for i in range(len(buckets)):
            bucket_counts.append(0)

        out_of_range_low_count = 0
        out_of_range_high_count = 0
        for var_idx, var_row in tqdm.tqdm(var_df.iterrows(), total=len(var_df), desc=f'Calculating Histogram Table for {target_col}...'):
            metric_val = var_row[target_col]
            if not pd.isna(metric_val):
                if metric_val >= bucket_lims[-1]:
                    metric_val = bucket_lims[-1]*0.999
                    out_of_range_high_count += 1
                elif metric_val < bucket_lims[0]:
                    metric_val = bucket_lims[0]
                    out_of_range_low_count += 1

                if mode == 'scalar':
                    bucket_idx = math.floor(metric_val/(bucket_lims[-1]/len(buckets)))
                if 'log' in mode:
                    log_n = int(mode.split('log')[1])
                    bucket_idx = math.floor(math.log(metric_val, log_n))
                try:
                    bucket_counts[bucket_idx] += 1
                except IndexError:
                    print(f'{bucket_lims}\tMetric:{target_col}\tValue:{metric_val}\tBucket:{bucket_idx}')
                    
        hist_row = pd.DataFrame({f'{target_col}':bucket_counts})
        hist_df = pd.concat([hist_df, hist_row], axis=1)
    
    return hist_df

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

def create_temporary_copies(f_path, n_copies):
    f_parts = os.path.basename(f_path).split('.')
    f_ext_list = [f_parts[-1]]
    f_parts = f_parts[:-1]
    while f_ext_list[0] in ['gz', 'tbi']:
        f_ext_list = [f_parts[-1]] + f_ext_list
        f_parts = f_parts[:-1]
        if len(f_parts) < 1:
            break
    f_ext = f'.{".".join(f_ext_list)}'

    f_name = os.path.basename(f_path).split(f_ext)[0]

    tmp_dir = tempfile.gettempdir()

    tmp_path_list = [os.path.join(tmp_dir, f'{f_name}.tmp.{n}.{f_ext}') for n in n_copies]
    for tmp_path in tmp_path_list:
        shutil.copy2(f_path, tmp_path)

    return tmp_path_list

###############
# main script #
###############
if __name__ == '__main__':
    tqdm.tqdm.pandas()

    # --- Step 1: define variables ---
    experiment = sys.argv[1]
    proj_dir = sys.argv[2]
    ref_dir = sys.argv[3]
    # target SNPs as a list of stings for example we are interested in C-to-T (C-to-U) so ['CT']
    target_snps = sys.argv[4].split(',')
    # (chrom, pos) in the format of the genomic reference
    target_edits = sys.argv[5]
    if target_edits != '':
        target_edits = [(e[0], int(e[1])) for e in [e_str.split(',') for e_str in target_edits.split(';')]]
    else:
        target_edits = [("",-1)]
    # name of the ref genome
    genome_name = os.path.split(ref_dir)[1]

    # full path to tsv output from GATK analysis
    tsv_path = os.path.join(proj_dir, '4_variants', f'{experiment}.variants.tsv.gz')
    # full path to gtf for the reference genome
    gtfgz_path = os.path.join(ref_dir, f'{genome_name}.sorted.gtf.gz')
    # full path to sample map
    sample_map_path = os.path.join(proj_dir, 'sample-map.csv')
    sample_map_df = pd.read_csv(sample_map_path, sep=',', compression='infer')

    # these must match the sample identifiers in the column titles of the variants tsv
    samples = list(sample_map_df['sample'])

    var_dir = os.path.join(proj_dir, '4_variants')
    out_dir = os.path.join(var_dir, 'raw-data')

    # output paths
    out_pq_path = os.path.join(var_dir, f'{experiment}.variants.annotated.parquet')
    analysis_out_dir = os.path.join(proj_dir, '4_variants', 'annotation-tsvs')
    os.makedirs(analysis_out_dir, exist_ok=True)

    var_df = import_VariantsToTable_tsv(tsv_path, samples=samples)

    # script

    #--- Step 2: import variant data and print relevant information to the terminal ---
    # convert nan genotypes to a string representation in var_df
    gt_cols = [f'{sample}.GT' for sample in samples]
    for col in gt_cols:
        var_df[col] = var_df[col].apply(lambda x: './.' if pd.isna(x) else x)

    # display dataframe info
    total_count = len(var_df)
    
    print(f'\nRaw Variants Total Entries:\t{total_count}')
    print(tabulate(var_df[['chrom', 'pos', 'ref', 'alt']+[f'{sample}.GT' for sample in samples]].head(), headers='keys', tablefmt='fancy_grid'))
    print('\n')

    #--- Step 3: determine variant strand and annotate with corresponding features ---
    init_dot_gt_counts = [(var_df[f'{sample}.GT'].value_counts()['./.'] if ('./.' in var_df[f'{sample}.GT'].value_counts().index) else 0) for sample in samples]

    # find the strand of the entry and hits for features at each pos
    print(f'Finding features and stranding variant hits, initial pass...')

    var_df[['strand', 'gtf_hits', 'warnings']] = var_df[['chrom', 'pos']].progress_apply(strand_id, gtfgz=pysam.TabixFile(gtfgz_path), axis=1, result_type='expand')

    # split rows with hits in both strands into one row for each strand
    ambig_strand_mask = var_df['strand'] == 'ambiguous'
    print(f'Splitting {len(var_df.loc[ambig_strand_mask])} ambiguous strands into sepparate rows for + and -...')

    new_row_list = []
    plus_hits_series = pd.Series(index=var_df.index)
    strand_list = []
    for var_idx, var_row in tqdm.tqdm(var_df.loc[ambig_strand_mask].iterrows(), total=len(var_df.loc[ambig_strand_mask]), desc='Fixing Ambiguous Strands'):
        plus_hits = []
        minus_hits = []
        for hit_dict in var_row['gtf_hits']:
            if hit_dict['strand'] == '+':
                plus_hits.append(hit_dict)
            else:
                minus_hits.append(hit_dict)
        
        minus_row = var_row.copy()
        
        var_df.at[var_idx, 'strand'] = '+'
        var_df.at[var_idx, 'gtf_hits'] = plus_hits

        minus_row['strand'] = '-'
        minus_row['gtf_hits'] = minus_hits
        new_row_list.append(minus_row)
    
    var_df = pd.concat([var_df, pd.DataFrame(new_row_list)], ignore_index=True)
    
    # reorder var df
    chrom_order = [f'chr{chr_num}' for chr_num in (list(range(1,23))+['M', 'X', 'Y'])]
    all_chroms = pd.Series(var_df['chrom'].unique())
    chrom_order = chrom_order + all_chroms[~all_chroms.isin(chrom_order)].tolist()
    
    var_df['chrom'] = pd.Categorical(var_df['chrom'], categories=chrom_order, ordered=True)
    var_df = var_df.sort_values(by=['chrom','pos'], ascending=True).reset_index(drop=True)
    
    print(f'New Total Entry Count: {len(var_df)}')
    print('')
    
    #--- Step 4: adjust ref and alt sequences based on strand ---
    # if the coding strand is assigned - adjust the sequence identity of the variant and sample genotypes
    print(f'Stranding Sequences...')
    seq_cols = ['ref', 'alt'] + [f'{sample}.GT' for sample in samples]
    var_df[seq_cols] = var_df[(['ref', 'alt', 'strand']+ [f'{sample}.GT' for sample in samples])].progress_apply(strand_sequences, axis=1, samples=samples, result_type='expand')
    print('')

    # display the number of uncalled genotypes for each sample and the new columns
    for i in range(len(samples)):
        print(f'Initial {samples[i]} ./. GTs:{init_dot_gt_counts[i]}')
    print('')

    # display stranding information for the dataset
    print(f'Total Hits: {len(var_df)}\n')

    value_counts = var_df['strand'].value_counts()
    strand_counts = pd.DataFrame({'strand':value_counts.index, 'count':value_counts.values}) 

    print(f'Strand Counts: {len(var_df)}')
    print(tabulate(strand_counts, headers='keys', tablefmt='fancy_grid'))
    print('')

    value_counts = pd.Series(zip(var_df['strand'], var_df['warnings'])).value_counts()
    strand_warning_counts = pd.DataFrame({'strand, warnings':value_counts.index, 'count':value_counts.values})
    
    print(f'Strands + Warnings: {len(var_df)}')
    print(tabulate(strand_warning_counts, headers='keys', tablefmt='fancy_grid'))
    print('')

    # --- Step 5: tabulate per sample nucleotides ---
    nt_results = []
    for sample in samples:
        new_nt_cols = [f'{sample}_ref', f'{sample}_all_1', f'{sample}_all_2']
        
        print(f'Tabulating Sequence Information for {sample}...')
        result = var_df[['ref', 'alt'] + [f'{sample}.GT' for sample in samples]].progress_apply(id_sample_nts, sample=sample, axis=1, result_type='expand', all_samples=samples)
        result.columns = new_nt_cols
        nt_results.append(result)
        print('')

    var_df = pd.concat([var_df] + nt_results, axis=1)

    # non ./. genotyped entries without nucleotide assignments will be displayed at the bottom of the output, the table should be empty
    sample_masks = []
    entry_count = 0
    for sample in samples:
        mask = ((var_df[f'{sample}.GT'] != './.') & (var_df[f'{sample}_ref'] == ''))
        sample_masks.append(mask)
        entry_count += len(var_df[mask])
    
    print(f'Non ./. Sample x Entries that are Unassigned ({entry_count}, {round((entry_count/(len(var_df)*len(samples)))*100,2)}%):')
    print('')

    # # --- Step 6: cross reference to clinvar ---
    # print('Cross Referencing to ClinVar...')
    # clinvar_cols = ['in_clinvar', 'clinvar_ids', 'clinvar_sig', 'clinvar_disease']
    # var_df[clinvar_cols] = var_df.progress_apply(check_in_clinvar, clinvar_vcf=pysam.VariantFile(clinvar_vcf_path, 'r'), samples=samples, axis=1, result_type='expand')
    # print('')

    # # print aggregated value counts for how many entries appear in ClinVar with a given significance and the head of the dataframe with new ClinVar columns
    # print(f'Total Entries in ClinVar:\t{len(var_df[var_df['in_clinvar']])}')
    # print(tabulate(pd.DataFrame(var_df.loc[var_df['in_clinvar'], 'clinvar_sig'].value_counts()), headers='keys', tablefmt='fancy_grid'))
    # print('')

    # --- Step 7: convert GTs to binary ---
    bin_results = []
    for sample in samples:
        gt_col = f'{sample}.GT'
        bin_gt_col = f'{sample}.binGT'
        input_cols = [gt_col, f'{sample}_ref']
        
        print(f'Converting GTs to Binary for {sample}...')
        result = var_df[input_cols].progress_apply(get_bin_gt, sample=sample, axis=1)
        result = pd.DataFrame({bin_gt_col:result})
        bin_results.append(result)

    var_df = pd.concat([var_df]+bin_results, axis=1)
    print('')

    # display aggregate binary GT counts
    print(f'Binary GT totals:')
    gt_counts_df = pd.DataFrame({'binary_gt':['1/1', '1/0', '0/1', '0/0', './.']})
    for sample in samples:
        gt_col = f'{sample}.binGT'
        simple_gt_col = pd.DataFrame(var_df[gt_col].apply(lambda x: x.replace('|', '/')))
        val_count_col = get_value_count_col(value_col_name=gt_col, count_index=list(gt_counts_df['binary_gt']), target_col=gt_col, df_slice=simple_gt_col)
        gt_counts_df = pd.concat([gt_counts_df, val_count_col], axis=1)
    print(tabulate(gt_counts_df, headers='keys', tablefmt='fancy_grid'))
    print('')

    # --- Step 8: calculate % ref values for each sample in all entries and %snp values for samples matching the target snp ---
    pct_ref_results = {}
    for sample in samples:
        ad_col = f'{sample}.AD'
        
        print(f'Calculating %ref for {sample}...')
        result = var_df[ad_col].progress_apply(calc_pct_ref)
        pct_ref_results[f'{sample}_pct_ref'] = result
        
        input_cols = [ad_col, f'{sample}_ref', f'{sample}_all_1', f'{sample}_all_2']
        
        print(f'Calculating %snp for {sample}...')
        for snp in target_snps:
            snp_mask = []
            snp_mask.append(var_df[f'{sample}_ref'] + var_df[f'{sample}_all_1'] == snp)
            snp_mask.append((var_df[f'{sample}_ref'] + var_df[f'{sample}_all_2']) == snp)
            snp_mask = reduce(lambda x,y:x|y,snp_mask)
            var_df.loc[snp_mask, f'{sample}_pct_snp'] = var_df[snp_mask].progress_apply(calc_pct_snp, sample=sample, target_snps=[snp], axis=1)

    var_df = pd.concat([var_df, pd.DataFrame(pct_ref_results)], axis=1)
    print('')

    # print total var_df entries and entries*samples and counts for entries where all or some samples have %_snp values
    total_count = len(var_df) * len(samples)
    pct_snp_masks = []
    pct_snp_count = 0
    for sample in samples:
        mask = var_df[f'{sample}_pct_snp'].apply(lambda x: not pd.isna(x))
        pct_snp_masks.append(mask)
        pct_snp_count += mask.value_counts()[True]

    all_pct_snp_mask = reduce(lambda x,y: (x&y), pct_snp_masks)
    any_pct_snp_mask = reduce(lambda x,y: (x|y), pct_snp_masks)

    print(f'Total Entries * Samples: {total_count}')
    print(f'Total Assigned pct_snp: {pct_snp_count}\n')
    print(f'Total Entries: {len(var_df)}')
    print(f'Entries with all samples pct_snp: {len(var_df[all_pct_snp_mask])}')
    print(f'Entries with not all samples pct_snp: {len(var_df[any_pct_snp_mask]) - len(var_df[all_pct_snp_mask])}\n')

    # print some var_df rows with all samples assigned %_snp and some var_df rows with only some samples assigned %_snp for visual inspection
    display_cols = ['chrom', 'pos', 'ref', 'alt'] + [f'{sample}_pct_snp' for sample in samples]
    print('All samples assigned pct_snp:')
    print(tabulate(var_df.loc[all_pct_snp_mask, display_cols].head(), headers='keys', tablefmt='fancy_grid'))
    print('')
    print('Some samples assigned pct_snp:')
    print(tabulate(var_df.loc[any_pct_snp_mask, display_cols].head(), headers='keys', tablefmt='fancy_grid'))
    print('')

    # #--- Step 9: create histogram tables for %_ref, %_snp and read depth ---
    # var_df = var_df.copy() # reset any fragmentation

    # target_cols = [f'{sample}_pct_ref' for sample in samples]
    # num_buckets = 10
    # bucket_lims = list(range(0,101,int(100/num_buckets)))
    # pct_hist_df = hist_table(mode='scalar', var_df=var_df, target_cols=target_cols, bucket_lims=bucket_lims)

    # target_cols = target_cols + [f'{sample}_pct_snp' for sample in samples]
    # pct_snp_hist_cols = pd.DataFrame()
    # for sample in samples:
    #     mask = var_df[f'{sample}_pct_snp'].apply(lambda x: not pd.isna(x))
    #     new_cols = var_df.loc[mask, [f'{sample}_pct_ref', f'{sample}_pct_snp']]
    #     pct_snp_hist_cols = pd.concat([pct_snp_hist_cols, new_cols], axis=1)
    # pct_snp_hist_df = hist_table(mode='scalar', var_df=pct_snp_hist_cols, target_cols=target_cols, bucket_lims=bucket_lims)

    # target_cols = [f'{sample}.DP' for sample in samples] + [f'{sample}.GQ' for sample in samples] + ['qual']
    # num_buckets = 21
    # bucket_lims = list(2**i for i in range(num_buckets+1))
    # dp_hist_df = hist_table(mode='log2', var_df=var_df, target_cols=target_cols, bucket_lims=bucket_lims)

    # print(tabulate(pct_hist_df, headers='keys', tablefmt='fancy_grid'))
    # print('\n')
    # print(tabulate(pct_snp_hist_df, headers='keys', tablefmt='fancy_grid'))
    # print('\n')
    # print(tabulate(dp_hist_df, headers='keys', tablefmt='fancy_grid'))
    # print('\n')

    #--- Step 10: write final outputs ---
    # output var_df as a parquet
    var_df.to_parquet(out_pq_path, compression='gzip')
    print(f'Annotated Variants Written to {out_pq_path}')

    # output various metrics tables as compressed tsvs
    export_dfs = {
        'bin-gt-counts':gt_counts_df, 
        'strand-counts':strand_counts,
        'strand-warning-counts':strand_warning_counts,
        # 'pct-hist':pct_hist_df,
        # 'pct-snp-hist':pct_snp_hist_df,
        # 'dp-hist':dp_hist_df
        }

    for suffix, export_df in export_dfs.items():
        out_tsv_path = os.path.join(analysis_out_dir, f'{experiment}.variants.{suffix}.tsv.gz')
        export_df.to_csv(out_tsv_path, sep='\t', index=False, compression='gzip')

    # compile metrics tables into an excel output
    excel_path = os.path.join(proj_dir, '4_variants', f'{experiment}.variants.summary-stats.xlsx')
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        for sheet_name, export_df in export_dfs.items():
            export_df.to_excel(writer, sheet_name=sheet_name, index=False)

    # print target edit %_ref and %_snp and output full var_df row for the target edit to a tsv file
    if target_edits[0][1] != -1:
        out_col_list = ['chrom', 'pos', 'ref', 'alt']
        
        for col_type in ['.DP', '_pct_ref', '_pct_snp']:
            for sample in samples:
                out_col_list.append(f'{sample}{col_type}')
        
        tgt_mask = [((var_df['chrom'] == e[0]) & (var_df['pos'] == e[1])) & ((var_df['ref'] + var_df['alt']).isin(target_snps)) for e in target_edits]
        tgt_mask = reduce(lambda x,y: x|y, tgt_mask)
        
        tgt_snp_df = var_df.drop(var_df[~tgt_mask].index)
        
        display_col_list = []
        for col_type in ['_pct_ref', '_pct_snp']:
            for condition in sample_map_df['condition'].unique():
                mean_col = f'{condition}_mean{col_type}'
                out_col_list.append(mean_col)
                display_col_list.append(mean_col)
                tgt_snp_df[mean_col] = tgt_snp_df.apply(calc_sample_mean, col_type=col_type, condition=condition, sample_map_df=sample_map_df, axis=1)

        print('Target Edit:')
        print(tabulate(tgt_snp_df[display_col_list], headers='keys', tablefmt='fancy_grid'))
        print('')

        out_tsv_path = os.path.join(proj_dir, '4_variants', f'{experiment}.variants.tgt-edit.tsv')
        tgt_snp_df[out_col_list].to_csv(out_tsv_path, sep='\t', index=False)
