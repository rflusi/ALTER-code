import pandas as pd

def get_quant_files(sample_map, proj_dir):
    test_sample = sample_map['sample'].iloc[0]
    samples = list(sample_map['sample'])
    if os.path.isfile(f'{proj_dir}/5_expression/salmon/padded-results/{test_sample}/quant.sf.gz'):
        return expand('{proj_dir}/5_expression/salmon/padded-results/{sample}/quant.sf.gz', proj_dir=proj_dir, sample=samples)
    else:
        return expand('{proj_dir}/5_expression/salmon/padded-results/{sample}/quant.sf', proj_dir=proj_dir, sample=samples)


rule deg_analysis: 
    input:
        padded_salmon_results = lambda wc: get_quant_files(
            sample_map=SAMPLE_MAP_DICT[wc.experiment], proj_dir=os.path.join(PROJ_DIR, EXP_DIR_NAMES_DICT[wc.experiment])
        )
    output:
        raw_count_tsv = f'{PROJ_DIR}/{{exp_dir}}/5_expression/deseq2/2_txi-counts/{{experiment}}.raw-counts.gene-level.tsv.gz'
    params:
        rmd = f'{SMK_DIR}/scripts/deseq2.Rmd',
        output_format = 'html_document',
        experiment = lambda wc: wc.experiment,
        exp_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}',
        out_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}/5_expression/deseq2',
        ref_dir = REF_DIR,
        ref_db_dir = REF_DB_DIR,
        deseq_fc_cutoff = DESEQ_FC_CUTOFF,
        deseq_padj_cutoff = DESEQ_PADJ_CUTOFF,
        go_padj_cutoff = GO_PADJ_CUTOFF,
        go_q_cutoff = GO_Q_CUTOFF,
    conda:
        f'{SMK_DIR}/env/mrnaseq_deseq.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.deseq.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] DeSeq2 Analysis for {params.experiment}"
            echo ""

            Rscript -e "rmarkdown::render('{params.rmd}', output_format = '{params.output_format}', output_dir = '{params.out_dir}')" \
                {params.experiment} \
                {params.exp_dir} \
                {params.ref_dir} \
                {params.ref_db_dir} \
                {params.deseq_fc_cutoff} \
                {params.deseq_padj_cutoff} \
                {params.go_padj_cutoff} \
                {params.go_q_cutoff}

            echo ""
            echo "[DONE] DeSeq2 Analysis for {params.experiment}"
            echo ""
        }} 2>&1 | tee {log}
        '''

rule mirtar_off_tgt_freq: 
    input:
        raw_count_tsv = f'{PROJ_DIR}/{{exp_dir}}/5_expression/deseq2/2_txi-counts/{{experiment}}.raw-counts.gene-level.tsv.gz',
        off_tgt_overlap_dfs = lambda wc: expand(
            '{proj_dir}/{exp_dir}/4_variants/hit-id/{experiment}.overlap.{snp}.parquet',
            proj_dir=PROJ_DIR, exp_dir=EXP_DIR_NAMES_DICT[wc.experiment], experiment=wc.experiment, snp=[SNP for SNP in ALL_TGT_SNPS]
        )
    output:
        mirtar_xlsx = f'{PROJ_DIR}/{{exp_dir}}/5_expression/deseq2/3_results/{{experiment}}.degs.mirna-freq.xlsx',
        off_tgt_xlsx = expand(
            '{proj_dir}/{{exp_dir}}/5_expression/deseq2/3_results/{{experiment}}.degs.off-tgt-freq.{snp}.xlsx',
            proj_dir=PROJ_DIR, snp=[SNP for SNP in ALL_TGT_SNPS]
        )
    params:
        py_script = f'{SMK_DIR}/scripts/deg-mirtar-off-tgt-freq.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}',
        var_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}/4_variants',
        ref_dir = REF_DIR,
        ref_db_dir = REF_DB_DIR,
        deseq_fc_cutoff = DESEQ_FC_CUTOFF,
        deseq_padj_cutoff = DESEQ_PADJ_CUTOFF,
        snps = ','.join([SNP for SNP in ALL_TGT_SNPS])
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.deg-mirtar-off-tgt-freq.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] Determining Frequency of miRNA Targets and Off-Targets in DEGs for {params.experiment}"
            echo ""

            python3 {params.py_script} \
                {params.experiment} \
                {params.exp_dir} \
                {params.var_dir} \
                {params.ref_dir} \
                {params.ref_db_dir} \
                {params.deseq_fc_cutoff} \
                {params.deseq_padj_cutoff} \
                {params.snps}

            echo ""
            echo "[DONE] Determining Frequency of miRNA Targets and Off-Targets in DEGs for {params.experiment}"
            echo ""
        }} 2>&1 | tee {log}
        '''