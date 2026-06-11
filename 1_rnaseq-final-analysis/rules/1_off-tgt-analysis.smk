rule annotate_var_table: 
    input:
        var_tsvgz = f'{PROJ_DIR}/{{exp_dir}}/4_variants/{{experiment}}.variants.tsv.gz'
    output:
        annotated_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/{{experiment}}.variants.annotated.parquet'
    params:
        py_script = f'{SMK_DIR}/scripts/annotate-var-table.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}',
        ref_dir = REF_DIR,
        tgt_snps = ','.join(ALL_TGT_SNPS),
        tgt_edits = lambda wc: ';'.join(edit for edit in TGT_EDIT_DICT[wc.experiment])
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.annotate.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] Variant Annotation for {params.experiment}"
            echo ""

            python3 {params.py_script} \
                {params.experiment} \
                {params.exp_dir} \
                {params.ref_dir} \
                {params.tgt_snps} \
                "{params.tgt_edits}"

            echo ""
            echo "[DONE] Variant Annotation for {params.experiment}"
            echo ""
        }} 2>&1 | tee {log}
        '''

rule id_off_tgt_hits:
    input:
        annotated_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/{{experiment}}.variants.annotated.parquet'
    output:
        overlap_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.overlap.{{snp}}.parquet',
        overlap_excel = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.overlap.{{snp}}.xlsx'
    params:
        py_script = f'{SMK_DIR}/scripts/off-tgt-hit-id.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}',
        var_dir = lambda wc: f'{PROJ_DIR}/{wc.exp_dir}/4_variants',
        tgt_snp = lambda wc: wc.snp,
        mutect_mode = 'False'
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}off-tgt-hits-{{snp}}.log'
    shell:
        '''
        {{
        echo ""
        echo "[RUN] {params.tgt_snp} Off-Target Hit Filtering for {params.experiment}"
        echo ""

        python3 {params.py_script} \
            {params.experiment} \
            {params.exp_dir} \
            {params.var_dir} \
            {params.tgt_snp} \
            {params.mutect_mode}

        echo ""
        echo "[DONE] {params.tgt_snp} Off-Target Hit Filtering for {params.experiment}"
        echo ""
        }} 2>&1 | tee {log}
        '''

rule align_off_tgt_hits:
    input:
        overlap_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.overlap.{{snp}}.parquet',
        transcripts_raw_tsv = f'{PROJ_DIR}/{{exp_dir}}/5_expression/salmon/analysis/salmon-transcript-raw-counts.tsv.gz',
        transcript_lookup_pkl = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align/transcript-id-to-contig.pkl.gz'
    output:
        align_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align/{{experiment}}.transcripts.{{snp}}.parquet'
    params:
        py_script = f'{SMK_DIR}/scripts/off-tgt-align-to-transc.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = f'{PROJ_DIR}/{{exp_dir}}',
        var_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants',
        ref_dir = REF_DIR,
        tgt_snp = lambda wc: wc.snp,
        align_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align',
        out_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align/{{snp}}'
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.align-off-tgts-{{snp}}.log'
    shell:
        '''
        {{
        echo ""
        echo "[RUN] {params.tgt_snp} Off-Target Hits Align to Transcripts for {params.experiment}"
        echo ""

        mkdir -p {params.align_dir}
        mkdir -p {params.out_dir}

        python3 {params.py_script} \
            {params.experiment} \
            {params.exp_dir} \
            {params.var_dir} \
            {params.ref_dir} \
            {params.tgt_snp}

        echo ""
        echo "[DONE] {params.tgt_snp} Off-Target Hits Align to Transcripts for {params.experiment}"
        echo ""
        }} 2>&1 | tee {log}
        '''

rule off_tgt_seq_consensus:
    input:
        align_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align/{{experiment}}.transcripts.{{snp}}.parquet'
    output:
        consensus_xlsx = f'{PROJ_DIR}/{{exp_dir}}/4_variants/{{experiment}}.transcripts.{{snp}}.o{{seq_consensus_offset}}-consensus.xlsx'
    params:
        py_script = f'{SMK_DIR}/scripts/off-tgt-seq-consensus.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = f'{PROJ_DIR}/{{exp_dir}}',
        ref_dir = REF_DIR,
        tgt_snp = lambda wc: wc.snp,
        seq_consensus_offset = lambda wc: wc.seq_consensus_offset,
        seq_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants/seq-context',
        out_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants/seq-context/{{snp}}-o{{seq_consensus_offset}}'
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.seq-consensus-{{snp}}-o{{seq_consensus_offset}}.log'
    shell:
        '''
        {{
        echo ""
        echo "[RUN] {params.tgt_snp} Off-Target Sequence Context Analysis (Offset={params.seq_consensus_offset}) for {params.experiment}"
        echo ""

        mkdir -p {params.out_dir}

        python3 {params.py_script} \
            {params.experiment} \
            {params.exp_dir} \
            {params.ref_dir} \
            {params.tgt_snp} \
            {params.seq_consensus_offset}

        echo ""
        echo "[DONE] {params.tgt_snp} Off-Target Sequence Context Analysis (Offset={params.seq_consensus_offset}) for {params.experiment}"
        echo ""
        }} 2>&1 | tee {log}
        '''

rule off_tgt_mirtarbase_freq:
    input:
        overlap_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.overlap.{{snp}}.parquet',
        gene_raw_tsv = f'{PROJ_DIR}/{{exp_dir}}/5_expression/salmon/analysis/salmon-gene-raw-counts.tsv.gz',
    output:
        out_xlsx = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.miRTar-freq.{{snp}}.xlsx',
        # single col with entrez_id_data for overlap_df, retains index for merging
        overlap_entrez_pq = f'{PROJ_DIR}/{{exp_dir}}/4_variants/hit-id/{{experiment}}.overlap.{{snp}}.entrez.parquet',
    params:
        py_script = f'{SMK_DIR}/scripts/off-tgt-miRTar-freq.py',
        experiment = lambda wc: wc.experiment,
        exp_dir = f'{PROJ_DIR}/{{exp_dir}}',
        var_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants',
        ref_db_dir = REF_DB_DIR,
        tgt_snp = lambda wc: wc.snp
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/{{experiment}}.mirtar-freq-{{snp}}.log'
    shell:
        '''
        {{
        echo ""
        echo "[RUN] {params.tgt_snp} Off-Target Gene Frequency in mirTarBase {params.experiment}"
        echo ""

        python3 {params.py_script} \
            {params.experiment} \
            {params.exp_dir} \
            {params.var_dir} \
            {params.ref_db_dir} \
            {params.tgt_snp}

        echo ""
        echo "[DONE] {params.tgt_snp} Off-Target Gene Frequency in mirTarBase {params.experiment}"
        echo ""
        }} 2>&1 | tee {log}
        '''

rule off_tgt_align_mirna:
    input:
        align_pqs = [
            f'{PROJ_DIR}/{EXP_DIR_NAMES_DICT[experiment]}/4_variants/transcript-align/{experiment}.transcripts.{{snp}}.parquet'
            for experiment in EXPERIMENT_NAMES
        ],
        mirbase_pq = f'{REF_DB_DIR}/miRBase/mirbase.parquet',
        mirbase_patterns = f'{REF_DB_DIR}/miRBase/miRBase.seed-patterns.pkl.gz',
        mirtar_tgt_sites_pq = f'{REF_DB_DIR}/miRtarBase/miRTarBase.tgt-sites.hsa.parquet',
        mirtar_patterns = f'{REF_DB_DIR}/miRtarBase/miRTarBase.tgt-sites.patterns.pkl.gz',
        transcripts_expressed_pq = f'{RESOURCE_DIR}/expr-transcripts/{'.'.join(EXPERIMENT_NAMES)}.transcripts.expressed.parquet',
        transcripts_random_pq = f'{RESOURCE_DIR}/expr-transcripts/{'.'.join(EXPERIMENT_NAMES)}.transcripts.expressed.random.parquet'
    output:
        out_xlsx_list = f'{AGG_RESULTS_DIR}/{'.'.join(EXPERIMENT_NAMES)}.mirna-align.{{snp}}.o{{mirna_offset}}.xlsx',
        positional_xlsx_list = f'{AGG_RESULTS_DIR}/{'.'.join(EXPERIMENT_NAMES)}.mirna-align.{{snp}}.o{{mirna_offset}}.pos.xlsx',
        mirna_align_pq_list = f'{RESOURCE_DIR}/transcript-overlap/{'.'.join(EXPERIMENT_NAMES)}.transcripts.{{snp}}.mirna.o{{mirna_offset}}.parquet',
        
    params:
        py_script = f'{SMK_DIR}/scripts/off-tgt-mirna-align.py',
        experiment_csv = ','.join(EXPERIMENT_NAMES),
        proj_dir = PROJ_DIR,
        exp_dir_csv = ','.join([os.path.join(PROJ_DIR, EXP_DIR_NAMES_DICT[e]) for e in EXPERIMENT_NAMES]),
        ref_dir = REF_DIR,
        ref_db_dir = REF_DB_DIR,
        tgt_snp = lambda wc: wc.snp,
        mirna_offset = lambda wc: wc.mirna_offset,
        resource_dir = RESOURCE_DIR,
        agg_res_dir = AGG_RESULTS_DIR
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{'.'.join(EXPERIMENT_NAMES)}.mirna-align-{{snp}}-o{{mirna_offset}}.log'
    shell:
        '''
        {{
        echo ""
        echo "[RUN] Align miRNAs to {params.tgt_snp} Off-Target Seq Context (Offset = {params.mirna_offset})"
        echo ""

        python3 {params.py_script} \
            {params.experiment_csv} \
            {params.proj_dir} \
            {params.exp_dir_csv} \
            {params.ref_dir} \
            {params.ref_db_dir} \
            {params.tgt_snp} \
            {params.mirna_offset} \
            {params.resource_dir} \
            {params.agg_res_dir}


        echo ""
        echo "[DONE] Align miRNAs to {params.tgt_snp} Off-Target Seq Context (Offset = {params.mirna_offset})"
        echo ""
        }} 2>&1 | tee {log}
        '''
