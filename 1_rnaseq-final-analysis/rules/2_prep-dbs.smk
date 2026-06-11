rule prep_mirtarbase_tgt_sites: 
    input:
        mirtar_tgt_sites_csv = f'{REF_DB_DIR}/miRtarBase/MicroRNA_Target_Sites.csv.gz'
    output:
        mirtar_tgt_sites_pq = f'{REF_DB_DIR}/miRtarBase/miRTarBase.tgt-sites.hsa.parquet',
        mirtar_patterns = f'{REF_DB_DIR}/miRtarBase/miRTarBase.tgt-sites.patterns.pkl.gz'
    params:
        py_script = f'{SMK_DIR}/scripts/prep-mirtarbase-tgt-site.py',
        ref_db_dir = REF_DB_DIR
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/prep-mirtarbase-tgt-site.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] miRTarBase Target Site Data Prep"
            echo ""

            python3 {params.py_script} \
                {params.ref_db_dir}

            echo ""
            echo "[DONE] miRTarBase Target Site Data Prep"
            echo ""
        }} 2>&1 | tee {log}
        '''

rule prep_mirbase: 
    input:
        mirbase_dat = f'{REF_DB_DIR}/miRBase/miRNA.dat.gz'
    output:
        mirbase_pq = f'{REF_DB_DIR}/miRBase/mirbase.parquet',
        mirbase_patterns = f'{REF_DB_DIR}/miRBase/miRBase.seed-patterns.pkl.gz'
    params:
        py_script = f'{SMK_DIR}/scripts/prep-mirbase.py',
        ref_db_dir = REF_DB_DIR
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/prep-mirbase.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] miRBase Data Prep"
            echo ""

            python3 {params.py_script} \
                {params.ref_db_dir}

            echo ""
            echo "[DONE] miRBase Data Prep"
            echo ""
        }} 2>&1 | tee {log}
        '''

rule prep_expressed_transcripts: 
    input:
        salmon_raw_transcripts = expand(
            '{proj_dir}/{exp_dir}/5_expression/salmon/analysis/salmon-transcript-raw-counts.tsv.gz',
            proj_dir=PROJ_DIR, exp_dir=list(EXP_DIR_NAMES_DICT.values())
        ),
        align_pqs = [
            f'{PROJ_DIR}/{EXP_DIR_NAMES_DICT[experiment]}/4_variants/transcript-align/{experiment}.transcripts.{snp}.parquet'
            for experiment in EXPERIMENT_NAMES
            for snp in ALL_TGT_SNPS
        ]
    output:
        transcripts_expressed_pq = f'{RESOURCE_DIR}/expr-transcripts/{'.'.join(EXPERIMENT_NAMES)}.transcripts.expressed.parquet',
        transcripts_random_pq = f'{RESOURCE_DIR}/expr-transcripts/{'.'.join(EXPERIMENT_NAMES)}.transcripts.expressed.random.parquet'
    params:
        py_script = f'{SMK_DIR}/scripts/prep-expressed-transcripts.py',
        resource_dir = RESOURCE_DIR,
        experiment_csv = ','.join(EXPERIMENT_NAMES),
        proj_dir = PROJ_DIR,
        exp_dir_csv = ','.join([os.path.join(PROJ_DIR, EXP_DIR_NAMES_DICT[e]) for e in EXPERIMENT_NAMES]),
        ref_dir = REF_DIR,
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/prep-expr-trans.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] Expressed Transcripts Background Data Prep"
            echo ""

            python3 {params.py_script} \
                {params.experiment_csv} \
                {params.proj_dir} \
                {params.exp_dir_csv} \
                {params.ref_dir} \
                {params.resource_dir}

            echo ""
            echo "[DONE] Expressed Transcripts Background Data Prep"
            echo ""
        }} 2>&1 | tee {log}
        '''

# could change to single common resource location
rule prep_transcript_contig_lookup: 
    input:
        transcripts_fa = f'{REF_DIR}/{REF_GENOME_NAME}.transcripts.filtered.fa.gz'
    output:
        transcript_lookup_pkl = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align/transcript-id-to-contig.pkl.gz'
    params:
        py_script = f'{SMK_DIR}/scripts/prep-transcript-contig-lookup.py',
        align_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants/transcript-align',
        exp_dir = f'{PROJ_DIR}/{{exp_dir}}',
        var_dir = f'{PROJ_DIR}/{{exp_dir}}/4_variants',
        ref_dir = REF_DIR
    conda:
        f'{SMK_DIR}/env/mrnaseq_py.yaml'
    log:
        f'logs/{{exp_dir}}/prep-transcript-lookup.log'
    shell:
        '''
        {{
            echo ""
            echo "[RUN] Preparing Transcript Contig Lookup"
            echo ""

            mkdir -p {params.align_dir}

            python3 {params.py_script} \
                {params.exp_dir} \
                {params.var_dir} \
                {params.ref_dir}

            echo ""
            echo "[DONE] Prepared Transcript Contig Lookup"
            echo ""
        }} 2>&1 | tee {log}
        '''