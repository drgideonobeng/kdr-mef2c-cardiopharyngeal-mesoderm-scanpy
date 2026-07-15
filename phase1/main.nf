#!/usr/bin/env nextflow
// =============================================================================
// phase1/main.nf  — Kdr;Mef2cAHF Phase 1: per-sample QC through clustering
// Python/scanpy port. Structure mirrors the R/Seurat phase1/main.nf 1:1:
// meta-map channels, per-step publishDir, sample-prefixed output filenames.
// =============================================================================

nextflow.enable.dsl = 2

def log_params() {
    log.info """
        ╔══════════════════════════════════════════╗
        ║   scRNA-seq — Phase 1 QC (scanpy port)   ║
        ╚══════════════════════════════════════════╝
        samples_yaml : ${params.samples_yaml}
        raw_dir      : ${params.raw_dir}
        outdir       : ${params.outdir}
        ── Step 03 QC thresholds (Global Defaults)─
        min_features : ${params.defaults.min_features}
        max_features : ${params.defaults.max_features}
        min_counts   : ${params.defaults.min_counts}
        max_counts   : ${params.defaults.max_counts}
        max_mt_pct   : ${params.defaults.max_mt_pct}
        ───────────────────────────────────────────
    """.stripIndent()
}

// Helper function to fetch sample-specific QC thresholds, falling back to
// the global default when no per-sample override exists.
def getQCHelper(sample_id, metric) {
    if (params.samples && params.samples.containsKey(sample_id) && params.samples[sample_id].containsKey(metric)) {
        return params.samples[sample_id][metric]
    }
    return params.defaults[metric]
}

// ── processes ─────────────────────────────────────────────────────────────────

process CREATE_ANNDATA_OBJECT {
    tag { "${meta.id}" }
    label 'medium'
    container params.container
    publishDir path: { "${params.outdir}/01_create_anndata_object" }, mode: 'copy'

    input:
    tuple val(meta), path(tenx_dir)

    output:
    tuple val(meta), path("${meta.id}_01_unfiltered.h5ad"), emit: adata_unfiltered

    script:
    """
    python ${projectDir}/bin/01_create_anndata_object.py \\
        --tenx_dir  ${tenx_dir} \\
        --sample_id ${meta.id} \\
        --genotype  ${meta.genotype} \\
        --out       ${meta.id}_01_unfiltered.h5ad
    """
}

process VISUALIZE_QC {
    tag { "${meta.id}" }
    label 'medium'
    container params.container
    publishDir path: { "${params.outdir}/02_visualize_qc" }, mode: 'copy'

    input:
    tuple val(meta), path(adata_in), val(min_feat), val(max_feat), val(min_cnt), val(max_cnt), val(max_mt)

    output:
    tuple val(meta), path("${meta.id}_02_with_qc.h5ad"),  emit: adata_with_qc
    tuple val(meta), path("${meta.id}_02_qc_plots.pdf"),  emit: qc_plots

    script:
    """
    python ${projectDir}/bin/02_visualize_qc.py \\
        --in_h5ad        ${adata_in} \\
        --sample_id      ${meta.id} \\
        --min_features   ${min_feat} \\
        --max_features   ${max_feat} \\
        --min_counts     ${min_cnt} \\
        --max_counts     ${max_cnt} \\
        --max_percent_mt ${max_mt} \\
        --out_h5ad       ${meta.id}_02_with_qc.h5ad \\
        --out_pdf        ${meta.id}_02_qc_plots.pdf
    """
}

process FILTER_CELLS {
    tag { "${meta.id}" }
    label 'medium'
    container params.container
    publishDir path: { "${params.outdir}/03_filter_cells" }, mode: 'copy'

    input:
    tuple val(meta), path(adata_in), val(min_feat), val(max_feat), val(min_cnt), val(max_cnt), val(max_mt)

    output:
    tuple val(meta), path("${meta.id}_03_filtered.h5ad"), emit: adata_filtered

    script:
    """
    python ${projectDir}/bin/03_filter_cells.py \\
        --in_h5ad        ${adata_in} \\
        --min_features   ${min_feat} \\
        --max_features   ${max_feat} \\
        --min_counts     ${min_cnt} \\
        --max_counts     ${max_cnt} \\
        --max_percent_mt ${max_mt} \\
        --out            ${meta.id}_03_filtered.h5ad
    """
}

process DETECT_DOUBLETS {
    tag { "${meta.id}" }
    label 'high'
    container params.container
    publishDir path: { "${params.outdir}/04_detect_doublets" }, mode: 'copy'

    input:
    tuple val(meta), path(adata_in)

    output:
    tuple val(meta), path("${meta.id}_04_singlets.h5ad"),        emit: adata_singlets
    tuple val(meta), path("${meta.id}_04_doublet_summary.csv"),  emit: doublet_summary
    tuple val(meta), path("${meta.id}_04_doublet_plots.pdf"),    emit: doublet_plots

    script:
    """
    python ${projectDir}/bin/04_detect_doublets.py \\
        --in_h5ad     ${adata_in} \\
        --sample_id   ${meta.id} \\
        --seed        ${params.seed} \\
        --out_h5ad    ${meta.id}_04_singlets.h5ad \\
        --out_summary ${meta.id}_04_doublet_summary.csv \\
        --out_pdf     ${meta.id}_04_doublet_plots.pdf
    """
}

process NORMALIZE_DATA {
    tag { "${meta.id}" }
    label 'high'
    container params.container
    publishDir path: { "${params.outdir}/05_normalize_data" }, mode: 'copy'

    input:
    tuple val(meta), path(adata_in)

    output:
    tuple val(meta), path("${meta.id}_05_normalized.h5ad"), emit: adata_normalized

    script:
    """
    python ${projectDir}/bin/05_normalize_data.py \\
        --in_h5ad ${adata_in} \\
        --n_hvgs  ${params.n_hvgs} \\
        --target_sum ${params.target_sum} \\
        --seed    ${params.seed} \\
        --out     ${meta.id}_05_normalized.h5ad
    """
}

process DIM_REDUCTION_AND_CLUSTER {
    tag { "${meta.id}" }
    label 'high'
    container params.container
    publishDir path: { "${params.outdir}/06_dim_reduction_and_cluster" }, mode: 'copy'

    input:
    tuple val(meta), path(adata_in)

    output:
    tuple val(meta), path("${meta.id}_06_phase1.h5ad"), emit: adata_phase1
    tuple val(meta), path("${meta.id}_06_elbow.pdf"),   emit: elbow
    tuple val(meta), path("${meta.id}_06_umap.pdf"),    emit: umap

    script:
    """
    python ${projectDir}/bin/06_dim_reduction_and_cluster.py \\
        --in_h5ad     ${adata_in} \\
        --n_pcs       ${params.n_pcs} \\
        --n_neighbors ${params.n_neighbors} \\
        --cluster_res ${params.cluster_res} \\
        --seed        ${params.seed} \\
        --out_rds     ${meta.id}_06_phase1.h5ad \\
        --out_elbow   ${meta.id}_06_elbow.pdf \\
        --out_umap    ${meta.id}_06_umap.pdf
    """
}

// ── workflow ──────────────────────────────────────────────────────────────────
workflow {
    log_params()

    def yamlSlurper = new groovy.yaml.YamlSlurper()
    def parsedYaml  = yamlSlurper.parse(file(params.samples_yaml))

    // Safe structure evaluation (handles both Maps and flat Lists)
    def sampleList = (parsedYaml instanceof Map && parsedYaml.containsKey('samples')) ?
                     parsedYaml.samples : parsedYaml

    ch_samples = Channel.fromList(sampleList)
        .map { s ->
            def meta = [id: s.id, genotype: s.genotype]

            // Explicit path (relative to raw_dir) or fall back to the
            // id-based convention "${raw_dir}/${id}/filtered_feature_bc_matrix".
            def tenx_path = s.containsKey('path') ?
                "${params.raw_dir}/${s.path}" :
                "${params.raw_dir}/${s.id}/filtered_feature_bc_matrix"

            [ meta, file(tenx_path) ]
        }

    CREATE_ANNDATA_OBJECT(ch_samples)

    ch_qc_inputs = CREATE_ANNDATA_OBJECT.out.adata_unfiltered.map { meta, h5ad ->
        [ meta, h5ad, getQCHelper(meta.id, 'min_features'), getQCHelper(meta.id, 'max_features'), getQCHelper(meta.id, 'min_counts'), getQCHelper(meta.id, 'max_counts'), getQCHelper(meta.id, 'max_mt_pct') ]
    }
    VISUALIZE_QC(ch_qc_inputs)

    ch_filter_inputs = VISUALIZE_QC.out.adata_with_qc.map { meta, h5ad ->
        [ meta, h5ad, getQCHelper(meta.id, 'min_features'), getQCHelper(meta.id, 'max_features'), getQCHelper(meta.id, 'min_counts'), getQCHelper(meta.id, 'max_counts'), getQCHelper(meta.id, 'max_mt_pct') ]
    }
    FILTER_CELLS(ch_filter_inputs)

    DETECT_DOUBLETS(FILTER_CELLS.out.adata_filtered)
    NORMALIZE_DATA(DETECT_DOUBLETS.out.adata_singlets)
    DIM_REDUCTION_AND_CLUSTER(NORMALIZE_DATA.out.adata_normalized)
}
