nextflow.enable.dsl=2

// Rewriter for C-Origami publishDir:
// turns "origami_chrN/<sample>/prediction/..." -> "prediction/..."
//def corigamiPublish = { name ->
//  def p = name.toString()
//  if( !p.startsWith('origami_') ) return null
//  def parts = p.tokenize('/')     // faster than split('/') and avoids regex
//  parts.size() >= 3 ? parts.drop(2).join('/') : null
//}

// ---------- safe fallbacks (CLI still overrides) ----------
params.input_pattern = params.input_pattern ?: './Input_BAMs/*.bam'
params.outdir        = params.outdir        ?: 'results'
params.chroms        = params.chroms        ?: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19'
def chroms_csv = (params.chroms as String)

// -----------------------------
// Processes (bodies same as yours, with cpus fixes)
// -----------------------------

process BamChrAppend {
  publishDir "${params.outdir}/${sample}/bam_prep", mode: 'copy'
  cpus 2

  input:
    tuple val(sample), path(bam)

  output:
    tuple val(sample), path("${sample}.chr.bam"), emit: prepped_bam_ch

  script:
  """
  set -euo pipefail
  module load ${params.moduleLoad}
  bash "${params.BAM_prep_script}" \
    -i "${bam}" \
    -o "${sample}.chr.bam" \
    -t ${task.cpus}

  [ -f "${sample}.chr.bam.bai" ] || samtools index -@ ${task.cpus} "${sample}.chr.bam"
  """
}

process MaxATACPrepare {
  publishDir "${params.outdir}/${sample}/maxatac", mode: 'copy'
  cpus 8

  input:
    tuple val(sample), path(chr_bam)

  output:
    tuple val(sample), path("prepare"), emit: prep_ch

  script:
  """
  set -euo pipefail
  mkdir -p prepare

  BAM_ABS=\$(readlink -f "${chr_bam}")
  OUT_DIR="\$(pwd)/prepare"

  mamba activate ${params.maxATAC_mamba_environment}
  module load ${params.moduleLoad}
  bash ${params.maxATAC_prepare_script} \
    -s ${sample} \
    -b "\$BAM_ABS" \
    -o "\$OUT_DIR" \
    -t ${task.cpus} \
    -C "${chroms_csv}"
  touch prepare/.done_prepare
  """
}

process MaxATACPredictChrom {
  publishDir "${params.outdir}/${sample}/maxatac/predict/${chrom}", mode: 'copy'
  cpus 4

  input:
    tuple val(sample), path(prep_dir)
    each chrom

  output:
    tuple val(sample), val(chrom), path(prep_dir), path("predict_${chrom}"), emit: pred_ch

  script:
  """
  set -euo pipefail
  mkdir -p predict_${chrom}
  
  PREP_ABS=\$(readlink -f "${prep_dir}")
  OUT_DIR="\$(pwd)/predict_${chrom}"

  mamba activate ${params.maxATAC_mamba_environment}
  module load ${params.moduleLoad}
  bash ${params.maxATAC_predict_script} \
    -s ${sample} \
    -c ${chrom} \
    -i "\$PREP_ABS" \
    -o "\$OUT_DIR" \
    -t ${task.cpus} \
    -d ${params.maxatac_data_dir}
  touch predict_${chrom}/.done_predict
  """
}

process MaxATACBenchmark {
  publishDir "${params.outdir}/${sample}/maxatac", mode: 'copy'
  cpus 2

  input:
    tuple val(sample), path(prep_dir), val(pred_dirs)    // <-- add prep_dir here

  output:
    tuple val(sample), path("benchmark"), emit: bench_ch

  script:
  def atac_name   = params.bench_atac_name ?: "${sample}_IS_slop20_RP20M_minmax01.bw"
  def tf_name     = params.benchmark_tf     ?: 'CTCF'
  def chroms_b    = params.benchmark_chroms ?: chroms_csv
  def threshold   = params.atac_threshold   ?: 0.30
  def bl_bed      = params.mm10_blacklist_bed
  def chrom_sizes = params.corigami_chrom_sizes
  """
  set -euo pipefail
  mkdir -p benchmark

  # Make prep dir absolute and build ATAC path from it
  PREP_ABS=\$(readlink -f "${prep_dir}")
  ATAC="\$PREP_ABS/${atac_name}"

  # Write one prediction directory per line (use Groovy-expanded list)
  cat > benchmark/${sample}.pred_dirs.txt <<'LIST'
${pred_dirs.collect{ it.toString() }.join('\n')}
LIST
  listfile="benchmark/${sample}.pred_dirs.txt"

  [[ -f "\$ATAC" ]] || { echo "[error] ATAC not found: \$ATAC"; exit 1; }

  mamba activate ${params.maxATAC_mamba_environment}
  module load ${params.moduleLoad}
  bash ${params.maxATAC_sampled_pred_benchmark_script} \
    -s ${sample} \
    -l "\$listfile" \
    -o benchmark \
    -a "\$ATAC" \
    -b ${bl_bed} \
    -c ${chrom_sizes} \
    -f ${tf_name} \
    -C "${chroms_b}" \
    -t ${threshold}

  touch benchmark/.done_benchmark
  """
}

process COrigamiPredictChrom {
  publishDir "${params.outdir}/${sample}/corigami", mode: 'copy'
//  publishDir "${params.outdir}/${sample}/corigami/${chrom}", mode: 'copy', pattern: "origami_${chrom}/**", saveAs: corigamiPublish 
  cpus 2

  input:
    tuple val(sample), val(chrom), path(prep_dir), path(maxatac_pred_dir)

  output:
    tuple val(sample), val(chrom), path("origami_${chrom}"), emit: origami_ch

  script:
  def atac_name = params.corigami_atac_name ?: (params.corigami_atac_suffix ? "${sample}${params.corigami_atac_suffix}" : "${sample}_IS_slop20_RP20M_minmax01.bw")
  def ctcf_name = params.corigami_ctcf_name ?: (params.corigami_ctcf_suffix ? "${sample}${params.corigami_ctcf_suffix}" : "${sample}_CTCF.bw")
  def cs        = params.corigami_chrom_sizes
  def seq_dir   = params.corigami_seq_dir
  def model     = params.corigami_model_ckpt
  def win       = params.corigami_win   ?: 2097152
  def step      = params.corigami_step  ?: win
  """
  set -euo pipefail
  mkdir -p origami_${chrom}

  PREP_ABS=\$(readlink -f "${prep_dir}")
  PRED_ABS=\$(readlink -f "${maxatac_pred_dir}")

  ATAC="\$PREP_ABS/${atac_name}"
  CTCF="\$PRED_ABS/${ctcf_name}"

  echo "[corigami] ATAC=\$ATAC"
  echo "[corigami] CTCF=\$CTCF"

  [[ -f "\$ATAC" ]] || { echo "[error] ATAC not found: \$ATAC"; exit 1; }
  [[ -f "\$CTCF" ]] || { echo "[error] CTCF not found: \$CTCF"; exit 1; }

  export SAMPLE="${sample}"
  export OUT="origami_${chrom}"
  export CS="${cs}"
  export SEQ_DIR="${seq_dir}"
  export MODEL="${model}"
  export ATAC="\$ATAC"
  export CTCF="\$CTCF"
  export WIN="${win}"
  export STEP="${step}"

  mamba activate ${params.Corigami_mamba_environment}
  module load ${params.moduleLoad}
  bash ${params.C_Origami_predict_script} ${chrom}

  touch origami_${chrom}/.done_corigami
  """
}

// -----------------------------
// Workflow (single unnamed entrypoint)
// -----------------------------
workflow {
  // Inputs -> (sample, bam)
  def raw_bams = Channel.fromPath(params.input_pattern, checkIfExists: true)
  def bam_only = raw_bams.filter { it.name.endsWith('.bam') && !it.name.endsWith('.bam.bai') }
  def bam_ch   = bam_only.map { bam -> tuple(bam.baseName, bam) }

  // 1) chr-append
  BamChrAppend(bam_ch)
  def prepped_bam_ch = BamChrAppend.out.prepped_bam_ch

  // 2) maxATAC prepare
  MaxATACPrepare(prepped_bam_ch)
  def prep_ch = MaxATACPrepare.out.prep_ch

  // 3) per-chrom maxATAC predict
  def chroms_list = (params.chroms as String).split(',').collect { it.trim() }.findAll { it }
  MaxATACPredictChrom(prep_ch, chroms_list)
  def pred_ch = MaxATACPredictChrom.out.pred_ch

  // Debug: peek at what the channel carries
  pred_ch.view { s, c, prep, pdir -> "pred_ch: s=${s} c=${c} prep=${prep} pred=${pdir}" }

  // 4) optional benchmark
  if (params.run_benchmark) {
    def pred_by_sample = pred_ch
        .map { s, c, prep, dir -> tuple(s, prep, dir) }  // (s, prep, pred)
        .groupTuple()                                     // (s, [prep...], [pred...])
        .map { s, preps, preds ->
            // sanity: all prep paths for a sample should be the same
            def uniq = preps.unique()
            if( uniq.size() != 1 ) {
              log.warn "Multiple prep dirs for sample ${s}; using first: ${uniq}"
            }
            tuple(s, uniq[0], preds)                      // (s, prep_dir, [predict_dirs])
        }
  
    // Optional: peek once to verify
    // pred_by_sample.take(5).view { s, prep, preds -> "bench_in: ${s}\n  prep=${prep}\n  n_preds=${preds.size()}" }
  
    MaxATACBenchmark(pred_by_sample)
  } 

  // 5) C-Origami per chrom
  COrigamiPredictChrom(pred_ch)
}
