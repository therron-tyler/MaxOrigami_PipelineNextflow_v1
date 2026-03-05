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

process MergeGroupBams {
  publishDir "${params.outdir}/${group}/bam_group", mode: 'copy'
  cpus 8

  input:
    tuple val(group), path(bams)

  output:
    tuple val(group), path("${group}.merged.bam"), emit: merged_group_bam_ch

  script:
  """
  set -euo pipefail
  module load ${params.moduleLoad}

  # Merge then sort+index
  samtools merge -@ ${task.cpus} -f ${group}.merged.unsorted.bam ${bams.join(' ')}
  samtools sort  -@ ${task.cpus} -o ${group}.merged.bam ${group}.merged.unsorted.bam
  rm -f ${group}.merged.unsorted.bam
  samtools index -@ ${task.cpus} ${group}.merged.bam
  """
}

process CountUsableReads {
  cpus 1

  input:
    tuple val(group), path(bam)

  output:
    tuple val(group), path(bam), stdout, emit: group_reads_ch

  script:
  """
  set -euo pipefail
  module load ${params.moduleLoad}

  samtools view -c -q ${params.ds_mapq} -F ${params.ds_exclude_flags} "${bam}"
  """
}

process DownsampleGroupBam {
  publishDir "${params.outdir}/${group}/bam_group", mode: 'copy'
  cpus 8

  input:
    tuple val(group), path(bam), val(n_reads), val(target_reads)

  output:
    tuple val(group), path("${group}.group_ds.bam"), emit: group_ds_bam_ch

  script:
  """
  set -euo pipefail
  module load ${params.moduleLoad}
  
  echo "[ds] group=${group} n_reads=${n_reads} target=${target_reads}"
  
  if [[ ${n_reads} -le ${target_reads} ]]; then
    cp "${bam}" "${group}.group_ds.bam"
  else
    prob=\$(python3 - <<PY
  n=${n_reads}
  t=${target_reads}
  print(t/float(n))
  PY
  )
  
    sarg=\$(python3 - <<PY
  seed=${params.ds_seed}
  p=float("${'$'}prob")
  print(f"{seed}.{int(p*1e6):06d}")
  PY
  )
  
    samtools view -@ ${task.cpus} -b -s "\$sarg" "${bam}" | \
      samtools sort -@ ${task.cpus} -o "${group}.group_ds.bam" -
  fi
  
  samtools index -@ ${task.cpus} "${group}.group_ds.bam"

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

process HiCTileVizCompare {
  publishDir "${params.outdir}/${params.hic_viz_outdir}/${sample_a}_vs_${sample_b}/${chrom}", mode: 'copy'
  cpus 16
  memory '30 GB'
  time '02:00:00'

  input:
    tuple val(sample_a), val(sample_b), val(chrom),
          path(origami_a), path(prep_a), path(pred_a),
          path(origami_b), path(prep_b), path(pred_b)

  output:
    path("viz_out"), emit: viz_out

  script:
  """
  set -euo pipefail
  module load ${params.moduleLoad}

  NPY_A="${origami_a}/npy"
  NPY_B="${origami_b}/npy"

  ATAC_A="${prep_a}/${sample_a}${params.hic_viz_atac_suffix_a}"
  ATAC_B="${prep_b}/${sample_b}${params.hic_viz_atac_suffix_b}"

  CTCF_A="${pred_a}/${sample_a}${params.hic_viz_ctcf_suffix}"
  CTCF_B="${pred_b}/${sample_b}${params.hic_viz_ctcf_suffix}"

  # Optional override for B (e.g. quantile matched)
  if [[ -n "${params.hic_viz_ctcf_override_b ?: ''}" ]]; then
    CTCF_B="${params.hic_viz_ctcf_override_b}"
  fi

  echo "[viz] NPY_A=\$NPY_A"
  echo "[viz] NPY_B=\$NPY_B"
  echo "[viz] ATAC_A=\$ATAC_A"
  echo "[viz] ATAC_B=\$ATAC_B"
  echo "[viz] CTCF_A=\$CTCF_A"
  echo "[viz] CTCF_B=\$CTCF_B"

  [[ -d "\$NPY_A" ]] || { echo "[error] missing \$NPY_A"; exit 1; }
  [[ -d "\$NPY_B" ]] || { echo "[error] missing \$NPY_B"; exit 1; }
  [[ -f "\$ATAC_A" ]] || { echo "[error] missing \$ATAC_A"; exit 1; }
  [[ -f "\$ATAC_B" ]] || { echo "[error] missing \$ATAC_B"; exit 1; }
  [[ -f "\$CTCF_A" ]] || { echo "[error] missing \$CTCF_A"; exit 1; }
  [[ -f "\$CTCF_B" ]] || { echo "[error] missing \$CTCF_B"; exit 1; }

  mkdir -p viz_out

  python3 ${params.hic_viz_script} \\
    --input "\$NPY_A" \\
    --input-b "\$NPY_B" \\
    --pattern "${chrom}_*.npy" \\
    --tile-width-bp 2097152 \\
    --cmap ${params.hic_viz_cmap} \\
    --label-a ${sample_a} --label-b ${sample_b} \\
    --ctcf-bw "\$CTCF_A" --atac-bw "\$ATAC_A" \\
    --ctcf-bw-b "\$CTCF_B" --atac-bw-b "\$ATAC_B" \\
    --tick-mb ${params.hic_viz_tick_mb} \\
    --workers ${task.cpus} --colorbar \\
    --output viz_out
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

    // --- read groupfile: Sample is filename like Ly6CHi_Old_WT1.chr.bam
  def group_map_ch = Channel
    .fromPath(params.groupfile)
    .splitCsv(header:true)
    .map { row ->
      def s = row.Sample.toString().replaceAll(/\.bam$/, '')  // -> Ly6CHi_Old_WT1.chr
      tuple(s, row.Group.toString())
    }

  // attach group to each BAM
  def bam_with_group = prepped_bam_ch
    .join(group_map_ch, by: 0)                      // (sample, bam, group)
    .map { sample, bam, group -> tuple(group, bam) } // (group, bam)

  // group -> [bam1,bam2,bam3,bam4]
  def grouped_bams = bam_with_group.groupTuple()

  // merge replicates
  MergeGroupBams(grouped_bams)
  def merged_group_bam_ch = MergeGroupBams.out.merged_group_bam_ch  // (group, group.merged.bam)

  // count usable reads per group bam
  CountUsableReads(merged_group_bam_ch)
  def group_reads_ch = CountUsableReads.out.group_reads_ch          // (group, bam, reads)

  // choose target: auto = min across groups unless params.group_downsample_reads set
  def min_group_reads_ch = group_reads_ch
    .map { g, bam, n -> n as long }
    .reduce { a, b -> Math.min(a as long, b as long) }

  def target_reads_ch = min_group_reads_ch.map { minN ->
    def tgt = (params.group_downsample_reads as long) > 0 ? Math.min(minN, params.group_downsample_reads as long) : minN
    log.info "[group-ds] min usable reads across groups = ${minN}; target = ${tgt}"
    return tgt
  }

  // attach target to each group bam

  def ds_in_ch = group_reads_ch
    .map { g, bam, n ->
      tuple(g, bam, (n.toString().trim() as long))
    }
    .combine(target_reads_ch)
    .map { g, bam, n_reads, tgt ->
      tuple(g, bam, (n_reads as long), (tgt as long))
    }

  DownsampleGroupBam(ds_in_ch)
  def group_ds_bam_ch = DownsampleGroupBam.out.group_ds_bam_ch        // (group, group.group_ds.bam)

  // 2) maxATAC prepare
  MaxATACPrepare(group_ds_bam_ch)
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
  def origami_ch = COrigamiPredictChrom.out.origami_ch

  if (params.run_hic_tile_viz) {
  
    // pred_ch emits:   (sample, chrom, prep_dir, pred_dir)
    // origami_ch emits:(sample, chrom, origami_dir)
  
    // join => (sample, chrom, prep_dir, pred_dir, origami_dir)
    def bundle = pred_ch
      .join(origami_ch, by: [0,1])
      .map { s, chrom, prep, pred, orig ->
        tuple(s, chrom, orig, prep, pred)     // (sample, chrom, orig, prep, pred)
      }

    bundle
      .filter { s, c, o, p, pr -> c == (params.hic_viz_chrom as String) }
      .view   { s, c, o, p, pr -> "bundle: s=$s c=$c orig=$o prep=$p pred=$pr" } 

    // pairs_ch emits: (pair_id, a, b)
    def pairs_ch = Channel
      .fromList(params.hic_viz_pairs as List)
      .map { p ->
        def toks = p.toString().split(':')
        def a = toks[0]
        def b = toks[1]
        tuple("${a}_vs_${b}", a, b)
      }
 
    pairs_ch.view { pid, a, b -> "pair: $pid  A=$a  B=$b" }
 
    // Key pairs by sample_a
    def a_keyed = pairs_ch.map { pid, a, b -> tuple(a, pid, a, b) }  // key = a
  
    // join bundle with a_keyed by sample => flattened tuple:
    // (sample, chrom, orig, prep, pred, pid, a, b)
    def A_side = bundle
      .join(a_keyed, by: 0)
      .map { s, chrom, orig, prep, pred, pid, a, b ->
        tuple(pid, a, b, chrom, orig, prep, pred)
      }
  
    // Key pairs by sample_b
    def b_keyed = pairs_ch.map { pid, a, b -> tuple(b, pid, a, b) }  // key = b
  
    def B_side = bundle
      .join(b_keyed, by: 0)
      .map { s, chrom, orig, prep, pred, pid, a, b ->
        tuple(pid, a, b, chrom, orig, prep, pred)
      }
  
    // Join A & B by (pair_id, chrom)
    // output flattened:
    // (pid, a, b, chrom, origA, prepA, predA, origB, prepB, predB)
    def compare_in = A_side
      .join(B_side, by: [0,3])
      .map { pid, chrom, a, b, origA, prepA, predA, a2, b2, origB, prepB, predB ->
    
        // optional sanity check
        if( a != a2 || b != b2 ) {
          log.warn "Pair mismatch for ${pid} ${chrom}: A_side=${a}:${b} B_side=${a2}:${b2}"
        }
    
        tuple(a, b, chrom, origA, prepA, predA, origB, prepB, predB)
      }
  
    // If you ONLY want one chrom, keep this filter
    if (params.hic_viz_chrom) {
      def chr = params.hic_viz_chrom as String
      compare_in = compare_in.filter { a, b, c, oa, pa, pra, ob, pb, prb -> c == chr }
    }
  
    HiCTileVizCompare(compare_in)
  }
}
