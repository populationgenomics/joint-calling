version 1.0


workflow SampleQcWorkflow {
  input {
    String mt_dir
    String work_bucket
    String local_tmp_dir
    File sample_map
  }

  call SampleQC {
    input:
      input_mt_dir = mt_dir,
      work_bucket = work_bucket,
      out_ht_dir = work_bucket + "/sample_qc.ht",
      local_tmp_dir = local_tmp_dir
  }

  call ComputeSex {
    input:
      input_mt_dir = mt_dir,
      work_bucket = work_bucket,
      out_ht_dir = work_bucket + "/sex.ht",
      local_tmp_dir = local_tmp_dir
  }

  call ComputeQcMt {
    input:
      input_mt = mt_dir,
      work_bucket = work_bucket,
      out_mt = work_bucket + "/qc.mt",
      local_tmp_dir = local_tmp_dir
  }

  call ComputeHardFilters {
    input:
      input_mt_dir = mt_dir,
      sample_map = sample_map,
      work_bucket = work_bucket,
      out_ht_dir = work_bucket + "/hard_filtered_samples.ht",
      sex_ht_dir = ComputeSex.sex_ht_dir,
      biallelic_qc_ht_dir = SampleQC.bi_alleleic_qc_ht_dir,
      local_tmp_dir = local_tmp_dir
  }

  call PCRelate {
    input:
      input_mt = ComputeQcMt.qc_mt,
      work_bucket = work_bucket,
      out_ht = work_bucket + "/relatedness.ht",
      local_tmp_dir = local_tmp_dir
  }

  output {
#    File sample_qc_ht_mrk = SampleQC.qc_ht_mrk
#    File sex_ht_mrk = ImputeSex.sex_ht_mrk
    String sample_qc_ht_dir = SampleQC.sample_qc_ht_dir
    String sex_ht_dir = ComputeSex.sex_ht_dir
    String hard_filters_ht_dir = ComputeHardFilters.hard_filters_ht_dir
  }
}

task SampleQC {
  input {
    String input_mt_dir
    String work_bucket
    String out_ht_dir
    String local_tmp_dir
  }

  output {
    String sample_qc_ht_dir = out_ht_dir
    String bi_alleleic_qc_ht_dir = basename(out_ht_dir, ".ht") + "_bi_allelic.ht"
#    File   qc_ht_mrk = out_ht_dir + "/_SUCCESS"
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{input_mt_dir} \
       --bucket ~{work_bucket} \
       --local-tmp-dir ~{local_tmp_dir} \
       --reuse \
       sample-qc \
       --out-ht ~{out_ht_dir} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}

task ComputeQcMt {
  input {
    String input_mt
    String work_bucket
    String out_mt
    String local_tmp_dir
  }

  output {
    String qc_mt = out_mt
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{input_mt} \
       --bucket ~{work_bucket} \
       --local-tmp-dir ~{local_tmp_dir} \
       --reuse \
       compute-qc-mt \
       --out-mt ~{out_mt} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}

task ComputeSex {
  input {
    String input_mt_dir
    String work_bucket
    String out_ht_dir
    String local_tmp_dir
  }

  output {
    String sex_ht_dir = out_ht_dir
#    File   sex_ht_mrk = out_ht_dir + "/_SUCCESS"
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{input_mt_dir} \
       --bucket ~{work_bucket} \
       --local-tmp-dir ~{local_tmp_dir} \
       --reuse \
       impute-sex \
       --out-ht ~{out_ht_dir} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}

task ComputeHardFilters {
  input {
    String input_mt_dir
    String biallelic_qc_ht_dir
    String sex_ht_dir
    String work_bucket
    String out_ht_dir
    String local_tmp_dir
    File sample_map
  }

  output {
    String hard_filters_ht_dir = out_ht_dir
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{input_mt_dir} \
       --bucket ~{work_bucket} \
       --local-tmp-dir ~{local_tmp_dir} \
       --reuse \
       compute-hard-filters \
       --sample-map ~{sample_map} \
       --out-ht ~{out_ht_dir} \
       --biallelic-qc-ht ~{biallelic_qc_ht_dir} \
       --sex-ht ~{sex_ht_dir} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}
task PCRelate {
  input {
    String input_mt
    String work_bucket
    String out_ht
    String local_tmp_dir
  }

  output {
    String relatedness_ht = out_ht
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{input_mt} \
       --bucket ~{work_bucket} \
       --local-tmp-dir ~{local_tmp_dir} \
       --reuse \
       run-pc-relate \
       --out-ht ~{out_ht} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}
