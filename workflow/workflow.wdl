version 1.0


import "SampleQC.wdl" as SampleQC


workflow CpgQC {
  input {
    File   sample_map
    String work_bucket
    String logs_dir_path
    String? mt_path
  }

  String mt = select_first([mt_path, work_bucket + "/genomes.mt"])

  call CombineGVCFs {
    input:
      sample_map = sample_map,
      out_mt_dir = mt,
      work_bucket = work_bucket,
      logs_dir_path = logs_dir_path
  }

  call SampleQC.SampleQcWorkflow as SampleQcWfl {
    input:
      mt_dir = CombineGVCFs.combined_mt_dir,
      sample_map = sample_map,
      work_bucket = work_bucket,
      logs_dir_path = logs_dir_path
  }

  output {
#    File sample_qc_ht_mrk = SampleQcWfl.sample_qc_ht_mrk
#    File sex_ht_mrk       = SampleQcWfl.sex_ht_mrk
    String sample_qc_ht_dir = SampleQcWfl.sample_qc_ht_dir
    String sex_ht_dir = SampleQcWfl.sex_ht_dir
    String hard_filters_ht_dir = SampleQcWfl.hard_filters_ht_dir
  }
}

task CombineGVCFs {
  input {
    File sample_map
    String out_mt_dir
    String work_bucket
    String logs_dir_path
  }
  output {
    String combined_mt_dir = out_mt_dir
#    File   combined_mt_mrk = out_mt_dir + "/_SUCCESS"
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{out_mt_dir} \
       --bucket ~{work_bucket} \
       --log-dir ~{logs_dir_path} \
       --reuse \
       combine-gvcfs \
       --sample-map ~{sample_map} \
  >>>
#  runtime {
#    memory: "10 GiB"
#    disks: "local-disk 50 HDD"
#    preemptible: 1
#    docker: "quay.io/biocontainers/hail:0.2.33--py37hc9558a2_0"
#  }
}
