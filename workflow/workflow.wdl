version 1.0


workflow CpgQC {
  input {
    File   sample_map
    String combined_mt_file
    String work_bucket
    String logs_dir_path
  }

  call CombineGVCFs {
    input:
      sample_map = sample_map,
      combined_mt_file = combined_mt_file,
      work_bucket = work_bucket,
      logs_dir_path = logs_dir_path
  }

  call SampleQC {
    input:
      combined_mt_file = CombineGVCFs.combined_mt_file,
      combined_mt_dir = CombineGVCFs.combined_mt_dir,
      work_bucket = work_bucket,
      logs_dir_path = logs_dir_path
  }

  output {
    String final_mt_dir = SampleQC.qced_mt_dir
    File final_mt_file = SampleQC.qced_mt_file
  }
}

task CombineGVCFs {
  input {
    File sample_map
    String combined_mt_file
    String work_bucket
    String logs_dir_path
  }
  output {
    String combined_mt_dir = "${combined_mt_file}"
    File combined_mt_file = "${combined_mt_file}" + "/_SUCCESS"
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{combined_mt_file} \
       --bucket ~{work_bucket} \
       --log-dir ~{logs_dir_path} \
       combine-gvcfs \
       --sample-map ~{sample_map} \
  >>>
}

task SampleQC {
  input {
    String combined_mt_dir
    File combined_mt_file
    String work_bucket
    String logs_dir_path
  }

  output {
    String qced_mt_dir = basename(combined_mt_dir, ".mt") + ".qc.mt"
    File qced_mt_file = basename(combined_mt_dir, ".mt") + ".qc.mt/_SUCCESS"
  }
  command <<<
    python /Users/vlad/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt ~{combined_mt_dir} \
       --bucket ~{work_bucket} \
       --log-dir ~{logs_dir_path} \
       sample-qc \
  >>>
}
