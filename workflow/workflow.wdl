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
      work_bucket = work_bucket,
      logs_dir_path = logs_dir_path
  }

  output {
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
    File combined_mt_file = "${combined_mt_file}"
  }
  command <<<
    python /Users/vsaveliev/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt-file ~{combined_mt_file} \
       --bucket ~{work_bucket} \
       --log-dir ~{logs_dir_path} \
       combine-gvcfs \
       --sample-map ~{sample_map} \
  >>>
}

task SampleQC {
  input {
    File combined_mt_file
    String work_bucket
    String logs_dir_path
  }

  output {
    File qced_mt_file = basename(combined_mt_file, ".mt") + ".qc.mt"
  }
  command <<<
    python /Users/vsaveliev/CPG/cpg_qc/scripts/cpg_qc.py \
       --mt-file ~{combined_mt_file} \
       --bucket ~{work_bucket} \
       --log-dir ~{logs_dir_path} \
       sample-qc \
  >>>
}
