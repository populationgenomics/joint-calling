python vqsr_batch.py \
	--combined_gvcf gs://cpg-fewgenomes-test/vqsr-testing/sites.vcf.bgz \
	--callset_name fewgenomes \
	--billing_project fewgenomes \
	--num_gvcfs 48 \
	--keep_scratch \
	--output_bucket gs://playground-au/vqsr
