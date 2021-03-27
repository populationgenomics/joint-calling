SAMPLE_MAP=gs://cpg-fewgenomes-temporary/joint-calling/100genomes-gvcf-gvcf-local.tsv
OUTPUT_BUCKET=gs://cpg-fewgenomes-temporary/joint-calling/test-run
CALLSET=fewgenomes

CMD="vqsr_batch.py\
 --sample-map $SAMPLE_MAP\
 --output_bucket $OUTPUT_BUCKET \
 --callset_name $CALLSET\
 --keep_scratch"

analysis-runner \
--dataset $CALLSET \
--access-level test \
--output-dir $OUTPUT_BUCKET \
--description "Joint-calling wofkflow" \
"$CMD"
