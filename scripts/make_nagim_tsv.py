"""
Prepare input tsv file for nagim to run with the workflow as follows:
batch_workflow.py --input-tsv gs://cpg-nagim-test/joint_calling/samples.tsv
"""

import os
import pandas as pd

os.system('gsutil ls \'gs://cpg-nagim-test-upload/gvcf/*.vcf.gz\' > gvcflist.txt')

datas = []
with open('gvcflist.txt') as in_f:
    for line in in_f:
        gvcf_path = line.strip()
        sample_name = os.path.basename(gvcf_path).replace('.hard-filtered.g.vcf.gz', '')
        datas.append({'s': sample_name, 'external_id': sample_name, 'gvcf': gvcf_path})

df = pd.DataFrame(datas)
df.to_csv('gs://cpg-nagim-test/joint_calling/samples.tsv', sep='\t', index=False)
