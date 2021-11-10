"""
Prepare input tsv file for nagim to run with the workflow as follows:
batch_workflow.py --input-tsv gs://cpg-nagim-test/joint_calling/samples.tsv
"""

import os
import pandas as pd

os.system('gsutil ls \'gs://cpg-nagim-test-upload/gvcf/*.vcf.gz\' > gvcflist.txt')

default_entry = {
    's': None,
    'external_id': None,
    'project': None,
    'continental_pop': '-',
    'subpop': '-',
    'gvcf': '-',
    'cram': '-',
    'crai': '-',
    'realign_cram': '-',
    'realign_crai': '-',
    'batch': '-',
    'operation': 'add',
    'flowcell_lane': '-',
    'library_id': '-',
    'platform': '-',
    'centre': '-',
    'r_contamination': None,
    'r_chimera': None,
    'r_duplication': None,
    'median_insert_size': None,
    'fam_id': '-',
    'mat_id': '-',
    'pat_id': '-',
    'sex': '-',
    'age': None,
}

datas = []
with open('gvcflist.txt') as in_f:
    for line in in_f:
        gvcf_path = line.strip()
        sample_name = os.path.basename(gvcf_path).replace('.hard-filtered.g.vcf.gz', '')
        entry = default_entry.copy()
        entry.update(
            {
                's': sample_name,
                'external_id': sample_name,
                'gvcf': gvcf_path,
                'project': 'nagim',
            }
        )
        datas.append(entry)

df = pd.DataFrame(datas)
out_path = 'gs://cpg-nagim-test/joint_calling/samples.tsv'
df.to_csv(out_path, sep='\t', index=False)
print(f'Written samples TSV to {out_path}')
