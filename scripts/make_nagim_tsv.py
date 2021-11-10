"""
Prepare input tsv file for nagim to run with the workflow as follows:
batch_workflow.py --input-tsv gs://cpg-nagim-test/joint_calling/samples.tsv
"""

import os
import pandas as pd
from joint_calling import utils

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)

run_cmd('gsutil ls \'gs://cpg-nagim-test-upload/gvcf/*.vcf.gz\' > gvcflist.txt')

default_entry = {
    's': None,
    'external_id': None,
    'project': None,
    'continental_pop': '-',
    'subpop': '-',
    'gvcf': '-',
    'topostproc_gvcf': '-',
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
        if not utils.file_exists(gvcf_path):
            print(f'gvcf doesnt exist for {sample_name}: {gvcf_path}')
            continue
            
        if not utils.file_exists(gvcf_path + '.tbi'):
            print(f'tbi doesnt exist for {sample_name}: {gvcf_path + ".tbi"}')
            continue
        
        print(f'adding sample {sample_name} with gvcf {gvcf_path}')
        entry = default_entry.copy()
        entry.update(
            {
                's': sample_name,
                'external_id': sample_name,
                'topostproc_gvcf': gvcf_path,
                'project': 'nagim',
            }
        )
        datas.append(entry)

df = pd.DataFrame(datas)
out_fname = 'samples.tsv'
out_fpath = f'gs://cpg-nagim-test/joint_calling/{out_fname}'
df.to_csv(out_fname, sep='\t', index=False)
run_cmd(f'gsutil cp {out_fname} {out_fpath}')
print(f'Written samples TSV to {out_fpath}')
