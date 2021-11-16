"""
Prepare input tsv file for nagim to run with the workflow as follows:
batch_workflow.py --input-tsv gs://cpg-nagim-test/joint_calling/samples.tsv
"""
import os
import pandas as pd


def run_cmd(cmd):
    print(cmd)
    os.system(cmd)


gvcf_list_fpath = 'gvcflist.txt'
tbi_list_fpath = 'tbilist.txt'
if not os.path.isfile(gvcf_list_fpath):
    run_cmd(
        f'gsutil ls \'gs://cpg-nagim-test-upload/gvcf/*.vcf.gz\' > {gvcf_list_fpath}'
    )
if not os.path.isfile(tbi_list_fpath):
    run_cmd(
        f'gsutil ls \'gs://cpg-nagim-test-upload/gvcf/*.vcf.gz.tbi\' > {tbi_list_fpath}'
    )

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
    'resequencing_label': '-',
    'primary_study': '-',
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

gvcf_by_sample = dict()
with open(gvcf_list_fpath) as gvcf_f:
    for line in gvcf_f:
        fname = line.strip()
        sample = os.path.basename(fname).replace('.hard-filtered.g.vcf.gz', '')
        gvcf_by_sample[sample] = fname

tbi_by_sample = dict()
with open(tbi_list_fpath) as tbi_f:
    for line in tbi_f:
        fname = line.strip()
        sample = os.path.basename(fname).replace('.hard-filtered.g.vcf.gz.tbi', '')
        if sample not in gvcf_by_sample:
            print(f'Found TBI without GVCF: {fname}')
        else:
            tbi_by_sample[sample] = fname

datas = []
for sname, gvcf_path in gvcf_by_sample.items():
    tbi_path = tbi_by_sample.get(sname)
    if not tbi_path:
        print(f'tbi doesnt exist for {gvcf_path}')
        continue

    print(f'adding sample {sname} with gvcf {gvcf_path}')
    entry = default_entry.copy()
    entry.update(
        {
            's': sname,
            'external_id': sname,
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
