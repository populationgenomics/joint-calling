#!/usr/bin/env python3

"""
Adds Batch job that drive Somalier for ancestry/pedigree
"""

import logging
import subprocess
from os.path import join
from typing import Optional, List, Tuple
import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job
from joint_calling import resources, utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def pedigree_checks(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    overwrite: bool,  # pylint: disable=unused-argument
    output_suffix: str,
    relatedness_bucket: str,
    web_bucket: str,
    web_url: str,
    tmp_bucket: str,
    ped_fpath: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Job, str, str, str]:
    """
    Add somalier and peddy based jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, a path to a fixed PED file if able to recover, and a path to a file
    with relatedness information for each sample pair
    """
    extract_jobs = []
    fp_file_by_sample = dict()

    def get_project_bucket(_proj):
        if _proj in ['syndip', 'giab']:
            proj_bucket = f'gs://cpg-reference/validation/{_proj}'
        else:
            proj_bucket = f'gs://cpg-{_proj}-{output_suffix}'
        return proj_bucket

    for sn, proj, gvcf_path in zip(samples_df.s, samples_df.project, samples_df.gvcf):
        proj_bucket = get_project_bucket(proj)
        fp_file_by_sample[sn] = join(proj_bucket, 'fingerprints', f'{sn}.somalier')
        if utils.can_reuse(fp_file_by_sample[sn], overwrite):
            extract_jobs.append(b.new_job(f'Somalier extract, {sn} [reuse]'))
        else:
            j = b.new_job(f'Somalier extract, {sn}')
            j.image(utils.SOMALIER_IMAGE)
            j.memory('standard')
            if gvcf_path.endswith('.bam'):
                j.cpu(4)
                j.storage(f'200G')
            elif gvcf_path.endswith('.cram'):
                j.cpu(4)
                j.storage(f'50G')
            else:
                j.cpu(2)
                j.storage(f'10G')
            if depends_on:
                j.depends_on(*depends_on)

            input_file = b.read_input_group(
                base=gvcf_path,
                index=gvcf_path + '.tbi',
            )

            sites = b.read_input(resources.SOMALIER_SITES)
            reference = b.read_input_group(
                base=resources.REF_FASTA,
                fai=resources.REF_FASTA + '.fai',
                dict=resources.REF_FASTA.replace('.fasta', '')
                .replace('.fna', '')
                .replace('.fa', '')
                + '.dict',
            )
            j.command(
                f"""set -ex
                
                somalier extract -d extracted/ --sites {sites} -f {reference.base} \\
                {input_file['base']}
                
                mv extracted/*.somalier {j.output_file}
                """
            )
            b.write_output(j.output_file, fp_file_by_sample[sn])
            extract_jobs.append(j)

    relate_j = b.new_job(f'Somalier relate')
    relate_j.image(utils.SOMALIER_IMAGE)
    relate_j.cpu(1)
    relate_j.memory('standard')  # ~ 4G/core ~ 4G
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    relate_j.storage(f'{1 + len(extract_jobs) // 4000 * 1}G')
    relate_j.depends_on(*extract_jobs)
    fp_files = [b.read_input(fp) for sn, fp in fp_file_by_sample.items()]

    if ped_fpath:
        ped_file = b.read_input(ped_fpath)
    else:
        ped_fpath = join(tmp_bucket, 'samples.ped')
        samples_df['Family.ID'] = samples_df['fam_id']
        samples_df['Individual.ID'] = samples_df['s']
        samples_df['Father.ID'] = samples_df['pat_id']
        samples_df['Mother.ID'] = samples_df['mat_id']
        samples_df['Sex'] = samples_df['sex']
        samples_df['Phenotype'] = 0
        samples_df[
            ['Family.ID', 'Individual.ID', 'Father.ID', 'Mother.ID', 'Sex', 'Phenotype']
        ].to_csv(
            ped_fpath,
            sep='\t',
            index=False,
        )
        ped_file = b.read_input(ped_fpath)

    relate_j.command(
        f"""set -e

        cat {ped_file} | grep -v Family.ID > samples.ped 

        somalier relate \\
        {' '.join(fp_files)} \\
        --ped samples.ped \\
        -o related \\
        --infer

        ls
        mv related.html {relate_j.output_html}
        mv related.pairs.tsv {relate_j.output_pairs}
        mv related.samples.tsv {relate_j.output_samples}
      """
    )

    # Copy somalier outputs to buckets
    sample_hash = utils.hash_sample_ids(samples_df['s'])
    prefix = join(relatedness_bucket, 'somalier', sample_hash, 'somalier')
    somalier_samples_path = f'{prefix}.samples.tsv'
    somalier_pairs_path = f'{prefix}.pairs.tsv'
    b.write_output(relate_j.output_samples, somalier_samples_path)
    b.write_output(relate_j.output_pairs, somalier_pairs_path)
    # Copy somalier HTML to the web bucket
    rel_path = join('loader', sample_hash, 'somalier.html')
    somalier_html_path = join(web_bucket, rel_path)
    somalier_html_url = f'{web_url}/{rel_path}'
    b.write_output(relate_j.output_html, somalier_html_path)

    check_j = b.new_job(f'Check relatedness and sex')
    check_j.image(utils.PEDDY_IMAGE)
    check_j.cpu(1)
    check_j.memory('standard')  # ~ 4G/core ~ 4G

    script_name = 'check_pedigree.py'
    try:
        script_path = (
            subprocess.check_output(f'which {script_name}', shell=True).decode().strip()
        )
    except subprocess.CalledProcessError:
        script_path = join(utils.SCRIPTS_DIR, script_name)

    with open(script_path) as f:
        script = f.read()
    check_j.command(
        f"""set -ex
cat <<EOT >> {script_name}
{script}
EOT
python {script_name} \
--somalier-samples {relate_j.output_samples} \
--somalier-pairs {relate_j.output_pairs} \
{('--somalier-html ' + somalier_html_url) if somalier_html_url else ''}
    """
    )
    check_j.command(f"""
# fix ped
cat {relate_j.output_samples} | cut -f1-6 | grep -v ^# > {check_j.fixed_ped}
ls $(dirname {check_j.fixed_ped})
cat {check_j.fixed_ped}
    """
    )
    fixed_ped_fpath = join(relatedness_bucket, 'samples.ped')
    b.write_output(check_j.fixed_ped, fixed_ped_fpath)

    check_j.depends_on(relate_j)
    return check_j, fixed_ped_fpath, somalier_samples_path, somalier_pairs_path
