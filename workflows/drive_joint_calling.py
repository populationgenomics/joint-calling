#!/usr/bin/env python

"""
Drive the joint calling workflow.

Use with analysis runner. Copy the script into the scripts folder of the analysis
repository (e.g. tob-wgs), and run:

(test)
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-test" \
    --description "joint calling test" \
    --access-level test \
    scripts/drive_joint_calling.py --is-test --callset tob-wgs

(prod)
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" \
    --description "joint calling prod" \
    --access-level standard \
    scripts/drive_joint_calling.py --callset tob-wgs \
        --version v0 --batch 0 --batch 1 --to temporary
"""

import subprocess
import click


def run_cmd(cmd):
    """Print the command and run"""
    print(cmd)
    subprocess.run(
        cmd,
        shell=True,
        check=False,
    )


@click.command()
@click.option('--is-test', 'is_test', is_flag=True)
@click.option('--callset', 'callset_name', type=str, required=True)
@click.option('--version', 'callset_version', type=str, required=True)
@click.option('--batch', 'callset_batches', type=str, multiple=True, required=True)
@click.option(
    '--from',
    'input_bucket_suffix',
    type=click.Choice(['main', 'test']),
    default='test',
    help='The bucket type to read from (default: test)',
)
@click.option(
    '--to',
    'output_bucket_suffix',
    type=click.Choice(['analysis', 'temporary']),
    default='temporary',
    help='The bucket type to write to (default: temporary)',
)
@click.option(
    '--skip-qc',
    'skip_qc',
    is_flag=True,
    help='Do not attempt to find QC CSV files. '
    'Get sample IDs from the GVCF file names.',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
def main(
    is_test,
    callset_name,
    callset_version,
    callset_batches,
    input_bucket_suffix,
    output_bucket_suffix,
    skip_qc,
    keep_scratch,
):
    """
    Runs the the joint-calling batch workflow script
    """
    batches_cmdl = ' '.join([f'--batch {b}' for b in callset_batches])
    run_cmd(
        'PYTHONPATH=$PWD/../joint-calling python ../joint-calling/workflows/batch_workflow.py '
        + f'--callset {callset_name} '
        + f'--version {callset_version} '
        + ('--batch 0 ' if is_test else f'{batches_cmdl} ')
        + '--from '
        + ('test ' if is_test else f'{input_bucket_suffix} ')
        + '--to '
        + ('temporary ' if is_test else f'{output_bucket_suffix} ')
        + ('--skip-qc ' if skip_qc else '')
        + ('--keep-scratch ' if keep_scratch else '')
        + f'--billing-project {callset_name} '
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
