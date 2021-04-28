#!/usr/bin/env python

"""
To run the batch_workflow.py with the analysis runner. Assuming 
the joint-calling repositroy is a submodule of the dataset repository.

# Test access level - reads from `test`, writes to `temporary`:
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-test" \
    --description "joint calling test" \
    --access-level test \
    joint-calling/workflows/drive_joint_calling.py \
    --access-level test \
    --callset tob-wgs

# Standard access level - reads from `main`, writes to `temporary`:
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" \
    --description "joint calling standard" \
    --access-level standard \
    joint-calling/workflows/drive_joint_calling.py \
    --access-level standard \
    --callset tob-wgs \
    --version v0 --batch 0 --batch 1

# Full access level - reads from `main`, writes matrix tables to `main`, 
# the rest to `temporary`:
$ analysis-runner \
    --dataset tob-wgs \
    --output-dir "gs://cpg-tob-wgs-main/joint-calling" \
    --description "joint calling full" \
    --access-level full \
    joint-calling/workflows/drive_joint_calling.py \
    --access-level full \
    --callset tob-wgs \
    --version v0 --batch 0 --batch 1
"""

import subprocess
import click
import logging

logger = logging.getLogger('combine_gvcfs')
logger.setLevel('INFO')


def run_cmd(cmd):
    """Print the command and run"""
    print(cmd)
    subprocess.run(
        cmd,
        shell=True,
        check=False,
    )


@click.command()
@click.option('--access-level', 'access_level', type=click.Choice(['test', 'standard', 'full']))
@click.option('--callset', 'callset_name', type=str, required=True)
@click.option('--version', 'callset_version', type=str, required=True)
@click.option('--batch', 'callset_batches', type=str, multiple=True)
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
    type=click.Choice(['main', 'temporary']),
    default='temporary',
    help='The bucket type to write matrix tables to (default: temporary)',
)
@click.option(
    '--skip-qc',
    'skip_qc',
    is_flag=True,
    help='Do not attempt to find QC CSV files. '
    'Get sample IDs from the GVCF file names.',
)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
def main(
    access_level,
    callset_name,
    callset_version,
    callset_batches,
    input_bucket_suffix,
    output_bucket_suffix,
    skip_qc,
    keep_scratch,
    overwrite,
):
    """
    Runs the the joint-calling batch workflow script
    """
    if not callset_batches:
        if access_level == 'test':
            callset_batches = ['0']
        else:
            raise click.BadParameter(
                'Please specify batch numbers to process with --batch'
            )
    
    run_cmd(
        'PYTHONPATH=$PWD/joint-calling python joint-calling/workflows/batch_workflow.py '
        + f'--callset {callset_name} '
        + f'--version {callset_version} '
        + ' '.join([f'--batch {b}' for b in callset_batches])
        + '--from '
        + ('test ' if is_test else f'{input_bucket_suffix} ')
        + '--to '
        + ('temporary ' if is_test else f'{output_bucket_suffix} ')
        + ('--skip-qc ' if skip_qc else '')
        + ('--keep-scratch ' if keep_scratch else '')
        + f'--billing-project {callset_name} '
        + f'{"--overwrite " if overwrite else ""}'
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
