from typing import List, Tuple

import os
import hailtop.batch as hb

DEFAULT_REFERENCE_PATH = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
)
SITES_PATH = "gs://cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/sites.hg38.vcf.gz"
SOMALIER_CONTAINER = "brentp/somalier:v0.2.13"

def validate_ped(pedfile):
    from ped_parser import FamilyParser

    with open(pedfile) as fp:
        parser = FamilyParser(
            family_info=fp, family_type='ped'
        )
        print(f"Successfully validated pedigree with {len(parser.familes)} familes.")
    return True

def check_pedfile_and_relatedness(pedfile, somalier_relatedness_file):
    from ped_parser import FamilyParser

    with open(pedfile) as fp:
        my_parser = FamilyParser(
            family_info=fp, family_type='ped'
        )
        # validate relatedness

    return False


def create_input_group_from_read(batch, read):
    if not isinstance(read, str):
        raise NotImplementedError(
            f"Somalier check doesn't support reads of type {type(read)}"
        )

    kwargs = {'base': read}
    if read.endswith('.cram'):
        # cram.crai
        kwargs['index'] = read + '.crai'
    elif read.endswith('.bam'):
        # bam.bai
        kwargs['index'] = read + '.bai'
    elif read.endswith('.fastq') or read.endswith('.fastq.gz'):
        # no extension for fastqs
        pass
    else:
        raise NotImplementedError("Couldn't recognise read with extension")

    return batch.read_input_group(**kwargs)


def create_somalier_extract_job(
    batch,
    basename,
    read_group,
    reference_group,
    sites,
    somalier_container=SOMALIER_CONTAINER,
):
    somalier_job = (
        batch.new_job(f'somalier-{basename}')
        .image(somalier_container)
        .storage('80Gi')
        .memory('16Gi')
    )

    output_file = somalier_job.out_somalier
    somalier_job.command(
        f"""
mkdir out
somalier extract -d out \\
-f {reference_group.base} \\
--sites {sites} \\
{read_group.cram}
# there only one output, so this should be fine
cp out/*.somalier {output_file} 
"""
    )

    return somalier_job, output_file


def main(
    output_dir: str,
    pedfile_path: str,
    read_paths: List[str],
    reference_path=DEFAULT_REFERENCE_PATH,
    sites_path=SITES_PATH,
    somalier_container=SOMALIER_CONTAINER,
    batch: hb.Batch = None,
    # checkpointing
    somalier_files=None,
) -> Tuple[hb.Batch, hb.Job]:
    if batch is None:
        batch = hb.Batch('somalier-pedigree-checks')

    pedfile = batch.read_input(pedfile_path)
    reference = batch.read_input_group(
        base=reference_path,
        dict=reference_path.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
        fai=reference_path + '.fai',
    )
    sites = batch.read_input(sites_path)

    # generate somalier outputs
    if somalier_files:
        somalier_files = map(batch.read_input, somalier_files)
    else:
        reads = [create_input_group_from_read(batch, read) for read in read_paths]

        somalier_files = []
        for read, read_path in zip(reads, read_paths):
            basename = os.splitext(os.path.basename(read_path))[0]
            _, somalier_file = create_somalier_extract_job(
                batch=batch,
                read_group=read,
                basename=basename,
                sites=sites,
                reference_group=reference,
                somalier_container=somalier_container,
            )
            batch.write_output(
                somalier_file,
                os.path.join(output_dir, 'somalier', basename + ".somalier"),
            )
            somalier_files.append(somalier_file)

    report_job = batch.new_job("combiner").image(somalier_container)

    samples_report = report_job.out_samples
    pairs_report = report_job.out_pairs

    report_job.command(
        f"""
somalier relate --ped {pedfile} {' '.join(somalier_files)} 
cp somalier.samples.tsv {samples_report}
cp somalier.pairs.tsv {pairs_report}
"""
    )

    batch.write_output(samples_report, os.path.join(output_dir, 'somalier.samples.tsv'))
    batch.write_output(pairs_report, os.path.join(output_dir, 'somalier.pairs.tsv'))

    last_job = report_job

    return batch, last_job


def na12878_test():
    output_bucket = (
        'cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/2021-06-16_attempt-2'
    )
    output_dir = 'gs://' + output_bucket
    ped_path = 'gs://cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/na12878.ped'
    crams = [
        'gs://cpg-seqr-testing/NA12878_trio/NA12878.cram',
        'gs://cpg-seqr-testing/NA12878_trio/NA12891.cram',
        'gs://cpg-seqr-testing/NA12878_trio/NA12892.cram',
    ]

    somalier_files = [
        "gs://cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/2021-06-16_attempt-2/somalier/NA12878.cram.somalier",
        "gs://cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/2021-06-16_attempt-2/somalier/NA12891.cram.somalier",
        "gs://cpg-seqr-testing/NA12878_trio/mfranklin-ancestry/2021-06-16_attempt-2/somalier/NA12892.cram.somalier",
    ]

    backend = hb.ServiceBackend(
        billing_project='seqr',
        bucket=output_bucket,
    )
    batch = hb.Batch('na12878-ped-inference', backend=backend)
    batch = main(
        batch=batch,
        output_dir=output_dir,
        pedfile_path=ped_path,
        read_paths=crams,
        somalier_files=somalier_files,
    )
    batch.run(dry_run=False)


if __name__ == "__main__":
    # na12878_test()
    load_ped('/Users/michael.franklin/Downloads/cpg_acute.ped')
