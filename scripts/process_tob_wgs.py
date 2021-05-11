"""
Process TOB-WGS files within upload bucket

"""

import os
import hailtop.batch as hb
from joint_calling.upload_processor import batch_move_files


def filter_files():
    """ Determine files that should be moved to main vs archive bucket """
    files_for_main = []
    files_for_archive = []
    prefix = 'gs://cpg-tob-wgs-upload/'
    with open('../test/data/tob_wgs_batch2.txt') as f:
        for file_path in f:
            file_path = file_path.strip()
            if file_path.endswith(
                ('.g.vcf.gz', '.g.vcf.gz.tbi', '.g.vcf.gz.md5', 'csv')
            ):
                files_for_main.append(file_path[len(prefix) :])
            else:
                files_for_archive.append(file_path[len(prefix) :])

    return files_for_main, files_for_archive


if __name__ == '__main__':

    # Process input file of sample names
    main_files, archive_files = filter_files()

    # Setting up inputs for batch_move_files
    upload_prefix = os.path.join('cpg-tob-wgs-upload')
    main_prefix = os.path.join('cpg-tob-wgs-main', 'gVCFs', 'batch2')
    docker_image = os.environ.get('DOCKER_IMAGE')
    key = os.environ.get('GSA_KEY')

    # Moving the files to the main bucket
    main_batch = hb.Batch(name='Move Files to Main')
    batch_move_files(
        main_batch,
        main_files,
        upload_prefix,
        main_prefix,
        docker_image,
        key,
    )

    main_batch.run()
    archive_prefix = os.path.join('cpg-tob-archive', 'cram', 'batch2')

    # Moving the files to the archive bucket
    archive_batch = hb.Batch(name='Move Files to Archive')
    batch_move_files(
        archive_batch,
        archive_files,
        upload_prefix,
        archive_prefix,
        docker_image,
        key,
    )
    archive_batch.run()
