from typing import Callable

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path
from hailtop.batch import Batch
from hailtop.batch.job import Job

from larcoh.hb.batch import setup_batch
from larcoh.providers.cpg.inputs import CpgInputProvider
from larcoh.providers.cpg.metamist import Metamist
from larcoh.targets import Cohort
from larcoh.types import SequencingType
from larcoh.utils import can_reuse, timestamp, slugify


class Pipeline:
    def __init__(self):
        self.output_version = get_config()['workflow']['output_version']
        self.analysis_prefix = to_path(
            dataset_path(f'joint-calling/{self.output_version}', category='analysis')
        )

        run_id = get_config()['workflow'].get('run_id', timestamp())
        analysis_dataset = get_config()['workflow']['dataset']
        name = get_config()['workflow'].get('name')
        description = get_config()['workflow'].get('description')
        name = name or description or analysis_dataset
        name = slugify(name)
        description = description or name
        description += f': run_id={run_id}'

        self.cohort = Cohort(
            analysis_dataset_name=analysis_dataset,
            name=name,
            sequencing_type=SequencingType.parse(
                get_config()['workflow']['sequencing_type']
            ),
        )

        CpgInputProvider(Metamist(self.cohort.analysis_dataset.name)).populate_cohort(
            cohort=self.cohort,
            dataset_names=get_config()['workflow'].get('datasets'),
            skip_samples=get_config()['workflow'].get('skip_samples'),
            only_samples=get_config()['workflow'].get('only_samples'),
            skip_datasets=get_config()['workflow'].get('skip_datasets'),
        )

        self.tmp_prefix = self.cohort.analysis_dataset.tmp_prefix() / run_id

        description += f' [{self.cohort.sequencing_type.value}]'
        if ds_set := set(d.name for d in self.cohort.get_datasets()):
            description += ' ' + ', '.join(sorted(ds_set))
        self.b = setup_batch(description=description)
        self.tasks = []

    def add_task(self, task: 'Task'):
        self.tasks.append(task)

    def create_task(
        self,
        name: str,
        work_fn: Callable[['Task'], list[Job]],
        inputs: list[Path] | None = None,
        outputs: list[Path] | None = None,
    ):
        self.add_task(
            Task(
                name=name,
                pipeline=self,
                work_fn=work_fn,
                inputs=inputs,
                outputs=outputs,
            )
        )

    def run(self):
        self.b.run()


class Task:
    def __init__(
        self,
        name: str,
        work_fn: Callable[['Task'], list[Job]],
        pipeline: 'Pipeline',
        inputs: list[Path] | None = None,
        outputs: list[Path] | None = None,
    ):
        self.name = name
        self.work_fn = work_fn
        self.inputs = inputs or []
        self.outputs = outputs or []
        self.pipeline = pipeline

    def run(self) -> list[Job]:
        overwrite = get_config()['workflow'].get('overwrite', False)
        if self.outputs and can_reuse(self.outputs, overwrite):
            return [self.pipeline.b.new_job(f'{self.name} [reuse]')]
        return self.work_fn(self)
