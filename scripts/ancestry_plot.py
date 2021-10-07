#!/usr/bin/env python

"""
Plot ancestry PCA analysis results
"""

import logging
from collections import Counter
from typing import List, Iterable, Optional

import click
import pandas as pd
import numpy as np
import hail as hl
from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.transform import factor_cmap, factor_mark
from bokeh.plotting import ColumnDataSource, figure
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
from bokeh.models import CategoricalColorMapper, HoverTool

from joint_calling import utils, resources
from joint_calling import _version


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--eigenvalues',
    'eigenvalues_path',
    required=True,
    callback=utils.get_validation_callback(ext='txt', must_exist=True),
)
@click.option(
    '--scores-ht',
    'scores_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--loadings-ht',
    'loadings_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--provided-pop-ht',
    'provided_pop_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option(
    '--inferred-pop-ht',
    'inferred_pop_ht_path',
    required=True,
    callback=utils.get_validation_callback(ext='ht', must_exist=True),
)
@click.option('--out-path-pattern', 'out_path_pattern')
@click.option(
    '--hail-billing',
    'hail_billing',
    required=True,
    help='Hail billing account ID.',
)
def main(  # pylint: disable=too-many-arguments,too-many-locals,missing-function-docstring
    eigenvalues_path: str,
    scores_ht_path: str,
    loadings_ht_path: str,
    provided_pop_ht_path: str,
    inferred_pop_ht_path: str,
    out_path_pattern: Optional[str],
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

    produce_plots(
        eigenvalues_path=eigenvalues_path,
        scores_ht_path=scores_ht_path,
        loadings_ht_path=loadings_ht_path,
        provided_pop_ht_path=provided_pop_ht_path,
        inferred_pop_ht_path=inferred_pop_ht_path,
        out_path_pattern=out_path_pattern,
    )


def produce_plots(
    eigenvalues_path: str,
    scores_ht_path: str,
    loadings_ht_path: str,
    provided_pop_ht_path: str,
    inferred_pop_ht_path: str,
    out_path_pattern: Optional[str] = None,
    number_of_pcs: Optional[int] = None,
):
    """
    Generate plots in png and html formats, write for each PC (of n_pcs) and
    scope ("study", "continental_pop", "subpop", plus for loadings) into
    file paths defined by `out_path_pattern`.
    """
    scores = hl.read_table(scores_ht_path)
    sample_names = scores.s.collect()
    provided_pop_ht = hl.read_table(provided_pop_ht_path)
    inferred_pop_ht = hl.read_table(inferred_pop_ht_path)
    scores = scores.annotate(
        continental_pop=hl.case()
        .when(
            provided_pop_ht[scores.s].continental_pop != '',
            provided_pop_ht[scores.s].continental_pop,
        )
        .default(
            inferred_pop_ht[scores.s].pop
            + ' ['
            + provided_pop_ht[scores.s].project
            + ', inferred]'
        ),
        subpop=hl.case()
        .when(provided_pop_ht[scores.s].subpop != '', provided_pop_ht[scores.s].subpop)
        .default(
            inferred_pop_ht[scores.s].pop
            + ' ['
            + provided_pop_ht[scores.s].project
            + ', inferred]'
        ),
        study=provided_pop_ht[scores.s].project,
    )

    eigenvalues = hl.import_table(
        eigenvalues_path, no_header=True, types={'f0': hl.tfloat}
    ).f0.collect()
    eigenvalues = pd.to_numeric(eigenvalues)
    variance = np.divide(eigenvalues[1:], float(eigenvalues.sum())) * 100
    variance = variance.round(2)
    number_of_pcs = number_of_pcs or len(eigenvalues) - 1

    plots = []

    studies = scores.study.collect()
    unique_studies = remove_duplicates(studies)
    bg_studies = [s for s in unique_studies if s == 'gnomad']
    fg_studies = [s for s in unique_studies if s != 'gnomad']

    plots.extend(
        _plot_study(
            number_of_pcs,
            variance,
            scores,
            studies,
            sample_names,
            bg_studies,
            fg_studies,
            out_path_pattern=out_path_pattern,
        )
    )

    plots.extend(
        _plot_continental_pop(
            number_of_pcs,
            variance,
            scores,
            studies,
            sample_names,
            bg_studies,
            fg_studies,
            out_path_pattern=out_path_pattern,
        )
    )

    plots.extend(
        _plot_subcontinental_pop(
            number_of_pcs,
            variance,
            scores,
            studies,
            sample_names,
            bg_studies,
            fg_studies,
            out_path_pattern=out_path_pattern,
        )
    )

    plots.extend(
        _plot_loadings(
            number_of_pcs, loadings_ht_path, out_path_pattern=out_path_pattern
        )
    )

    return plots


def _plot_study(
    number_of_pcs,
    variance,
    scores,
    studies,
    sample_names,
    bg_studies,
    fg_studies,
    out_path_pattern=None,
):
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    cntr: Counter = Counter(studies)
    labels = [f'{x} ({cntr[x]})' for x in studies]
    unique_labels = list(Counter(labels).keys())

    plots = []
    for i in range(number_of_pcs - 1):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title='Study',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=scores.scores[pc1].collect(),
                y=scores.scores[pc2].collect(),
                label=labels,
                samples=sample_names,
                study=studies,
            )
        )
        plot.scatter(
            'x',
            'y',
            alpha=0.5,
            marker=factor_mark(
                'study',
                ['cross'] * len(bg_studies) + ['circle'] * len(fg_studies),
                bg_studies + fg_studies,
            ),
            source=source,
            size=4,
            color=factor_cmap('label', ['#1b9e77', '#d95f02'], unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)
        if out_path_pattern:
            plot_filename = out_path_pattern.format(scope='study', pci=pc2, ext='png')
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')
            html = file_html(plot, CDN, 'my plot')
            plot_filename_html = out_path_pattern.format(
                scope='study', pci=pc2, ext='html'
            )
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def _plot_continental_pop(
    number_of_pcs,
    variance,
    scores,
    studies,
    sample_names,
    bg_studies,
    fg_studies,
    out_path_pattern=None,
):
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    labels = scores.continental_pop.collect()
    cntr = Counter(labels)
    labels = [f'{x} ({cntr[x]})' for x in labels]
    unique_labels = list(Counter(labels).keys())

    plots = []
    for i in range(number_of_pcs - 1):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title='Continental Population',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=scores.scores[pc1].collect(),
                y=scores.scores[pc2].collect(),
                label=labels,
                samples=sample_names,
                study=studies,
            )
        )
        plot.scatter(
            'x',
            'y',
            alpha=0.6,
            marker=factor_mark(
                'study',
                ['cross'] * len(bg_studies) + ['circle'] * len(fg_studies),
                bg_studies + fg_studies,
            ),
            source=source,
            size=2,
            color=factor_cmap('label', turbo(len(unique_labels)), unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)
        if out_path_pattern:
            plot_filename = out_path_pattern.format(
                scope='continental_pop', pci=pc2, ext='png'
            )
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')
            html = file_html(plot, CDN, 'my plot')
            plot_filename_html = out_path_pattern.format(
                scope='continental_pop', pci=pc2, ext='html'
            )
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def _plot_subcontinental_pop(
    number_of_pcs,
    variance,
    scores,
    studies,
    sample_names,
    bg_studies,
    fg_studies,
    out_path_pattern=None,
):
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    labels = scores.subpop.collect()
    cntr = Counter(labels)
    labels = [f'{x} ({cntr[x]})' for x in labels]
    unique_labels = list(Counter(labels).keys())

    plots = []
    for i in range(number_of_pcs - 1):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title='Subpopulation',
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]})%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
        )
        source = ColumnDataSource(
            dict(
                x=scores.scores[pc1].collect(),
                y=scores.scores[pc2].collect(),
                label=labels,
                samples=sample_names,
                study=studies,
            )
        )
        plot.scatter(
            'x',
            'y',
            alpha=0.6,
            marker=factor_mark(
                'study',
                ['cross'] * len(bg_studies) + ['circle'] * len(fg_studies),
                bg_studies + fg_studies,
            ),
            source=source,
            size=2,
            color=factor_cmap('label', turbo(len(unique_labels)), unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)
        if out_path_pattern:
            plot_filename = out_path_pattern.format(scope='subpop', pci=pc2, ext='png')
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')
            html = file_html(plot, CDN, 'my plot')
            plot_filename_html = out_path_pattern.format(
                scope='subpop', pci=pc2, ext='html'
            )
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def _plot_loadings(number_of_pcs, loadings_ht_path, out_path_pattern=None):
    loadings_ht = hl.read_table(loadings_ht_path)
    plots = []
    gtf_ht = hl.experimental.import_gtf(
        resources.GENCODE_GTF,
        reference_genome='GRCh38',
        skip_invalid_contigs=True,
        min_partitions=12,
        force_bgz=True,
    )
    for i in range(number_of_pcs - 1):
        pc = i + 1
        plot = manhattan_loadings(
            iteration=i,
            gtf=gtf_ht,
            loadings=loadings_ht,
            title='Loadings of PC ' + str(pc),
            collect_all=True,
        )
        plots.append(plot)
        if out_path_pattern:
            plot_filename = out_path_pattern.format(scope='loadings', pci=pc, ext='png')
            with hl.hadoop_open(plot_filename, 'wb') as f:
                get_screenshot_as_png(plot).save(f, format='PNG')
            html = file_html(plot, CDN, 'my plot')
            plot_filename_html = out_path_pattern.format(
                scope='loadings', pci=pc, ext='html'
            )
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def manhattan_loadings(
    iteration,
    gtf,
    loadings,
    title=None,
    size=4,
    hover_fields=None,
    collect_all=False,
    n_divisions=500,
):
    """modify hail manhattan plot"""
    palette = [
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf',
    ]

    # add gene names, p-values, and locus info
    loadings = loadings.annotate(gene_names=gtf[loadings.locus].gene_name)
    pvals = hl.abs(loadings.loadings[iteration])
    locus = loadings.locus

    if hover_fields is None:
        hover_fields = {}

    hover_fields['locus'] = hl.str(locus)
    hover_fields['gene'] = hl.str(loadings.gene_names)

    source_pd = (
        hl.plot.plots._collect_scatter_plot_data(  # pylint: disable=protected-access
            ('_global_locus', locus.global_position()),
            ('_pval', pvals),
            fields=hover_fields,
            n_divisions=None if collect_all else n_divisions,
        )
    )
    source_pd['p_value'] = source_pd['_pval']
    source_pd['_contig'] = [locus.split(':')[0] for locus in source_pd['locus']]

    observed_contigs = set(source_pd['_contig'])
    ref = locus.dtype.reference_genome
    observed_contigs = [
        contig for contig in ref.contigs.copy() if contig in observed_contigs
    ]

    contig_ticks = [
        ref._contig_global_position(contig)  # pylint: disable=protected-access
        + ref.contig_length(contig) // 2
        for contig in observed_contigs
    ]
    color_mapper = CategoricalColorMapper(
        factors=ref.contigs, palette=palette[:2] * int((len(ref.contigs) + 1) / 2)
    )

    p = figure(
        title=title, x_axis_label='Chromosome', y_axis_label='Loadings', width=1000
    )
    (
        p,
        _,
        legend,
        _,
        _,
        _,
    ) = hl.plot.plots._get_scatter_plot_elements(  # pylint: disable=protected-access
        p,
        source_pd,
        x_col='_global_locus',
        y_col='_pval',
        label_cols=['_contig'],
        colors={'_contig': color_mapper},
        size=size,
    )
    legend.visible = False
    p.xaxis.ticker = contig_ticks
    p.xaxis.major_label_overrides = dict(zip(contig_ticks, observed_contigs))
    p.select_one(HoverTool).tooltips = [
        t for t in p.select_one(HoverTool).tooltips if not t[0].startswith('_')
    ]

    return p


def remove_duplicates(x: Iterable) -> List:
    """
    Removes duplicates from a list, keeps order
    """
    return list(dict.fromkeys(x))


if __name__ == '__main__':
    main()  # pylint: disable=E1120
