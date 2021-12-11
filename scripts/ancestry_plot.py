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
    '--meta-tsv',
    'meta_tsv_path',
    callback=utils.get_validation_callback(ext='tsv', must_exist=True),
    required=True,
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
    meta_tsv_path: str,
    out_path_pattern: Optional[str],
    hail_billing: str,  # pylint: disable=unused-argument
):
    local_tmp_dir = utils.init_hail(__file__)

    meta_ht = utils.parse_input_metadata(meta_tsv_path, local_tmp_dir)

    produce_plots(
        eigenvalues_path=eigenvalues_path,
        scores_ht_path=scores_ht_path,
        loadings_ht_path=loadings_ht_path,
        inferred_pop_ht_path=inferred_pop_ht_path,
        meta_ht=meta_ht,
        out_path_pattern=out_path_pattern,
    )


def key_by_external_id(ht, sample_map_ht):
    """
    Assuming ht.s is a CPG id, replaces it with external ID, using map in sample_map_ht
    """
    ht = ht.annotate(internal_id=ht.s).key_by('internal_id')
    ht = (
        ht.annotate(
            s=hl.if_else(
                hl.is_defined(sample_map_ht[ht.internal_id]),
                sample_map_ht[ht.internal_id].external_id,
                ht.internal_id,
            )
        )
        .key_by('s')
        .drop('internal_id')
    )
    return ht


def produce_plots(
    eigenvalues_path: str,
    scores_ht_path: str,
    loadings_ht_path: str,
    inferred_pop_ht_path: str,
    meta_ht: hl.Table,
    out_path_pattern: Optional[str] = None,
    number_of_pcs: Optional[int] = None,
):
    """
    Generate plots in HTML format, write for each PC (of n_pcs) and
    scope ("study", "continental_pop", "subpop", plus for loadings) into
    file paths defined by `out_path_pattern`.
    """
    scores_ht = hl.read_table(scores_ht_path)
    inferred_pop_ht = hl.read_table(inferred_pop_ht_path)

    sample_map_ht = meta_ht.select('external_id')
    scores_ht = key_by_external_id(scores_ht, sample_map_ht).cache()
    provided_pop_ht = key_by_external_id(meta_ht, sample_map_ht).cache()
    inferred_pop_ht = key_by_external_id(inferred_pop_ht, sample_map_ht).cache()

    scores_ht = scores_ht.annotate(
        continental_pop=hl.case()
        .when(
            provided_pop_ht[scores_ht.s].continental_pop != '',
            provided_pop_ht[scores_ht.s].project
            + ' ['
            + provided_pop_ht[scores_ht.s].continental_pop
            + ']',
        )
        .default(
            provided_pop_ht[scores_ht.s].project
            + ' [inferred '
            + inferred_pop_ht[scores_ht.s].pop
            + ']'
        ),
        subpop=hl.case()
        .when(
            provided_pop_ht[scores_ht.s].subpop != '',
            provided_pop_ht[scores_ht.s].project
            + ' ['
            + provided_pop_ht[scores_ht.s].subpop
            + ']',
        )
        .default(
            provided_pop_ht[scores_ht.s].project
            + ' [inferred '
            + inferred_pop_ht[scores_ht.s].pop
            + ']'
        ),
        study=provided_pop_ht[scores_ht.s].project,
    ).cache()

    eigenvalues_ht = hl.import_table(
        eigenvalues_path, no_header=True, types={'f0': hl.tfloat}
    ).f0.collect()
    eigenvalues_ht = pd.to_numeric(eigenvalues_ht)
    variance = np.divide(eigenvalues_ht[1:], float(eigenvalues_ht.sum())) * 100
    variance = variance.round(2)
    number_of_pcs = number_of_pcs or len(eigenvalues_ht) - 1

    plots = []

    sample_names = scores_ht.s.collect()
    studies = scores_ht.study.collect()
    unique_studies = remove_duplicates(studies)
    bg_studies = [s for s in unique_studies if s == 'gnomad']
    fg_studies = [s for s in unique_studies if s != 'gnomad']

    plots.extend(
        _plot_study(
            number_of_pcs,
            variance,
            scores_ht,
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
            scores_ht,
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
            scores_ht,
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
            width=1000,
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
            width=1000,
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
            size=4,
            color=factor_cmap('label', turbo(len(unique_labels)), unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)
        if out_path_pattern:
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
            width=1000,
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
            size=4,
            color=factor_cmap('label', turbo(len(unique_labels)), unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)
        if out_path_pattern:
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
