#!/usr/bin/env python

"""
Plot ancestry PCA analysis results
"""

import logging
import click
import pandas as pd
import hail as hl
from bokeh.io.export import get_screenshot_as_png
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.transform import factor_cmap
from bokeh.plotting import ColumnDataSource, figure
from bokeh.palettes import turbo  # pylint: disable=no-name-in-module
from bokeh.models import CategoricalColorMapper, HoverTool

from joint_calling import utils
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
    '--assigned-pop-ht',
    'assigned_pop_ht_path',
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
    assigned_pop_ht_path: str,
    out_path_pattern: str,
    hail_billing: str,  # pylint: disable=unused-argument
):
    utils.init_hail(__file__)

    produce_plots(
        eigenvalues_path=eigenvalues_path,
        scores_ht_path=scores_ht_path,
        loadings_ht_path=loadings_ht_path,
        assigned_pop_ht_path=assigned_pop_ht_path,
        out_path_pattern=out_path_pattern,
    )


def produce_plots(
    eigenvalues_path: str,
    scores_ht_path: str,
    loadings_ht_path: str,
    assigned_pop_ht_path: str,
    out_path_pattern: str,
):
    """
    Generate plots in png and html formats, write for each PC (of n_pcs) and
    scope ("study", "continental_pop", "subpop", plus for loadings) into
    file paths defined by `out_path_pattern`.
    """
    assigned_pop_ht = hl.read_table(assigned_pop_ht_path)
    scores = hl.read_table(scores_ht_path)
    scores = scores.annotate(study=assigned_pop_ht[scores.s].project)
    sample_names = scores.s.collect()
    labels = scores.study.collect()
    study = list(set(labels))
    tooltips = [('labels', '@label'), ('samples', '@samples')]
    eigenvalues = hl.import_table(
        eigenvalues_path, no_header=True, types={'f0': hl.tfloat}
    )
    eigenvalues = eigenvalues.to_pandas()
    eigenvalues.columns = ['eigenvalue']
    eigenvalues = pd.to_numeric(eigenvalues.eigenvalue)[1:]
    variance = eigenvalues.divide(float(eigenvalues.sum())) * 100
    variance = variance.round(2)

    # Get number of PCs
    number_of_pcs = len(eigenvalues)

    # plot by study
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
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=4,
            color=factor_cmap('label', ['#1b9e77', '#d95f02'], study),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plot_filename = out_path_pattern.format(scope='study', pci=pc2, ext='png')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = out_path_pattern.format(scope='study', pci=pc2, ext='html')
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)

    # plot by continental population
    scores = scores.annotate(continental_pop=assigned_pop_ht[scores.s].continental_pop)
    labels = scores.continental_pop.collect()
    # Change TOB-WGS 'none' values to 'TOB-WGS'
    labels = ['TOB-NFE' if x is None else x for x in labels]
    continental_population = list(set(labels))
    tooltips = [('labels', '@label'), ('samples', '@samples')]

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
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=4,
            color=factor_cmap(
                'label', turbo(len(continental_population)), continental_population
            ),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
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

    # plot by subpopulation
    scores = scores.annotate(subpop=assigned_pop_ht[scores.s].subpop)
    labels = scores.subpop.collect()
    labels = ['TOB-NFE' if x is None else x for x in labels]
    sub_population = list(set(labels))
    tooltips = [('labels', '@label'), ('samples', '@samples')]

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
            )
        )
        plot.circle(
            'x',
            'y',
            alpha=0.5,
            source=source,
            size=4,
            color=factor_cmap('label', turbo(len(sub_population)), sub_population),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plot_filename = out_path_pattern.format(scope='subpop', pci=pc2, ext='png')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = out_path_pattern.format(
            scope='subpop', pci=pc2, ext='html'
        )
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)

    # Plot loadings
    loadings_ht = hl.read_table(loadings_ht_path)
    for i in range(0, number_of_pcs - 1):
        pc = i + 1
        plot = manhattan_loadings(
            pvals=hl.abs(loadings_ht.loadings[i]),
            locus=loadings_ht.locus,
            title='Loadings of PC ' + str(pc),
            collect_all=True,
        )
        plot_filename = out_path_pattern.format(scope='loadings', pci=pc, ext='png')
        with hl.hadoop_open(plot_filename, 'wb') as f:
            get_screenshot_as_png(plot).save(f, format='PNG')
        html = file_html(plot, CDN, 'my plot')
        plot_filename_html = out_path_pattern.format(
            scope='loadings', pci=pc, ext='html'
        )
        with hl.hadoop_open(plot_filename_html, 'w') as f:
            f.write(html)


def manhattan_loadings(
    pvals,
    locus=None,
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
    if locus is None:
        locus = pvals._indices.source.locus  # pylint: disable=protected-access

    ref = locus.dtype.reference_genome

    if hover_fields is None:
        hover_fields = {}

    hover_fields['locus'] = hl.str(locus)

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


if __name__ == '__main__':
    main()  # pylint: disable=E1120
