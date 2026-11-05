"""Publication-style plots for miPyRNA differential expression results.

Provides volcano, MA, heatmap, Venn and PCA plots. All functions return a
matplotlib figure; use :func:`savePlot` to write it to disk. Colours follow the
Okabe-Ito colourblind-safe palette and heatmaps use perceptually uniform maps.

:Author: Naveen Duhan
"""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib import patches
from matplotlib.colors import to_rgba
from sklearn.decomposition import PCA

# Order is (down, not significant, up).
DEG_COLORS = ('#D55E00', '#999999', '#0072B2')

_UP, _DOWN = DEG_COLORS[2], DEG_COLORS[0]

_RC: dict[Any, Any] = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'legend.frameon': False,
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
}

_VENN_COLORS = ['#E69F00', '#56B4E9', '#009E73', '#CC79A7']


def savePlot(fig, path, dpi=300):
    """
    Save a figure for publication. Use a .pdf or .svg path for vector output.

    :param fig: Figure returned by one of the plot functions.
    :param path: Output file path.
    :param dpi: Resolution for raster formats. Defaults to 300.
    """

    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')


def _offset_spines(ax, offset=6):
    ax.spines['left'].set_position(('outward', offset))
    ax.spines['bottom'].set_position(('outward', offset))


def _deg_state(lfc, significant, threshold):
    # 0 = down, 1 = not significant, 2 = up
    return np.where(significant & (lfc >= threshold), 2, np.where(significant & (lfc <= -threshold), 0, 1))


def _scatter_by_state(ax, x, y, state, color, dotsize, markerType, alpha):
    names = ['Down', 'Not significant', 'Up']
    for i in (1, 0, 2):
        mask = state == i
        ax.scatter(x[mask], y[mask], c=color[i], s=dotsize, marker=markerType, alpha=alpha,
                   linewidths=0, rasterized=True, label=f"{names[i]} (n={int(mask.sum())})")


def plotVolcano(degDF=None, comp=None, FOLD=2, pValue=0.05, color=DEG_COLORS, dim=(8, 5), font=14, dotsize=10, markerType='o', alpha=0.5):
    """
    Plot a volcano plot of log2 fold change against -log10 P value.

    :param degDF: All mature_id expression file.
    :param comp: Sample comparison.
    :param FOLD: Fold change threshold. Defaults to 2.
    :param pValue: P value threshold. Defaults to 0.05.
    :param color: Colors for (down, not significant, up).
    :param dim: Dimensions of the plot. Defaults to (8, 5).
    :param font: Axis label font size. Defaults to 14.
    :param dotsize: Dotsize on plot. Defaults to 10.
    :param markerType: Marker shape. Defaults to 'o'.
    :param alpha: Transparency of points. Defaults to 0.5.
    :returns: matplotlib figure.
    :raises ValueError: If the comparison columns are missing or hold no data.
    """

    if degDF is None or comp is None:
        raise ValueError("degDF and comp are required")

    PVAL = f"pvalue({comp})"
    LFC = f"logFC({comp})"

    missing = [c for c in (PVAL, LFC) if c not in degDF.columns]
    if missing:
        raise ValueError(f"Columns not found for comparison {comp!r}: {missing}")

    data = degDF[[LFC, PVAL]].dropna()
    if data.empty:
        raise ValueError(f"No data to plot for comparison {comp!r}")

    lfc = data[LFC].to_numpy()
    pval = data[PVAL].to_numpy()
    y = -np.log10(np.clip(pval, np.finfo(float).tiny, None))
    threshold = np.log2(FOLD)
    state = _deg_state(lfc, pval <= pValue, threshold)

    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=dim, constrained_layout=True)
        _scatter_by_state(ax, lfc, y, state, color, dotsize, markerType, alpha)

        ax.axvline(-threshold, color='#555555', linestyle='--', linewidth=0.8)
        ax.axvline(threshold, color='#555555', linestyle='--', linewidth=0.8)
        ax.axhline(-np.log10(pValue), color='#555555', linestyle='--', linewidth=0.8)

        ax.set_xlabel(r'$\log_{2}$(fold change)', fontsize=font)
        ax.set_ylabel(r'$-\log_{10}$(P value)', fontsize=font)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=3, markerscale=2)
        _offset_spines(ax)

    return fig


def plotMA(degDF=None, countDF=None, comp=None, FOLD=2, FDR=0.05, color=DEG_COLORS, font=14, dim=(8, 5), dotsize=8, markerType='o', alpha=0.5):
    """
    Plot an MA plot of log2 fold change against log2 mean count.

    :param degDF: All mature_id expression file.
    :param countDF: Raw counts file with one column per sample.
    :param comp: Sample comparison, written as 'group1-group2'.
    :param FOLD: Fold change threshold. Defaults to 2.
    :param FDR: FDR threshold. Defaults to 0.05.
    :param color: Colors for (down, not significant, up).
    :param font: Axis label font size. Defaults to 14.
    :param dim: Dimensions of the plot. Defaults to (8, 5).
    :param dotsize: Dotsize on plot. Defaults to 8.
    :param markerType: Marker shape. Defaults to 'o'.
    :param alpha: Transparency of points. Defaults to 0.5.
    :returns: matplotlib figure.
    :raises ValueError: If the comparison columns are missing or hold no data.
    """

    if degDF is None or countDF is None or comp is None:
        raise ValueError("degDF, countDF and comp are required")

    LFC = f"logFC({comp})"
    PVAL = f"FDR({comp})"

    missing = [c for c in (PVAL, LFC) if c not in degDF.columns]
    if missing:
        raise ValueError(f"Columns not found for comparison {comp!r}: {missing}")

    group1, group2 = comp.split("-")[:2]
    mean1 = countDF.filter(regex=group1, axis=1).mean(axis=1)
    mean2 = countDF.filter(regex=group2, axis=1).mean(axis=1)

    final = pd.concat([degDF[[LFC, PVAL]], mean1.rename('first'), mean2.rename('second')], axis=1).fillna(0)
    final['mean'] = final[['first', 'second']].mean(axis=1)
    final = final.loc[final['mean'] > 0]

    if final.empty:
        raise ValueError(f"No data to plot for comparison {comp!r}")

    lfc = final[LFC].to_numpy()
    x = np.log2(final['mean'].to_numpy())
    state = _deg_state(lfc, final[PVAL].to_numpy() <= FDR, np.log2(FOLD))

    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=dim, constrained_layout=True)
        _scatter_by_state(ax, x, lfc, state, color, dotsize, markerType, alpha)

        ax.axhline(0, color='#555555', linestyle='--', linewidth=0.8)

        ax.set_xlabel(r'$\log_{2}$(mean count)', fontsize=font)
        ax.set_ylabel(r'$\log_{2}$(fold change)', fontsize=font)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=3, markerscale=2)
        _offset_spines(ax)

    return fig


def plotHeatmap(degDF=None, combinations=None, num=50, figdim=(12, 10), extraColumns=False, type='counts'):
    """
    Plot a heatmap of the top and bottom mature_ids by fold change or counts.

    Fold change heatmaps use a diverging map centred at zero. Count heatmaps
    show log2(count + 1) on a perceptually uniform map.

    :param degDF: All mature_id expression file or Counts file.
    :param combinations: All sample combinations.
    :param num: Total number of mature_ids to plot. Default is 50 (25 up and 25 down).
    :param figdim: Figure dimensions.
    :param extraColumns: Drop the 'Name' and 'Description' columns first.
    :param type: Heatmap to create on counts/degs.
    :returns: Tuple of (figure, axes).
    :raises ValueError: If type is not 'degs' or 'counts'.
    """

    if degDF is None or combinations is None:
        raise ValueError("degDF and combinations are required")

    if extraColumns:
        degDF = degDF.drop(columns=['Name', 'Description'])

    degDF = degDF.set_index('mature_id')

    if type.lower() == 'degs':
        degDF = degDF.fillna(0)
        com2 = [f"logFC({i})" for i in combinations]
        topmature_id = degDF.nlargest(int(num / 2), com2)
        botmature_id = degDF.nsmallest(int(num / 2), com2)
        fin = pd.concat([topmature_id, botmature_id])[com2]
        fin.columns = list(combinations)
        heat: dict[str, Any] = {'cmap': 'RdBu_r', 'center': 0}
        label = r'$\log_{2}$(fold change)'

    elif type.lower() == 'counts':
        degDF = degDF.replace(0, np.nan).dropna(axis=1, how="all")
        topmature_id = degDF.nlargest(int(num / 2), degDF.columns)
        botmature_id = degDF.nsmallest(int(num / 2), degDF.columns)
        fin = np.log2(pd.concat([topmature_id, botmature_id])[degDF.columns] + 1)
        heat = {'cmap': 'viridis'}
        label = r'$\log_{2}$(count + 1)'

    else:
        raise ValueError(f"type must be 'degs' or 'counts', got {type!r}")

    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=figdim, constrained_layout=True)
        sns.heatmap(fin, ax=ax, cbar=True, cbar_kws={'shrink': .30, 'label': label}, **heat)
        ax.set_ylabel('')

    return fig, ax


def _insertEllipse(ax, x, y, w, h, a, fillcolor):
    ax.add_patch(patches.Ellipse(xy=(x, y), width=w, height=h, angle=a, linewidth=2, color=fillcolor))


def _insertText(ax, x, y, text, fontsize=None, col="black", ha="center", va="center", fontweight=600):
    ax.text(x, y, text, horizontalalignment=ha, verticalalignment=va, fontsize=fontsize, fontweight=fontweight, color=col)


def _mature_idrateCollection(data=None):
    if data is None:
        raise ValueError("data is required")

    N = len(data)
    sets = [set(d) for d in data]
    union = set().union(*sets)
    collections = {}
    for n in range(1, 2 ** N):
        key = format(n, f'0{N}b')
        value = set(union)
        for i in range(N):
            value = value & sets[i] if key[i] == '1' else value - sets[i]
        collections[key] = value
    return collections


# ellipses (x, y, width, height, angle), label centres in key order, comparison names (x, y, ha, va)
_VENN_LAYOUT = {
    2: {
        'ellipses': [(0.375, 0.3, 0.5, 0.5, 0.0), (0.625, 0.3, 0.5, 0.5, 0.0)],
        'points': [(0.74, 0.30), (0.26, 0.30), (0.50, 0.30)],
        'names': [(0.20, 0.56, 'right', 'bottom'), (0.80, 0.56, 'left', 'bottom')],
    },
    3: {
        'ellipses': [(0.333, 0.633, 0.5, 0.5, 0.0), (0.666, 0.633, 0.5, 0.5, 0.0), (0.500, 0.310, 0.5, 0.5, 0.0)],
        'points': [(0.50, 0.27), (0.73, 0.65), (0.61, 0.46), (0.27, 0.65), (0.39, 0.46), (0.50, 0.65), (0.50, 0.51)],
        'names': [(0.15, 0.87, 'right', 'bottom'), (0.85, 0.87, 'left', 'bottom'), (0.50, 0.02, 'left', 'top')],
    },
    4: {
        'ellipses': [(0.350, 0.400, 0.72, 0.45, 140.0), (0.450, 0.500, 0.72, 0.45, 140.0),
                     (0.544, 0.500, 0.72, 0.45, 40.0), (0.644, 0.400, 0.72, 0.45, 40.0)],
        'points': [(0.85, 0.42), (0.66, 0.72), (0.77, 0.59), (0.32, 0.72), (0.69, 0.30), (0.50, 0.66), (0.65, 0.50), (0.14, 0.42),
                   (0.50, 0.17), (0.30, 0.30), (0.38, 0.25), (0.23, 0.59), (0.63, 0.26), (0.35, 0.50), (0.50, 0.38)],
        'names': [(0.13, 0.18, 'right', 'center'), (0.18, 0.83, 'right', 'bottom'),
                  (0.82, 0.83, 'left', 'bottom'), (0.87, 0.18, 'left', 'top')],
    },
}


def plotVenn(DEGFile=None, FOLD=2, comparisons=None, degLabel="", fontsize=14, figsize=(12, 12), dpi=300):
    """
    Plot a Venn diagram of differentially expressed mature_ids across 2-4 comparisons.

    :param DEGFile: Filtered deg excel file containing samples sheet-wise.
    :param FOLD: Fold change threshold. Defaults to 2.
    :param comparisons: Comparison list with 2 to 4 entries.
    :param degLabel: "total" for one count per region, otherwise up (blue) and down (orange) counts.
    :param fontsize: Font size. Defaults to 14.
    :param figsize: Figure size. Defaults to (12, 12).
    :param dpi: Figure DPI resolution. Defaults to 300.
    :returns: matplotlib figure.
    :raises ValueError: If the number of comparisons is not between 2 and 4.
    """

    if DEGFile is None or comparisons is None:
        raise ValueError("DEGFile and comparisons are required")

    n = len(comparisons)
    if n not in _VENN_LAYOUT:
        raise ValueError(f"Venn diagrams need 2 to 4 comparisons, got {n}")

    total, up, down = [], [], []
    threshold = np.log2(FOLD)

    for com in comparisons:
        df = pd.read_excel(DEGFile, sheet_name=com)
        lfc = df[f"logFC({com})"]
        total.append(df['mature_id'])
        up.append(df.loc[lfc >= threshold, 'mature_id'])
        down.append(df.loc[lfc <= -threshold, 'mature_id'])

    layout = _VENN_LAYOUT[n]
    keys = [format(i, f'0{n}b') for i in range(1, 2 ** n)]

    with plt.rc_context(_RC):
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_axis_off()
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)

        for e, c in zip(layout['ellipses'], _VENN_COLORS):
            _insertEllipse(ax, e[0], e[1], e[2], e[3], e[4], to_rgba(c, 0.4))

        if degLabel == "total":
            counts = {k: str(len(v)) for k, v in _mature_idrateCollection(total).items()}
            for (x, y), k in zip(layout['points'], keys):
                _insertText(ax, x, y, counts[k], fontsize=fontsize)
        else:
            up_counts = {k: str(len(v)) for k, v in _mature_idrateCollection(up).items()}
            down_counts = {k: str(len(v)) for k, v in _mature_idrateCollection(down).items()}
            for (x, y), k in zip(layout['points'], keys):
                _insertText(ax, x, y, up_counts[k], col=_UP, fontsize=fontsize)
                _insertText(ax, x, y - 0.03, down_counts[k], col=_DOWN, fontsize=fontsize)
            _insertText(ax, 0.02, 0.04, 'Up', col=_UP, fontsize=fontsize, ha='left')
            _insertText(ax, 0.02, 0.01, 'Down', col=_DOWN, fontsize=fontsize, ha='left')

        for name, (x, y, ha, va) in zip(comparisons, layout['names']):
            _insertText(ax, x, y, name, fontsize=fontsize, ha=ha, va=va)

    return fig


def pcaPlot(ncountdf=None, legends=False, fontsize=14, figsize=(12, 12), dpi=300):
    """
    Plot the first two principal components of normalized counts.

    Axis labels give the percentage of variance explained by each component.

    :param ncountdf: Normalized counts file with mature_ids as rows and samples as columns.
    :param legends: Show a legend of sample names. Defaults to False.
    :param fontsize: Font size. Defaults to 14.
    :param figsize: Figure size. Defaults to (12, 12).
    :param dpi: Figure DPI resolution. Defaults to 300.
    :returns: matplotlib figure.
    """

    if ncountdf is None:
        raise ValueError("ncountdf is required")

    ndf = ncountdf.T

    pca = PCA(n_components=2)
    pca_df = pd.DataFrame(pca.fit_transform(ndf), columns=['PC1', 'PC2'], index=ndf.index)
    var1, var2 = pca.explained_variance_ratio_ * 100

    colors = sns.color_palette('colorblind', n_colors=len(pca_df.index))

    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)

        for i, sample in enumerate(pca_df.index):
            ax.scatter(pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i], color=colors[i], s=60,
                       edgecolors='white', linewidths=0.5, label=sample)

        texts = [ax.text(pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i], sample, ha='right', va='bottom')
                 for i, sample in enumerate(pca_df.index)]
        adjust_text(texts, ax=ax, arrowprops={'arrowstyle': '-', 'color': 'k', 'lw': 0.5})

        if legends:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Samples', ncol=len(pca_df.index) // 15 + 1)

        ax.set_xlabel(f'PC1 ({var1:.1f}% variance)', fontsize=fontsize)
        ax.set_ylabel(f'PC2 ({var2:.1f}% variance)', fontsize=fontsize)
        _offset_spines(ax)

    return fig
