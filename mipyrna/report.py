#!/usr/bin/env python3

"""Report Module

This module writes a PDF summary and a metadata JSON for a miPyRNA run.

Features:
    - Run summary, input samples and the settings used
    - Per-sample read statistics from trimming, collapsing and alignment
    - Known miRNAs by database and the top novel miRNA candidates
    - Differential expression summary, top miRNAs and volcano plots
    - Target and enrichment counts
    - Figures: read lengths, database and score summaries, PCA, sample correlation, volcano, MA and heatmap
    - Links to every output directory and file, software versions and model checksum

Configuration:
    Input is the report directory and the metadata dictionary returned by the
    workflow. Sections whose data is missing are left out. Output is
    ``miPyRNA_report.pdf`` and ``workflow_metadata.json`` in the report directory.

Dependencies:
    - reportlab: PDF layout
    - matplotlib: Figures
    - pandas: Table handling

Classes / Functions / Exceptions:
    - generate_report: Write the PDF report and metadata files.

:Author: Naveen Duhan
"""

import collections
import gzip
import io
import itertools
import json
import os
import re
from datetime import datetime
from importlib import metadata as importlib_metadata
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log='report')

NAVY = colors.HexColor("#1f3a5f")
LIGHT = colors.HexColor("#eef2f7")
GRID = colors.HexColor("#c5ccd6")

STYLES = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=STYLES["BodyText"], fontSize=9.5, leading=13.5, spaceAfter=6)
CELL = ParagraphStyle("cell", parent=STYLES["BodyText"], fontSize=8, leading=10)
CELLB = ParagraphStyle("cellb", parent=CELL, fontName="Helvetica-Bold", textColor=colors.white)
H1 = ParagraphStyle("h1", parent=STYLES["Heading1"], fontSize=14, textColor=NAVY, spaceBefore=14, spaceAfter=6)
H2 = ParagraphStyle("h2", parent=STYLES["Heading2"], fontSize=11, textColor=NAVY, spaceBefore=8, spaceAfter=4)
CAPTION = ParagraphStyle("caption", parent=BODY, fontSize=8.5, leading=11, textColor=colors.HexColor("#444444"))
TITLE = ParagraphStyle("title", parent=STYLES["Title"], fontSize=24, textColor=NAVY, alignment=0, spaceAfter=4)
SUB = ParagraphStyle("sub", parent=BODY, fontSize=11, textColor=colors.HexColor("#555555"))

SIG_FDR = 0.05
SIG_LOGFC = 1.0

_figures = itertools.count(1)


def _version():
    """Get the installed miPyRNA version.

    Returns:
        str: Version string, or ``'unknown'`` if the package is not installed.
    """
    try:
        return importlib_metadata.version("mipyrna")
    except importlib_metadata.PackageNotFoundError:
        return "unknown"


def _read_table(path):
    """Read a result table, choosing the reader from the file extension.

    Args:
        path (str): Path to a ``.tsv``, ``.txt``, ``.csv`` or ``.xlsx`` file, may be None.

    Returns:
        pandas.DataFrame or None: The table, or None if the file is missing or unreadable.
    """
    if not path or not os.path.exists(path):
        return None
    try:
        if path.endswith(".xlsx"):
            return pd.read_excel(path)
        return pd.read_csv(path, sep="," if path.endswith(".csv") else "\t")
    except Exception as exc:
        log.warning(f"Could not read {path}: {exc}")
        return None


def _num(value):
    """Format a number with thousands separators.

    Args:
        value: An integer or float, or None.

    Returns:
        str: The formatted number, or ``'-'`` for None.
    """
    if value is None:
        return "-"
    return f"{value:,}" if isinstance(value, (int, np.integer)) else f"{value:,.2f}"


def _table(rows, widths=None, header=True):
    """Build a styled table that wraps long cells.

    Args:
        rows (list): List of rows, the first being the header when ``header`` is True.
        widths (list, optional): Column widths in points. Defaults to equal widths across the page.
        header (bool): Whether the first row is a header. Defaults to True.

    Returns:
        reportlab.platypus.Table: The table.
    """
    cells = [
        [Paragraph(str(c) if str(c).startswith("<a ") else escape(str(c)), CELLB if header and i == 0 else CELL) for c in row]
        for i, row in enumerate(rows)
    ]
    widths = widths or [17 * cm / len(rows[0])] * len(rows[0])
    table = Table(cells, colWidths=widths, repeatRows=1 if header else 0)
    style = [
        ("GRID", (0, 0), (-1, -1), 0.4, GRID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1), [colors.white, LIGHT]),
    ]
    if header:
        style.append(("BACKGROUND", (0, 0), (-1, 0), NAVY))
    table.setStyle(TableStyle(style))
    return table


def _frame(df, columns, limit=10, digits=3):
    """Turn the first rows of a DataFrame into table rows.

    Args:
        df (pandas.DataFrame): Source table.
        columns (list): Preferred columns; those not in ``df`` are skipped.
        limit (int): Maximum number of rows. Defaults to 10.
        digits (int): Decimals for floats. Defaults to 3.

    Returns:
        list: Header row followed by the data rows.
    """
    cols = [c for c in columns if c in df.columns]
    rows = [cols]
    for _i, r in df[cols].head(limit).iterrows():
        rows.append([f"{v:.{digits}g}" if isinstance(v, float) else v for v in r])
    return rows


def _figure(fig, width, caption):
    """Render a matplotlib figure with a numbered caption.

    Args:
        fig (matplotlib.figure.Figure): The figure.
        width (float): Image width in points.
        caption (str): Caption text, without the figure number.

    Returns:
        reportlab.platypus.KeepTogether: The image and its caption.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    buf.seek(0)
    w, h = ImageReader(buf).getSize()
    buf.seek(0)
    return KeepTogether([
        Image(buf, width=width, height=width * h / w),
        Paragraph(f"Figure {next(_figures)}. {caption}", CAPTION),
    ])


def _link(path, base):
    """Make a clickable link, relative to the report directory so the run folder can be moved.

    Args:
        path (str): File or directory to link to.
        base (str): Directory the link is relative to.

    Returns:
        str: Markup for a link, or the plain relative path if ``path`` does not exist.
    """
    rel = os.path.relpath(path, base)
    if not os.path.exists(path):
        return escape(rel)
    return f'<a href="{escape(rel)}" color="#1a56b0">{escape(rel)}</a>'


def _log_stats(path, pattern):
    """Extract numbers from a tool log with a regular expression.

    Args:
        path (str): Log file path.
        pattern (str): Regular expression with one group per number.

    Returns:
        tuple or None: The captured groups, or None if the file or match is missing.
    """
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as fh:
        match = re.search(pattern, fh.read())
    return match.groups() if match else None


def _collapsed_stats(path):
    """Count the unique sequences, total reads and read lengths of a collapsed FASTA.

    Args:
        path (str): Collapsed FASTA whose headers end in ``_x<count>``.

    Returns:
        tuple or None: ``(unique sequences, total reads, Counter of read length to reads)``,
        or None if the file is missing.
    """
    if not path or not os.path.exists(path):
        return None
    unique = total = 0
    lengths = collections.Counter()
    count = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                count = int(line.rsplit("_x", 1)[-1])
                unique += 1
                total += count
            else:
                lengths[len(line.strip())] += count
    return unique, total, lengths


def _read_length_figure(lengths):
    """Plot the read length distribution of each sample.

    Args:
        lengths (dict): Sample name to Counter of read length to reads.

    Returns:
        matplotlib.figure.Figure: Share of reads at each length from 18 to 30 nt.
    """
    fig = Figure(figsize=(6.5, 3))
    ax = fig.subplots()
    xs = list(range(18, 31))
    for name, counter in lengths.items():
        total = sum(counter.values()) or 1
        ax.plot(xs, [100 * counter.get(x, 0) / total for x in xs], marker="o", ms=3, lw=1, label=name)
    ax.set_xlabel("Read length (nt)")
    ax.set_ylabel("Reads (%)")
    ax.set_xticks(xs)
    ax.legend(fontsize=6, ncol=2, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def _style(fig):
    """Remove the top and right frame lines and tighten the layout.

    Args:
        fig (matplotlib.figure.Figure): The figure.
    """
    for ax in fig.axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()


def _logcpm(counts_df):
    """Convert raw counts to log2 counts per million.

    Args:
        counts_df (pandas.DataFrame): miRNA ID column followed by one count column per sample.

    Returns:
        pandas.DataFrame: log2(CPM + 1) indexed by miRNA ID.
    """
    mat = counts_df.set_index(counts_df.columns[0])
    return np.log2(mat / mat.sum() * 1e6 + 1)


def _discovery_figure(known, novel, summary):
    """Plot known miRNAs per database, mature length and model score.

    Args:
        known (pandas.DataFrame or None): Known miRNA table.
        novel (pandas.DataFrame or None): Novel miRNA table.
        summary (pandas.DataFrame or None): Known miRNAs per database.

    Returns:
        matplotlib.figure.Figure or None: The figure, or None if there is nothing to plot.
    """
    sets = [(df, label, color) for df, label, color in ((known, "Known", "#1f3a5f"), (novel, "Novel", "#c0392b")) if df is not None and len(df)]
    if not sets:
        return None
    has_db = summary is not None and len(summary) > 0
    fig = Figure(figsize=(7.5, 2.6))
    axes = list(np.atleast_1d(fig.subplots(1, 2 + int(has_db))))
    if has_db and summary is not None:
        ax = axes.pop(0)
        ax.barh(summary["database"], summary["miRNAs"], color="#1f3a5f")
        ax.invert_yaxis()
        ax.set_xlabel("Known miRNAs")
    length, score = axes
    for df, label, color in sets:
        if "mature_length" in df.columns:
            length.hist(df["mature_length"].dropna(), bins=range(18, 29), alpha=0.7, color=color, label=label)
        if "score" in df.columns:
            score.hist(df["score"].dropna(), bins=20, alpha=0.7, color=color, label=label)
    length.set_xlabel("Mature length (nt)")
    length.set_ylabel("miRNAs")
    length.legend(fontsize=7, frameon=False)
    score.set_xlabel("Model score")
    _style(fig)
    return fig


def _pca_figure(counts_df, groups):
    """Plot a PCA of the samples and their correlation.

    Args:
        counts_df (pandas.DataFrame): miRNA ID column followed by one count column per sample.
        groups (dict): Sample name to group.

    Returns:
        matplotlib.figure.Figure or None: The figure, or None with fewer than three samples or miRNAs.
    """
    logc = _logcpm(counts_df)
    logc = logc[logc.var(axis=1) > 0]
    if logc.shape[1] < 3 or len(logc) < 3:
        return None
    x = logc.T.to_numpy()
    x = x - x.mean(axis=0)
    u, s, _vt = np.linalg.svd(x, full_matrices=False)
    pcs = u * s
    var = 100 * s**2 / (s**2).sum()
    palette = ["#1f3a5f", "#c0392b", "#2e8b57", "#d4a017", "#6a3d9a", "#7f7f7f"]
    names = list(logc.columns)
    fig = Figure(figsize=(7.5, 3.2))
    pca, heat = fig.subplots(1, 2)
    for i, group in enumerate(sorted({groups.get(n, "-") for n in names})):
        idx = [j for j, n in enumerate(names) if groups.get(n, "-") == group]
        pca.scatter(pcs[idx, 0], pcs[idx, 1], s=25, color=palette[i % len(palette)], label=group)
    pca.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    pca.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    pca.legend(fontsize=7, frameon=False)
    corr = logc.corr().to_numpy()
    im = heat.imshow(corr, cmap="viridis", vmin=corr.min(), vmax=1)
    heat.set_xticks(range(len(names)))
    heat.set_xticklabels(names, rotation=90, fontsize=6)
    heat.set_yticks(range(len(names)))
    heat.set_yticklabels(names, fontsize=6)
    fig.colorbar(im, ax=heat, fraction=0.046)
    _style(fig)
    return fig


def _deg_figure(deg, comparison):
    """Plot the volcano and MA plots of one comparison.

    Args:
        deg (pandas.DataFrame): DESeq2 table with ``baseMean``, ``logFC`` and ``FDR`` columns suffixed with the comparison.
        comparison (str): Comparison name, for example ``'ago1-Col0'``.

    Returns:
        matplotlib.figure.Figure: Volcano plot and MA plot with significant miRNAs highlighted.
    """
    d = deg[[f"baseMean({comparison})", f"logFC({comparison})", f"FDR({comparison})"]].dropna()
    base, lfc, fdr = (d.iloc[:, i].to_numpy() for i in range(3))
    sig = (fdr < SIG_FDR) & (np.abs(lfc) >= SIG_LOGFC)
    fig = Figure(figsize=(7, 3))
    volcano, ma = fig.subplots(1, 2)
    y = -np.log10(np.clip(fdr, 1e-300, None))
    volcano.scatter(lfc[~sig], y[~sig], s=6, c="#9aa5b4")
    volcano.scatter(lfc[sig], y[sig], s=8, c="#c0392b")
    volcano.axhline(-np.log10(SIG_FDR), ls="--", lw=0.6, c="grey")
    volcano.set_xlabel("log2 fold change")
    volcano.set_ylabel("-log10 FDR")
    mean = np.log10(base + 1)
    ma.scatter(mean[~sig], lfc[~sig], s=6, c="#9aa5b4")
    ma.scatter(mean[sig], lfc[sig], s=8, c="#c0392b")
    ma.axhline(0, ls="--", lw=0.6, c="grey")
    ma.set_xlabel("log10 mean count")
    ma.set_ylabel("log2 fold change")
    _style(fig)
    return fig


def _heatmap_figure(counts_df, names, groups):
    """Plot a heatmap of row-scaled expression for selected miRNAs.

    Args:
        counts_df (pandas.DataFrame): miRNA ID column followed by one count column per sample.
        names (list): miRNA IDs to show.
        groups (dict): Sample name to group, used to order the columns.

    Returns:
        matplotlib.figure.Figure: Heatmap of z-scored log2 CPM.
    """
    logc = _logcpm(counts_df)
    cols = sorted(logc.columns, key=lambda c: groups.get(c, ""))
    z = logc.loc[[n for n in names if n in logc.index], cols]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    fig = Figure(figsize=(5.5, 0.2 * len(z) + 1.4))
    ax = fig.subplots()
    im = ax.imshow(z.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=6)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(z.index, fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.04, label="z-score")
    fig.tight_layout()
    return fig


def _decorate(canvas, doc):
    """Draw the page footer.

    Args:
        canvas (reportlab.pdfgen.canvas.Canvas): Page canvas.
        doc (reportlab.platypus.SimpleDocTemplate): The document.
    """
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#777777"))
    canvas.drawString(2 * cm, 1.2 * cm, f"miPyRNA {_version()}")
    canvas.drawRightString(A4[0] - 2 * cm, 1.2 * cm, f"Page {doc.page}")
    canvas.restoreState()


def _overview(meta, params, counts):
    """Build the title block and the run summary.

    Args:
        meta (dict): Workflow metadata.
        params (dict): Run settings.
        counts (dict): Counts of known, novel, DE, target and enrichment records.

    Returns:
        list: Report flowables.
    """
    samples = meta.get("samples", {})
    groups = sorted({v[0] for v in samples.values()})
    elapsed = meta.get("elapsed_seconds")
    rows = [
        ["Species", f"{params.get('species', '-')} ({params.get('species_type', '-')})"],
        ["Samples", f"{len(samples)} in {len(groups)} group(s): {', '.join(groups) or '-'}"],
        ["Comparisons", ", ".join(meta.get("comparisons", [])) or "-"],
        ["Known miRNAs", _num(counts["known"])],
        ["Novel miRNAs", _num(counts["novel"])],
        ["Significantly changed miRNAs", f"{_num(counts['deg'])} (FDR < {SIG_FDR}, |log2FC| >= {SIG_LOGFC:g})" if counts["deg"] is not None else "-"],
        ["miRNA-target interactions", _num(counts["targets"])],
        ["Enriched terms", _num(counts["enrichment"])],
        ["Run time", f"{elapsed / 3600:.2f} h" if elapsed else "-"],
    ]
    return [
        Paragraph("miPyRNA analysis report", TITLE),
        Paragraph(f"Generated {datetime.now():%Y-%m-%d %H:%M} with miPyRNA {_version()}", SUB),
        Spacer(1, 10),
        Paragraph("Summary", H1),
        _table(rows, [5.5 * cm, 11.5 * cm], header=False),
    ]


def _genome_stats(path):
    """Count the sequences and bases of a genome FASTA.

    Args:
        path (str): FASTA path, plain or gzipped, may be None.

    Returns:
        tuple or None: ``(sequences, bases)``, or None if the file is missing.
    """
    if not path or not os.path.exists(path):
        return None
    opener = gzip.open if path.endswith(".gz") else open
    seqs = bases = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                seqs += 1
            else:
                bases += len(line.strip())
    return seqs, bases


def _inputs(meta, params):
    """Build the input data, samples and settings sections.

    Args:
        meta (dict): Workflow metadata.
        params (dict): Run settings.

    Returns:
        list: Report flowables.
    """
    genome = params.get("genome")
    stats = _genome_stats(genome)
    trimmer = params.get("trimmer")
    adapter = params.get("adapter")
    trimming = "skipped" if not trimmer else f"{trimmer}, " + (f"adapter {adapter}" if adapter else "adapter detected automatically")
    databases = {"plants": "miRBase, PmiREN", "animals": "miRBase, MirGeneDB"}.get(params.get("species_type"), "miRBase")
    rows = [
        ["Species", f"{params.get('species', '-')} ({params.get('species_type', '-')})"],
        ["Reference genome", genome or "-"],
        ["Genome size", f"{_num(stats[1])} bp in {_num(stats[0])} sequences" if stats else "-"],
        ["Sample sheet", params.get("input_file", "-")],
        ["Reads directory", params.get("input_path", "-")],
        ["Library type", "paired-end" if params.get("paired") else "single-end"],
        ["Trimming", trimming],
        ["Known miRNA databases", databases],
    ]
    out: list = [Paragraph("Input data", H1), _table(rows, [4.5 * cm, 12.5 * cm], header=False)]
    samples = meta.get("samples", {})
    if samples:
        rows = [["Sample", "Group", "Input file(s)"]]
        rows += [[k, v[0], ", ".join(os.path.basename(f) for f in v[1:])] for k, v in samples.items()]
        out += [Paragraph("Samples", H2), _table(rows, [4.5 * cm, 4 * cm, 8.5 * cm])]
    shown = {"species", "species_type", "genome", "input_file", "input_path", "paired", "trimmer", "adapter"}
    rest = [[k, "-" if v is None else str(v)] for k, v in params.items() if k not in shown]
    if rest:
        out += [Paragraph("Settings", H2), _table([["Setting", "Value"]] + rest, [5.5 * cm, 11.5 * cm])]
    return out


def _reads(meta, dirs):
    """Build the read statistics section from the trimming, collapse and alignment outputs.

    Args:
        meta (dict): Workflow metadata.
        dirs (dict): Stage name to directory.

    Returns:
        list: Report flowables.
    """
    samples = meta.get("samples", {})
    if not samples:
        return []
    rows = [["Sample", "Raw reads", "After trimming", "Unique sequences", "Aligned reads (%)"]]
    lengths = {}
    for name in samples:
        trim = _log_stats(
            os.path.join(dirs.get("quality", ""), "trimmomatic_results", "logs", f"trimmomatic_{name}.err"),
            r"Input Reads: (\d+) Surviving: (\d+) \(([\d.]+)%\)",
        )
        collapsed = _collapsed_stats(os.path.join(dirs.get("alignment", ""), "collapsed_reads", f"{name}_collapsed_reads.txt"))
        aligned = _log_stats(
            os.path.join(dirs.get("alignment", ""), "bowtie_results", "logs", f"bowtie_{name}.err"),
            r"reads processed: (\d+)\s+# reads with at least one alignment: (\d+) \(([\d.]+)%\)",
        )
        if collapsed:
            lengths[name] = collapsed[2]
        rows.append([
            name,
            _num(int(trim[0])) if trim else "-",
            f"{_num(int(trim[1]))} ({trim[2]}%)" if trim else "-",
            _num(collapsed[0]) if collapsed else "-",
            f"{_num(int(aligned[1]))} ({aligned[2]}%)" if aligned else "-",
        ])
    out = [
        Paragraph("Read processing", H1),
        Paragraph(
            "Reads were trimmed of the 3' adapter, collapsed to unique sequences with their read counts, and aligned to "
            "the genome with Bowtie. Aligned reads are counted over the unique sequences, so the percentage is the share "
            "of unique sequences that mapped, not the share of all reads. Raw and trimmed read counts are available "
            "for Trimmomatic only.",
            BODY,
        ),
        _table(rows, [3.6 * cm, 3 * cm, 3.8 * cm, 3 * cm, 3.6 * cm]),
    ]
    if lengths:
        out += [
            Spacer(1, 8),
            _figure(_read_length_figure(lengths), 14 * cm, "Read length distribution of the collapsed reads, weighted by read count."),
        ]
    return out


def _discovery(meta, params, base, known, novel):
    """Build the known and novel miRNA section.

    Args:
        meta (dict): Workflow metadata.
        params (dict): Run settings.
        base (str): Report directory, used for links.
        known (pandas.DataFrame or None): Known miRNA table.
        novel (pandas.DataFrame or None): Novel miRNA table.

    Returns:
        list: Report flowables.

    """
    out: list = [
        Paragraph("miRNA discovery", H1),
        Paragraph(
            "Aligned reads were clustered into candidate loci. Precursors around each locus were folded with ViennaRNA "
            f"and scored with the bundled {params.get('species_type', '')} neural network model. A candidate is kept if "
            f"its score is at least {params.get('score_cutoff', '-')}, its randfold p-value is at most "
            f"{params.get('randfold_pvalue_cutoff', '-')} and it has at least {params.get('min_mature_reads', '-')} "
            "mature reads. Candidates matching a known miRNA (BLAST identity above 94% over more than 19 nt) are reported "
            "as known and the rest as novel. Novel miRNAs are computational predictions and should be validated.",
            BODY,
        ),
    ]
    summary = _read_table(meta.get("database_summary"))
    if summary is not None and len(summary):
        out += [Paragraph("Known miRNAs by database", H2), _table([list(summary.columns)] + summary.values.tolist(), [8 * cm, 4 * cm])]
        out.append(Paragraph("Shared means found in more than one database.", CAPTION))
    fig = _discovery_figure(known, novel, summary)
    if fig is not None:
        out += [Spacer(1, 6), _figure(fig, 16 * cm, "Known miRNAs per database, mature miRNA length and model score of known and novel miRNAs.")]
    stars = [(label, int((df["star_reads"] >= max(1, params.get("min_star_reads") or 1)).sum()), len(df))
             for label, df in (("Known", known), ("Novel", novel)) if df is not None and len(df) and "star_reads" in df.columns]
    if stars:
        rule = "Candidates without star reads were removed." if params.get("require_star") else "Star reads were not required, so candidates without them can be present."
        out += [
            Paragraph("Star read support", H2),
            Paragraph(
                "The star strand is the other arm of the miRNA duplex, and reads from it are evidence of Dicer processing. "
                + " ".join(f"{n:,} of {total:,} {label.lower()} miRNAs have star reads." for label, n, total in stars)
                + f" {rule}",
                BODY,
            ),
        ]
    if known is not None and len(known):
        out += [Paragraph("Known miRNAs (top 10 by score)", H2)]
        top = known.sort_values("score", ascending=False) if "score" in known.columns else known
        out.append(_table(_frame(top, ["Novel_ID", "Known_ID", "Database", "chrom", "strand", "mature_reads", "star_reads", "score"]), None))
    if novel is not None and len(novel): 
        out += [Paragraph("Novel miRNAs (top 10 by score)", H2)]
        top = novel.sort_values("score", ascending=False) if "score" in novel.columns else novel
        out.append(_table(_frame(top, ["Novel_ID", "chrom", "strand", "mature_length", "mature_reads", "star_reads", "score"]), None))
    if meta.get("known_file") or meta.get("novel_file"):
        out.append(Paragraph(
            "Full tables: " + ", ".join(_link(meta[k], base) for k in ("known_file", "novel_file") if meta.get(k)) + ".",
            BODY,
        ))
    return out


def _expression(meta, base, counts_df, deg):
    """Build the quantification and differential expression section.

    Args:
        meta (dict): Workflow metadata.
        base (str): Report directory, used for links.
        counts_df (pandas.DataFrame or None): Raw read count table.
        deg (pandas.DataFrame or None): DESeq2 results.

    Returns:
        list: Report flowables.
    """
    out = []
    if counts_df is not None:
        sample_cols = counts_df.columns[1:]
        rows = [["Sample", "Counted reads", "miRNAs with reads"]]
        rows += [[c, _num(int(counts_df[c].sum())), _num(int((counts_df[c] > 0).sum()))] for c in sample_cols]
        out += [
            Paragraph("Quantification", H1),
            Paragraph(
                f"Collapsed reads were aligned to the {len(counts_df)} known and novel mature miRNAs and counted per sample.",
                BODY,
            ),
            _table(rows, [6 * cm, 5 * cm, 5 * cm]),
        ]
        groups = {k: v[0] for k, v in meta.get("samples", {}).items()}
        fig = _pca_figure(counts_df, groups)
        if fig is not None:
            out += [Spacer(1, 6), _figure(fig, 15 * cm, "Principal component analysis of log2 counts per million (left) and Pearson correlation between samples (right).")]
    if deg is None:
        return out
    comparisons = [c[4:-1] for c in deg.columns if c.startswith("FDR(")]
    out += [
        Paragraph("Differential expression", H1),
        Paragraph(
            f"Differential expression was tested with DESeq2. A miRNA is called significant at FDR below {SIG_FDR} and "
            f"an absolute log2 fold change of at least {SIG_LOGFC:g}. The fold change is the first group relative to the second.",
            BODY,
        ),
    ]
    rows = [["Comparison", "Tested", "Up", "Down"]]
    for comp in comparisons:
        d = deg[[f"logFC({comp})", f"FDR({comp})"]].dropna()
        fdr, lfc = d.iloc[:, 1], d.iloc[:, 0]
        rows.append([comp, _num(len(d)), _num(int(((fdr < SIG_FDR) & (lfc >= SIG_LOGFC)).sum())), _num(int(((fdr < SIG_FDR) & (lfc <= -SIG_LOGFC)).sum()))])
    out.append(_table(rows, [6 * cm, 3.5 * cm, 3.5 * cm, 3.5 * cm]))
    for comp in comparisons:
        out.append(Paragraph(f"Top miRNAs: {escape(comp)}", H2))
        d = deg[[deg.columns[0], f"baseMean({comp})", f"logFC({comp})", f"FDR({comp})"]].dropna().sort_values(f"FDR({comp})")
        d.columns = ["miRNA", "baseMean", "log2FC", "FDR"]
        out.append(_table(_frame(d, list(d.columns)), None))
        out += [Spacer(1, 6), _figure(_deg_figure(deg, comp), 14 * cm, f"Volcano plot (left) and MA plot (right) of {escape(comp)}. Red points are significant.")]
        if counts_df is not None and len(d) > 1:
            groups = {k: v[0] for k, v in meta.get("samples", {}).items()}
            heat = _heatmap_figure(counts_df, d["miRNA"].head(25).tolist(), groups)
            out += [Spacer(1, 6), _figure(heat, 10 * cm, f"Expression of the 25 miRNAs with the lowest FDR in {escape(comp)}, scaled per miRNA.")]
    if meta.get("deg_file"):
        out.append(Paragraph(f"Full table: {_link(meta['deg_file'], base)}.", BODY))
    return out


def _targets(meta, base, counts):
    """Build the target prediction and enrichment section.

    Args:
        meta (dict): Workflow metadata.
        base (str): Report directory, used for links.
        counts (dict): Counts of known, novel, DE, target and enrichment records.

    Returns:
        list: Report flowables.
    """
    if not meta.get("targets_file"):
        return []
    out = [
        Paragraph("Targets and enrichment", H1),
        Paragraph(
            f"miRanda predicted {_num(counts['targets'])} miRNA-target interactions"
            + (f", and g:Profiler found {_num(counts['enrichment'])} enriched terms among the target genes" if meta.get("enrichment_file") else "")
            + ". Predicted targets are computational and are not evidence of regulation.",
            BODY,
        ),
        Paragraph("Files: " + ", ".join(_link(meta[k], base) for k in ("targets_file", "enrichment_file") if meta.get(k)) + ".", BODY),
    ]
    return out


def _files(meta, base):
    """Build the table of output directories and files with links.

    Args:
        meta (dict): Workflow metadata.
        base (str): Report directory, used for links.

    Returns:
        list: Report flowables.
    """
    labels = {
        "quality": "Quality control and trimming",
        "alignment": "Collapsed reads, Bowtie index and alignments",
        "discovery": "Known and novel miRNAs, FASTA files, BLAST results",
        "quantification": "Read counts per miRNA",
        "differential": "DESeq2 results",
        "targets": "miRanda target predictions",
        "annotation": "g:Profiler enrichment",
        "report": "This report and the run metadata",
    }
    rows = [["Directory", "Contents"]]
    rows += [[_link(path, base), labels[k]] for k, path in meta.get("stage_dirs", {}).items() if k in labels]
    if meta.get("outdir"):
        rows.append([_link(os.path.join(meta["outdir"], "mipyrna.log"), base), "Log of the whole run"])
    files = {
        "known_file": "Known miRNAs",
        "novel_file": "Novel miRNAs",
        "counts_file": "Raw read counts",
        "deg_file": "DESeq2 results",
        "targets_file": "miRanda targets",
        "enrichment_file": "g:Profiler enrichment",
    }
    frows = [["File", "Contents"]] + [[_link(meta[k], base), label] for k, label in files.items() if meta.get(k)]
    return [
        Paragraph("Output files", H1),
        Paragraph("Links are relative to this report, so they work when the whole run folder is copied together.", BODY),
        _table(rows, [9 * cm, 8 * cm]),
        Spacer(1, 8),
        _table(frows, [9 * cm, 8 * cm]),
    ]


def _software(meta):
    """Build the software versions and model section.

    Args:
        meta (dict): Workflow metadata with ``versions`` and ``model``.

    Returns:
        list: Report flowables.
    """
    versions = meta.get("versions")
    if not versions:
        return []
    rows = [["Component", "Version"], ["miPyRNA", _version()]] + [[k, "not installed" if v is None else v] for k, v in versions.items()]
    out = [Paragraph("Software", H1), _table(rows, [5.5 * cm, 11.5 * cm])]
    if meta.get("model"):
        out += [
            Spacer(1, 6),
            Paragraph("Bundled model (SHA-256)", H2),
            _table([["Model", "Checksum"]] + [[k, v] for k, v in meta["model"].items()], [4 * cm, 13 * cm]),
        ]
    return out


def generate_report(outdir=".", metadata=None):
    """Write the PDF analysis report and the workflow metadata JSON.

    Args:
        outdir (str): Report directory. Defaults to the current directory.
        metadata (dict, optional): Workflow metadata. Used keys are ``outdir``,
            ``stage_dirs``, ``samples``, ``comparisons``, ``parameters``, ``versions``,
            ``model``, ``elapsed_seconds``, ``database_summary`` and the
            ``known_file``, ``novel_file``, ``counts_file``, ``deg_file``,
            ``targets_file`` and ``enrichment_file`` paths.

    Returns:
        tuple: ``(report_path, metadata_path)``.
    """
    global _figures
    _figures = itertools.count(1)
    os.makedirs(outdir, exist_ok=True)
    report_path = os.path.join(outdir, "miPyRNA_report.pdf")
    meta_path = os.path.join(outdir, "workflow_metadata.json")
    meta = metadata or {}

    with open(meta_path, "w") as jh:
        json.dump(meta, jh, indent=2)

    known = _read_table(meta.get("known_file"))
    novel = _read_table(meta.get("novel_file"))
    counts_df = _read_table(meta.get("counts_file"))
    deg = _read_table(meta.get("deg_file"))
    targets = _read_table(meta.get("targets_file"))
    enrichment = _read_table(meta.get("enrichment_file"))

    n_deg = None
    if deg is not None:
        sig = pd.Series(False, index=deg.index)
        for c in [c[4:-1] for c in deg.columns if c.startswith("FDR(")]:
            sig |= (deg[f"FDR({c})"] < SIG_FDR) & (deg[f"logFC({c})"].abs() >= SIG_LOGFC)
        n_deg = int(sig.sum())
    counts = {
        "known": len(known) if known is not None else 0,
        "novel": len(novel) if novel is not None else 0,
        "deg": n_deg,
        "targets": len(targets) if targets is not None else 0,
        "enrichment": len(enrichment) if enrichment is not None else 0,
    }
    params = meta.get("parameters", {})
    dirs = meta.get("stage_dirs", {})

    story = (
        _overview(meta, params, counts)
        + _inputs(meta, params)
        + _reads(meta, dirs)
        + _discovery(meta, params, outdir, known, novel)
        + _expression(meta, outdir, counts_df, deg)
        + _targets(meta, outdir, counts)
        + _files(meta, outdir)
        + _software(meta)
    )
    doc = SimpleDocTemplate(
        report_path,
        pagesize=A4,
        leftMargin=2 * cm,
        rightMargin=2 * cm,
        topMargin=2 * cm,
        bottomMargin=2 * cm,
        title="miPyRNA analysis report",
        author="miPyRNA",
    )
    doc.build(story, onFirstPage=_decorate, onLaterPages=_decorate)
    log.info(f"Report written to {report_path}")
    return report_path, meta_path
