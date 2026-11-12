"""Differential Expression Module

This module runs differential expression analysis of miRNAs from raw read
counts using R packages through rpy2, and filters the results.

Features:
    - DESeq2 analysis for all or per-comparison subsets of samples
    - edgeR analysis with or without replicates
    - Filtering of DEGs by FDR and fold change
    - Stacked bar plot of up- and down-regulated miRNAs per comparison

Configuration:
    Inputs are a count DataFrame, a target table with a ``sample`` column
    and a list of comparisons written as ``"A-B"``.

Dependencies:
    - R with DESeq2, edgeR and limma
    - rpy2: Python to R bridge
    - pandas, numpy, matplotlib

Classes / Functions / Exceptions:
    - runDESeq2: Differential expression with DESeq2.
    - run_edgeR: Differential expression with edgeR.
    - degFilter: Filter DEGs and summarize them.

:Author: Naveen Duhan
"""

import logging
from functools import cache
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log="diff")


@cache
def _rpy2():
    """Lazily initialize rpy2 so importing this module does not boot embedded R."""
    from rpy2 import robjects
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    from rpy2.robjects import Formula, numpy2ri, pandas2ri
    from rpy2.robjects.conversion import get_conversion, localconverter, py2rpy
    from rpy2.robjects.packages import importr

    rpy2_logger.setLevel(logging.ERROR)
    r_env: Any = robjects.r
    to_df = r_env('function(x) data.frame(x)')
    return robjects, Formula, numpy2ri, pandas2ri, get_conversion, localconverter, py2rpy, importr, r_env, to_df


def _split_contrast(co, levels=None):
    """Split a comparison into its two group names and display label.

    Args:
        co (str or tuple): Comparison written as ``"A-B"``, ``"A_vs_B"`` or ``(A, B)``.
        levels (iterable, optional): Group names present in the design table, used to
            resolve group names that contain hyphens.

    Returns:
        tuple: ``(group_a, group_b, label)``.

    Raises:
        ValueError: If the comparison is ambiguous or does not match two groups.
    """
    if isinstance(co, (tuple, list)) and len(co) == 2:
        c1, c2 = str(co[0]), str(co[1])
        return c1, c2, f"{c1}-{c2}"
    label = str(co)
    if levels:
        known = {str(lv) for lv in levels}
        matches = [
            (label[:idx], label[idx + 1 :])
            for idx, ch in enumerate(label)
            if ch == "-" and label[:idx] in known and label[idx + 1 :] in known
        ]
        if len(matches) == 1:
            return matches[0][0], matches[0][1], label
        if "_vs_" in label:
            left, right = label.split("_vs_", 1)
            if left in known and right in known:
                return left, right, label
        if len(matches) > 1:
            raise ValueError(f"Ambiguous contrast {label!r} for levels {sorted(known)}")
    if "_vs_" in label:
        c1, c2 = label.split("_vs_", 1)
        return c1, c2, label
    if label.count("-") == 1:
        c1, c2 = label.split("-", 1)
        return c1, c2, label
    raise ValueError(f"Cannot resolve contrast {label!r} against design levels")


def runDESeq2(countDF=None, targetFile=None, design='sample',combination=None,  mirna_column='mature_id',  mmg =False, subset=False, lib=None):
    """Run DESeq2 differential expression analysis on raw read counts.

    Args:
        countDF (pandas.DataFrame): Raw counts with a miRNA ID column and one column per sample.
        targetFile (pandas.DataFrame): Target table with a ``sample`` column.
        design (str): Design column used in the formula and contrasts. Defaults to ``'sample'``.
        combination (list): Comparisons written as ``"A-B"``.
        mirna_column (str): Name of the miRNA ID column. Defaults to ``'mature_id'``.
        mmg (bool): Unused. Defaults to False.
        subset (bool): Fit a separate model per comparison using only its samples. Defaults to False.
        lib (str, optional): R library path containing DESeq2.

    Returns:
        pandas.DataFrame: miRNA IDs with ``baseMean``, ``logFC``, ``lfcSE``,
        ``stat``, ``pvalue`` and ``FDR`` columns for each comparison, suffixed
        with ``(comparison)``.

    Raises:
        ValueError: If ``countDF``, ``targetFile`` or ``combination`` is None.
        Exception: If DESeq2 cannot be imported.
    """
    if countDF is None or targetFile is None or combination is None:
        raise ValueError("countDF, targetFile and combination are required")

    (
        robjects,
        Formula,
        numpy2ri,
        pandas2ri,
        get_conversion,
        localconverter,
        py2rpy,
        importr,
        _r,
        to_dataframe,
    ) = _rpy2()

    try:
        deseq = importr('DESeq2') if lib is None else importr('DESeq2', lib_loc=lib)
    except Exception:
        log.error("DESeq2 installation was not found ")
        raise

    countDF = countDF.copy()
    mirna_id = countDF[[mirna_column]].values
    countDF.set_index(mirna_column, inplace=True)
    levels = set(targetFile[design].astype(str))

    deseq_results=pd.DataFrame(mirna_id,columns=[mirna_column])
    if not len(combination):
        return deseq_results

    # subset: one iteration per combination; otherwise all count data as a single set
    loopSamples = np.array(range(len(combination))) if subset else np.array([0])
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:

        if subset:

            c1, c2, _ = _split_contrast(combination[j], levels)

            subTF = targetFile[targetFile[design].astype(str).isin([c1, c2])]
            subDF = countDF.loc[:, [col for col in subTF.index if col in countDF.columns]]

            with localconverter(get_conversion() + pandas2ri.converter + numpy2ri.converter):

                count_matrix = py2rpy(subDF)

                design_matrix = py2rpy(subTF)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix, design= design_formula)

            comb = [combination[j]] # put the combination in comb for constrast

        else:

            with localconverter(get_conversion() + pandas2ri.converter + numpy2ri.converter):

                count_matrix = py2rpy(countDF)

                design_matrix = py2rpy(targetFile)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix,  design= design_formula)

            comb=combination
    
        dds1 = deseq.DESeq(dds, quiet=True) # run deseq2


        # Itrate through all the given combination for DEG comparision 

        for co in comb:

            c1, c2, co_label = _split_contrast(co, levels)

            R_contrast = robjects.vectors.StrVector(np.array([design,c1,c2]))

            result = deseq.results(dds1, contrast=R_contrast)

            result=to_dataframe(result)

            with localconverter(get_conversion() + pandas2ri.converter + numpy2ri.converter):

                result = robjects.conversion.rpy2py(result)

            result=pd.DataFrame(result)

            result['padj']=result['padj'].replace(np.nan,1)

            result['log2FoldChange'] = result['log2FoldChange'].replace(np.nan, 0)
            result.columns = ['baseMean','logFC','lfcSE','stat','pvalue','FDR']
            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+co_label+")" for s in result.columns]  # Add combination names to the column names 

            deseq_results.reset_index(drop=True, inplace=True)

            deseq_results = pd.concat([deseq_results,result],axis=1)

    return deseq_results

def run_edgeR(countDF=None, targetFile=None, combination=None,  mirna_column='mature_id',  mmg=False, subset=False, replicate=True, bcv= 0.4, lib=None):
    """Run edgeR differential expression analysis on raw read counts.

    Args:
        countDF (pandas.DataFrame): Raw counts with a miRNA ID column and one column per sample.
        targetFile (pandas.DataFrame): Target table with a ``sample`` column.
        combination (list): Comparisons written as ``"A-B"``.
        mirna_column (str): Name of the miRNA ID column. Defaults to ``'mature_id'``.
        mmg (bool): Unused. Defaults to False.
        subset (bool): Fit a separate model per comparison using only its samples. Defaults to False.
        replicate (bool): False if there are no replicates, in which case ``bcv`` is used. Defaults to True.
        bcv (float): Biological coefficient of variation used without replicates. Defaults to 0.4.
        lib (str, optional): R library path containing edgeR and limma.

    Returns:
        pandas.DataFrame: miRNA IDs with ``logFC``, ``logCPM``, ``LR``,
        ``pvalue`` and ``FDR`` columns for each comparison, suffixed with
        ``(comparison)``.

    Raises:
        ValueError: If ``countDF``, ``targetFile`` or ``combination`` is None.
        Exception: If edgeR or limma cannot be imported.
    """
    if countDF is None or targetFile is None or combination is None:
        raise ValueError("countDF, targetFile and combination are required")

    (
        robjects,
        _Formula,
        numpy2ri,
        pandas2ri,
        get_conversion,
        localconverter,
        py2rpy,
        importr,
        r,
        to_dataframe,
    ) = _rpy2()

    try:
        if lib is None:
            edgeR = importr('edgeR')
            importr('limma')
        else:
            edgeR = importr('edgeR', lib_loc=lib)
            importr('limma', lib_loc=lib)
            
    except Exception:
        log.error("edgeR installation not found. Please install edgeR")
        raise

    countDF = countDF.copy()
    mirna_id = countDF[[mirna_column]].values
    countDF.set_index(mirna_column, inplace=True)
    all_ids = countDF.index
    levels = set(targetFile['sample'].astype(str))

    edgeR_results=pd.DataFrame(mirna_id,columns=[mirna_column])
    if not len(combination):
        return edgeR_results

    # subset: one iteration per combination; otherwise all count data as a single set
    loopSamples = np.array(range(len(combination))) if subset else np.array([0])
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:

        if subset:

            c1, c2, _ = _split_contrast(combination[j], levels)

            subTF = targetFile[targetFile['sample'].astype(str).isin([c1, c2])]
            active_df = countDF.loc[:, [col for col in subTF.index if col in countDF.columns]]
            groups = subTF['sample']
            comb = [combination[j]]

        else:

            active_df = countDF[countDF.sum(axis=1)>2]
            groups = targetFile['sample']
            comb = list(combination)

        with localconverter(get_conversion() + pandas2ri.converter + numpy2ri.converter):
            count_matrix = py2rpy(active_df)
            group = py2rpy(groups)

        dds = edgeR.DGEList(counts=count_matrix,group=group)
        dds = edgeR.calcNormFactors(dds)
        r.assign('group', group)
        r.assign('dds', dds)
        design_mat = r('model.matrix(~0+dds$samples$group,data=dds$samples)')
        level_names = [str(lv) for lv in r('levels(dds$samples$group)')]

        if replicate:
            dds = edgeR.estimateGLMCommonDisp(dds, design_mat)
            dds = edgeR.estimateGLMTrendedDisp(dds, design_mat)
            dds = edgeR.estimateGLMTagwiseDisp(dds, design_mat)
            fit = edgeR.glmFit(dds, design_mat)
        else:
            fit = edgeR.glmFit(dds, design_mat, dispersion=float(bcv**2))

        for co in comb:
            c1, c2, co_label = _split_contrast(co, levels)
            contrast_vec = robjects.FloatVector(
                [1.0 if lv == c1 else (-1.0 if lv == c2 else 0.0) for lv in level_names]
            )
            lrt = edgeR.glmLRT(fit, contrast=contrast_vec)
            deg = edgeR.topTags(lrt, n=active_df.shape[0], sort_by="none")
            result = to_dataframe(deg)

            with localconverter(robjects.default_converter + pandas2ri.converter + numpy2ri.converter):
                result = robjects.conversion.rpy2py(result)

            result = pd.DataFrame(result, index=active_df.index).reindex(all_ids)
            result.columns = ['logFC','logCPM','LR','pvalue','FDR']
            result['logFC'] = result['logFC'].fillna(0.0)
            result['FDR'] = result['FDR'].fillna(1.0)
            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+co_label+")" for s in result.columns]  # Add combination names to the column names 

            edgeR_results.reset_index(drop=True, inplace=True)

            edgeR_results = pd.concat([edgeR_results,result],axis=1)
    
    return edgeR_results

def degFilter(degDF=None, CompareList=None, FDR=0.05, FOLD=2, plot=True, figsize=(10,6),text_size=14, replicate=True):
    """Filter differential expression results by fold change and FDR.

    Args:
        degDF (pandas.DataFrame): Results of all comparisons with a ``mature_id`` column.
        CompareList (list): Comparison names.
        FDR (float): FDR cutoff, applied only when ``replicate`` is True. Defaults to 0.05.
        FOLD (float): Fold change cutoff; its log2 is used. Defaults to 2.
        plot (bool): Draw a bar plot of DEGs per comparison. Defaults to True.
        figsize (tuple): Unused. Defaults to (10, 6).
        text_size (int): Font size of the plot text. Defaults to 14.
        replicate (bool): True if the data have replicates and an FDR is available. Defaults to True.

    Returns:
        dict: ``summary`` (DEG counts per comparison), ``filtered``,
        ``filteredup`` and ``filtereddown`` (DataFrames per comparison) and
        ``plot`` (matplotlib figure or None).

    Raises:
        ValueError: If ``degDF`` or ``CompareList`` is None.
    """
    if degDF is None or CompareList is None:
        raise ValueError("degDF and CompareList are required")

    fig = None
    Up = []
    Down = []
    Total = []
    DEGs = {}
    Ups = {}
    Downs = {}
    # summary = pd.DataFrame()

    degDF = degDF.set_index('mature_id')
        
    for c in CompareList:

        suffix = f"({c})"
        dk = degDF.loc[:, [col for col in degDF.columns if col.endswith(suffix)]]

        FDRR = "FDR("+c+")"

        LFC = "logFC("+c+")"

        if replicate:

            fdr = dk[dk[FDRR]<=FDR].dropna()

            upDF = fdr[fdr[LFC]>=np.log2(FOLD)]

            downDF = fdr[fdr[LFC]<=-np.log2(FOLD)]
        else:
            upDF = dk[dk[LFC]>=np.log2(FOLD)]

            downDF = dk[dk[LFC]<=-np.log2(FOLD)]

        Up.append(upDF.shape[0])

        Down.append(downDF.shape[0])

        Total.append(upDF.shape[0]+downDF.shape[0])
        
        upDF =upDF.reset_index()
        downDF =downDF.reset_index()
        final = pd.concat([upDF, downDF], axis=0)
      
        DEGs[c] = final
        Ups[c] = upDF
        Downs[c] = downDF

    summary =pd.DataFrame({"Comparisons": CompareList, "Total_DEGs": Total, "Up_DEGs": Up, "Down_DEGs": Down})

    if plot== True:

        labels= summary['Comparisons'].values.tolist()

        hg = 15 if len(labels) > 10 else len(labels)
            

        updata= summary['Up_DEGs'].values.tolist()
        downdata = summary['Down_DEGs'].values.tolist()
        my_range=list(range(1,len(summary.index)+1))
        fig, ax = plt.subplots(figsize=(hg,hg),dpi=300)
        ax.barh(labels,updata, height=0.5, color='mediumseagreen')
        ax.barh(labels,downdata, left=updata, height=0.5, color='salmon')
        plt.xticks(fontsize=text_size, weight='bold')
        plt.yticks(fontsize=text_size, weight='bold')
        plt.xlabel("Number of Genes", fontsize=text_size, weight='bold', labelpad=20)
        plt.ylabel("Comparisons", fontsize=text_size, weight='bold', labelpad=20)
        plt.legend(['Up-regulated', 'Down-regulated'],  ncol =2, loc='center', fontsize=6, bbox_to_anchor=(0.5, 1.1))


        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_bounds(0, len(my_range))
        # add some space between the axis and the plot
        ax.spines['left'].set_position(('outward', 8))
        ax.spines['bottom'].set_position(('outward', 5))
        # plt.savefig('deg.png', dpi=300, bbox_inches='tight')
        if replicate:
            plt.title(f'Filter DEGs (Fold:{FOLD} and FDR:{FDR})', loc='center', fontsize=text_size, weight='bold', pad=20)
        else:
            plt.title(f'Filter DEGs (Fold:{FOLD} )', loc='center',fontsize=12, weight='bold', pad=50)
        fig.tight_layout()
            
    return {'summary': summary, "filtered": DEGs,"filteredup":Ups, "filtereddown":Downs, "plot": fig}

       
