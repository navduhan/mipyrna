#!/usr/bin/env python3

'''
:Title: This module finds differentially expressed genes from raw read counts

:Created: October 22, 2021

:Author: Naveen Duhan
'''

import logging
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula, numpy2ri, pandas2ri
from rpy2.robjects.conversion import get_conversion, localconverter, py2rpy
from rpy2.robjects.packages import importr

from mipyrna.logger import MiPyRNALogger

rpy2_logger.setLevel(logging.ERROR)

log = MiPyRNALogger(mode='a', log="diff")

r: Any = robjects.r

to_dataframe = r('function(x) data.frame(x)')

numpy2ri.activate()


def runDESeq2(countDF=None, targetFile=None, design='sample',combination=None,  mirna_column='mature_id',  mmg =False, subset=False, lib=None):
    """
    This function is a wrapper to DESeq2 package in R for differeantial expression analysis from raw read counts.

    :param countFile: Raw read count file. 
    
    :param targetFile: Tab-delimited target file with replication and sample name

    :param design: [description]. Defaults to None.

    :param combination: Comparison list contaning samples to compare.

    :param gene_column: First column in raw read count file. Defaults to 'Gene'.

    :param mmg: True if raw read counts are from multimapped gene groups.

    :param subset: If runDESeq2 subset raw read count according to comparison.

    :param lib: library path of DESeq2 to use.

    :returns: DataFrame

    :rtype: A datafram containing all gene differantial expression in all combinations.
        

    """
    if countDF is None or targetFile is None or combination is None:
        raise ValueError("countDF, targetFile and combination are required")

    try:
        deseq = importr('DESeq2') if lib is None else importr('DESeq2', lib_loc=lib)
    except Exception:
        log.error("DESeq2 installation was not found ")
        raise


    mirna_id = countDF[[mirna_column]].values
    countDF.set_index(mirna_column, inplace=True)

    deseq_results=pd.DataFrame(mirna_id,columns=[mirna_column])

    # subset: one iteration per combination; otherwise all count data as a single set
    loopSamples = np.array(range(len(combination))) if subset else np.array([1])
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:

        if subset:

            c1,c2 = combination[j].split("-")

            subDF = countDF.filter(regex=f'{c1}|{c2}') # Subset the count data for one combination

            subTF = targetFile[targetFile['sample'].str.contains(f'{c1}|{c2}')] # Subset the target data for one combination

            with localconverter(get_conversion() + pandas2ri.converter):

                count_matrix = py2rpy(subDF)

                design_matrix = py2rpy(subTF)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix, design= design_formula)

            comb = [combination[j]] # put the combination in comb for constrast

        else:

            with localconverter(get_conversion() + pandas2ri.converter):

                count_matrix = py2rpy(countDF)

                design_matrix = py2rpy(targetFile)

            designFormula="~ "+design

            design_formula = Formula(designFormula)

            dds=deseq.DESeqDataSetFromMatrix(countData=count_matrix, colData=design_matrix,  design= design_formula)

            comb=combination
    
        dds1 = deseq.DESeq(dds, quiet=True) # run deseq2


        # Itrate through all the given combination for DEG comparision 

        for co in comb:

            c1,c2=co.split("-")

            R_contrast = robjects.vectors.StrVector(np.array([design,c1,c2]))

            result = deseq.results(dds1, contrast=R_contrast)

            result=to_dataframe(result)

            with localconverter(get_conversion() + pandas2ri.converter):

                result = robjects.conversion.rpy2py(result)

            result=pd.DataFrame(result)

            result['padj']=result['padj'].replace(np.nan,1)

            result['log2FoldChange'] = result['log2FoldChange'].replace(np.nan, 0)
            result.columns = ['baseMean','logFC','lfcSE','stat','pvalue','FDR']
            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+co+")" for s in result.columns]  # Add combination names to the column names 

            deseq_results.reset_index(drop=True, inplace=True)

            deseq_results = pd.concat([deseq_results,result],axis=1)

    return deseq_results

def run_edgeR(countDF=None, targetFile=None, combination=None,  mirna_column='mature_id',  mmg=False, subset=False, replicate=True, bcv= 0.4, lib=None):
    """
    This function is a wrapper to edgeR package in R for differeantial expression analysis from raw read counts.

    :param countFile: Raw read count file. 
    
    :param targetFile: Tab-delimited target file with replication and sample name

    :param design: [description]. Defaults to None.

    :param combination: Comparison list contaning samples to compare.

    :param gene_column: First column in raw read count file. Defaults to 'Gene'.

    :param mmg: True if raw read counts are from multimapped gene groups.

    :param subset: If runDESeq2 subset raw read count according to comparison.

    :param replicate: False if there are no replicates.

    :param bcv:  Biological coefficient of variation if there are no replicate.

    :param lib: library path of DESeq2 to use.

    :returns: DataFrame

    :rtype: A datafram containing all gene differantial expression in all combinations.
    """
    if countDF is None or targetFile is None or combination is None:
        raise ValueError("countDF, targetFile and combination are required")

    try:
        if lib is None:
            edgeR = importr('edgeR')
            limma = importr('limma')
        else:
            edgeR = importr('edgeR', lib_loc=lib)
            limma = importr('limma', lib_loc=lib)
            
    except Exception:
        log.error("edgeR installation not found. Please install edgeR")
        raise

    mirna_id = countDF[[mirna_column]].values
    countDF.set_index(mirna_column, inplace=True)

    edgeR_results=pd.DataFrame(mirna_id,columns=[mirna_column])

    # subset: one iteration per combination; otherwise all count data as a single set
    loopSamples = np.array(range(len(combination))) if subset else np.array([1])
    
    # Now loopthrough combination in loopSamples 

    for j in loopSamples:
        

        if subset:

            c1,c2 = combination[j].split("-")

            subTF = targetFile[targetFile['sample'].str.contains(f'{c1}|{c2}')] # Subset the target data for one combination

            groups = subTF['sample']

            with localconverter(get_conversion() + pandas2ri.converter):

                count_matrix = py2rpy(countDF)

                group = py2rpy(groups)

            dds = edgeR.DGEList(counts=count_matrix,group=group)
    
            dds = edgeR.calcNormFactors(dds)
            
            r.assign('group', group)

            r.assign('dds', dds)

            design = r('model.matrix(~0+dds$samples$group,data=dds$samples)')

            design = pd.DataFrame(design)

            col= r('levels(dds$samples$group)')

            design.columns = list(col)

            design =design.to_records(index=False)

            cont= robjects.vectors.StrVector([combination[j]])
           
            contrasts = limma.makeContrasts(contrasts=cont,levels=design)

        else:

            counts = np.asarray(countDF,dtype=int)

            cpm = (counts * 1e6) / counts.sum(axis=0) 

            cpm =pd.DataFrame(data=cpm,index=countDF.index,columns=countDF.columns)

            cpm = countDF[cpm.sum(axis=1)>1]

            countDF = countDF[countDF.sum(axis=1)>2]

            groups = targetFile['sample']

            with localconverter(get_conversion() + pandas2ri.converter):

                count_matrix = py2rpy(countDF)

                group = py2rpy(groups)

            dds = edgeR.DGEList(counts=count_matrix,group=group)
    
            dds = edgeR.calcNormFactors(dds)
            
            r.assign('group', group)

            r.assign('dds', dds)

            design = r('model.matrix(~0+dds$samples$group,data=dds$samples)')

            design = pd.DataFrame(design)

            col= r('levels(dds$samples$group)')

            design.columns = list(col)

            design =design.to_records(index=False)
            
            cont= robjects.vectors.StrVector(combination)
           
            contrasts = limma.makeContrasts(contrasts=cont,levels=design)
        
        # Itrate through all the given combination for DEG comparision 
    

    # counts = np.asarray(countDF,dtype=int)
        if replicate:
            dds = edgeR.estimateGLMCommonDisp(dds, design)

            dds = edgeR.estimateGLMTrendedDisp(dds, design)

            dds = edgeR.estimateGLMTagwiseDisp(dds, design)

            fit = edgeR.glmFit(dds, design)
        else:

            fit = edgeR.glmFit(dds, design, dispersion=float(bcv**2))

        for i in range(len(contrasts.T)):

            lrt = edgeR.glmLRT(fit, contrast=contrasts.T[i])
            
            deg = edgeR.topTags(lrt,countDF.shape[0] )
            
            result=to_dataframe(deg)

            with localconverter(robjects.default_converter + pandas2ri.converter):

                result = robjects.conversion.rpy2py(result)
            
            result=pd.DataFrame(result)

            result.columns = ['logFC','logCPM','LR','pvalue','FDR']

            result.reset_index(drop=True, inplace=True)

            result.columns=[s+"("+combination[i]+")" for s in result.columns]  # Add combination names to the column names 

            edgeR_results.reset_index(drop=True, inplace=True)

            edgeR_results = pd.concat([edgeR_results,result],axis=1)
    
    return edgeR_results

def degFilter(degDF=None, CompareList=None, FDR=0.05, FOLD=2, plot=True, figsize=(10,6),text_size=14, replicate=True):
    """
    This function filter all gene expression file based on given FOLD and FDR

    :param degDF: A datafram containing all gene differantial expression in all combinations.

    :param CompareList: A list of all the sample comparison.

    :param FDR: False Discovery Rate for filtering DEGs. Defaults to 0.05.

    :param FOLD: Fold change value. The log2 of the value will be calculated. Defaults to 2.

    :param plot: True if want to plot DEGs per sample on barplot. Defaults to True.
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

        dk = degDF.filter(regex=c, axis=1)

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

       
