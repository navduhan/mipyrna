
import numpy as np
import pandas as pd
import collections
import pyfaidx as pf
from tqdm import tqdm
from bx.intervals.cluster import ClusterTree


class Cluster_reads():

    def __init__(self,  cluster_distance=10, min_size=3,min_length=20, max_length=24, key='align_id'):


        self.cluster_distance = cluster_distance
        self.min_size = min_size
        self.key = key
        self.min = min_length
        self.max = max_length

    def get_clusters(self, samples):
        output = {}
        for key, sample in tqdm(samples.items(), desc="Processing samples"):
            filtered_reads = sample[2][(sample[2].length<=self.max) & (sample[2].length>=self.min)]
            input_file =filtered_reads.copy()
            if self.key in input_file.columns:
                input_file.set_index(self.key, inplace=True)

            cluster_trees = collections.defaultdict(lambda:
                ClusterTree(self.cluster_distance, self.min_size))
            for i, row in input_file.iterrows():
                chrom = row['name']
                cluster_trees[chrom].insert(row.start, row.end, row.name)
                clustertrees = dict(cluster_trees)
            cgroups = []
            icount = 1
            for chrom, cltree in clustertrees.items():
            
                for start, end, ids in cltree.getregions():
                    c = input_file.loc[ids].copy()
                    c['cl_start'] = start
                    c['cl_end'] = end
                    c['cluster'] = icount
                    c = c.groupby(['strand']).filter(lambda x: len(x) > 1)
                    cgroups.append(c)
                    icount +=1
            if len(cgroups) > 0:
                df = pd.concat(cgroups)
            else:
                df = pd.DataFrame()

        #     clusts = df.groupby(['name', 'cluster', 'cl_start', 'cl_end', 'strand']) \
		#    .agg(reads='sum', length=pd.NamedAgg(column='length', aggfunc='max')) \
		#    .reset_index() \
		#    .rename(columns={'cl_start': 'start', 'cl_end': 'end'})
            clusts = df.groupby(['name','cluster','cl_start','cl_end','strand'])\
                        .agg({'reads':np.sum,'length':np.max})\
                        .reset_index()\
                        .rename(columns={'cl_start':'start','cl_end':'end'})

            clusts['clust_size'] = clusts.end-clusts.start
            output[key]= [sample[0], sample[1], df, clusts]
            
        return output
