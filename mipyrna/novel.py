import RNA
import re
import os
import sys
import shutil
import pandas as pd
import numpy as np
import forgi.graph.bulge_graph as fgb
from  mipyrna.reads import Read_process
from keras.models import load_model
import pkg_resources
from mipyrna import utility as mu
from mipyrna.features import Precursor_Features
from mipyrna.logger import MiPyRNALogger
from mipyrna.mirbase import MirBase

log = MiPyRNALogger(mode='a', log='novel')


class Novel_miRNA:

    def __init__(self,  genome=None, read_cutoff=50, species=None, species_type='plants', outdir="."):

        self.genome = genome
        self.readsCutoff = read_cutoff
        self.species = species
        self.species_type = species_type
        self.outdir = outdir

        return

    def get_precursor_structure(self, seq):

        try:
            x = RNA.fold(seq)
        except Exception as e:
            print (e)
            return
        return x

    def get_bulge(self, seq):
        
            seq = seq.replace('T','U')
            pstructure, pmfe = RNA.fold(seq)

            bulge_graph = fgb.BulgeGraph.from_dotbracket(pstructure, seq)
            
            return bulge_graph

    def close_match(self,  preseq, seq):

        string_to_match = re.compile('|'.join(seq[:i]+'.'+seq[i+1:] for i in range(len(seq))))
        m = string_to_match.findall(preseq)
        if len(m)==0: return -1
        idx = preseq.find(m[0])
        return idx

    def find_subseq(self, seq, string):

        for i in range(16,4,-4):
            c = seq.find(string[:i])
            if c != -1: return c
        return self.close_match(seq, string)

    def check_mature(self, precursor, mature):

            bg = self.get_bulge(precursor)
            if bg is None:
                return
            start = self.find_subseq(precursor, mature)+1
            end = start + len(mature)
            loops = list(bg.hloop_iterator())
            if len(loops)==0:
                return 'too many loops'
            l = list(bg.define_range_iterator('h0'))[0]
            if (start<l[0] and end>l[0]) or (start<l[1] and end>l[1]):
                return 'mature in loop'
            p = list(bg.stem_bp_iterator('s0'))[0]
            
            if start == 0:
                return 'mature not found'
            if end>p[1]+2:
                return '3p mature outside'
            if p[0]>start+2:
                return 'mature unpaired'
            return 'ok'
        
        
    def find_star_sequence(self, precursor, mature):
        """Estimate the star sequence from a given mature and precursor."""

        bg = self.get_bulge(precursor)
        start = self.find_subseq(precursor, mature)+1
        end = start + len(mature)
        
        for p in bg.adjacent_stem_pairs_iterator():
            pp = [p[0], p[2]]
            for s in pp:
                xx = list(bg.stem_bp_iterator(s))
                for i in xx:
                    if start == i[0]:
                       
                        star = precursor[i[1]-len(mature)+2:i[1]+2]
                        return star
                    elif start == i[1]:
                       
                        star = precursor[i[0]-len(mature)+1:i[0]+2]
                        return star
        return

    def predict_precursor(self, df):
        positives = df['features'].values.tolist()
        file = pkg_resources.resource_filename('mipyrna', "data/plants_F5_100_train.h5")
        # file='data/plants_F5_100_train.h5'
        myModel=load_model(file,compile=False)
        samples, labels = mu.preprocess(positives, np.array([1,0]), 128) 
        Samples=samples.reshape(samples.shape[0], 1, samples.shape[1],1)
        prob=myModel.predict(Samples)
        prob = prob.tolist()
        df= pd.DataFrame(prob, columns=['miRNA', 'Non-miRNA'])
        
        return df

    def get_putative_precursor(self, filtered_cluster=None,  other=None,  seqlength=20, loop=20, precursor_step= 10 ):
        
        global mature_status
        putative_precursor = []
        
        x = filtered_cluster.iloc[0]
        ref_coords = (x['name'], x.start-10, x.end+10, x.strand)
        # first get a sequence from reference genome based on coordinates
        refseq =  Read_process(reference=self.genome).get_sequence(ref_coords)
        # check whether any sequencing reads is in sequence above
        mature_mirna = Read_process(reference=self.genome).get_consensus_read(ref_precursor=refseq, aligned_cluster=filtered_cluster)
        # for i, x in filtered_cluster.iterrows():
        coords = (x['name'], x.start, x.start, x.strand)
        chrom, start, end, strand = coords
        mature_mirna =x.seq
    
        if mature_mirna != None:
            seqlen = len(mature_mirna)
        else:
            seqlen = seqlength

        # get 5' precursor


        for i in range(10, 60, precursor_step):
        
            start5 = start - 1 
            end5 = start + 2 * seqlen-1 + loop + i
            coords = [chrom,start5,end5,strand]
        
            precursor_seq = Read_process(reference=self.genome).get_sequence(coords)
            
            if precursor_seq == None or 'N' in precursor_seq:
                continue
            try:
                struct,sc = self.get_precursor_structure(precursor_seq)
            
                mature_status = self.check_mature(precursor_seq, mature_mirna)
                
                if mature_status is None:
                    continue
            except Exception:
                continue
            
            k = Precursor_Features(precursor_seq)

            m = k.generate_features()

            if m == None:
                    continue
                
            maturecounts = filtered_cluster.reads.sum()
            putative_precursor.append({'chrom':chrom, 'mature':mature_mirna, 'mature_length':len(mature_mirna), 'precursor':precursor_seq, 'precursor_length':len(precursor_seq),'start':start5,'end':end5,'mature_start':start,
                    'strand':strand, 'struct':struct,'mfe':sc,  'mature_check': mature_status, 'features': m, 'mature_reads': maturecounts, 'position': '5-prime'})

        for i in range(10,60, precursor_step):

            start3 = start - (loop + seqlen + i)
            end3 = end + seqlen + 1 
            coords = [chrom,start3,end3,strand]
            
            precursor_seq = Read_process(reference=self.genome).get_sequence(coords)
            
            if precursor_seq == None or 'N' in precursor_seq:
                continue
            try:
                struct,sc = self.get_precursor_structure(precursor_seq)
            
                mature_status = self.check_mature(precursor_seq, mature_mirna)
                
                if mature_status is None:
                    continue
            except Exception:
                continue
            

            k = Precursor_Features(sequence=precursor_seq)

            m = k.generate_features()

            if m == None:

                    continue
            maturecounts = filtered_cluster.reads.sum()
            
            putative_precursor.append({'chrom':chrom, 'mature':mature_mirna, 'mature_length':len(mature_mirna), 'precursor':precursor_seq, 'precursor_length':len(precursor_seq),'start':start3,'end':end3,'mature_start':start,
                    'strand':strand, 'struct':struct,'mfe':sc,  'mature_check': mature_status, 'features': m, 'mature_reads': maturecounts,  'position': '3-prime'})
        putative_precursor = pd.DataFrame(putative_precursor)
        
        if len(putative_precursor)>0:
            P = putative_precursor.iloc[0].copy()
            
        else:
            return
        
        star = self.find_star_sequence(P.precursor, P.mature)
        starcounts = 0
        if other is not None and star != None:
            s = self.find_subseq(P.precursor, star)

            ss = P.start+s; se = ss+len(star)

            sreads = other[(other.start>=ss-2) & (other.end<=se+3)]
            starcounts = sreads.reads.sum()
            
        putative_precursor['star_reads'] = starcounts
        putative_precursor['star'] = star
        putative_precursor['cluster'] = x.cluster
        
        return putative_precursor
    
    def filter_miRNA(self, df):
        # Group by the 'mature' column and aggregate other columns with desired functions
        grouped = df.groupby('mature').agg({'chrom': 'first', 'mature_length': 'first', 'precursor': 'first',
                                            'precursor_length': 'first', 'start': 'first', 'end': 'first',
                                            'mature_start': 'first', 'strand': 'first', 'struct': 'first',
                                            'mfe': 'first', 'mature_check': 'first', 'mature_reads': 'first',
                                            'position': 'first', 'star_reads': 'first', 'star': 'first',
                                            'cluster': 'first', 'score': 'max'}).reset_index()

        # Create a dictionary from the grouped DataFrame
        final = {row['mature']: row.tolist()[1:] for _, row in grouped.iterrows()}

        final_miRNAs = pd.DataFrame.from_dict(final, orient='index', columns=grouped.columns[1:]).reset_index()
        final_miRNAs.columns = ['mature'] + list(grouped.columns[1:])

        try:
            final_miRNAs['Novel_ID'] = [f"{self.species}_mipyrna_{final_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                                        for i in range(len(final_miRNAs))]
        except Exception:
            pass

        return final_miRNAs
    
    def _write_fasta(self, df, output, type='known', seq='mature'):
        mature_seq = []
        with open(output,'w') as fp:
            for i, r in df.iterrows():
                if seq=='mature':
                    if type=='known':
                        if r['Known_ID'] in  mature_seq:
                            continue
                        else:
                            fp.write(f">{r['Known_ID']}\n{r['mature']}\n")
                            mature_seq.append(r['Known_ID'])
                    else:
                        fp.write(f">{r['Novel_ID']}\n{r['mature']}\n")
                if seq == 'hairpin':
                    if type=='known':
                        if r['Known_ID'] in  mature_seq:
                            continue
                        else:
                            fp.write(f">{r['Known_ID']}\n{r['precursor']}\n")
                            mature_seq.append(r['Known_ID'])
                    else:
                        fp.write(f">{r['Novel_ID']}\n{r['precursor']}\n")
        return output
    
    def predict_known_miRNAs(self, novel_df):

        blast = mu.make_directory(os.path.join(self.outdir, "blast_results"))
        novel_mature = os.path.join(self.outdir, "all_sample_pooled_mature.fa")

        mir = MirBase(species=self.species, speciesClass=self.species_type, outdir= blast)

        mirbase_data = mir.get_species()

        execPATH = shutil.which('makeblastdb')

        if execPATH is None:

            log.error("ncbi-blast not found in path")

            sys.exit(1)
        else:

            blastdb_command = f"makeblastdb -in {mirbase_data['mature']} -dbtype 'nucl'"

        mu.local_run(program='makeblastdb', command=blastdb_command, outdir=blast, message=f"building database for {mirbase_data['mature']}")
        
        execPATH = shutil.which('blastn')

        if execPATH is None:

            log.error("ncbi-blast+ not found in path")

            sys.exit(1)
        else:
            
            blast_out = os.path.join(blast, "blast_mature_vs_miRbase.txt")
            blastn_command = f"blastn -db {mirbase_data['mature']} -query {novel_mature} -evalue 100 -out {blast_out} -outfmt 6 -word_size 4"
            mu.local_run(program='blastn', command=blastn_command, outdir=blast, message=f"similarity search against miRbase")

        blast_results = pd.read_csv(blast_out, sep="\t", names=['Novel_ID', 'Known_ID', 'pident', 'length', 'mismatch', 'gapopen', 'qstart','qend', 'sstart', 'send', 'evalue', 'bitscore'])

        filterd_results = blast_results[(blast_results.qstart<2) & (blast_results.length>19) &(blast_results.pident>94)].drop_duplicates('Novel_ID') 

        known_miRNAs_dataframe = filterd_results[['Novel_ID', 'Known_ID']]

        known_miRNAs = known_miRNAs_dataframe.merge(novel_df, on=['Novel_ID'],  indicator=True)

        known_miRNAs =known_miRNAs.drop("_merge", axis=1)

        novel_miRNAs = pd.merge(known_miRNAs['Novel_ID'],novel_df, on=['Novel_ID'], how='outer', indicator=True).query('_merge == "right_only"').drop("_merge", axis=1)

        try:
            known_miRNAs['raw_ID'] = known_miRNAs['Novel_ID']
            known_miRNAs['Novel_ID'] = [f"{self.species}_miPyRNA_{known_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                                        for i in range(len(known_miRNAs))]
        except Exception:
            pass
        
        try:
            novel_miRNAs['raw_ID'] = novel_miRNAs['Novel_ID']
            novel_miRNAs['Novel_ID'] = [f"{self.species}_miPyRNA_{novel_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                                        for i in range(len(novel_miRNAs))]
        except Exception:
            pass

        known_miRNAs.to_excel(f"{self.outdir}/known_miRNAs.xlsx", index=False)
        novel_miRNAs.to_excel(f"{self.outdir}/novel_miRNAs.xlsx", index=False)
        
        novel_fasta = self._write_fasta(df=novel_miRNAs, output=os.path.join(f"{self.outdir}/novel_mature.fa"), type='novel')
        known_fasta = self._write_fasta(df=known_miRNAs, output=os.path.join(f"{self.outdir}/known_mature.fa"), type='known')

        return known_miRNAs, novel_miRNAs,  known_fasta, novel_fasta


    def get_novel_miRNA(self,samples=None, cluster_cutoff = 28):
       
        outputs = []
        
        mirna_temp = mu.make_directory(os.path.join(self.outdir, "temp"))
        
        for key, sample in samples.items():
            
            reads = [] #stores reads associated with each cluster
            N = []

            clusts = sample[3][sample[3].reads>=self.readsCutoff]
            for i,c in clusts.iterrows():
                df = sample[2][sample[2].cluster==c.cluster]
                df = df.sort_values('reads',ascending=False)
                if c.clust_size<cluster_cutoff:
                    df['mature'] = True
                    reads.append(df)
                    N.append(self.get_putative_precursor(filtered_cluster=df))
                else:
                    anchor = df.iloc[0]
                    st = anchor.start
                    end = anchor.end
                    mm = df.loc[(abs(df.start-st)<=3) & (abs(df.end-end)<=5)].copy()
                    if mm.reads.sum() <= self.readsCutoff:
                        continue
                    mm['mature'] = True
                    reads.append(mm)
                    other = df.loc[~df.index.isin(mm.index)].copy()
                    other['mature'] = False
                    reads.append(other)
                    
                    N.append(self.get_putative_precursor(filtered_cluster=mm, other=other))
                
            final_data = pd.concat(N)
            
            final_reads = pd.concat(reads)

            pre_prob = self.predict_precursor(final_data)
        
            final_data['score'] = pre_prob['miRNA']

            outputs.append(final_data)

            final_data.to_csv(f"{mirna_temp}/{sample[2]}_raw_miRNA.txt", sep="\t", index=False)
            
        final = pd.concat(outputs)

        final_miRNA = self.filter_miRNA(final)

        final_miRNA.to_csv(os.path.join(self.outdir, "all_sample_pooled_miRNAs.txt"), sep="\t", index=False)
        final_data = final_miRNA[(final_miRNA['mature_check']=='ok') & (final_miRNA['score']>=0.5)]
        final_data.to_csv(os.path.join(self.outdir, "all_sample_filtered_miRNAs.txt"), sep="\t", index=False)
        mature_out = os.path.join(self.outdir,"all_sample_pooled_mature.fa")
        hairpin_out = os.path.join(self.outdir,"all_sample_pooled_hairpin.fa")
        
        self._write_fasta(final_miRNA,mature_out,type='raw',seq='mature')
        self._write_fasta(final_miRNA,hairpin_out,type='raw',seq='hairpin')

        known_miRNAs, novel_miRNAs, known_fasta, novel_fasta = self.predict_known_miRNAs(novel_df=final_miRNA)
        
        return known_miRNAs, novel_miRNAs, known_fasta, novel_fasta
    

    
