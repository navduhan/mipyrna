#!/usr/bin/python

"""
Title: Feature generation class for Novel miRNA precursor prediction
Author: Naveen Duhan
Date: August 23, 2022
"""

# import statements
import RNA
import forgi.graph.bulge_graph as fgb



class Precursor_Features:

    def __init__(self, sequence=None):
        """Feature Generation Class

        Args:
            sequence (_type_, optional): _description_. Defaults to None.
        """

        self.seq = sequence.replace('T', 'U')

        try:

            data = self.get_buldge_graph()

            self.buldge_graph = data['bulge_graph']

            self.structure = data['structure']

            self.mfe = data['mfe']
        except Exception:
            pass

        return

    def get_buldge_graph(self):

        pstructure, pmfe = RNA.fold(self.seq)

        bulge_graph = fgb.BulgeGraph.from_dotbracket(pstructure, self.seq)

        return {'structure': pstructure, 'mfe': pmfe, 'bulge_graph': bulge_graph}

    def find_biggest_stem(self):

        big_stem = (-1, 'a')

        for stem in self.buldge_graph.stem_iterator():

            if self.buldge_graph.stem_length(stem) > big_stem[0]:

                big_stem = (self.buldge_graph.stem_length(stem), stem)

        return big_stem[0]

    def find_stem_pairs(self):

        stem_pairs = []

        try:

            for stem in self.buldge_graph.stem_iterator():

                for p in self.buldge_graph.stem_bp_iterator(stem):

                    stem_pairs.append(
                        (self.buldge_graph.seq[p[0]-1], self.buldge_graph.seq[p[1]-1]))
        except Exception:
            stem_pairs = []

        return {'stem_pairs': stem_pairs, 'stem_length': len(stem_pairs)}

    def find_matched_stem(self):

        base_pairs = {'A': 'T', 'T': 'A','U': 'A', 'A': 'U', 'G': 'C', 'C': 'G', 'R': 'R','K': 'K','Y': 'Y', 'S': 'S', 'M': 'M', 'W': 'W', 'H': 'H', 'B': 'B', 'V': 'V', 'D': 'D', 'N': 'N'}
        matched_stem = []

        for pair in self.find_stem_pairs()['stem_pairs']:

            if str(pair[1]) == base_pairs[str(pair[0])]:

                matched_stem.append(True)

            else:

                matched_stem.append(False)

        return {'matched_stem': matched_stem, 'mismatched_stem': matched_stem.count(False)}

    def find_loops(self):

        try:
            h0seq = self.buldge_graph.get_define_seq_str('h0')[0]

            loops = len(list(self.buldge_graph.hloop_iterator()))

            loop_length = self.buldge_graph.get_bulge_dimensions('h0')[0]

            loop_gc = sum(1 for x in h0seq if x in ['G', 'C'])/len(h0seq)*100

        except:

            h0seq = ''

            loops = 0

            loop_length = 0

            loop_gc = 0

        return {'loop_gc': loop_gc, 'total_loops': loops, 'loop_length': loop_length}

    def find_bulges(self):

        bulges = []

        for bulge in self.buldge_graph.iloop_iterator():

            bulges.append(self.buldge_graph.get_bulge_dimensions(bulge))

        total_bulges = len(bulges)

        try:
            longest_bulge = max(max(zip(*bulges)))
        except:
            longest_bulge = 0

        matched_bulges = [True if i[0] == i[1] else False for i in bulges]

        symmetric_bulge = matched_bulges.count(True)

        asymmetric_bulge = matched_bulges.count(False)

        return {'bulges': bulges, 'total': total_bulges, 'longest_bulge': longest_bulge, 'symmetric_bulge': symmetric_bulge, 'asymmetric_bulge': asymmetric_bulge}

    def calculate_GC(self):

        GC_count = sum(1 for x in self.seq if x in ['G', 'C'])

        if len(self.seq) == 0:

            return 0

        return float(GC_count) / len(self.seq) * 100

    def generate_nucleotide_composition(self):

        nuc = ['A', 'C', 'U', 'G']

        nuc_comp = {}

        for n in nuc:

            nuc_comp[n] = round(float(self.seq.count(n)/len(self.seq))*100, 2)

        return nuc_comp

    def generate_dinucleotide_composition(self):

        N = len(self.seq)

        nuc = ['A', 'C', 'U', 'G']

        dpc = {}

        for i in nuc:

            for j in nuc:

                dp = i + j

                dp_count = round(float(self.seq.count(dp)) / (N - 1) * 100, 2)

                dpc[dp] = dp_count

        return dpc

    def generate_trinucleotide_composition(self):

        N = len(self.seq)

        nuc = ['A', 'C', 'U', 'G']

        tpc = {}

        for i in nuc:

            for j in nuc:

                for k in nuc:

                    tp = i + j+k

                    tp_count = round(
                        float(self.seq.count(tp)) / (N - 2) * 100, 2)

                    tpc[tp] = tp_count
        return tpc

    def generate_triplets(self):

        structure_string = ['(((', '((.', '(..', '(.(',
                            '.((', '.(.', '..(', '...']

        nuc = ['A', 'G', 'U', 'C']

        structure_triplets = {}

        for i in nuc:

            for j in structure_string:

                structure_triplets[i+j] = 0

        self.structure = self.structure.replace(')', '(')

        l = len(self.seq)-len(self.seq) % 3

        for i in range(0, l, 3):

            n = self.seq[i+1]+self.structure[i:i+3]

            if n in structure_triplets:

                structure_triplets[n] += 1

        return structure_triplets

    def generate_features(self):

        nfeatures = {}

        stems = self.find_stem_pairs()

        if stems['stem_length'] <= 1:

            return

        loops = self.find_loops()
        bulges = self.find_bulges()
        mathed_stem = self.find_matched_stem()

        # Assemble all the features in one dictionary
        nfeatures['length'] = len(self.seq)
        nfeatures['mfe'] = self.mfe
        nfeatures['GC'] = self.calculate_GC()
        nfeatures['loops'] = loops['total_loops']
        nfeatures['loop_GC'] = loops['loop_gc']
        nfeatures['loop_length'] = loops['loop_length']
        nfeatures['stem_length'] = stems['stem_length']
        nfeatures['longest_stem'] = self.find_biggest_stem()
        nfeatures['bulges'] = bulges['total']
        nfeatures['longest_bulge'] = bulges['longest_bulge']
        nfeatures['symmetric_bulge'] = bulges['symmetric_bulge']
        nfeatures['asymmetric_bulge'] = bulges['asymmetric_bulge']
        nfeatures['mismatched_stem'] = mathed_stem['mismatched_stem']
        nfeatures.update(self.generate_nucleotide_composition())
        nfeatures.update(self.generate_dinucleotide_composition())
        nfeatures.update(self.generate_trinucleotide_composition())
        nfeatures.update(self.generate_triplets())

        return nfeatures
