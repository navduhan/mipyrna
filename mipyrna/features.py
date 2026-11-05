"""Structure and sequence features of miRNA precursors.

Defines Precursor_Features, which folds a precursor with ViennaRNA and builds
the feature dictionary scored by the bundled model. Must stay identical to the
training code so features match what the models were trained on.

:Author: Naveen Duhan
"""

from __future__ import annotations

import forgi.graph.bulge_graph as fgb
import RNA

# The nucleotide/structure-triplet encoding crosses the middle base of each 3-nt
# window with the dot-bracket pattern of that window: 4 bases x 8 patterns = 32
# columns, in the order generate_triplets() emits them.
TRIPLET_BASES = ('A', 'G', 'U', 'C')
TRIPLET_STRUCTURES = ('(((', '((.', '(..', '(.(',
                      '.((', '.(.', '..(', '...')



class Precursor_Features:

    def __init__(self, sequence: str | None = None):
        """Feature Generation Class

        Args:
            sequence (str | None, optional): RNA sequence. Defaults to None.
        """

        if sequence is None:
            raise ValueError("a sequence is required")

        self.seq = sequence.replace('T', 'U')

        self.drop_reason: str | None = None

        self.fallback_note: str | None = None

        self.buldge_graph: fgb.BulgeGraph | None = None
        self.structure: str | None = None
        self.mfe: float | None = None

        try:

            data = self.get_buldge_graph()

            self.buldge_graph = data['bulge_graph']

            self.structure = data['structure']

            self.mfe = data['mfe']
        except (ValueError, TypeError, RuntimeError, IndexError, KeyError) as exc:
            # Keep the reason so callers can report why a sequence was skipped.
            self.buldge_graph = None
            self.drop_reason = f"structure prediction failed: {type(exc).__name__}: {exc}"


    def _require_bulge_graph(self) -> fgb.BulgeGraph:
        """Return the bulge graph, raising if structure prediction failed."""
        if self.buldge_graph is None:
            raise RuntimeError("structure prediction did not produce a bulge graph")

        return self.buldge_graph


    def get_buldge_graph(self):
        """Fold the sequence and return its structure, MFE and bulge graph."""

        pstructure, pmfe = RNA.fold(self.seq)

        bulge_graph = fgb.BulgeGraph.from_dotbracket(pstructure, self.seq)

        return {'structure': pstructure, 'mfe': pmfe, 'bulge_graph': bulge_graph}

    def find_biggest_stem(self):
        """Return the length of the longest stem."""

        graph = self._require_bulge_graph()

        big_stem = (-1, 'a')

        for stem in graph.stem_iterator():

            if graph.stem_length(stem) > big_stem[0]:

                big_stem = (graph.stem_length(stem), stem)

        return big_stem[0]

    def find_stem_pairs(self):
        """Return the base pairs of all stems and their count."""

        # forgi >= 2.0 returns bg.seq as a 1-based Sequence object, and
        # stem_bp_iterator() yields 1-based positions. Convert once to a plain
        # string so p[0] - 1 addresses the base at position p[0].
        graph = self._require_bulge_graph()

        sequence = str(graph.seq)

        stem_pairs = []

        for stem in graph.stem_iterator():

            for p in graph.stem_bp_iterator(stem):

                stem_pairs.append((sequence[p[0] - 1], sequence[p[1] - 1]))

        return {'stem_pairs': stem_pairs, 'stem_length': len(stem_pairs)}

    def find_matched_stem(self):
        """Flag each stem pair as Watson-Crick matched or not, and count the mismatches."""

        base_pairs = {'A': 'U', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'R','K': 'K','Y': 'Y', 'S': 'S', 'M': 'M', 'W': 'W', 'H': 'H', 'B': 'B', 'V': 'V', 'D': 'D', 'N': 'N'}
        matched_stem = []

        for pair in self.find_stem_pairs()['stem_pairs']:

            if str(pair[1]) == base_pairs[str(pair[0])]:

                matched_stem.append(True)

            else:

                matched_stem.append(False)

        return {'matched_stem': matched_stem, 'mismatched_stem': matched_stem.count(False)}

    def find_loops(self):
        """Return the hairpin loop count, length and GC percentage; zero-filled if unavailable."""

        try:
            graph = self._require_bulge_graph()

            h0seq = graph.get_define_seq_str('h0')[0]

            loops = len(list(graph.hloop_iterator()))

            loop_length = graph.get_bulge_dimensions('h0')[0]

            loop_gc = sum(1 for x in h0seq if x in ['G', 'C'])/len(h0seq)*100

        except (ValueError, TypeError, RuntimeError, IndexError, KeyError, ZeroDivisionError) as exc:

            # Zero-filled rather than fatal, but recorded so it is not silent.
            self.fallback_note = f"loop features zero-filled: {type(exc).__name__}: {exc}"

            h0seq = ''

            loops = 0

            loop_length = 0

            loop_gc = 0

        return {'loop_gc': loop_gc, 'total_loops': loops, 'loop_length': loop_length}

    def find_bulges(self):
        """Return interior loop (bulge) counts, the longest, and symmetric and asymmetric counts."""

        graph = self._require_bulge_graph()

        bulges = []

        for bulge in graph.iloop_iterator():

            bulges.append(graph.get_bulge_dimensions(bulge))

        total_bulges = len(bulges)

        try:
            longest_bulge = max(max(zip(*bulges)))
        except (ValueError, TypeError):

            longest_bulge = 0

        matched_bulges = [i[0] == i[1] for i in bulges]

        symmetric_bulge = matched_bulges.count(True)

        asymmetric_bulge = matched_bulges.count(False)

        return {'bulges': bulges, 'total': total_bulges, 'longest_bulge': longest_bulge, 'symmetric_bulge': symmetric_bulge, 'asymmetric_bulge': asymmetric_bulge}

    def calculate_GC(self):
        """Return the GC percentage of the sequence."""

        GC_count = sum(1 for x in self.seq if x in ['G', 'C'])

        if len(self.seq) == 0:

            return 0

        return float(GC_count) / len(self.seq) * 100

    def generate_nucleotide_composition(self):
        """Return the percentage of each nucleotide."""

        nuc = ['A', 'C', 'U', 'G']

        nuc_comp = {}

        for n in nuc:

            nuc_comp[n] = round(float(self.seq.count(n)/len(self.seq))*100, 2)

        return nuc_comp

    def generate_dinucleotide_composition(self):
        """Return the percentage of each dinucleotide."""

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
        """Return the percentage of each trinucleotide."""

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
        """Count middle-base and dot-bracket triplets (4 bases x 8 structure patterns)."""

        structure_string = list(TRIPLET_STRUCTURES)

        nuc = list(TRIPLET_BASES)

        structure_triplets = {}

        for i in nuc:

            for j in structure_string:

                structure_triplets[i+j] = 0

        if self.structure is None:
            raise RuntimeError("structure prediction did not produce a dot-bracket string")

        self.structure = self.structure.replace(')', '(')

        l = len(self.seq)-len(self.seq) % 3

        for i in range(0, l, 3):

            n = self.seq[i+1]+self.structure[i:i+3]

            if n in structure_triplets:

                structure_triplets[n] += 1

        return structure_triplets

    def generate_features(self):
        """Return the full feature dictionary, or None if the precursor has no usable structure."""

        nfeatures = {}

        if self.buldge_graph is None:

            return

        stems = self.find_stem_pairs()

        if stems['stem_length'] <= 1:

            self.drop_reason = 'fewer than 2 stem pairs'

            return

        loops = self.find_loops()
        bulges = self.find_bulges()
        mathed_stem = self.find_matched_stem()

        # Assemble all the features in one dictionary
        nfeatures['length'] = len(self.seq)
        nfeatures['mfe'] = self.mfe
        # Named GC_content so the GC dinucleotide below cannot overwrite it.
        nfeatures['GC_content'] = self.calculate_GC()
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
