"""Known miRNA Scan Module

This module detects known miRNAs directly from reads aligned to the known
precursors, independently of the discovery path. A mature is reported as
detected when enough reads start at its annotated 5-prime end, and the same
fold, minimum free energy, mature-position check and star estimate are added as
for novel candidates.

Features:
    - Resolution of every known mature to a precursor, from the miRBase hairpin
      set or from a folded window around its genomic locus
    - Alignment of the collapsed reads to the resolved precursors
    - Per-arm read counting with the modal 5-prime offset and its consistency
    - Structure, minimum free energy, mature-position check and star sequence
    - A ``known_scan.tsv`` table and a merge into the known miRNA table

Configuration:
    The constructor takes the known matures, the collapsed samples, the genome
    and the bowtie index of the genome, and the read threshold, 5-prime
    tolerance and unique-read switch. Results go to the output directory.

Dependencies:
    - Bowtie 1 through ``mipyrna.aligner``
    - ViennaRNA and forgi for folding through ``mipyrna.novel``
    - SAM parsing through ``mipyrna.reads``
    - pandas

Classes / Functions / Exceptions:
    - read_sources: Read the database of each known mature.
    - precursor_name_variants: Names a mature may have as a precursor.
    - merge_known_calls: Merge alignment calls into the known miRNA table.
    - drop_detected_novel: Remove novel candidates that the scan detected as known.
    - KnownScan: Resolve precursors, align reads, count arms and write the table.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os
import re
from collections import Counter, defaultdict

import pandas as pd

from mipyrna import model_store as ms
from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.features import Precursor_Features
from mipyrna.logger import MiPyRNALogger
from mipyrna.novel import (
    LOCAL_EXTENSION,
    describe_precursor,
    find_subseq,
    fold_hairpins,
)
from mipyrna.reads import Read_process

log = MiPyRNALogger(mode="a", log="known_scan")

# A read may extend this far past the annotated 3-prime end, matching the
# downstream window of the quantifier.
THREE_PRIME_SLACK = 5
# Share of the reads that must start at the modal 5-prime position.
MIN_CONSISTENCY = 0.5


def read_sources(path):
    """Read the database of every known mature from a two-column table.

    Args:
        path (str): Table with ``Known_ID`` and ``Database`` columns.

    Returns:
        dict: miRNA ID to comma-separated database names, or an empty dictionary
        when the table is missing.
    """

    if not path or not os.path.exists(path):
        return {}
    table = pd.read_csv(path, sep="\t")
    return dict(zip(table["Known_ID"].astype(str), table["Database"].astype(str)))


def precursor_name_variants(mature_id):
    """List the precursor names a mature ID may have.

    A miRBase mature such as ``ath-miR159a-3p`` belongs to ``ath-MIR159a``; a
    PmiREN mature such as ``Ath-miR159a`` may be written with either case.

    Args:
        mature_id (str): Mature miRNA ID.

    Returns:
        list: Candidate precursor names, most likely first.
    """

    base = re.sub(r"\.\d+$", "", str(mature_id))
    base = re.sub(r"-[35]p$", "", base)
    names = []
    for name in (base, base.replace("-miR", "-MIR"), base.replace("-MIR", "-miR")):
        for variant in (name, name[:1].upper() + name[1:], name[:1].lower() + name[1:]):
            if variant not in names:
                names.append(variant)
    return names


def _fold_window(reader, chrom, start, end, strand, mature, max_extension=LOCAL_EXTENSION):
    """Fold the genomic window around an annotated mature and pick its hairpin.

    The hairpins are folded by :func:`mipyrna.novel.fold_hairpins`. A hairpin
    whose mature passes the structure check wins over one that does not; between
    equals, the lower minimum free energy wins.

    Args:
        reader (Read_process): Reference reader.
        chrom (str): Chromosome.
        start (int): Leftmost genomic coordinate of the mature.
        end (int): Rightmost genomic coordinate of the mature (exclusive).
        strand (str): ``'+'`` or ``'-'``.
        mature (str): Mature sequence in the orientation of the strand.
        max_extension (int): Nucleotides folded on each side of the mature.

    Returns:
        dict or None: Precursor, structure and coordinates, or None if the
        window cannot be folded.
    """

    best = None
    for candidate in fold_hairpins(reader, chrom, start, end, strand, mature, max_extension):
        offset = candidate["hairpin_start"]
        candidate["mature"] = candidate["precursor"][offset : offset + len(mature)]
        candidate["hairpin"] = "genome"
        candidate["hairpin_source"] = "genome"
        if _better(candidate, best):
            best = candidate
    return best


def _better(candidate, best):
    """Say whether a folded candidate is a better precursor than the current best.

    A hairpin whose mature passes the structure check wins over one that does
    not; between equals, the lower minimum free energy wins.

    Args:
        candidate (dict): Folded candidate.
        best (dict): Current best candidate.

    Returns:
        bool: True when the candidate should replace the best.
    """

    if best is None:
        return True
    if (candidate["mature_check"] == "ok") != (best["mature_check"] == "ok"):
        return candidate["mature_check"] == "ok"
    return float(candidate["mfe"]) < float(best["mfe"])


def _db_precursor(mature_id, mature, precursor):
    """Build a known entry from a database precursor.

    Args:
        mature_id (str): Mature miRNA ID.
        mature (str): Mature sequence.
        precursor (str): Precursor sequence that holds the mature.

    Returns:
        dict or None: Precursor, structure and position, or None when the mature
        is not found in the precursor.
    """

    pos = precursor.find(mature)
    if pos < 0:
        pos = find_subseq(precursor, mature)
    if pos < 0:
        return None
    found = precursor[pos : pos + len(mature)]
    try:
        entry = describe_precursor(precursor, found, start=pos)
    except (TypeError, ValueError, RuntimeError) as exc:
        log.warning(f"Could not fold the precursor of {mature_id}: {exc}")
        return None
    entry.update(
        mature=found,
        mature_length=len(found),
        precursor=precursor,
        precursor_length=len(precursor),
        hairpin_start=pos,
        hairpin_source="miRBase",
    )
    return entry


class KnownScan:
    """Detect known miRNAs from reads aligned to known precursors."""

    def __init__(
        self,
        matures=None,
        samples=None,
        genome=None,
        index_dir=None,
        hairpin_fasta=None,
        outdir=".",
        sources=None,
        min_reads=10,
        tolerance=2,
        unique=False,
        cpu=8,
        mem=20,
        jobs=1,
        species_type="plants",
        rfam=None,
    ):
        """Initialize the scan.

        Args:
            matures (iterable): ``(Known_ID, mature)`` pairs, the known reference.
            samples (dict): Sample name to ``[group, replicate, collapsed reads path]``.
            genome (str): Reference genome FASTA.
            index_dir (str): Directory holding the bowtie index of the genome.
            hairpin_fasta (str, optional): miRBase hairpin FASTA of the species.
            outdir (str): Output directory.
            sources (dict, optional): miRNA ID to comma-separated database names.
            min_reads (int): Minimum reads on the precursor arm. Defaults to 10.
            tolerance (int): Allowed offset of a read 5-prime end from the annotated start. Defaults to 2.
            unique (bool): Count a read only when it maps to one precursor. Defaults to False.
            cpu (int): Threads for alignment. Defaults to 8.
            mem (int): Memory in GB for alignment. Defaults to 20.
            jobs (int): Samples aligned at once. Defaults to 1.
            species_type (str): ``'plants'`` or ``'animals'``. Defaults to ``'plants'``.
            rfam (dict, optional): Sample name to a dictionary of read name to the set of
                Rfam sequences the read matched. Defaults to None.
        """

        self.matures = list(matures or [])
        self.samples = samples or {}
        self.genome = genome
        self.index_dir = index_dir or ""
        self.hairpin_fasta = hairpin_fasta
        self.outdir = outdir
        self.sources = sources or {}
        self.min_reads = min_reads
        self.tolerance = tolerance
        self.unique = unique
        self.cpu = cpu
        self.mem = mem
        self.jobs = jobs
        self.species_type = species_type
        self.rfam = rfam or {}
        self.reference = os.path.join(outdir, "known_hairpins.fa")

    def _genome_hits(self, records, name):
        """Align sequences to the genome and return their loci.

        The queries are aligned to the genome index of the alignment stage with
        :class:`mipyrna.aligner.Bowtie_Aligner`, and the SAM output is read with
        :class:`mipyrna.reads.Read_process`.

        Args:
            records (list): ``(name, sequence)`` pairs.
            name (str): Prefix for the query, alignment and log files.

        Returns:
            dict: Query name to a list of ``(chrom, start, end, strand)`` tuples.
        """

        if not records:
            return {}
        query = mu.write_fasta(records, os.path.join(self.outdir, f"{name}_queries.fa"))
        config = str(ir.files("mipyrna").joinpath("param/bowtie_query.ini"))
        query_dir = mu.make_directory(os.path.join(self.outdir, f"known_{name}"))
        index = os.path.join(query_dir, "bowtie_index")
        if not os.path.exists(index):
            os.symlink(os.path.join(self.index_dir, "bowtie_index"), index)
        aln = Bowtie_Aligner(
            ref_genome=self.genome, configFile=config, outdir=query_dir, slurm=False
        )
        with mu.without_array():
            sam = aln.run_alignment(
                samplesDict={"query": [name, "0", query]},
                fileType="fasta",
                outType="SAM",
                cpu=self.cpu,
                mem=self.mem,
                jobs=1,
            )["query"][2]
        aligned = Read_process().aligned_reads(
            {"query": [name, "0", sam]}, alignType="SAM", jobs=1
        )
        hits = defaultdict(list)
        for row in aligned["query"][2].itertuples():
            hits[row.read_name].append((row.name, int(row.start), int(row.end), row.strand))
        for loci in hits.values():
            loci.sort()
        return hits

    def _database_precursors(self, matures):
        """Resolve matures against the miRBase hairpins and place them on the genome.

        Args:
            matures (list): ``(Known_ID, mature)`` pairs.

        Returns:
            tuple: Resolved entries and the matures that had no usable hairpin.
        """

        hairpins = (
            dict(mu.read_fasta(self.hairpin_fasta))
            if self.hairpin_fasta and os.path.exists(self.hairpin_fasta)
            else {}
        )
        entries = {}
        pending = []
        for mid, mature in matures:
            precursor = None
            for variant in precursor_name_variants(mid):
                if variant in hairpins:
                    precursor = hairpins[variant]
                    break
            if precursor is None:
                pending.append((mid, mature))
                continue
            entry = _db_precursor(mid, mature, precursor)
            if entry is None:
                pending.append((mid, mature))
                continue
            entry["hairpin"] = next(
                (v for v in precursor_name_variants(mid) if v in hairpins), mid
            )
            entries[mid] = entry

        queries = {}
        paralog_hairpins = defaultdict(list)
        for mid, entry in entries.items():
            queries[entry["hairpin"]] = entry["precursor"]
            prefix = f"{mid.split('-', 1)[0].lower()}-" if "-" in mid else ""
            mseq = entry["mature"]
            for hname, hseq in hairpins.items():
                if hname == entry["hairpin"] or (prefix and not hname.lower().startswith(prefix)):
                    continue
                pos = hseq.find(mseq)
                if pos >= 0:
                    queries[hname] = hseq
                    paralog_hairpins[mid].append((hname, pos))
        loci = self._genome_hits(list(queries.items()), "precursor")
        self._all_loci = {}
        for mid, entry in entries.items():
            hits = loci.get(entry["hairpin"], [])
            extra_hits = []
            for hname, p_offset in paralog_hairpins.get(mid, ()):
                for hchrom, hstart, hend, hstrand in loci.get(hname, []):
                    hmature = hstart + p_offset if hstrand == "+" else hend - p_offset - entry["mature_length"]
                    extra_hits.append((hchrom, hstrand, hmature, hstart, hend))
            if not hits:
                self._all_loci[mid] = extra_hits
                entry.update(chrom="", start="", end="", mature_start="", strand="")
                continue
            offset = entry["hairpin_start"]
            resolved_hits = []
            for hchrom, hstart, hend, hstrand in hits:
                hmature = hstart + offset if hstrand == "+" else hend - offset - entry["mature_length"]
                resolved_hits.append((hchrom, hstrand, hmature, hstart, hend))
            self._all_loci[mid] = resolved_hits + extra_hits
            chrom, strand, mature_start, start, end = resolved_hits[0]
            entry.update(
                chrom=chrom,
                start=start,
                end=end,
                mature_start=mature_start,
                strand=strand,
            )
        return entries, pending

    def _genome_precursors(self, matures):
        """Resolve matures through their genomic locus.

        Args:
            matures (list): ``(Known_ID, mature)`` pairs.

        Returns:
            dict: miRNA ID to the folded precursor entry.
        """

        hits = self._genome_hits(matures, "mature")
        reader = Read_process(reference=self.genome)
        entries = {}
        for mid, mature in matures:
            best = None
            for chrom, start, end, strand in hits.get(mid, [])[:20]:
                try:
                    entry = _fold_window(reader, chrom, start, end, strand, mature)
                except Exception as exc:
                    log.warning(f"Could not place {mid} at {chrom}:{start}: {exc}")
                    continue
                if entry is None:
                    continue
                if _better(entry, best):
                    best = entry
            if best is None:
                log.warning(f"No precursor found for {mid}")
                continue
            entries[mid] = best
        return entries

    def resolve(self):
        """Resolve every known mature to a folded precursor.

        Returns:
            pandas.DataFrame: One row per resolved mature with its precursor,
            structure, free energy, mature-position check and coordinates.
        """

        entries, pending = self._database_precursors(self.matures)
        log.info(
            f"Known scan: {len(entries)} precursors from miRBase, {len(pending)} from the genome"
        )
        entries.update(self._genome_precursors(pending))

        rows = []
        for mid, _mature in self.matures:
            entry = entries.get(mid)
            if entry is None:
                continue
            row = {"Known_ID": mid, "Database": self.sources.get(mid, "")}
            row.update(entry)
            rows.append(row)
        table = pd.DataFrame(rows)
        if len(table):
            table["hairpin_start"] = table["hairpin_start"].astype(int)
        log.info(f"Known scan: {len(table)} of {len(self.matures)} matures resolved")
        return table

    def _align(self, table):
        """Align the collapsed reads to the resolved precursors.

        Args:
            table (pandas.DataFrame): Resolved precursors.

        Returns:
            dict: Sample name to ``[group, replicate, SAM path]``.
        """

        mu.write_fasta(
            [(mid, seq) for mid, seq in zip(table["Known_ID"], table["precursor"])],
            self.reference,
        )
        config = str(ir.files("mipyrna").joinpath("param/bowtie_mirna.ini"))
        align_dir = mu.make_directory(os.path.join(self.outdir, "known_align"))
        aln = Bowtie_Aligner(
            ref_genome=self.reference, configFile=config, outdir=align_dir, slurm=False
        )
        aln.build_index()
        return aln.run_alignment(
            samplesDict=self.samples, outType="SAM", cpu=self.cpu, mem=self.mem, jobs=self.jobs
        )

    def _count(self, table, aligned, star_windows):
        """Count the reads that start on the mature and star arms.

        Args:
            table (pandas.DataFrame): Resolved precursors.
            aligned (dict): Sample name to the aligned reads DataFrame.
            star_windows (dict): miRNA ID to the star position and length.

        Returns:
            pandas.DataFrame: Per mature and per sample counts, the modal
            5-prime offset, its consistency, the modal read length and the call.
        """

        windows = {
            row.Known_ID: (int(row.hairpin_start), len(str(row.mature)))
            for row in table.itertuples()
        }
        peak = defaultdict(lambda: defaultdict(float))
        unique = defaultdict(lambda: defaultdict(float))
        star = defaultdict(lambda: defaultdict(float))
        rfam_copies = defaultdict(lambda: defaultdict(float))
        rfam_ids = defaultdict(set)
        starts = defaultdict(Counter)
        lengths = defaultdict(Counter)

        for sample, frame in aligned.items():
            per_read = defaultdict(dict)
            per_star = defaultdict(dict)
            for row in frame.itertuples():
                mid = row.name
                window = windows.get(mid)
                if window is None:
                    continue
                start, end = int(row.start), int(row.end)
                copies = int(row.reads)
                star_window = star_windows.get(mid)
                if star_window is not None:
                    sposition, slength = star_window
                    if (
                        abs(start - sposition) <= self.tolerance
                        and end <= sposition + slength + THREE_PRIME_SLACK
                    ):
                        per_star[row.read_name][mid] = copies
                position, length = window
                offset = start - position
                if abs(offset) > self.tolerance:
                    continue
                if end > position + length + THREE_PRIME_SLACK:
                    continue
                per_read[row.read_name][mid] = (offset, end - start, copies)

            for mids in per_star.values():
                for mid, copies in mids.items():
                    star[sample][mid] += copies

            for read_name, mids in per_read.items():
                hit = self.rfam.get(sample, {}).get(read_name)
                for mid, (offset, length, copies) in mids.items():
                    peak[sample][mid] += copies
                    starts[mid][offset] += copies
                    lengths[mid][length] += copies
                    if len(mids) == 1:
                        unique[sample][mid] += copies
                    if hit:
                        rfam_copies[sample][mid] += copies
                        rfam_ids[mid].update(hit)

        rows = []
        samples = list(aligned)
        for row in table.itertuples():
            mid = row.Known_ID
            row_counts = {sample: peak[sample].get(mid, 0.0) for sample in samples}
            row_unique = {sample: unique[sample].get(mid, 0.0) for sample in samples}
            total = sum(row_counts.values())
            total_unique = sum(row_unique.values())
            observed = starts.get(mid, Counter())
            if observed:
                modal, modal_count = min(
                    observed.items(), key=lambda item: (-item[1], abs(item[0]), item[0])
                )
            else:
                modal, modal_count = 0, 0
            consistency = modal_count / total if total else 0.0
            length_mode = (
                min(lengths[mid].items(), key=lambda item: (-item[1], item[0]))[0]
                if lengths.get(mid)
                else 0
            )
            basis = total_unique if self.unique else total
            detected = (
                basis >= self.min_reads
                and abs(modal) <= self.tolerance
                and consistency >= MIN_CONSISTENCY
            )
            record = {
                "Known_ID": mid,
                "hairpin": row.hairpin,
                "hairpin_source": row.hairpin_source,
                "Database": row.Database,
                "precursor": row.precursor,
                "precursor_length": row.precursor_length,
                "mature": row.mature,
                "mature_length": row.mature_length,
                "chrom": row.chrom,
                "start": row.start,
                "end": row.end,
                "mature_start": row.mature_start,
                "strand": row.strand,
                "struct": row.struct,
                "mfe": row.mfe,
                "mature_check": row.mature_check,
                "position": row.position,
                "star": row.star,
                "total_reads": total,
                "unique_reads": total_unique,
                "max_reads": max(row_counts.values()) if row_counts else 0.0,
                "star_reads": sum(star[sample].get(mid, 0.0) for sample in samples),
                "rfam_reads": sum(rfam_copies[sample].get(mid, 0.0) for sample in samples),
                "Rfam": ",".join(sorted(rfam_ids.get(mid, ()))) or "-",
                "modal_5p_offset": modal,
                "five_prime_consistency": consistency,
                "length_mode": length_mode,
                "detected": detected,
            }
            for sample in samples:
                record[f"{sample}_reads"] = row_counts[sample]
            rows.append(record)
        return pd.DataFrame(rows)

    def scan(self):
        """Resolve precursors, align the reads, count them and write the table.

        Returns:
            pandas.DataFrame: The table written to ``known_scan.tsv``.
        """

        log.info(f"Known scan: resolving {len(self.matures)} known matures")
        table = self.resolve()
        if len(table) == 0:
            log.warning("Known scan: no precursor could be resolved")
            table["rfam_reads"] = pd.Series(dtype=float)
            table["Rfam"] = pd.Series(dtype=str)
            table.to_csv(os.path.join(self.outdir, "known_scan.tsv"), sep="\t", index=False)
            return table

        aligned = self._align(table)
        star_windows = {}
        for known_id, precursor, star in zip(table["Known_ID"], table["precursor"], table["star"]):
            if isinstance(star, str) and star:
                position = str(precursor).find(star)
                if position >= 0:
                    star_windows[known_id] = (position, len(star))

        frames = Read_process().aligned_reads(aligned, alignType="SAM", jobs=self.jobs)
        counts = self._count(table, {name: value[2] for name, value in frames.items()}, star_windows)
        counts = self._score_structure(counts)
        counts.attrs["all_loci"] = dict(getattr(self, "_all_loci", {}))
        counts.to_csv(os.path.join(self.outdir, "known_scan.tsv"), sep="\t", index=False)
        log.info(
            f"Known scan: {int(counts['detected'].sum())} of {len(counts)} known matures called"
        )
        return counts

    def _score_structure(self, counts):
        """Score the resolved precursors with the bundled model.

        Args:
            counts (pandas.DataFrame): Scan table.

        Returns:
            pandas.DataFrame: The table with ``ml_score``, ``score`` and ``features`` columns.
        """

        features = []
        for precursor in counts["precursor"]:
            try:
                features.append(Precursor_Features(str(precursor)).generate_features())
            except Exception:
                features.append(None)
        valid = [i for i, feature in enumerate(features) if feature is not None]
        counts["features"] = [feature if feature is not None else "" for feature in features]
        counts["ml_score"] = float("nan")
        if valid:
            try:
                predictions = ms.predict([features[i] for i in valid], self.species_type)
                counts.loc[counts.index[valid], "ml_score"] = predictions
            except Exception as exc:
                log.warning(f"Known scan: model scoring skipped: {exc}")
        counts["score"] = counts["ml_score"]
        return counts


def _to_int(value):
    """Convert a coordinate value to int, or return None when missing or empty."""
    if value is None or isinstance(value, bool) or pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def merge_known_calls(known_df, scan_df, species):
    """Merge the alignment calls into the known miRNA table.

    A locus found by both paths keeps the discovery row and is marked ``both``.
    A locus found only by alignment is appended, and a detected name already
    carried by an existing row, as an alias, same mature sequence or on the same
    arm of the same precursor, is counted as a sibling instead of becoming a
    duplicate row.

    Args:
        known_df (pandas.DataFrame): Known miRNAs of the discovery path.
        scan_df (pandas.DataFrame): Table written by :meth:`KnownScan.scan`.
        species (str): miRBase species code, used to build IDs of new rows.

    Returns:
        tuple: Merged table, number of added rows, number of loci found by both
        paths and number of detected names already carried by an existing row.
    """

    columns = list(known_df.columns) if isinstance(known_df, pd.DataFrame) else []
    merged = known_df.copy() if isinstance(known_df, pd.DataFrame) else pd.DataFrame(columns=columns)
    detected = scan_df[scan_df["detected"].astype(bool)] if len(scan_df) else scan_df
    lookup = {row.Known_ID: row for row in detected.itertuples()} if len(detected) else {}

    merged["Detection"] = "discovery"
    merged["alignment_reads"] = 0.0
    if len(merged):
        for index, known_id in zip(merged.index, merged["Known_ID"]):
            call = lookup.get(known_id)
            if call is not None:
                merged.at[index, "Detection"] = "both"
                merged.at[index, "alignment_reads"] = float(call.total_reads)

    both = sum((merged["Detection"] == "both").tolist()) if len(merged) else 0
    primary = set(merged["Known_ID"].astype(str)) if len(merged) else set()
    aliases = {}
    if "Known_IDs" in merged.columns and len(merged):
        for index, names in zip(merged.index, merged["Known_IDs"].astype(str)):
            for name in names.split(";"):
                if name and name not in aliases:
                    aliases[name] = index

    has_arm = "position" in merged.columns and "position" in detected.columns
    by_precursor_arm = {}
    if len(merged) and "precursor" in merged.columns:
        for index, row in merged.iterrows():
            arm = str(row["position"]) if has_arm else ""
            by_precursor_arm[(str(row["precursor"]), arm)] = index
    by_mature = (
        {str(row["mature"]): index for index, row in merged.iterrows()}
        if len(merged) and "mature" in merged.columns
        else {}
    )

    added = []
    added_by_mature = {}
    added_by_precursor_arm = {}
    added_ranks = {}
    siblings = 0
    counted = set()
    for call in detected.itertuples() if len(detected) else []:
        if call.Known_ID in primary or call.Known_ID in counted:
            continue
        counted.add(call.Known_ID)
        if call.Known_ID in aliases:
            siblings += 1
            continue
        call_arm = str(getattr(call, "position", "")) if has_arm else ""
        precursor_key = (str(call.precursor), call_arm)
        sibling = by_mature.get(str(call.mature))
        if sibling is None:
            sibling = by_precursor_arm.get(precursor_key)
        if sibling is not None:
            if "Known_IDs" in merged.columns:
                names = str(merged.at[sibling, "Known_IDs"])
                if call.Known_ID not in names.split(";"):
                    merged.at[sibling, "Known_IDs"] = f"{names};{call.Known_ID}"
            siblings += 1
            continue

        rank = (
            float(getattr(call, "unique_reads", 0.0)),
            float(getattr(call, "star_reads", 0.0)),
            float(getattr(call, "total_reads", 0.0)),
        )
        staged_idx = added_by_mature.get(str(call.mature))
        if staged_idx is None:
            staged_idx = added_by_precursor_arm.get(precursor_key)
        if staged_idx is not None:
            staged = added[staged_idx]
            existing_ids = staged["Known_IDs"].split(";")
            if call.Known_ID not in existing_ids:
                existing_ids.append(call.Known_ID)
            if rank > added_ranks[staged_idx]:
                existing_ids = [call.Known_ID, *[i for i in existing_ids if i != call.Known_ID]]
                identifier = f"{species}_miPyRNA_{call.chrom}_miR_{len(merged) + staged_idx + 1}"
                staged.update(
                    Novel_ID=identifier,
                    Known_ID=call.Known_ID,
                    chrom=call.chrom,
                    mature=call.mature,
                    mature_length=call.mature_length,
                    precursor=call.precursor,
                    precursor_length=call.precursor_length,
                    start=call.start,
                    end=call.end,
                    mature_start=call.mature_start,
                    strand=call.strand,
                    struct=call.struct,
                    mfe=call.mfe,
                    mature_check=call.mature_check,
                    features=call.features,
                    mature_reads=call.total_reads,
                    position=call.position,
                    star_reads=call.star_reads,
                    five_prime_consistency=float(
                        getattr(call, "five_prime_consistency", float("nan"))
                    ),
                    star=call.star,
                    rfam_reads=float(getattr(call, "rfam_reads", 0.0)),
                    Rfam=str(getattr(call, "Rfam", "-")),
                    score=call.score,
                    ml_score=call.ml_score,
                    Database=call.Database,
                    raw_ID=call.Known_ID,
                    alignment_reads=float(call.total_reads),
                )
                added_ranks[staged_idx] = rank
                added_by_precursor_arm[precursor_key] = staged_idx
            staged["Known_IDs"] = ";".join(existing_ids)
            siblings += 1
            continue

        staged_idx = len(added)
        number = len(merged) + staged_idx + 1
        identifier = f"{species}_miPyRNA_{call.chrom}_miR_{number}"
        added.append(
            {
                "Novel_ID": identifier,
                "Known_ID": call.Known_ID,
                "Known_IDs": call.Known_ID,
                "chrom": call.chrom,
                "mature": call.mature,
                "mature_length": call.mature_length,
                "precursor": call.precursor,
                "precursor_length": call.precursor_length,
                "start": call.start,
                "end": call.end,
                "mature_start": call.mature_start,
                "strand": call.strand,
                "struct": call.struct,
                "mfe": call.mfe,
                "mature_check": call.mature_check,
                "features": call.features,
                "mature_reads": call.total_reads,
                "position": call.position,
                "star_reads": call.star_reads,
                "inconsistent_fraction": float("nan"),
                "five_prime_consistency": float(
                    getattr(call, "five_prime_consistency", float("nan"))
                ),
                "star": call.star,
                "randfold_pvalue": float("nan"),
                "cluster": float("nan"),
                "rfam_reads": float(getattr(call, "rfam_reads", 0.0)),
                "Rfam": str(getattr(call, "Rfam", "-")),
                "score": call.score,
                "ml_score": call.ml_score,
                "heuristic_score": float("nan"),
                "Database": call.Database,
                "raw_ID": call.Known_ID,
                "Detection": "alignment",
                "alignment_reads": float(call.total_reads),
            }
        )
        added_by_mature[str(call.mature)] = staged_idx
        added_by_precursor_arm[precursor_key] = staged_idx
        added_ranks[staged_idx] = rank

    if added:
        merged = pd.concat([merged, pd.DataFrame(added)], ignore_index=True)
    return merged, len(added), both, siblings


def drop_detected_novel(novel_df, scan_df, window=5):
    """Remove novel candidates that the scan detected as a known miRNA.

    A candidate is dropped when it sits at the same locus as a detected known
    mature or overlaps a detected known precursor hairpin, or when it shares the
    precursor sequence of one.

    Args:
        novel_df (pandas.DataFrame): Novel miRNAs of the discovery path.
        scan_df (pandas.DataFrame): Table written by :meth:`KnownScan.scan`.
        window (int): Nucleotides of allowed difference between mature starts. Defaults to 5.

    Returns:
        tuple: Filtered novel table and the number of removed candidates.
    """

    if not isinstance(novel_df, pd.DataFrame) or len(novel_df) == 0 or len(scan_df) == 0:
        return novel_df, 0
    detected = scan_df[scan_df["detected"].astype(bool)]
    if len(detected) == 0:
        return novel_df, 0

    loci = defaultdict(list)
    hairpin_intervals = defaultdict(list)
    extra_loci = scan_df.attrs.get("all_loci", {}) if hasattr(scan_df, "attrs") else {}
    for row in detected.itertuples():
        chrom = str(getattr(row, "chrom", "") or "").strip()
        strand = str(getattr(row, "strand", "") or "").strip()
        mstart = _to_int(getattr(row, "mature_start", None))
        hstart = _to_int(getattr(row, "start", None))
        hend = _to_int(getattr(row, "end", None))
        if chrom and strand and mstart is not None:
            loci[(chrom, strand)].append(mstart)
        if chrom and strand and hstart is not None and hend is not None:
            hairpin_intervals[(chrom, strand)].append((hstart, hend))
        for echrom, estrand, emstart, ehstart, ehend in extra_loci.get(getattr(row, "Known_ID", ""), ()):
            loci[(str(echrom), str(estrand))].append(int(emstart))
            hairpin_intervals[(str(echrom), str(estrand))].append((int(ehstart), int(ehend)))
    precursors = (
        set(detected["precursor"].dropna().astype(str))
        if "precursor" in detected.columns
        else set()
    )

    keep = []
    removed = 0
    for index, row in novel_df.iterrows():
        chrom = str(row.get("chrom", "") or "").strip()
        strand = str(row.get("strand", "") or "").strip()
        start = _to_int(row.get("mature_start"))
        nstart = _to_int(row.get("start"))
        nend = _to_int(row.get("end"))
        precursor = str(row.get("precursor", ""))
        hits = loci.get((chrom, strand), [])
        near = start is not None and any(abs(start - position) <= window for position in hits)
        overlaps_hairpin = (
            nstart is not None
            and nend is not None
            and any(a <= nend and b >= nstart for a, b in hairpin_intervals.get((chrom, strand), []))
        )
        if near or overlaps_hairpin or (precursor and precursor in precursors):
            removed += 1
            continue
        keep.append(index)
    return novel_df.loc[keep].reset_index(drop=True), removed
