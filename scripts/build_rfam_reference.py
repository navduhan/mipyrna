"""Prepare the labelled Rfam bundle offline; never called by the analysis workflow.

Run from the repository root with PYTHONPATH=.: python scripts/build_rfam_reference.py
    --cache /path/to/downloads --outdir /path/to/new/bundle --jobs 4
"""

import argparse
import gzip
import json
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import SeqIO

from mipyrna import utility as mu

RELEASE = "15.1"
BASE = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{RELEASE}/"
# Pin the metadata and its schema so the RNA-type selection cannot silently change.
METADATA = {
    "database_files/family.txt.gz": "5662343d24b6f8dab566977fabebeb9df4bbec47fb8991400749fc642f03df69",
    "database_files/family.sql": "f4c4b1bd0d3fc5d99ee2678cda6e5996800a7de77e2319bd597816ace45e34b0",
}


def fetch_file(name, cache):
    """Reuse a completed download; publish new downloads only after success."""
    path = cache / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = path.with_name(path.name + ".part")
        for attempt in range(3):
            try:
                mu.download_file(BASE + name, str(temporary))
                temporary.replace(path)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
    return path


def read_families(path):
    """Read family metadata; family FASTAs contain seed plus full-region records."""
    with gzip.open(path, "rt", encoding="latin1") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    return [
        {"accession": row[0], "name": row[1], "type": row[18], "expected_records": int(row[14]) + int(row[15])}
        for row in sorted(rows)
    ]


def add_family(family, path, sequences):
    """Merge exact sequences while retaining every non-miRNA family membership."""
    label = f"{family['accession']}:{family['name']}"
    count = 0
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).upper().replace("U", "T")
            if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
                raise ValueError(f"Invalid full-family sequence in {path}: {record.id}")
            sequences.setdefault(sequence, set()).add(label)
            count += 1
    if count == 0 and family["expected_records"]:
        raise ValueError(f"{family['accession']}: the downloaded family FASTA is empty")
    return count


def build_reference(cache, outdir, jobs=4):
    """Download non-miRNA full families and write a deterministic, checksummed archive."""
    outdir.mkdir(parents=True, exist_ok=False)
    cache.mkdir(parents=True, exist_ok=True)
    sources = []
    for name, expected in METADATA.items():
        path = fetch_file(name, cache)
        if mu.file_sha256(path) != expected:
            raise ValueError(f"Pinned Rfam metadata checksum mismatch: {path}")
        sources.append((name, path.stat().st_size, expected))
    families = read_families(cache / "database_files/family.txt.gz")
    selected = [f for f in families if "miRNA" not in f["type"].replace(";", " ").split()]
    print(f"Rfam {RELEASE}: {len(selected)} non-miRNA families; {len(families) - len(selected)} miRNA families excluded", flush=True)

    def download(family):
        return fetch_file(f"fasta_files/{family['accession']}.fa.gz", cache)

    sequences = {}
    records = 0
    count_differences = {}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for index, (family, path) in enumerate(zip(selected, pool.map(download, selected), strict=True), 1):
            count = add_family(family, path, sequences)
            family["source_records"] = count
            records += count
            if count != family["expected_records"]:
                count_differences[family["accession"]] = {"metadata": family["expected_records"], "fasta": count}
            sources.append((f"fasta_files/{path.name}", path.stat().st_size, mu.file_sha256(path)))
            if index % 100 == 0 or index == len(selected):
                print(f"Prepared {index}/{len(selected)} families; {records:,} source records, {len(sequences):,} unique sequences", flush=True)

    fasta = outdir / "ncRNA_rfam.fa"
    mu.write_fasta((
        (f"rfam_{index}|{'|'.join(sorted(sequences[sequence]))}", sequence)
        for index, sequence in enumerate(sorted(sequences), 1)
    ), str(fasta))
    with open(outdir / "rfam_families.tsv", "w") as handle:
        handle.write("accession\tname\ttype\tincluded\tmetadata_records\tsource_records\n")
        for family in families:
            included = "miRNA" not in family["type"].replace(";", " ").split()
            handle.write(f"{family['accession']}\t{family['name']}\t{family['type']}\t{int(included)}\t{family['expected_records']}\t{family.get('source_records', '-')}\n")
    for name in ("README", "COPYING"):
        source = fetch_file(name, cache)
        (outdir / name).write_bytes(source.read_bytes())
        sources.append((name, source.stat().st_size, mu.file_sha256(source)))
    with open(outdir / "rfam_sources.tsv", "w") as handle:
        handle.write("url\tbytes\tsha256\n")
        for name, size, checksum in sorted(sources):
            handle.write(f"{BASE + name}\t{size}\t{checksum}\n")
    manifest = {
        "release": RELEASE,
        "base_url": BASE,
        "policy": "all non-miRNA full-family sequences; exact sequence deduplication; all retained family memberships",
        "header_format": "rfam_<number>|RFxxxxx:family_name[|RFyyyyy:family_name...]",
        "included_families": len(selected),
        "excluded_miRNA_families": len(families) - len(selected),
        "source_records": records,
        "unique_sequences": len(sequences),
        "multi_family_sequences": sum(len(labels) > 1 for labels in sequences.values()),
        "metadata_count_differences": count_differences,
        "fasta_sha256": mu.file_sha256(fasta),
        "families_sha256": mu.file_sha256(outdir / "rfam_families.tsv"),
        "sources_sha256": mu.file_sha256(outdir / "rfam_sources.tsv"),
        "builder_sha256": mu.file_sha256(__file__),
    }
    with open(outdir / "rfam_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    archive = outdir / "ncRNA_rfam.tar.xz"
    with tarfile.open(archive, "w:xz", preset=6) as tar:
        for name in ("ncRNA_rfam.fa", "rfam_families.tsv", "rfam_sources.tsv", "rfam_manifest.json", "README", "COPYING"):
            info = tar.gettarinfo(str(outdir / name), arcname=name)
            info.uid = info.gid = info.mtime = 0
            info.uname = info.gname = ""
            info.mode = 0o644
            with open(outdir / name, "rb") as handle:
                tar.addfile(info, handle)
    checksum = mu.file_sha256(archive)
    (outdir / "ncRNA_rfam.tar.xz.sha256").write_text(f"{checksum}  {archive.name}\n")
    print(json.dumps(manifest, indent=2), flush=True)
    print(f"Bundle: {archive.stat().st_size:,} bytes; SHA-256 {checksum}", flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    build_reference(args.cache, args.outdir, args.jobs)
