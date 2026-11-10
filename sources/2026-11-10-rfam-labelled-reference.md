# Labelled Rfam reference options

Verified: 10 November 2026. Initial research and the subsequently approved
preparation/validation are recorded below. No GSE13605 benchmark was launched.

## Selection criteria

Official, redistributable Rfam sequence data with stable family accessions,
human-readable family names and RNA types; suitable for the existing Bowtie
read-screening path. The reference must distinguish miRNAs from other ncRNAs.

## Verified sources

- Release README: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/README
- Family FASTAs: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/
- Family metadata: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/family.txt.gz
- Metadata schema: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/family.sql
- Seed alignments: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.seed.gz
- Stockholm annotation definitions: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/USERMAN
- Redistribution terms: https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/COPYING
- Database documentation: https://docs.rfam.org/en/latest/database.html

The current official release is 15.1 (January 2026), with 4,227 families.
Rfam is distributed under CC0, permitting redistribution in a software bundle.
Use the versioned URLs above, not the mutable CURRENT directory, for a build.

## Findings

### Full family FASTAs plus family metadata: preferred for exact read screening

The FTP directory supplies a separate compressed FASTA for each family, named
`RFxxxxx.fa.gz`. The original sequence headers need not include the family name:
membership is supplied by the filename. Join that accession to `family.txt.gz`
using its schema to recover `rfam_id`, description and RNA type.

Examples verified against the downloaded metadata:

| Accession | Family name | Type |
| --- | --- | --- |
| RF00001 | 5S_rRNA | Gene; rRNA; |
| RF00002 | 5_8S_rRNA | Gene; rRNA; |
| RF00003 | U1 | Gene; snRNA; splicing; |
| RF00005 | tRNA | Gene; tRNA; |
| RF00027 | let-7 | Gene; miRNA; |

The small full-family file `RF00002.fa.gz` was downloaded and checked. Its first
header is `>L78065.1/3758-3910`. A derived bundle could prefix this header with
`RF00002|5_8S_rRNA|` and report `RF00002:5_8S_rRNA` in the existing Rfam column.

Metadata contains 1,598 miRNA families and 2,629 other families. Never use the
unfiltered complete Rfam collection as a contaminant-removal reference, because
it includes genuine miRNA precursors. Exclude by the metadata type, not a
family-name heuristic. The precise inclusion policy beyond excluding miRNAs
remains a user decision: all other types, or a narrower contaminant list.

Summing the rounded sizes in the official FASTA directory gives approximately
423 MiB compressed for all families, or 395 MiB excluding miRNA families. The
rRNA/tRNA/snRNA subset, including snoRNAs, is approximately 316 MiB. These are
rough source-download sizes, not final bundle sizes or pipeline benchmark
measurements. Only the small samples were downloaded during initial research;
the subsequent full preparation is recorded below.

### Seed alignments: compact but less sequence coverage

`Rfam.seed.gz` is 5,928,432 bytes compressed. It contains annotated Stockholm
alignments with family accession (`AC`), name (`ID`), description (`DE`), type
(`TP`) and curated representative sequences. A labelled FASTA can be generated
by removing alignment gaps and normalising U to T.

This is a convenient compact option, but exact short-read screening against
seed sequences alone can miss variants present in the full family FASTAs.
It is not an equivalent drop-in replacement for full sequence coverage.

## Original bundle preparation proposal (subsequently implemented below)

1. Pin release 15.1 and the desired RNA-type inclusion policy.
2. Download selected full-family FASTAs and metadata from that release.
3. Normalise sequences, deduplicate exact sequences and retain every family
   membership when an identical sequence occurs in several families.
4. Bundle one compressed FASTA with labelled identifiers, a small family lookup
   table if needed, the licence and a provenance/checksum manifest.
5. Update the existing Rfam module to display family labels. Keep preparation
   separate from pipeline execution; no runtime network calls are required.

The current one-hit Bowtie screen would still report representative matching
families, not necessarily every possible family for an ambiguous short read.
Identical-sequence deduplication must not silently discard family memberships.
Exact matches provide an ncRNA alert, not a covariance-model family assignment.

The original unlabelled bundle has no verified release provenance in this
investigation. Do not join it to current coordinates and assume all labels are
correct; rebuilding from a pinned, labelled source is more reproducible.

## Reproducible retrieval and checksums

Only these small files were downloaded into the approved temporary research
directory; they were not added to the package data:

```sh
curl -fL https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/family.sql -o family.sql
curl -fL https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/family.txt.gz -o family.txt.gz
curl -fL https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.seed.gz -o Rfam.seed.gz
curl -fL https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/fasta_files/RF00002.fa.gz -o RF00002.fa.gz
shasum -a 256 family.sql family.txt.gz Rfam.seed.gz RF00002.fa.gz
```

| File | SHA-256 |
| --- | --- |
| family.sql | f4c4b1bd0d3fc5d99ee2678cda6e5996800a7de77e2319bd597816ace45e34b0 |
| family.txt.gz | 5662343d24b6f8dab566977fabebeb9df4bbec47fb8991400749fc642f03df69 |
| Rfam.seed.gz | 41f014f4ab5628620935f089b1bb77987ebc162fda7c72bff368a32c9cccabac |
| RF00002.fa.gz | bf9dc7aac7cc1b52a7ff117c6174ec2f13bd6a09c03a65c57915d03360257e71 |

The family SQL dump declares latin1. Read its tab-delimited data accordingly;
do not assume UTF-8. Family accession, name, description and type are columns
1, 2, 4 and 19, respectively, in this release's schema.

## Approved preparation follow-up

The user approved all non-miRNA full families, with exact-sequence deduplication
and all retained family memberships. `scripts/build_rfam_reference.py` performs
this preparation separately from the workflow. It pins the family metadata and
schema checksums, checks complete gzip streams, and records source checksums.
The archive carries a labelled FASTA, family table, source manifest, build
manifest, upstream README and CC0 licence. No runtime downloads are introduced.

Family FASTAs include both seed and full-region records. Their counts are not
always identical to the metadata's `num_seed + num_full`: RF00004 (U2), for
example, has 71,722 FASTA records versus 71,720 in the metadata. These are
upstream files from the same pinned release, not a reason to truncate the FASTA.
The builder records metadata and actual counts separately and lists differences
in the manifest. A truncated gzip stream remains a hard error.

The label-aware pipeline passed 15 tests on iris in
`$B/rfam_label_verify.NBmslL`. One test uses real Bowtie with eight matching
reference entries, verifying that one representative hit is retained rather
than suppressed. That check also identified that Bowtie rejects `--strata`
with `-k 1`; the redundant option was removed. All matches are exact (`-v 0`),
so every reported alignment already belongs to the same mismatch stratum.

These are preparation and focused software checks, not full GSE13605 benchmark
results. Preparation is isolated at `$B/rfam15_1_build`; live HPC source files,
the shared environment and published results are unchanged.

## Full build and content audit on iris

The completed build is at `$B/rfam15_1_build/bundle_v2`. Its independent
streaming audit passed, checking every FASTA record for strict sequence order
(therefore no duplicates), normalised bases, valid included-family labels and
the absence of excluded-family labels. All 2,629 selected families are present.

| Item | Result |
| --- | --- |
| Included non-miRNA families | 2,629 |
| Excluded miRNA families | 1,598 |
| Downloaded FASTA records | 9,402,960 |
| Unique normalised sequences | 5,763,899 |
| Sequences shared by multiple retained families | 82 |
| Maximum FASTA identifier length | 86 characters |
| Families with metadata/FASTA count differences | 182 |
| Initial gzip bundle size | 124,929,615 bytes (119.14 MiB) |
| Initial gzip bundle SHA-256 | d788fd276ce4fb222677b8f3063d7411af33b7e1e86e6aeac27c5c9cea296b05 |
| FASTA SHA-256 | c6d6f7858259f6a2c0d72b09a62f000100e61b617dc1d534d1c1f8561d228066 |

Family-table and source-manifest checksums also passed. The detailed discrepancies
are recorded, without discarding source records. Because the gzip archive exceeds
the normal GitHub per-file limit, the final builder uses standard-library XZ
compression. An intermediate gzip/XZ comparison verified identical decompressed
tar bytes. No sequences or memberships were removed to reduce size.

## Final shipped archive

The final builder was rerun from the cached versioned sources into
`$B/rfam15_1_build/bundle_final`. Its sequence, family-table and source-manifest
checksums and all content statistics match the independently audited build.

- Package archive: `mipyrna/data/ncRNA_rfam.tar.xz`
- Archive size: 80,693,572 bytes (76.96 MiB)
- Archive SHA-256: `d60005af30d7d1efd4ad7538b56021be641bee75d187115553e50fe89252d85f`
- Final builder SHA-256: `abfdd0504d2065742e04a941f17ac875d8686632a9afc09abac6d8aaf593e38f`
- Archive checksum sidecar: `mipyrna/data/ncRNA_rfam.tar.xz.sha256`
- Members: `ncRNA_rfam.fa`, `rfam_families.tsv`, `rfam_sources.tsv`,
  `rfam_manifest.json`, upstream `README` and `COPYING`.

The old unlabelled `.tar.gz` package resource was replaced, not retained as a
second copy. XZ needs no additional runtime dependency. The runtime verifies the
extracted FASTA checksum and expands shared sequence labels before read-level
counts are deduplicated. Resume rejects pre-labelled-reference discovery caches.

All 15 focused tests passed for the XZ version on iris at
`$B/rfam_xz_verify.keaaZY`, including real Bowtie screening. A separate test
exercised the runtime reader on the actual final 77 MiB archive: full extraction,
FASTA checksum verification, manifest loading and empty-hit handling all passed;
alignment itself was mocked for that full-archive extraction check.

The final archive was copied locally and its checksum verified before replacing
the package resource. Live HPC pipeline code, the shared environment and all
published results remain unchanged. Full GSE13605 comparisons remain pending.
