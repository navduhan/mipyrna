"""Environment checks and run provenance for miPyRNA.

:Author: Naveen Duhan
"""

import shutil
import sys
from importlib import metadata

REQUIRED_TOOLS = ('bowtie', 'bowtie-build', 'faSomeRecords', 'makeblastdb', 'blastn', 'miranda')
OPTIONAL_TOOLS = ('fastqc', 'trimmomatic', 'flexbar', 'trim_galore')

# Versions the bundled models were trained and validated with.
PINNED = {'tensorflow': '2.20.0', 'keras': '3.11.3', 'numpy': '2.4.6',
          'ViennaRNA': '2.7.2', 'forgi': '2.2.3'}


def package_versions():
    """Return the installed version of each pinned library, or None if missing."""
    versions: dict[str, str | None] = {'python': sys.version.split()[0]}
    for name in PINNED:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def provenance(species_type):
    """Return library versions and the model checksum for the run report."""
    from mipyrna import model_store as ms

    return {'versions': package_versions(), 'model': ms.model_info(species_type)}


def run_check():
    """Print the state of tools, libraries and bundled models; return an exit code."""
    failed = False

    print('External tools')
    for tool in REQUIRED_TOOLS + OPTIONAL_TOOLS:
        path = shutil.which(tool)
        required = tool in REQUIRED_TOOLS
        status = path or ('MISSING' if required else 'missing (optional)')
        print(f'  {tool:15s} {status}')
        failed |= path is None and required

    print('Python libraries')
    for name, version in package_versions().items():
        expected = PINNED.get(name)
        note = '' if expected is None or version == expected else f'  (validated with {expected})'
        print(f'  {name:15s} {version}{note}')
        failed |= version is None

    print('Bundled models')
    try:
        from mipyrna import model_store as ms

        for species_type in ('plants', 'animals'):
            ms.load_model_for(species_type)
            print(f'  {species_type:15s} ok')
    except Exception as exc:
        print(f'  FAILED: {exc}')
        failed = True

    return 1 if failed else 0
