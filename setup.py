import setuptools 
from pathlib import Path
import tarfile
import os

def install_package_data():
    package_data_dir = 'mipyrna/data'
    tarfile_path = os.path.join(package_data_dir, 'ncRNA_rfam.tar.gz')
    extract_dir = package_data_dir
    
    with tarfile.open(tarfile_path, 'r:gz') as tar:
        tar.extractall(extract_dir)

        
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
setuptools.setup(
    name='miPyRNA',
    version='0.1',
    description ="miPyRNA a python based Small RNA Sequencing analysis package",
    url="https://bioinfo.usu.edu/pyseqrna",
    author="Naveen Duhan",
    author_email = "naveen.duhan@usu.edu",
    license='MIT',
    packages = setuptools.find_packages(),
    include_package_data=True,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
            'console_scripts': [
                    'mipyrna = mipyrna.__main__:main',
            ]
    },
#     install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
install_requires = ['certifi','cffi','charset-normalizer','cycler','dill','fonttools','future',
'idna','importlib-metadata','Jinja2','kiwisolver','Markdown','MarkupSafe','matplotlib','numpy',
'packaging','pandas','patsy','Pillow','psutil','pycparser','pyfastx','pyparsing','pysam',
'python-dateutil','pytz','pytz-deprecation-shim','requests','scipy','seaborn','six',
'statsmodels','tzdata','tzlocal','urllib3','waiting','zipp','wheel','openpyxl'
],

 classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',


        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    package_data={'': ['param/*.ini', 'data/*']},
    # python_requires='==3.8'
)


install_package_data()