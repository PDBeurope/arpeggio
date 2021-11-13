from setuptools import setup, find_namespace_packages
import arpeggio

setup(
    name="arpeggio",
    version=arpeggio.__version__,
    description="Arpeggio calculates interatomic contacts based on the rules defined in CREDO.",
    project_urls={
        "Source code": "https://github.com/PDBeurope/arpeggio",
        "Documentation": "https://github.com/PDBeurope/arpeggio",
        "Author's repository": "https://github.com/harryjubb/arpeggio",
        "Paper": "https://doi.org/10.1016/j.jmb.2016.12.004",
        "Web server": "https://biosig.unimelb.edu.au/arpeggioweb/",
    },
    author="Harry Jubb",
    author_email="harry.jubb@sanger.ac.uk",
    license="GNU General Public License v3.0",
    keywords="arpeggio PDB contacts bioinformatics mmCIF protein ligand interactions CREDO",
    packages=find_namespace_packages(),
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=["pdbecif>=1.5", "numpy", "biopython"],
    python_requires='>=3.7',
    extras_require={
        "tests": ["pytest", "pytest-cov"]
    },
    entry_points={
        "console_scripts": ["arpeggio=arpeggio.scripts.process_protein_cli:main"]
    },
)
