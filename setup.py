from setuptools import setup

setup(
    name='Arpeggio',
    version='1.1',
    description='Arpeggio calculates interatomic contacts based on the rules defined in CREDO.',
    project_urls={
        'Source code': 'https://github.com/harryjubb/arpeggio',
        'Documentation': 'https://github.com/harryjubb/arpeggio',
        'Paper': 'https://doi.org/10.1016/j.jmb.2016.12.004',
        'Web server': 'http://biosig.unimelb.edu.au/arpeggioweb/'
    },
    author='Harry Jubb',
    author_email='harry.jubb@sanger.ac.uk',
    license='GNU General Public License v3.0',
    keywords='arpeggio PDB mmCIF protein ligand interactions CREDO',
    packages=['arpeggio'],
    scripts=['bin/arpeggio',
             'bin/run_and_show.sh'],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=['pdbecif', 'numpy', 'biopython'],
    tests_require=['pytest'])
