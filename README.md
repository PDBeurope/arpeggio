# OpenCREDO

Outline
--------

OpenCREDO calculates interatomic contacts based on the rules defined in CREDO. The program will be freely available and require only Open Source dependencies.

Dependencies
------------

OpenCREDO is written in Python and currently has the following dependencies:

### Dependencies

- Python (v2.7)
- Numpy
- BioPython (>= v1.60)
- OpenBabel (with Python bindings)

### Recommended
- PyMOL (for visualising contacts)

OpenCREDO may work with earlier versions of BioPython, however these haven't been tested. It is recommended that each dependency be the latest version.

Running
-------

`python opencredo.py pdb [options]`

Use `python opencredo.py -h` for available options.

Output
------

TBA.