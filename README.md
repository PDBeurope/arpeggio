# Arpeggio

Outline
--------

Arpeggio calculates interatomic contacts based on the rules defined in [CREDO](http://marid.bioc.cam.ac.uk/credo). The program will be freely available and require only Open Source dependencies.

Dependencies
------------

Arpeggio is written in Python and currently has the following dependencies:

### Dependencies

- Python (v2.7)
- Numpy
- BioPython (>= v1.60)
- OpenBabel (with Python bindings)

### Recommended
- PyMOL (for visualising contacts)

Arpeggio may work with earlier versions of BioPython, however these haven't been tested. It is recommended that each dependency be the latest version.

Running
-------

`python arpeggio.py pdb [options]`

Use `python arpeggio.py -h` for available options.

Arpeggio doesn't do any checking of your PDB structure, other than what BioPython does by default. Alternate locations and missing density are not explicitly accounted for and may result in anomalous results. Please use with caution.

Output Files
------------

### `*.ari`

Atom-aromatic ring interactions.

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Atom   | `<chain_id>/<res_num><ins_code (stripped)><atom_name>` | Uniquely identifiable atom |
| Ring ID | integer | Internal number used to identify the aromatic ring |
| Ring Centroid | list | 3D coordinates of the centre of the ring |