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
| Atom   | string `<chain_id>/<res_num><ins_code (stripped)>/<atom_name>` | Uniquely identifies an atom |
| Ring ID | integer | Internal number used to identify the aromatic ring |
| Ring centroid | list | 3D coordinates of the centre of the ring |
| Interaction type | list | Type(s) of interaction this atom/ring are making |

### `*.ri`

Aromatic ring-aromatic ring interactions.

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Ring 1 ID | integer | Internal number used to identify the first aromatic ring |
| Ring 1 Residue   | string `<chain_id>/<res_num>` | Uniquely identifies the first ring's residue |
| Ring 1 centroid | list | 3D coordinates of the centre of the first ring |
| Ring 2 ID | integer | Internal number used to identify the second aromatic ring |
| Ring 2 Residue   | string `<chain_id>/<res_num>` | Uniquely identifies an the second ring's residue |
| Ring 2 centroid | list | 3D coordinates of the centre of the second ring |
| Inter or intra residue | string from (`INTER`, `INTRA_RESIDUE`) | States whether this ring-ring interaction is within the same residue (e.g. within a small molecule ligand), or between two different residues |
| Interacting entities | string from (`INTER`, `INTRA_NON_SELECTION`, `INTRA_SELECTION`) | Distinguishes how this interacting ring pair relates to the selected atoms: see below |

### `*.rings`

Aromatic rings found in the structure

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Ring ID | integer | Internal number used to identify the ring |
| Ring Residue   | string `<chain_id>/<res_num>` | Uniquely identifies the ring's residue |
| Ring centroid | list | 3D coordinates of the centre of the ring |

### `*.atomtypes`

Atom types for all of the atoms for which interactions are calculated for. This includes the selected atoms, and the atoms that those atoms interact with.

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Atom   | string `<chain_id>/<res_num><ins_code (stripped)>/<atom_name>` | Uniquely identifies an atom |
| Atom types | list | All the atom types that this atom possesses |

### `*.contacts`

Pairwise contacts between two individual atoms.

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Atom 1 | string `<chain_id>/<res_num><ins_code (stripped)>/<atom_name>` | Uniquely identifies the first atom in this contact |
| Atom 2 | string `<chain_id>/<res_num><ins_code (stripped)>/<atom_name>` | Uniquely identifies the second atom in this contact |
| Clash | boolean::integer | Denotes if the covalent radii of the two atoms are clashing, i.e. steric clash |
| Covalent | boolean::integer | Denotes if the two atoms appear to be covalently bonded |
| VdW Clash | boolean::integer | Denotes if the van der Waals radii of the two atoms are clashing |
| VdW | boolean::integer | Denotes if the van der Waals radii of the two atoms are interacting |
| Proximal | boolean::integer | Denotes the two atoms being > the VdW interaction distance, but with in 5 Angstroms of each other |
| Hydrogen Bond | boolean::integer | Denotes if the atoms form a hydrogen bond |
| Weak Hydrogen Bond | boolean::integer | Denotes if the atoms form a weak hydrogen bond |
| Halogen Bond | boolean::integer | Denotes if the atoms form a halogen bond |
| Ionic | boolean::integer | Denotes if the atoms may interact via charges |
| Metal Complex | boolean::integer | Denotes if the atoms are part of a metal complex |
| Aromatic | boolean::integer | Denotes two aromatic ring atoms interacting |
| Hydrophobic | boolean::integer | Denotes two hydrophobic atoms interacting |
| Carbonyl | boolean::integer | Denotes a carbonyl-carbon:carbonyl-carbon interaction |
| Polar | boolean::integer | Less strict hydrogen bonding (without angle terms) |
| Weak Polar | boolean::integer | Less strict weak hydrogen bonding (without angle terms) |
| Interacting entities | string from (`INTER`, `INTRA_NON_SELECTION`, `INTRA_SELECTION`, `SELECTION_WATER`, `NON_SELECTION_WATER`, `WATER_WATER`) | Distinguishes how this atom pair relates to the selected atoms: see below |

**Clash, Covalent, VdW Clash, VdW and Proximal interactions are mutually exclusive: Other interactions can occur simulataneously.**

Entity interactions:

- `INTER`: Between an atom from the user's selection and a non-selected atom
- `INTRA_SELECTION`: Between two atoms both in the user's selection
- `INTRA_NON_SELECTION`: Between two atoms that are both not in the user's selection
- `SELECTION_WATER`: Between an atom in the user's selection and a water molecule
- `NON_SELECTION_WATER`: Between an atom that is not in the user's selection and a water molecule
- `WATER_WATER`: Between two water molecules

### `*.sift`

Interaction fingerprints for individual atoms. These are binary (i.e., on/off) indications of an atom's interaction, not counts.

| Column | Datatype | Description |
| ------ | -------- | ----------- |
| Atom | string `<chain_id>/<res_num><ins_code (stripped)>/<atom_name>` | Uniquely identifies the atom |
| Clash | boolean::integer | Denotes if the atom is involved in a steric clash |
| Covalent | boolean::integer | Denotes if the atom appears to be covalently bonded |
| VdW Clash | boolean::integer | Denotes if the van der Waals radius of the atom is clashing with one or more other atoms |
| VdW | boolean::integer | Denotes if the van der Waals radius of the the atom is interacting with one or more other atoms |
| Proximal | boolean::integer | Denotes if the atom is > the VdW interaction distance, but with in 5 Angstroms of other atom(s) |
| Hydrogen Bond | boolean::integer | Denotes if the atom forms a hydrogen bond |
| Weak Hydrogen Bond | boolean::integer | Denotes if the atom forms a weak hydrogen bond |
| Halogen Bond | boolean::integer | Denotes if the atom forms a halogen bond |
| Ionic | boolean::integer | Denotes if the atom may interact via charges |
| Metal Complex | boolean::integer | Denotes if the atom is part of a metal complex |
| Aromatic | boolean::integer | Denotes an aromatic ring atom interacting with another aromatic ring atom |
| Hydrophobic | boolean::integer | Denotes hydrophobic interaction |
| Carbonyl | boolean::integer | Denotes a carbonyl-carbon:carbonyl-carbon interaction |
| Polar | boolean::integer | Less strict hydrogen bonding (without angle terms) |
| Weak Polar | boolean::integer | Less strict weak hydrogen bonding (without angle terms) |

### `*.specific.sift`

Interaction fingerprints for individual atoms. These are binary (i.e., on/off) indications of an atom's interaction, not counts.

The columns match the `*.sift` files, but the first 15 columns (after the atom identifier) denote only interactions between the selection made by the user, and non-selection atoms; the second 15 columns indicate interactions made within the selection made by the user; and the third 15 columns indicate interactions made with water only.