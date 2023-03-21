
# PDBe Arpeggio

This repository contains a refactored and expanded version of the [arpeggio](https://github.com/harryjubb/arpeggio) tool, modified and maintained by the [PDBe team](https://pdbe.org) to use as part of their weekly PDB release process.

`pdbe-arpeggio` is a python library, that can be used to calculate interatomic contacts in a protein based on the rules defined in [CREDO](https://pubmed.ncbi.nlm.nih.gov/19207418/). This library only supports protein structures in mmCIF format.

If you make use of `pdbe-arpeggio`, please cite the following article:

Harry C Jubb, Alicia P Higueruelo, Bernardo Ochoa-Montaño, Will R Pitt, David B Ascher, Tom L Blundell,
[Arpeggio: A Web Server for Calculating and Visualising Interatomic Interactions in Protein Structures](https://doi.org/10.1016/j.jmb.2016.12.004).
Journal of Molecular Biology,
Volume 429, Issue 3,
2017,
Pages 365-371,
ISSN 0022-2836,

# Disclaimer

This is a refactored version of the original [arpeggio](https://github.com/harryjubb/arpeggio) by Harry Jubb. Main changes are:

* Python 3 support
* Modular architecture to make `arpeggio` PIP installable.
* Support for mmCIF format to process also larger structures.
* Results in JSON format.

`pdbe-arpeggio` doesn't do any processing of input protein structure, other than what BioPython does by default. Alternate locations and missing density are not explicitly accounted for and may result in anomalous results. Please use with caution.

# Getting Started

## Web Interface

If you would like to run original version of arpeggio on a small number of individual structures, the easiest way to get started is to use the [web interface](http://biosig.unimelb.edu.au/arpeggioweb/).


## Installation

The easiest way to set up arpeggio is using [Conda](https://docs.conda.io/en/latest/). Create a conda environment and install arpeggio dependencies:

```bash
conda create -n arpeggio-env python=3.9 -c conda-forge gemmi openbabel biopython
conda activate arpeggio-env
pip install pdbe-arpeggio
```

## Dependencies

* [Open Babel](https://openbabel.org) >= 3.0.0
* [Biopython](https://biopython.org) >= 1.80
* [gemmi](https://gemmi.readthedocs.io/en/latest/index.html) >= 0.5.8


## Running arpeggio

`pdbe-arpeggio 1tqn_h.cif [options]`

e.g. 

`pdbe-arpeggio -s /A/508/ -o out/ 1tqn_h.cif`

Use `pdbe-arpeggio -h` for available options.



## Output

All the identified interactions are written in a JSON file.

### Interactions

#### atom-atom interactions

| Key | Interaction | Description |
| --- | ----------- | ----------- |
| clash | Clash | Denotes if the atom is involved in a steric clash. |
| covalent | Covalent | Denotes if the atom appears to be covalently bonded. |
| vdw_clash| VdW Clash | Denotes if the van der Waals radius of the atom is clashing with one or more other atoms. |
| vdw | VdW | Denotes if the van der Waals radius of the the atom is interacting with one or more other atoms. |
| proximal | Proximal | Denotes if the atom is > the VdW interaction distance, but within 5 Angstroms of other atom(s). |
| hbond | Hydrogen Bond | Denotes if the atom forms a hydrogen bond. |
| weak_hbond | Weak Hydrogen Bond  | Denotes if the atom forms a weak hydrogen bond. |
| xbond | Halogen Bond | Denotes if the atom forms a halogen bond. |
| ionic | Ionic | Denotes if the atom may interact via charges. |
| metal | Metal Complex | Denotes if the atom is part of a metal complex. |
| aromatic | Aromatic | Denotes an aromatic ring atom interacting with another aromatic ring atom. |
| hydrophobic | Hydrophobic | Denotes hydrophobic interaction. |
| carbonyl | Carbonyl | Denotes a carbonyl-carbon:carbonyl-carbon interaction. |
| polar | Polar | Less strict hydrogen bonding (without angle terms). |
| weak_polar | Weak Polar | Less strict weak hydrogen bonding (without angle terms).

#### atom-plane interactions

| Key | Interaction | Description |
| --- | ----------- | ----------- |
| CARBONPI | Carbon-PI  | Weakly electropositive carbon atom - $\Pi$ interactions [[ref]](https://doi.org/10.1016/j.bmc.2007.09.023) |
| CATIONPI | Cation-PI | Cation - $\Pi$ interactions [[ref]](https://doi.org/10.1002/prot.20417) |
| DONORPI | Donor-PI | Hydrogen Bond donor - $\Pi$ interactions [[ref]](https://doi.org/10.1016/j.bmc.2007.09.023) |
| HALOGENPI | Halogen-PI | Halogen bond donors - $\Pi$ [[ref]](https://doi.org/10.1073/pnas.0407607101) |
| METSULPHURPI | Sulphur-PI | Methionine sulphur - $\Pi$ ring interactions [[ref]](https://dx.doi.org/10.1074/jbc.M112.374504) |

#### plane-plane interactions

Follows nomenclature established by [Chakrabarti and Bhattacharyya (2007)](https://doi.org/10.1016/j.pbiomolbio.2007.03.016)

#### group-group/plane interactions

| Key | Interaction | Description |
| --- | ----------- | ----------- |
| AMIDEAMIDE | amide - amide  |  [[ref]](https://doi.org/10.1002/jcc.21212) |
| AMIDERING | amide - ring | [[ref]](https://doi.org/10.1016/0014-5793(86)80730-X) |

### Interacting entities

| Key  | Meaning |
| ---- |---------|
| INTER | Between an atom from the user's selection and a non-selected atom |
| INTRA_SELECTION | Between two atoms both in the user's selection |
| INTRA_NON_SELECTION | Between two atoms that are both not in the user's selection |
| SELECTION_WATER | Between an atom in the user's selection and a water molecule |
| NON_SELECTION_WATER | Between an atom that is not in the user's selection and a |water molecule
| WATER_WATER | Between two water molecules |

### Examples

#### atom-atom interaction

```json
   {
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "CB",
            "auth_seq_id": 313,
            "label_comp_id": "VAL",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "contact": [
            "proximal",
            "hydrophobic"
        ],
        "distance": 4.02,
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "CBB",
            "auth_seq_id": 508,
            "label_comp_id": "HEM",
            "label_comp_type": "B",
            "pdbx_PDB_ins_code": " "
        },
        "interacting_entities": "INTER",
        "type": "atom-atom"
    },
```

#### atom-plane interaction

```json
    {
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "O",
            "auth_seq_id": 523,
            "label_comp_id": "HOH",
            "label_comp_type": "W",
            "pdbx_PDB_ins_code": " "
        },
        "contact": [
            "DONORPI"
        ],
        "distance": 3.9,
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "C1A,C2A,C3A,C4A,NA",
            "auth_seq_id": 508,
            "label_comp_id": "HEM",
            "label_comp_type": "B",
            "pdbx_PDB_ins_code": " "
        },
        "interacting_entities": "INTER",
        "type": "atom-plane"
    },
```

#### plane-plane interaction

```json
{
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "C1B,C2B,C3B,C4B,NB",
            "auth_seq_id": 508,
            "label_comp_id": "HEM",
            "label_comp_type": "B",
            "pdbx_PDB_ins_code": " "
        },
        "contact": [
            "FT",
            "ET"
        ],
        "distance": 4.72,
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "CD1,CD2,CE1,CE2,CG,CZ",
            "auth_seq_id": 435,
            "label_comp_id": "PHE",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "interacting_entities": "INTER",
        "type": "plane-plane"
    },
```

#### group-group interaction

```json
    {
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "C,CA,N,O",
            "auth_seq_id": 308,
            "label_comp_id": "GLU",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "contact": [
            "AMIDEAMIDE"
        ],
        "distance": 4.29,
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "C,CA,N,O",
            "auth_seq_id": 310,
            "label_comp_id": "THR",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "interacting_entities": "INTRA_BINDING_SITE",
        "type": "group-group"
    },
```
