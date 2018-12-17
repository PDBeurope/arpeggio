import os

import numpy
from Bio import PDB
from Bio.PDB.StructureBuilder import StructureBuilder
from mmCif.mmcifIO import MMCIF2Dict
from openbabel import OBBitVec, OBElementTable, OBMol, OBResidueData

"""
mmCIF reading module to crate internal representation of proteins in
biopython and openbabel.
"""

# region common


def _get_res_id(atom_sites, i):
    return int(atom_sites['pdbe_label_seq_id'][i]
               if 'pdbe_label_seq_id' in atom_sites
               else atom_sites['auth_seq_id'][i])


def _get_ins_code(atom_sites, i):
    return ' ' if atom_sites['pdbx_PDB_ins_code'][i] == '?' else atom_sites['pdbx_PDB_ins_code'][i]


def _format_formal_charge(atom_sites, i):
    if 'pdbx_formal_charge' not in atom_sites:
        return 0

    charge = atom_sites['pdbx_formal_charge'][i]
    return 0 if charge in ('?', '.') else int(charge)


def _trim_models(atom_site):
    """Trim all the other models but the first one from the mmcif structure
    They are not used in the calculation and causes problems in Openbabel.

    Args:
        atom_site (dict): Parsed _atom_site category

    Returns:
        dict: first model from the _atom_site
    """
    output = {}
    length = sum(x == '1' for x in atom_site['pdbx_PDB_model_num'])

    for k, v in atom_site.items():
        output[k] = v[:length]

    return output
# endregion common


# region openbabel
def read_mmcif_to_openbabel(path):
    """Read in mmcif structure to to the openbabel.

    Args:
        path (str): path to the mmcif file.

    Raises:
        IOError: In case file does not exist

    Returns:
        openbabel.OBMol: openbabel structure
    """
    if not os.path.isfile(path):
        raise IOError('File {} not found'.format(path))

    parsed = MMCIF2Dict().parse(path)
    perceived_atom_site = list(parsed.values())[0]['_atom_site']
    _atom_site = _trim_models(perceived_atom_site)
    mol = _parse_atom_site_openbabel(_atom_site)

    return mol


def _parse_atom_site_openbabel(atom_sites):
    """Parse _atom_site record to OBMolecule

    Args:
        atom_sites (dict of str): Parsed mmcif structure of the input file.

    Returns:
        [OBMol]: openbabel representation of the protein structure.
    """
    table = OBElementTable()
    last_res_id = None
    last_res_name = None
    last_chain_id = None
    chain_num = 0
    res = None

    mol = OBMol()
    mol.SetChainsPerceived()
    mol.BeginModify()

    for i in range(len(atom_sites['id'])):
        current_res_id = _get_res_id(atom_sites, i)
        ins_code = _get_ins_code(atom_sites, i)

        if last_chain_id != atom_sites['auth_asym_id'][i]:
            chain_num += 1
            last_chain_id = atom_sites['auth_asym_id'][i]

        if current_res_id != last_res_id or atom_sites['auth_asym_id'][i] != last_chain_id or atom_sites['label_comp_id'][i] != last_res_name:
            last_res_id = current_res_id
            last_res_name = atom_sites['label_comp_id'][i]

            res = mol.NewResidue()
            res.SetChainNum(chain_num)
            res.SetNum(str(last_res_id))
            res.SetName(atom_sites['label_comp_id'][i])
            res.SetInsertionCode(ins_code)

        _init_openbabel_atom(table, mol, res, atom_sites, i)

    resdat = OBResidueData()
    resdat.AssignBonds(mol, OBBitVec())
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()
    mol.EndModify()

    return mol


def _init_openbabel_atom(table, mol, res, atom_sites, i):
    """Initialize openbabel atom

    Args:
        table (OBElementTable): Element table to translate element type
            to element numbers.
        mol (OBMol): Molecule the atom will be added to.
        res (OBResidue): Residue the atom will be added to.
        atom_sites (dict of str): Parsed mmcif structure of the input file.
        i (int): Pointer to the atom under question.

    Returns:
        [OBAtom]: openbabel Atom representation.
    """
    atom = mol.NewAtom(int(atom_sites['id'][i]))

    atom.SetVector(
        float(atom_sites['Cartn_x'][i]),
        float(atom_sites['Cartn_y'][i]),
        float(atom_sites['Cartn_z'][i])
    )
    atomic_num = table.GetAtomicNum(atom_sites['type_symbol'][i])
    atom.SetAtomicNum(atomic_num)
    atom.SetFormalCharge(_format_formal_charge(atom_sites, i))

    res.AddAtom(atom)
    res.SetHetAtom(atom, atom_sites['group_PDB'][i] == 'HETATM')
    res.SetSerialNum(atom, int(atom_sites['id'][i]))

    #_set_ob_occupancy(float(atom_sites['occupancy'][i]), atom)
    return atom


# def _set_ob_occupancy(occupancy, atom):
#     fp = OBPairFloatingPoint()
#     fp.SetAttribute('_atom_site_occupancy')
#     fp.SetValue(occupancy)
#     fp.SetOrigin('mmcif')
#     atom.SetData(fp)

# endregion


# region Biopython


def read_mmcif_to_biopython(path):
    """Read in mmcif protein structure and report its Biopython structure

    Args:
        path (str): Path to the mmcif file.

    Raises:
        ValueError: In case _atom_site table is not present in the file.

    Returns:
        Bio.PDB.Structure.Structure: BioPython PDB structure
    """
    if not os.path.isfile(path):
        raise IOError('File {} not found'.format(path))

    structure_builder = StructureBuilder()
    parsed = MMCIF2Dict().parse(path)

    file_name = os.path.basename(path).split('.')[0]
    structure_id = next((x for x in parsed), file_name).lower()
    structure_builder.init_structure(structure_id)

    try:
        perceived_atom_site = list(parsed.values())[0]['_atom_site']
        _atom_site = _trim_models(perceived_atom_site)
        _parse_atom_site_biopython(_atom_site, structure_builder)
    except KeyError:
        raise ValueError('The cif file does not contain _atom_site record')

    return structure_builder.get_structure()


def _get_hetero_flag(field, resn):
    """Converts PDB group_PDB field into the BioPython flag

    Args:
        field (str): value of the group_PDB field
        resn (str): Residue name as defined by the label_comp_id field.

    Returns:
        str: PDB heteroatom flag as defined by BioPython
    """
    if field == 'HETATM':
        if resn == 'HOH' or resn == 'WAT':
            return 'W'
        else:
            return 'H'
    else:
        return ' '


def _init_biopython_atom(builder, atom_sites, i):
    """Initializes a single atom in the PDB structure

    Args:
        builder (Bio.PDB.StructureBuilder.StructureBuilder): Bipython
            structure building object.
        atom_sites (dict of str): _atom_site dictionary of the mmcif file.
        i (int): Position within the _atom_site record.
    """
    x = float(atom_sites['Cartn_x'][i])
    y = float(atom_sites['Cartn_y'][i])
    z = float(atom_sites['Cartn_z'][i])
    coord = numpy.array((x, y, z), "f")

    try:
        builder.init_atom(atom_sites['label_atom_id'][i],
                          coord,
                          float(atom_sites['B_iso_or_equiv'][i]),
                          float(atom_sites['occupancy'][i]),
                          ' ' if atom_sites['label_alt_id'][i] == '.' else atom_sites['label_alt_id'][i],
                          atom_sites['label_atom_id'][i],
                          int(atom_sites['id'][i]),
                          atom_sites['type_symbol'][i].upper())
    except PDB.PDBExceptions.PDBConstructionException:
        builder.init_residue(atom_sites['label_comp_id'][i],
                             _get_hetero_flag(atom_sites, i),
                             _get_res_id(atom_sites, i),
                             _get_ins_code(atom_sites, i))
        builder.init_atom(atom_sites['label_atom_id'][i],
                          coord,
                          float(atom_sites['B_iso_or_equiv'][i]),
                          float(atom_sites['occupancy'][i]),
                          ' ' if atom_sites['label_alt_id'][i] == '.' else atom_sites['label_alt_id'][i],
                          atom_sites['label_atom_id'][i],
                          int(atom_sites['id'][i]),
                          atom_sites['type_symbol'][i].upper())


def _parse_atom_site_biopython(atom_sites, builder):
    """Parse mmcif atom site list to the BioPython atoms reprensentation.

    Args:
        atom_sites (dict of str): _atom_site dictionary of the mmcif file.
        builder (Bio.PDB.StructureBuilder.StructureBuilder): Bipython
            structure building object.
    """
    last_model = None
    last_chain_id = None
    last_res_name = None
    last_res_id = None

    for i in range(len(atom_sites['id'])):
        res_id = _get_res_id(atom_sites, i)
        hetero_flag = _get_hetero_flag(atom_sites['group_PDB'][i], atom_sites['label_comp_id'][i])
        ins_code = _get_ins_code(atom_sites, i)

        if last_model != atom_sites['pdbx_PDB_model_num'][i]:  # init model
            last_model = atom_sites['pdbx_PDB_model_num'][i]
            builder.init_model(int(last_model) - 1)

        if i == 0:
            builder.init_seg('    ')  # some biopython magic. https://github.com/biopython/biopython/blob/master/Bio/PDB/PDBParser.py#L211

        if last_chain_id != atom_sites['auth_asym_id'][i]:  # init chain
            last_chain_id = atom_sites['auth_asym_id'][i]
            builder.init_chain(last_chain_id)

        if last_res_id != res_id or last_chain_id != atom_sites['auth_asym_id'][i] or last_res_name != atom_sites['label_comp_id'][i]:
            last_res_id = res_id
            last_res_name = atom_sites['label_comp_id'][i]

            builder.init_residue(atom_sites['label_comp_id'][i],
                                 hetero_flag,
                                 res_id,
                                 ins_code)

        _init_biopython_atom(builder, atom_sites, i)
# endregion
