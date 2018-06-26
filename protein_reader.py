import os
import numpy

from mmCif.mmcifIO import MMCIF2Dict
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB import PDBIO, PDBParser
"""
mmCIF reading module.
"""


def read_protein_from_file(path):
    """Read in mmcif protein structure.

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
        _atom_site = list(parsed.values())[0]['_atom_site']
        _parse_atom_site(_atom_site, structure_builder)
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


def _init_atom(builder, atom_sites, i):
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

    builder.init_atom(atom_sites['label_atom_id'][i],
                      coord,
                      float(atom_sites['B_iso_or_equiv'][i]),
                      float(atom_sites['occupancy'][i]),
                      ' ' if atom_sites['label_alt_id'][i] == '.' else atom_sites['label_alt_id'][i],
                      atom_sites['label_atom_id'][i],
                      int(atom_sites['id'][i]),
                      atom_sites['type_symbol'][i])


def _parse_atom_site(atom_sites, builder):
    """Parse mmcif atom site list to the BioPython atoms reprensentation.

    Args:
        atom_sites (dict of str): _atom_site dictionary of the mmcif file.
        builder (Bio.PDB.StructureBuilder.StructureBuilder): Bipython
            structure building object.
    """
    current_model = None
    current_chain_id = None
    current_res_id = None

    for i in range(len(atom_sites['id'])):
        res_id = int(atom_sites['pdbe_label_seq_id'][i]
                     if 'pdbe_label_seq_id' in atom_sites
                     else atom_sites['label_seq_id'][i])
        hetero_flag = _get_hetero_flag(atom_sites['group_PDB'][i], atom_sites['label_comp_id'][i])
        ins_code = ' ' if atom_sites['pdbx_PDB_ins_code'][i] == '?' else atom_sites['pdbx_PDB_ins_code'][i]

        if current_model != atom_sites['pdbx_PDB_model_num'][i]:  # init model
            current_model = atom_sites['pdbx_PDB_model_num'][i]
            builder.init_model(int(current_model) - 1)

        if i == 0:
            builder.init_seg('    ')

        if current_chain_id != atom_sites['label_asym_id'][i]:  # init chain
            current_chain_id = atom_sites['label_asym_id'][i]
            builder.init_chain(current_chain_id)

        if current_res_id != res_id:
            current_res_id = res_id
            builder.init_residue(atom_sites['label_comp_id'][i],
                                 hetero_flag,
                                 res_id,
                                 ins_code)

        _init_atom(builder, atom_sites, i)
