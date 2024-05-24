"""
mmCIF reading module to crate internal representation of proteins in
biopython and openbabel.
"""

import os

import numpy
from Bio import PDB
from Bio.PDB.StructureBuilder import StructureBuilder
import gemmi
from openbabel import openbabel as ob
from arpeggio.core import config 

# region common


def _get_res_id(atom_sites, i):
    return int(
        atom_sites["pdbe_label_seq_id"][i]
        if "pdbe_label_seq_id" in atom_sites
        else atom_sites["auth_seq_id"][i]
    )


def _get_b_factor(atom_sites, i):
    return float(
        atom_sites["B_iso_or_equiv"][i] if "B_iso_or_equiv" in atom_sites else 20.0
    )


def _get_ins_code(atom_sites, i):
    return (
        " "
        if not atom_sites["pdbx_PDB_ins_code"][i]
        else atom_sites["pdbx_PDB_ins_code"][i]
    )


def _format_formal_charge(atom_sites, i):
    if "pdbx_formal_charge" not in atom_sites:
        return 0

    charge = atom_sites["pdbx_formal_charge"][i]
    return 0 if not charge else int(charge)


# endregion common

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
        raise IOError("File {} not found".format(path))
    
    st = gemmi.read_structure(path, merge_chain_parts=True)
    
    # Trim all the other models but the first one from the mmcif structure
    # They are not used in the calculation and causes problems in Openbabel.
    del st[1:]
    groups = gemmi.MmcifOutputGroups(True)
    groups.group_pdb = True
    cif_block = st.make_mmcif_block(groups)
    mol = _parse_atom_site_openbabel(cif_block)

    return mol


def _parse_atom_site_openbabel(cif_block):
    """Parse _atom_site record to OBMolecule

    Args:
        cif_block (gemmi.cif.Block): Parsed mmcif file.

    Returns:
        [OBMol]: openbabel representation of the protein structure.
    """
    atom_site = cif_block.get_mmcif_category("_atom_site.")

    last_res_id = None
    last_res_name = None
    last_chain_id = None
    chain_num = 0
    res = None

    mol = ob.OBMol()
    mol.SetChainsPerceived()
    mol.BeginModify()

    for i in range(len(atom_site["id"])):
        current_res_id = _get_res_id(atom_site, i)
        ins_code = _get_ins_code(atom_site, i)

        if last_chain_id != atom_site["auth_asym_id"][i]:
            chain_num += 1
            last_chain_id = atom_site["auth_asym_id"][i]

        if (
            current_res_id != last_res_id
            or atom_site["auth_asym_id"][i] != last_chain_id
            or atom_site["label_comp_id"][i] != last_res_name
        ):

            last_res_id = current_res_id
            last_res_name = atom_site["label_comp_id"][i]

            res = mol.NewResidue()
            res.SetChainNum(chain_num)
            res.SetNum(str(last_res_id))
            res.SetName(atom_site["label_comp_id"][i])
            res.SetInsertionCode(ins_code)

        _init_openbabel_atom(mol, res, atom_site, i)

    resdat = ob.OBResidueData()
    resdat.AssignBonds(mol)
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()


    if "_struct_conn." in cif_block.get_mmcif_category_names():
        parse_struct_conn_bonds(mol, cif_block)

    mol.EndModify()


    return mol


def parse_struct_conn_bonds(mol, cif_block):
    """Add bonds from _struct_conn to the ob molecule.

    Args:
        mol (ob.OBMol): OpenBabel molecule object
        cif_block (gemmi.cif.Block): Parsed mmcif file Block object.
    """
    
    struct_conn = cif_block.get_mmcif_category("_struct_conn.")
    pivot = list(struct_conn.keys())[0]
    
    for i in range(len(struct_conn[pivot])):
        __process_struct_conn(mol, cif_block, i)


def __process_struct_conn(mol, cif_block, index):
    """Process a single entry in the mmcif_dict

    Args:
        mol (ob.OBMol): Openbabel molecule object
        cif_block (gemmi.cif.Block): Parsed mmcif Block object.
        index (int): Position in the struct_conn category

    Returns:
        int: Id of the atom in the openbabel molecule
    """

    def find_atom_id(atom_site, residue):
        for i in range(0, len(atom_site["id"])):
            if (
                atom_site["auth_seq_id"][i] == residue[1]
                and atom_site["auth_asym_id"][i] == residue[0]
                and atom_site["label_atom_id"][i] == residue[2]
            ):

                return int(atom_site["id"][i])

        return 0

    struct_conn = cif_block.get_mmcif_category("_struct_conn.")
    atom_sites = cif_block.get_mmcif_category("_atom_site.")

    res_a = (
        struct_conn["ptnr1_auth_asym_id"][index],
        struct_conn["ptnr1_auth_seq_id"][index],
        struct_conn["ptnr1_label_atom_id"][index],
    )

    res_b = (
        struct_conn["ptnr2_auth_asym_id"][index],
        struct_conn["ptnr2_auth_seq_id"][index],
        struct_conn["ptnr2_label_atom_id"][index],
    )

    atom_id_a = find_atom_id(atom_sites, res_a)
    atom_id_b = find_atom_id(atom_sites, res_b)

    if atom_id_a != 0 and atom_id_b != 0:
        __add_bond_to_openbabel(mol, atom_id_a, atom_id_b)


def __add_bond_to_openbabel(mol, atom_id_a, atom_id_b):
    """Create bond in the openbabel molecule if it does not exist.

    Args:
        mol (OBMol): Openbabel mol representation
        atom_id_a (int): First atom id.
        atom_id_b (int): Neighbout atom id.
    """

    pivot = mol.GetAtomById(atom_id_a)
    for neighbor in ob.OBAtomAtomIter(pivot):
        if neighbor.GetId() == atom_id_b:
            return

    b = mol.NewBond()
    b.SetBegin(pivot)
    b.SetEnd(mol.GetAtomById(atom_id_b))
    b.SetBondOrder(1)


def _init_openbabel_atom(mol, res, atom_sites, i):
    """Initialize openbabel atom

    Args:
        mol (ob.OBMol): Molecule the atom will be added to.
        res (OBResidue): Residue the atom will be added to.
        atom_sites (dict of str): Parsed mmcif structure of the input file.
        i (int): Pointer to the atom under question.

    Returns:
        OBAtom: openbabel Atom representation.
    """
    atom = mol.NewAtom(int(atom_sites["id"][i]))

    atom.SetVector(
        float(atom_sites["Cartn_x"][i]),
        float(atom_sites["Cartn_y"][i]),
        float(atom_sites["Cartn_z"][i]),
    )

    type_symbol = atom_sites["type_symbol"][i]
    element = f'{type_symbol[0]}{type_symbol[1].lower()}' if len(type_symbol) == 2 else type_symbol

    atomic_num = ob.GetAtomicNum(element)
    atom.SetAtomicNum(atomic_num)
    atom.SetFormalCharge(_format_formal_charge(atom_sites, i))
    res.AddAtom(atom)
    res.SetHetAtom(atom, atom_sites["group_PDB"][i] == "HETATM")
    res.SetSerialNum(atom, int(atom_sites["id"][i]))

    return atom


# end region


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
        raise IOError("File {} not found".format(path))

    st = gemmi.read_structure(path, merge_chain_parts=True)
    del st[1:]
    groups = gemmi.MmcifOutputGroups(True)
    groups.group_pdb = True
    cif_block = st.make_mmcif_block(groups)
    structure_builder = StructureBuilder()

    file_name = os.path.basename(path).split(".")[0]
    
    if '_entry.' in cif_block.get_mmcif_category_names():
        structure_id = cif_block.find_value('_entry.id').lower()
    else:
        structure_id = file_name.split('_')[0].lower()
    
    structure_builder.init_structure(structure_id)

    try:
        _atom_site = cif_block.get_mmcif_category("_atom_site.")
        _parse_atom_site_biopython(_atom_site, structure_builder)
    except KeyError:
        raise ValueError("The cif file does not contain _atom_site record")

    return structure_builder.get_structure()


def _get_hetero_flag(field, resn):
    """Converts PDB group_PDB field into the BioPython flag

    Args:
        field (str): value of the group_PDB field
        resn (str): Residue name as defined by the label_comp_id field.

    Returns:
        str: PDB heteroatom flag as defined by BioPython
    """
    if field == "HETATM":
        if resn in ("HOH", "WAT"):
            return "W"
        return "H"
    return " "


def _init_biopython_atom(builder, atom_sites, i):
    """Initializes a single atom in the PDB structure

    Args:
        builder (Bio.PDB.StructureBuilder.StructureBuilder): Bipython
            structure building object.
        atom_sites (dict of str): _atom_site dictionary of the mmcif file.
        i (int): Position within the _atom_site record.
    """
    x = float(atom_sites["Cartn_x"][i])
    y = float(atom_sites["Cartn_y"][i])
    z = float(atom_sites["Cartn_z"][i])
    coord = numpy.array((x, y, z), "f")

    try:
        builder.init_atom(
            atom_sites["label_atom_id"][i],
            coord,
            _get_b_factor(atom_sites, i),
            float(atom_sites["occupancy"][i]),
            " "
            if not atom_sites["label_alt_id"][i]
            else atom_sites["label_alt_id"][i],
            atom_sites["label_atom_id"][i],
            int(atom_sites["id"][i]),
            atom_sites["type_symbol"][i].upper(),
        )
    except PDB.PDBExceptions.PDBConstructionException:
        builder.init_residue(
            atom_sites["label_comp_id"][i],
            _get_hetero_flag(atom_sites, i),
            _get_res_id(atom_sites, i),
            _get_ins_code(atom_sites, i),
        )
        builder.init_atom(
            atom_sites["label_atom_id"][i],
            coord,
            _get_b_factor(atom_sites, i),
            float(atom_sites["occupancy"][i]),
            " "
            if not atom_sites["label_alt_id"][i]
            else atom_sites["label_alt_id"][i],
            atom_sites["label_atom_id"][i],
            int(atom_sites["id"][i]),
            atom_sites["type_symbol"][i].upper(),
        )


def _parse_atom_site_biopython(atom_sites, builder):
    """Parse mmcif atom site list to the BioPython atoms representation.

    Args:
        atom_sites (dict of str): _atom_site dictionary of the mmcif file.
        builder (Bio.PDB.StructureBuilder.StructureBuilder): Biopython
            structure building object.
    """
    last_model = None
    last_ins_code = None
    last_chain_id = None
    last_res_name = None
    last_res_id = None

    for i in range(len(atom_sites["id"])):
        res_id = _get_res_id(atom_sites, i)
        hetero_flag = _get_hetero_flag(
            atom_sites["group_PDB"][i], atom_sites["label_comp_id"][i]
        )
        ins_code = _get_ins_code(atom_sites, i)

        if last_model != atom_sites["pdbx_PDB_model_num"][i]:  # init model
            last_model = atom_sites["pdbx_PDB_model_num"][i]
            builder.init_model(int(last_model) - 1)

        if i == 0:
            builder.init_seg(
                "    "
            )  # some biopython magic. https://github.com/biopython/biopython/blob/master/Bio/PDB/PDBParser.py#L211

        if last_chain_id != atom_sites["auth_asym_id"][i]:  # init chain
            last_chain_id = atom_sites["auth_asym_id"][i]
            builder.init_chain(last_chain_id)

        if (
            last_res_id != res_id
            or last_chain_id != atom_sites["auth_asym_id"][i]
            or last_res_name != atom_sites["label_comp_id"][i]
            or last_ins_code != ins_code
        ):
            last_res_id = res_id
            last_res_name = atom_sites["label_comp_id"][i]
            last_ins_code = ins_code

            builder.init_residue(
                atom_sites["label_comp_id"][i], hetero_flag, res_id, ins_code
            )

        _init_biopython_atom(builder, atom_sites, i)



def get_component_types(path):

    if not os.path.isfile(path):
        raise IOError("File {} not found".format(path))
    
    cif_block = gemmi.cif.read(path).sole_block()
    if '_chem_comp.' in cif_block.get_mmcif_category_names():
        chem_comp = cif_block.get_mmcif_category('_chem_comp.')
        component_types = {}
        for i in range(len(chem_comp['id'])):
            chem_comp_type = chem_comp['type'][i]
            if chem_comp_type:
                comp_type  = config.ComponentType.from_chem_comp_type(chem_comp_type.upper())
                if not comp_type == config.ComponentType(5).name:
                    component_types[chem_comp['id'][i]] = comp_type
                else:
                    if chem_comp['id'][i].upper() == 'HOH' or chem_comp['id'][i].upper() == 'DOD':
                        component_types[chem_comp['id'][i]] = config.ComponentType(8).name
                    else:
                        component_types[chem_comp['id'][i]] = config.ComponentType(7).name
            else:
                component_types[chem_comp['id'][i]] = config.ComponentType(9).name

        return component_types
    
    else:
        raise ValueError('Missing _chem_comp. category in mmcif')

    

# endregion
