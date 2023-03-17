"""Arpeggio helper methods
"""

import collections
import logging
import os
import platform

import numpy as np
from openbabel import openbabel as ob
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from arpeggio.core import SelectionError, SiftMatchError, config

# region general


def setup_filetype(filepath):
    """Sets up argument type to be consumedd by openbabel

    Args:
        filepath (str): Path to the source file

    Returns:
        str: file type
    """
    filetype = 'pdb'

    if filepath.endswith('.pdb') or filepath.endswith('.ent'):
        filetype = 'pdb'
    elif filepath.endswith('.mmcif') or filepath.endswith('.cif'):
        filetype = 'mmcif'

    return filetype


def rename_output_file(original_filename, new_extension):
    """

    Args:
        original_filename (str): [description]
        new_extension (str): [description]

    Returns:
        str: new output file name
    """
    original_extension = os.path.splitext(original_filename)[-1]

    if original_extension:
        return os.path.splitext(original_filename)[0] + new_extension

    return original_filename + new_extension


def max_mem_usage():
    """Returns maximum memory usage of the program thus far, in megabytes, as a string.

    Returns:
        str: Program message
    """

    try:
        import resource
        base = 1024.0 if platform.system() == 'Linux' else 1048576.0
        return str(round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / base, 2)) + ' MB'
    except Exception as err:
        return 'Resource usage information not available {}'.format(str(err))
# endregion

# region bonds identification


def is_hbond(donor, acceptor, vdw_comp_factor):
    """Determines if atoms share Hbond.

    Args:
        donor (OBAtom): First atom
        acceptor (OBAtom): Second atom
        vdw_comp_factor (float): Compensation factor for VdW radii
            dependent interaction types.

    Returns:
        int: 1/0 binary information whether or not the atoms share
            Hbond.
    """
    for hydrogen_coord in donor.h_coords:
        h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)

        if h_dist <= config.VDW_RADII['H'] + acceptor.vdw_radius + vdw_comp_factor:
            if get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= config.CONTACT_TYPES['hbond']['angle rad']:
                return 1

    return 0


def is_weak_hbond(donor, acceptor, vdw_comp_factor):
    """Determines if atoms share a weak Hbond.

    Args:
        donor (OBAtom): First atom
        acceptor (OBAtom): Second atom
        vdw_comp_factor (float): Compensation factor for VdW radii
            dependent interaction types.

    Returns:
        int: 1/0 binary information whether or not the atoms share
            weak Hbond.
    """
    for hydrogen_coord in donor.h_coords:
        h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)

        if h_dist <= config.VDW_RADII['H'] + acceptor.vdw_radius + vdw_comp_factor:
            if get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= config.CONTACT_TYPES['weak hbond']['angle rad']:
                return 1

    return 0


def is_halogen_weak_hbond(donor, halogen, ob_mol, vdw_comp_factor, bio_to_ob, ob_to_bio):
    """Determines if atoms share a weak halogen bond.

    org. desc: `nbr` WILL BE A BIOPYTHON ATOM ... HOPEFULLY
    Args:
        donor (OBAtom): Donor atom
        halogen (OBAtom): Halogen atom
        ob_mol (OBMol): [description]
        vdw_comp_factor (float): [description]
        bio_to_ob (dict of int): biopython to openbabel mapping
        ob_to_bio (dict of int): openbabel to biopython mapping

    Returns:
        int: 1/0 binary information whether or not the atoms share
            weak weak Hbond.
    """

    # entry 5t45 obabel does not identify bond between
    # SER 246 A (OG) and BEF 903 A (F1) which makes neighbour to be
    # None and then it crashes
    neighbour = get_single_bond_neighbour(ob_mol.GetAtomById(bio_to_ob[halogen]))
    if neighbour is None:
        return 0

    nbr = ob_to_bio[neighbour.GetId()]

    for hydrogen_coord in donor.h_coords:

        h_dist = np.linalg.norm(halogen.coord - hydrogen_coord)

        if h_dist <= config.VDW_RADII['H'] + halogen.vdw_radius + vdw_comp_factor:

            if config.CONTACT_TYPES['weak hbond']['cx angle min rad'] <= get_angle(nbr.coord, halogen.coord, hydrogen_coord) <= config.CONTACT_TYPES['weak hbond']['cx angle max rad']:

                return 1

    return 0


def is_xbond(donor, acceptor, ob_mol, bio_to_ob, ob_to_bio):
    """Determines if atoms share a weak halogen bond.

    Args:
        donor (OBAtom): Donor atom
        acceptor (OBAtom): Halogen atom
        ob_mol (OBMol): [description]

    Returns:
        int: 1/0 binary information whether or not the atoms share
            xbond.
    """

    # `nbr` WILL BE A BIOPYTHON ATOM
    # ... HOPEFULLY
    nbr = ob_to_bio[get_single_bond_neighbour(ob_mol.GetAtomById(bio_to_ob[donor])).GetId()]
    theta = get_angle(nbr.coord, donor.coord, acceptor.coord)

    if theta >= config.CONTACT_TYPES['xbond']['angle theta 1 rad']:
        return 1

    return 0


def update_atom_sift(atom, addition, contact_type='INTER'):
    """[summary]

    Args:
        atom ([type]): [description]
        addition ([type]): [description]
        contact_type (str, optional): Defaults to 'INTER'. [description]
    """

    atom.sift = [x or y for x, y in zip(atom.sift, addition)]

    if contact_type == 'INTER':
        atom.sift_inter_only = [x or y for x, y in zip(atom.sift_inter_only, addition)]

    if 'INTRA' in contact_type:
        atom.sift_intra_only = [x or y for x, y in zip(atom.sift_intra_only, addition)]

    if 'WATER' in contact_type:
        atom.sift_water_only = [x or y for x, y in zip(atom.sift_water_only, addition)]


def update_atom_fsift(atom, addition, contact_type='INTER'):
    """[summary]

    Args:
        atom ([type]): [description]
        addition ([type]): [description]
        contact_type (str, optional): Defaults to 'INTER'. [description]
    """

    atom.actual_fsift = [x or y for x, y in zip(atom.actual_fsift, addition)]

    if contact_type == 'INTER':
        atom.actual_fsift_inter_only = [x or y for x, y in zip(atom.actual_fsift_inter_only, addition)]

    if 'INTRA' in contact_type:
        atom.actual_fsift_intra_only = [x or y for x, y in zip(atom.actual_fsift_intra_only, addition)]

    if 'WATER' in contact_type:
        atom.actual_fsift_water_only = [x or y for x, y in zip(atom.actual_fsift_water_only, addition)]


def update_atom_integer_sift(atom, addition, contact_type='INTER'):
    """[summary]

    Args:
        atom ([type]): [description]
        addition ([type]): [description]
        contact_type (str, optional): Defaults to 'INTER'. [description]
    """

    atom.integer_sift = [x + y for x, y in zip(atom.sift, addition)]

    if contact_type == 'INTER':
        atom.integer_sift_inter_only = [x + y for x, y in zip(atom.sift_inter_only, addition)]

    if 'INTRA' in contact_type:
        atom.integer_sift_intra_only = [x + y for x, y in zip(atom.sift_intra_only, addition)]

    if 'WATER' in contact_type:
        atom.integer_sift_water_only = [x + y for x, y in zip(atom.sift_water_only, addition)]


def sift_xnor(sift1, sift2):
    """[summary]

    Args:
        sift1 ([type]): [description]
        sift2 ([type]): [description]

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """'''
    '''

    out = []

    for x, y in zip(sift1, sift2):

        if x and not y:  # TF
            out.append(0)

        elif not x and not y:  # FF
            out.append(1)

        elif x and y:  # TT
            out.append(1)

        elif not x and y:  # FT
            out.append(0)

        else:
            raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

    return out


def sift_match_base3(sift1, sift2):
    """[summary]
    0 = UNMATCHED
    1 = MATCHED
    2 = MATCH NOT POSSIBLE

    Assuming that sift1 is the potential SIFt, and sift2 is the actual SIFt.
    Args:
        sift1 ([type]): [description]
        sift2 ([type]): [description]

    Raises:
        SiftMatchError: [description]
        ValueError: [description]

    Returns:
        [type]: [description]
    """

    out = []

    for x, y in zip(sift1, sift2):

        if x and not y:  # TF
            out.append(0)  # UNMATCHED

        elif not x and not y:  # FF
            out.append(2)  # MATCH NOT POSSIBLE

        elif x and y:  # TT
            out.append(1)  # MATCHED

        elif not x and y:  # FT
            raise SiftMatchError

        else:
            raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

    return out


def human_sift_match(sift_match, feature_sift=config.FEATURE_SIFT):
    """Takes a base-3 SIFt indicating contact matched-ness and converts it to a human readable form.

    Args:
        sift_match ([type]): [description]
        feature_sift ([type], optional): Defaults to config.FEATURE_SIFT. [description]

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """

    terms = []

    for e, fp in enumerate(sift_match):

        if fp == 2:  # MATCH NOT POSSIBLE
            continue

        elif fp == 1:
            terms.append('Matched {}'.format(feature_sift[e]))

        elif fp == 0:
            terms.append('Unmatched {}'.format(feature_sift[e]))

        else:
            raise ValueError

    terms.sort()
    return ':'.join(terms)
# endregion


def int2(x):
    '''
    Return integer from base 2 number.

    Can accept a list/tuple of 0s and 1s.
    '''

    if isinstance(x, collections.Iterable):
        x = ''.join([str(k) for k in x])
    else:
        x = str(x)

    return int(x, 2)


def is_digit(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def int3(x):
    '''
    Return integer from base 3 number.

    Can accept a list/tuple of 0s, 1s and 2s.
    '''

    if isinstance(x, collections.Iterable):
        x = ''.join([str(k) for k in x])
    else:
        x = str(x)

    return int(x, 3)


def selection_parser(selection_list, atom_list):
    '''
    Selection syntax:

    /<chain_id>/<res_num>[<ins_code>]/<atom_name>

    HET:<het_id>

    Other formats will be rejected, for now.
    You can omit fields as long as the number of `/` is correct.

    Selections are additive, so if you select chain A with /A//,
    adding /A/91/C23 won't make any difference.
    '''

    final_atom_list = set([])

    for selection in selection_list:

        # selection_dict = {
        #    'chain': None,
        #    'residue_number': None,
        #    'atom_name': None
        # }

        chain = None
        residue_number = None
        insertion_code = ' '
        # residue_range = None # TODO
        atom_name = None

        current_atom_list = atom_list[:]

        original_selection = selection
        selection = selection.strip()

        if selection.startswith('RESNAME:'):

            selection = selection.replace('RESNAME:', '').strip()

            # RESNAMES ARE MAX LENGTH 3
            if len(selection) > 3:
                raise SelectionError(original_selection)

            current_atom_list = [x for x in current_atom_list if x.get_parent().resname.strip() == selection]

            for selected_atom in current_atom_list:
                final_atom_list.add(selected_atom)

        # SELECT ALL ORGANIC SMALL-MOLECULE LIGANDS
        elif selection.startswith('LIGANDS'):

            current_atom_list = [x for x in current_atom_list if
                                 not x.get_parent().is_polypeptide  # MUST NOT BE POLYPEPTIDE
                                 and len(x.get_parent().child_list) >= 5  # MIN NUMBER OF ATOMS
                                 and len(x.get_parent().child_list) <= 100  # MAX NUMBER OF ATOMS
                                 and 'C' in set([y.element for y in x.get_parent().child_list])  # MUST CONTAIN CARBON
                                 and x.get_parent().resname.strip().upper() not in config.COMMON_SOLVENTS  # MUST NOT BE COMMON SOLVENT
                                 and x.get_parent().resname.strip().upper() not in config.STANDARD_NUCLEOTIDES  # MUST NOT BE NUCLEOTIDE
                                 and not x.get_parent().resname.startswith('+')  # MUST NOT BE MODIFIED NUCLEOTIDE
                                 ]

            for selected_atom in current_atom_list:
                final_atom_list.add(selected_atom)

        elif selection.startswith('/'):

            selection = selection.lstrip('/').split('/')

            # print selection

            if len(selection) != 3:
                raise SelectionError(original_selection)

            # CHAIN
            if selection[0]:
                chain = selection[0]

            # RESIDUE AND INS CODE
            if selection[1]:

                if is_digit(selection[1]):

                    # JUST THE RESNUM
                    residue_number = int(selection[1])

                elif selection[1].isalnum():

                    # CHECK FOR VALID RESNUM+INSCODE
                    if selection[1][-1].isalpha() and is_digit(selection[1][:-1]):

                        residue_number = int(selection[1][:-1])
                        insertion_code = selection[1][-1]

                    else:
                        raise SelectionError(original_selection)

                else:
                    raise SelectionError(original_selection)

            # ATOM NAME
            if selection[2]:

                if not selection[2].isalnum() and "'" not in selection[2]:
                    raise SelectionError(original_selection)
                else:
                    atom_name = selection[2]

            # NOW MAKE THE SELECTION WITH BIOPYTHON
            if chain:
                current_atom_list = [x for x in current_atom_list if x.get_parent().get_parent().id == chain]

            if residue_number:
                current_atom_list = [x for x in current_atom_list if x.get_parent().id[1] == residue_number and x.get_parent().id[2] == insertion_code]

            if atom_name:
                current_atom_list = [x for x in current_atom_list if x.name == atom_name]

            for selected_atom in current_atom_list:
                final_atom_list.add(selected_atom)

        else:
            raise SelectionError(original_selection)

    final_atom_list = list(final_atom_list)

    if not final_atom_list:
        logging.error('Selection was empty.')
        raise SelectionError('entity not found')

    return list(final_atom_list)



def make_pymol_json(entity):
    '''
    Feed me a BioPython atom or BioPython residue.

    See `http://pymol.sourceforge.net/newman/user/S0220commands.html`.

    chain-identifier/resi-identifier/name-identifier
    chain-identifier/resi-identifier/
    '''

    if isinstance(entity, Atom):
        chain = entity.get_parent().get_parent()
        residue = entity.get_parent()
        atom_name = entity.name

        return {
            'label_comp_id': residue.resname,
            'auth_seq_id': residue.id[1],
            'auth_asym_id': chain.id,
            'auth_atom_id': atom_name,
            'pdbx_PDB_ins_code': residue.id[2]
        }

    elif isinstance(entity, Residue):
        chain = entity.get_parent()

        return {
            'label_comp_id': entity.resname,
            'auth_seq_id': entity.id[1],
            'auth_asym_id': chain.id,
            'pdbx_PDB_ins_code': entity.id[2]
        }

    else:
        raise TypeError('Cannot make a json object from non-Atom/Residue object.')


def make_pymol_string(entity):
    """Feed me a BioPython atom or BioPython residue.

    See `http://pymol.sourceforge.net/newman/user/S0220commands.html`.

    chain-identifier/resi-identifier/name-identifier
    chain-identifier/resi-identifier/

    Args:
        entity ([type]): [description]

    Raises:
        TypeError: [description]

    Returns:
        str: String representation of the entity
    """

    if isinstance(entity, Atom):

        chain = entity.get_parent().get_parent()
        residue = entity.get_parent()
        atom_name = entity.name

    elif isinstance(entity, Residue):
        chain = entity.get_parent()
        residue = entity
        atom_name = ''

    else:
        raise TypeError('Cannot make a PyMOL string from a non-Atom or Residue object.')

    res_num = residue.id[1]

    # ADD INSERTION CODE IF NEED BE
    if residue.id[2] != ' ':
        res_num = str(res_num) + residue.id[2]

    macro = '{}/{}/{}'.format(chain.id,
                              res_num,
                              atom_name)

    return macro


def get_single_bond_neighbour(ob_atom):
    """[summary]

    Args:
        ob_atom ([type]): [description]

    Returns:
        [type]: [description]
    """'''
    '''

    for bond in ob.OBAtomBondIter(ob_atom):

        if not (bond.GetBondOrder()==1 and not bond.IsAromatic()):
            continue

        current_neighbour = bond.GetNbrAtom(ob_atom)

        if current_neighbour.GetAtomicNum() == 1:
            continue

        return current_neighbour

    return None


def group_angle(group, point_coords, degrees=False, signed=False):
    """
    Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`

    `group` should be a dict with Numpy array 'center' and 'normal' attributes.
    `point_coords` should be a Numpy array.
    """

    cosangle = np.dot(group['normal'], point_coords) / (np.linalg.norm(group['normal']) * np.linalg.norm(point_coords))

    # GET THE ANGLE AS RADIANS
    rad = np.arccos(cosangle)

    if not degrees:
        return rad

    # CONVERT RADIANS INTO DEGREES
    # CONVERT INTO A SIGNED ANGLE
    if signed:
        rad = rad - np.pi if rad > np.pi / 2 else rad

    # RETURN DEGREES
    return rad * 180 / np.pi


def group_group_angle(group, group2, degrees=False, signed=False):
    """Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`

    `group` and `group2` should be a dict with Numpy array 'center' and 'normal' attributes.
    Args:
        group ([type]): [description]
        group2 ([type]): [description]
        degrees (bool, optional): Defaults to False. [description]
        signed (bool, optional): Defaults to False. [description]

    Returns:
        [type]: [description]
    """

    cosangle = np.dot(group['normal'], group2['normal']) / (np.linalg.norm(group['normal']) * np.linalg.norm(group2['normal']))

    # GET THE ANGLE AS RADIANS
    rad = np.arccos(cosangle)

    if not degrees:
        return rad

    # CONVERT RADIANS INTO DEGREES
    else:

        # CONVERT INTO A SIGNED ANGLE
        if signed:
            rad = rad - np.pi if rad > np.pi / 2 else rad

        # RETURN DEGREES
        return rad * 180 / np.pi


def get_angle(point_a, point_b, point_c):
    """Get the angle between three points in 3D space.
    Points should be supplied in Numpy array format.

    http://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space

    Args:
        point_a ([type]): [description]
        point_b ([type]): [description]
        point_c ([type]): [description]

    Returns:
        [type]: [description]
    """
    # In pseudo-code, the vector BA (call it v1) is:
    # v1 = {A.x - B.x, A.y - B.y, A.z - B.z}
    v1 = point_a - point_b

    # Similarly the vector BC (call it v2) is:
    # v2 = {C.x - B.x, C.y - B.y, C.z - B.z}
    v2 = point_c - point_b

    # The dot product of v1 and v2 is a function of the cosine of the angle between them
    # (it's scaled by the product of their magnitudes). So first normalize v1 and v2:

    # v1mag = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
    v1_mag = np.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])

    # v1norm = {v1.x / v1mag, v1.y / v1mag, v1.z / v1mag}
    v1_norm = np.array([v1[0] / v1_mag, v1[1] / v1_mag, v1[2] / v1_mag])

    # v2mag = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
    v2_mag = np.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])

    # v2norm = {v2.x / v2mag, v2.y / v2mag, v2.z / v2mag}
    v2_norm = np.array([v2[0] / v2_mag, v2[1] / v2_mag, v2[2] / v2_mag])

    # Then calculate the dot product:

    # res = v1norm.x * v2norm.x + v1norm.y * v2norm.y + v1norm.z * v2norm.z
    res = v1_norm[0] * v2_norm[0] + v1_norm[1] * v2_norm[1] + v1_norm[2] * v2_norm[2]

    # And finally, recover the angle:
    angle = np.arccos(res)

    if np.isnan(angle):
        logging.warning(f'Angle for <{point_a}, {point_b}, {point_c}> was NaN, setting to PI.')
        angle = np.pi

    return angle


def get_residue_name(entity):
    """
    Returns the name of residue

    Args:
        entity (Atom/Residue): Atom or Residue objects from BioPython
    
    Returns:
        str: Name of the residue
    """

    if isinstance(entity, Atom):
        residue = entity.get_parent()
        return residue.get_resname()
    
    elif isinstance(entity, Residue):
        return entity.get_resname()
    
    else:
        raise TypeError('Cannot return Residue from from non-Atom/non-Residue object.')
    
        






