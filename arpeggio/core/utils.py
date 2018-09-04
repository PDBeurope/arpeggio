import os
import sys
import collections
import logging
import platform

from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
import openbabel as ob

from arpeggio.core import SelectionError, config


def rename_output_file(original_filename, new_extension):

    original_extension = os.path.splitext(original_filename)[-1]

    if original_extension:
        return os.path.splitext(original_filename)[0] + new_extension
    else:
        return original_filename + new_extension


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

                if selection[1].isdigit():

                    # JUST THE RESNUM
                    residue_number = int(selection[1])

                elif selection[1].isalnum():

                    # CHECK FOR VALID RESNUM+INSCODE
                    if selection[1][-1].isalpha() and selection[1][:-1].isdigit():

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

    if len(final_atom_list) == 0:
        logging.error('Selection was empty.')
        sys.exit(1)

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

    else:
        raise TypeError('Cannot make a json object from non-Atom object.')

    return {
        'res_name': residue.resname,
        'auth_seq_id': residue.id[1],
        'auth_asym_id': chain.id,
        'auth_atom_id': atom_name,
        'pdbx_PDB_ins_code': residue.id[2]
    }


def make_pymol_string(entity):
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
    '''
    '''

    for bond in ob.OBAtomBondIter(ob_atom):

        if not bond.IsSingle():
            continue

        current_neighbour = bond.GetNbrAtom(ob_atom)

        if current_neighbour.IsHydrogen():
            continue

        return current_neighbour

    return None


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
