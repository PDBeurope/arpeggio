#!/usr/bin/python

##################################################################################################################################
#                                                                                                                                #
# CREDO-BASED DEFINITIONS:                                                                                                       #
#                                                                                                                                #
# https://bitbucket.org/harryjubb/credovi/src/ac69222542134bb86d26a24561e4f740466e2b0e/credovi/config/credo.json?at=default      #
# https://bitbucket.org/harryjubb/credovi/src/ac69222542134bb86d26a24561e4f740466e2b0e/credovi/structbio/structure.py?at=default #
#                                                                                                                                #
##################################################################################################################################

###########
# IMPORTS #
###########

import argparse
import collections
import logging
import operator
import sys
from collections import OrderedDict
from functools import reduce
from os import path
import platform

import numpy as np
import openbabel as ob
from Bio.PDB import NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Residue import Residue

import arpeggio.core.config as config
from arpeggio.core import HydrogenError, OBBioMatchError, AtomSerialError, SiftMatchError, SelectionError
from arpeggio.core import protein_reader

if __name__ == '__main__':

    # RING-RING INTERACTIONS
    # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`
    # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/sql/populate.sql?at=default`
    # `http://marid.bioc.cam.ac.uk/credo/about`
    with open(rename_output_file(filename, '.ri'), 'w') as fo:

        if args.headers:

            fo.write('{}\n'.format('\t'.join(
                ['ring_bgn_id',
                 'ring_bgn_residue',
                 'ring_bgn_centroid',
                 'ring_end_id',
                 'ring_end_residue',
                 'ring_end_centroid',
                 'interaction_type',
                 'residue_interaction',
                 'contact_type'
                 ]
            )))

        for ring in s.rings:

            ring_key = ring
            ring = s.rings[ring]

            for ring2 in s.rings:

                ring_key2 = ring2
                ring2 = s.rings[ring2]

                # CHECK THAT THE RINGS ARE INVOLVED WITH THE SELECTION OR BINDING SITE
                if ring_key not in selection_plus_ring_ids or ring_key2 not in selection_plus_ring_ids:
                    continue

                # NO SELFIES
                if ring_key == ring_key2:
                    continue

                # CHECK IF INTERACTION IS WITHIN SAME RESIDUE
                intra_residue = False

                if ring['residue'] == ring2['residue']:
                    intra_residue = True

                # DETERMINE CONTACT TYPE
                contact_type = ''

                if ring_key not in selection_ring_ids and ring_key2 not in selection_ring_ids:
                    contact_type = 'INTRA_NON_SELECTION'

                if ring_key in selection_plus_ring_ids and ring_key2 in selection_plus_ring_ids:
                    contact_type = 'INTRA_BINDING_SITE'

                if ring_key in selection_ring_ids and ring_key2 in selection_ring_ids:
                    contact_type = 'INTRA_SELECTION'

                if (ring_key in selection_ring_ids and ring_key2 not in selection_ring_ids) or (ring_key2 in selection_ring_ids and ring_key not in selection_ring_ids):
                    contact_type = 'INTER'

                # DETERMINE RING-RING DISTANCE
                distance = np.linalg.norm(ring['center'] - ring2['center'])

                if distance > config.CONTACT_TYPES['aromatic']['centroid_distance']:
                    continue

                theta_point = ring['center'] - ring2['center']
                # iota_point = ring2['center'] - ring['center']

                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                dihedral = abs(group_group_angle(ring, ring2, True, True))
                theta = abs(group_angle(ring, theta_point, True, True))

                # logging.info('Dihedral = {}     Theta = {}'.format(dihedral, theta))

                int_type = ''

                if dihedral <= 30.0 and theta <= 30.0:
                    int_type = 'FF'
                elif dihedral <= 30.0 and theta <= 60.0:
                    int_type = 'OF'
                elif dihedral <= 30.0 and theta <= 90.0:
                    int_type = 'EE'

                elif dihedral > 30.0 and dihedral <= 60.0 and theta <= 30.0:
                    int_type = 'FT'
                elif dihedral > 30.0 and dihedral <= 60.0 and theta <= 60.0:
                    int_type = 'OT'
                elif dihedral > 30.0 and dihedral <= 60.0 and theta <= 90.0:
                    int_type = 'ET'

                elif dihedral > 60.0 and dihedral <= 90.0 and theta <= 30.0:
                    int_type = 'FE'
                elif dihedral > 60.0 and dihedral <= 90.0 and theta <= 60.0:
                    int_type = 'OE'
                elif dihedral > 60.0 and dihedral <= 90.0 and theta <= 90.0:
                    int_type = 'EF'

                # DON'T COUNT INTRA-RESIDUE EDGE-TO-EDGE RING INTERACTIONS
                # TO AVOID INTRA-HETEROCYCLE INTERACTIONS
                # POTENTIALLY A BUG IF TWO RINGS ARE SEPARATED BY A LONG ALIPHATIC CHAIN
                # AND MAKE A GENUINE INTRA EE INTERACTION, BUT ASSUMING THIS IS RARE
                if intra_residue and int_type == 'EE':
                    continue

                # OUTPUT INTRA/INTER RESIDUE AS TEXT RATHER THAN BOOLEAN
                intra_residue_text = 'INTER_RESIDUE'

                if intra_residue:
                    intra_residue_text = 'INTRA_RESIDUE'

                # UPDATE RESIDUE RING INTEGER SIFTS

                # RING-RING INTERACTING SIFT
                # MUST BE INTER-RESIDUE
                # MUST BE OF CONTACT TYPE INTER, I.E. BETWEEN SELECTION AND NON-SELECTION
                # ALL 9 RING INTERACTION TYPES
                # ADDED ONLY FOR THE FIRST RING, AS THE RECIPROCAL INTERACTION
                # SHOULD BE COVERED FOR THE OTHER RING
                if contact_type == 'INTER' and not intra_residue:

                    for k, i_type in enumerate(('FF', 'OF', 'EE', 'FT', 'OT', 'ET', 'FE', 'OE', 'EF')):

                        if int_type == i_type:
                            ring['residue'].ring_ring_inter_integer_sift[k] = ring['residue'].ring_ring_inter_integer_sift[k] + 1

                # WRITE RING INTERACTION TO FILE
                output = [
                    ring['ring_id'],
                    make_pymol_string(ring['residue']),
                    list(ring['center']),
                    ring2['ring_id'],
                    make_pymol_string(ring2['residue']),
                    list(ring2['center']),
                    int_type,
                    intra_residue_text,
                    contact_type
                ]

                fo.write('{}\n'.format('\t'.join([str(x) for x in output])))

    # RINGS AND ATOM-RING INTERACTIONS
    with open(rename_output_file(filename, '.ari'), 'w') as fo, open(rename_output_file(filename, '.rings'), 'w') as ring_fo:

        if args.headers:

            fo.write('{}\n'.format('\t'.join(
                ['atom',
                 'ring_id',
                 'ring_residue',
                 'ring_centroid',
                 'interactions',
                 'residue_interactions',
                 'contact_type'
                 ]
            )))

            ring_fo.write('{}\n'.format('\t'.join(
                ['ring_id',
                 'ring_residue',
                 'ring_centroid'
                 ]
            )))

        for ring in s.rings:

            ring_key = ring
            ring = s.rings[ring]

            # CHECK RING IS INVOLVED IN THE SELECTION OR BINDING SITE
            if ring_key not in selection_plus_ring_ids:
                continue

            for atom in ns.search(ring['center'], config.CONTACT_TYPES['aromatic']['met_sulphur_aromatic_distance']):

                # IGNORE ANY HYDROGENS FOR THESE CONTACTS
                # IF HYDROGENS ARE PRESENT IN THE BIOPYTHON STRUCTURE FOR ANY REASON
                if atom.element.strip() == 'H':
                    continue

                # CHECK THAT THE ATOM IS INVOLVED IN THE SELECTION OR BINDING SITE
                if atom not in selection_plus:
                    continue

                # GET DISTANCE AND CHECK IF FAR ENOUGH
                distance = np.linalg.norm(atom.coord - ring['center'])

                # NO AROMATIC ATOM-RING INTERACTIONS
                if 'aromatic' in atom.atom_types:
                    continue

                # CHECK IF INTRA-RESIDUE
                intra_residue = False

                if ring['residue'] == atom.get_parent():
                    intra_residue = True

                # DETERMINE CONTACT TYPE
                contact_type = ''

                if ring_key not in selection_ring_ids and atom not in selection_set:
                    contact_type = 'INTRA_NON_SELECTION'

                if ring_key in selection_plus_ring_ids and atom in selection_plus:
                    contact_type = 'INTRA_BINDING_SITE'

                if ring_key in selection_ring_ids and atom in selection_set:
                    contact_type = 'INTRA_SELECTION'

                if (ring_key in selection_ring_ids and atom not in selection_set) or (atom in selection_set and ring_key not in selection_ring_ids):
                    contact_type = 'INTER'

                # DETERMINE INTERACTIONS
                potential_interactions = set([])

                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                theta = abs(group_angle(ring, ring['center'] - atom.coord, True, True))  # CHECK IF `atom.coord` or `ring['center'] - atom.coord`

                if distance <= config.CONTACT_TYPES['aromatic']['atom_aromatic_distance'] and theta <= 30.0:

                    if atom.element == 'C' and 'weak hbond donor' in atom.atom_types:
                        potential_interactions.add('CARBONPI')

                    if 'pos ionisable' in atom.atom_types:
                        potential_interactions.add('CATIONPI')

                    if 'hbond donor' in atom.atom_types:
                        potential_interactions.add('DONORPI')

                    if 'xbond donor' in atom.atom_types:
                        potential_interactions.add('HALOGENPI')

                if distance <= config.CONTACT_TYPES['aromatic']['met_sulphur_aromatic_distance']:

                    if atom.get_parent().resname == 'MET' and atom.element == 'S':
                        potential_interactions.add('METSULPHURPI')

                if not potential_interactions:
                    continue

                interaction_type = list(potential_interactions)[0]

                # OUTPUT INTRA/INTER RESIDUE AS TEXT RATHER THAN BOOLEAN
                intra_residue_text = 'INTER_RESIDUE'

                if intra_residue:
                    intra_residue_text = 'INTRA_RESIDUE'

                # logging.info('Atom: <{}>   Ring: <{}>  Theta = {}'.format(atom.get_full_id(), ring['ring_id'], theta))

                # RESIDUE RING-ATOM SIFT
                if contact_type == 'INTER' and not intra_residue:

                    for k, i_type in enumerate(('CARBONPI', 'CATIONPI', 'DONORPI', 'HALOGENPI', 'METSULPHURPI')):

                        for potential_interaction in potential_interactions:

                            if potential_interaction == i_type:

                                ring['residue'].ring_atom_inter_integer_sift[k] = ring['residue'].ring_atom_inter_integer_sift[k] + 1
                                atom.get_parent().atom_ring_inter_integer_sift[k] = atom.get_parent().atom_ring_inter_integer_sift[k] + 1

                                if atom.get_parent() in polypeptide_residues:

                                    if atom.name in config.MAINCHAIN_ATOMS:
                                        atom.get_parent().mc_atom_ring_inter_integer_sift[k] = atom.get_parent().mc_atom_ring_inter_integer_sift[k] + 1

                                    else:
                                        atom.get_parent().sc_atom_ring_inter_integer_sift[k] = atom.get_parent().sc_atom_ring_inter_integer_sift[k] + 1

                # WRITE ATOM-RING INTERACTION TO FILE
                output = [
                    make_pymol_string(atom),
                    ring['ring_id'],
                    make_pymol_string(ring['residue']),
                    list(ring['center']),
                    sorted(list(potential_interactions)),
                    intra_residue_text,
                    contact_type
                ]

                fo.write('{}\n'.format('\t'.join([str(x) for x in output])))

            # WRITE RING OUT TO RING FILE
            output = [
                ring['ring_id'],
                make_pymol_string(ring['residue']),
                list(ring['center'])
            ]

            ring_fo.write('{}\n'.format('\t'.join([str(x) for x in output])))

    # AMIDE-RING INTERACTIONS
    with open(rename_output_file(filename, '.amri'), 'w') as fo:

        if args.headers:

            fo.write('{}\n'.format('\t'.join(
                ['amide_id',
                 'residue',
                 'amide_centroid',
                 'ring_id',
                 'ring_residue',
                 'ring_centroid',
                 'interaction_type',
                 'residue_interactions',
                 'contact_type'
                 ]
            )))

        for amide in s.amides:

            amide_key = amide
            amide = s.amides[amide]

            # CHECK AMIDE IS INVOLVED IN THE SELECTION OR BINDING SITE
            if amide_key not in selection_plus_amide_ids:
                continue

            for ring in s.rings:

                ring_key = ring
                ring = s.rings[ring]

                # CHECK RING IS INVOLVED WITH THE SELECTION OR BINDING SITE
                if ring_key not in selection_plus_ring_ids:
                    continue

                # CHECK IF INTERACTION IS WITHIN SAME RESIDUE
                intra_residue = False

                if amide['residue'] == ring['residue']:
                    intra_residue = True

                # OUTPUT INTRA/INTER RESIDUE AS TEXT RATHER THAN BOOLEAN
                intra_residue_text = 'INTER_RESIDUE'

                if intra_residue:
                    intra_residue_text = 'INTRA_RESIDUE'

                # DETERMINE CONTACT TYPE
                contact_type = ''

                if amide_key not in selection_amide_ids and ring_key not in selection_ring_ids:
                    contact_type = 'INTRA_NON_SELECTION'

                if amide_key in selection_plus_amide_ids and ring_key in selection_plus_ring_ids:
                    contact_type = 'INTRA_BINDING_SITE'

                if amide_key in selection_amide_ids and ring_key in selection_ring_ids:
                    contact_type = 'INTRA_SELECTION'

                if (amide_key in selection_amide_ids and ring_key not in selection_ring_ids) or (ring_key in selection_ring_ids and amide_key not in selection_amide_ids):
                    contact_type = 'INTER'

                # DETERMINE AMIDE-RING DISTANCE
                distance = np.linalg.norm(amide['center'] - ring['center'])

                if distance > config.CONTACT_TYPES['amide']['centroid_distance']:
                    continue

                theta_point = amide['center'] - ring['center']

                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                dihedral = abs(group_group_angle(amide, ring, True, True))
                theta = abs(group_angle(amide, theta_point, True, True))

                # FACE-ON ORIENTATION ONLY
                if dihedral > 30.0 or theta > 30.0:
                    continue

                # IF IT'S SURVIVED SO FAR...(!)

                int_type = 'AMIDERING'

                # SIFTs
                if contact_type == 'INTER' and not intra_residue:
                    amide['residue'].amide_ring_inter_integer_sift[0] = amide['residue'].amide_ring_inter_integer_sift[0] + 1
                    ring['residue'].ring_amide_inter_integer_sift[0] = ring['residue'].ring_amide_inter_integer_sift[0] + 1

                # WRITE TO FILE
                output = [
                    amide['amide_id'],
                    make_pymol_string(amide['residue']),
                    list(amide['center']),
                    ring['ring_id'],
                    make_pymol_string(ring['residue']),
                    list(ring['center']),
                    int_type,
                    intra_residue_text,
                    contact_type,
                    # dihedral,
                    # theta,
                    # list(theta_point)
                ]

                fo.write('{}\n'.format('\t'.join([str(x) for x in output])))

    # AMIDE-AMIDE INTERACTIONS
    with open(rename_output_file(filename, '.amam'), 'w') as fo:

        if args.headers:
            fo.write('{}\n'.format('\t'.join(
                ['amide_bgn_id',
                 'amide_bgn_residue',
                 'amide_bgn_centroid',
                 'amide_end_id',
                 'amide_end_residue',
                 'amide_end_centroid',
                 'interaction_type',
                 'residue_interactions',
                 'contact_type'
                 ]
            )))

        for amide in s.amides:

            amide_key = amide
            amide = s.amides[amide]

            # CHECK AMIDE IS INVOLVED IN THE SELECTION OR BINDING SITE
            if amide_key not in selection_plus_amide_ids:
                continue

            for amide2 in s.amides:

                amide_key2 = amide2
                amide2 = s.amides[amide2]

                # NO SELFIES
                if amide_key == amide_key2:
                    continue

                # CHECK RING IS INVOLVED WITH THE SELECTION OR BINDING SITE
                if amide_key2 not in selection_plus_amide_ids:
                    continue

                # CHECK IF INTERACTION IS WITHIN SAME RESIDUE
                intra_residue = False

                if amide['residue'] == amide2['residue']:
                    intra_residue = True

                # OUTPUT INTRA/INTER RESIDUE AS TEXT RATHER THAN BOOLEAN
                intra_residue_text = 'INTER_RESIDUE'

                if intra_residue:
                    intra_residue_text = 'INTRA_RESIDUE'

                # DETERMINE CONTACT TYPE
                contact_type = ''

                if amide_key not in selection_amide_ids and amide_key2 not in selection_amide_ids:
                    contact_type = 'INTRA_NON_SELECTION'

                if amide_key in selection_plus_amide_ids and amide_key2 in selection_plus_amide_ids:
                    contact_type = 'INTRA_BINDING_SITE'

                if amide_key in selection_amide_ids and amide_key2 in selection_amide_ids:
                    contact_type = 'INTRA_SELECTION'

                if (amide_key in selection_amide_ids and amide_key2 not in selection_amide_ids) or (amide_key2 in selection_amide_ids and amide_key not in selection_amide_ids):
                    contact_type = 'INTER'

                # DETERMINE AMIDE-RING DISTANCE
                distance = np.linalg.norm(amide['center'] - amide2['center'])

                if distance > config.CONTACT_TYPES['amide']['centroid_distance']:
                    continue

                theta_point = amide['center'] - amide2['center']

                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                dihedral = abs(group_group_angle(amide, amide2, True, True))
                theta = abs(group_angle(amide, theta_point, True, True))

                # FACE-ON ORIENTATION ONLY
                if dihedral > 30.0 or theta > 30.0:
                    continue

                # IF IT'S SURVIVED SO FAR...(!)

                int_type = 'AMIDEAMIDE'

                # SIFT
                if contact_type == 'INTER' and not intra_residue:
                    amide['residue'].amide_amide_inter_integer_sift[0] = amide['residue'].amide_amide_inter_integer_sift[0] + 1

                # WRITE TO FILE
                output = [
                    amide['amide_id'],
                    make_pymol_string(amide['residue']),
                    list(amide['center']),
                    amide2['amide_id'],
                    make_pymol_string(amide2['residue']),
                    list(amide2['center']),
                    int_type,
                    intra_residue_text,
                    contact_type,
                    # dihedral,
                    # theta,
                    # list(theta_point)
                ]

                fo.write('{}\n'.format('\t'.join([str(x) for x in output])))

    logging.info('Program End. Maximum memory usage was {}.'.format(max_mem_usage()))
