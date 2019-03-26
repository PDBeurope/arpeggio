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
