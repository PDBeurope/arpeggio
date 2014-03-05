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
#import math
import operator
import resource
import sys

from collections import OrderedDict

import numpy as np

#from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch

import openbabel as ob

#############
# CONSTANTS #
#############

from config import DEFAULT_SIFT, ATOM_TYPES, CONTACT_TYPES, VDW_RADII, METALS, HALOGENS, CONTACT_TYPES_DIST_MAX, FEATURE_SIFT

###########
# CLASSES #
###########

class HydrogenError(Exception):
    
    def __init__(self):
        logging.error('Please remove all hydrogens from the structure then re-run.')

class OBBioMatchError(Exception):
    
    def __init__(self):
        logging.error('An OpenBabel atom could not be matched to a BioPython counterpart.')

#############
# FUNCTIONS #
#############

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
    
    Other formats will be rejected, for now.
    You can omit fields as long as the number of `/` is correct.
    
    Selections are additive, so if you select chain A with /A//,
    adding /A/91/C23 won't make any difference.
    '''
    
    final_atom_list = set([])
    
    for selection in selection_list:
        
        #selection_dict = {
        #    'chain': None,
        #    'residue_number': None,
        #    'atom_name': None
        #}
        
        chain = None
        residue_number = None
        insertion_code = ' '
        #residue_range = None # TODO
        atom_name = None
        
        current_atom_list = atom_list[:]
        
        original_selection = selection
        selection = selection.strip()
        
        if selection.startswith('/'):
            
            selection = selection.lstrip('/').split('/')
            
            #print selection
            
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
                
                if not selection[2].isalnum():
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
        sys.exit()
    
    return list(final_atom_list)

def make_pymol_string(atom):
    '''
    Feed me a BioPython atom.
    
    See `http://pymol.sourceforge.net/newman/user/S0220commands.html`.
    
    chain-identifier/resi-identifier/name-identifier
    '''
    
    chain = atom.get_parent().get_parent()
    residue = atom.get_parent()
    
    res_num = residue.id[1]
    
    # ADD INSERTION CODE IF NEED BE
    if residue.id[2] != ' ':
        res_num = str(res_num) + residue.id[2]
    
    macro =  '{}/{}/{}'.format(chain.id,
                               res_num,
                               atom.name)
    
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
    '''
    Returns maximum memory usage of the program thus far, in megabytes, as a string.
    '''
    
    return str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.0) + ' MB'

def get_angle(point_a, point_b, point_c):
    '''
    Get the angle between three points in 3D space.
    Points should be supplied in Numpy array format.
    
    http://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space
    '''
    
    #In pseudo-code, the vector BA (call it v1) is:
    #v1 = {A.x - B.x, A.y - B.y, A.z - B.z}
    v1 = point_a - point_b
    
    #Similarly the vector BC (call it v2) is:
    #v2 = {C.x - B.x, C.y - B.y, C.z - B.z}
    v2 = point_c - point_b
    
    #The dot product of v1 and v2 is a function of the cosine of the angle between them
    #(it's scaled by the product of their magnitudes). So first normalize v1 and v2:
    
    #v1mag = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
    v1_mag = np.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    
    #v1norm = {v1.x / v1mag, v1.y / v1mag, v1.z / v1mag}
    v1_norm = np.array([v1[0] / v1_mag, v1[1] / v1_mag, v1[2] / v1_mag])
    
    #v2mag = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
    v2_mag = np.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    
    #v2norm = {v2.x / v2mag, v2.y / v2mag, v2.z / v2mag}
    v2_norm = np.array([v2[0] / v2_mag, v2[1] / v2_mag, v2[2] / v2_mag])
    
    #Then calculate the dot product:
    
    #res = v1norm.x * v2norm.x + v1norm.y * v2norm.y + v1norm.z * v2norm.z
    res = v1_norm[0] * v2_norm[0] + v1_norm[1] * v2_norm[1] + v1_norm[2] * v2_norm[2]
    
    #And finally, recover the angle:
    angle = np.arccos(res)
    
    if np.isnan(angle):
        logging.warn('Angle for <{}, {}, {}> was NaN, setting to Pi.'.format(point_a, point_b, point_c))
        angle = np.pi
    
    return angle

def ring_angle(ring, point_coords, degrees=False, signed=False):
    '''
    Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`
    
    `ring` should be a dict with Numpy array 'center' and 'normal' attributes.
    `point_coords` should be a Numpy array.
    '''
    
    cosangle = np.dot(ring['normal'], point_coords) / (np.linalg.norm(ring['normal']) * np.linalg.norm(point_coords))

    # GET THE ANGLE AS RADIANS
    rad = np.arccos(cosangle)

    if not degrees: return rad

    # CONVERT RADIANS INTO DEGREES
    else:

        # CONVERT INTO A SIGNED ANGLE
        if signed: rad = rad -np.pi if rad > np.pi / 2 else rad

        # RETURN DEGREES
        return rad * 180 / np.pi
    
def ring_ring_angle(ring, ring2, degrees=False, signed=False):
    '''
    Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`
    
    `ring` and `ring2` should be a dict with Numpy array 'center' and 'normal' attributes.
    '''
    
    cosangle = np.dot(ring['normal'], ring2['normal']) / (np.linalg.norm(ring['normal']) * np.linalg.norm(ring2['normal']))

    # GET THE ANGLE AS RADIANS
    rad = np.arccos(cosangle)

    if not degrees: return rad

    # CONVERT RADIANS INTO DEGREES
    else:

        # CONVERT INTO A SIGNED ANGLE
        if signed: rad = rad -np.pi if rad > np.pi / 2 else rad

        # RETURN DEGREES
        return rad * 180 / np.pi

## GOLDEN SECTION SPIRAL
## THANK YOU BOSCOH!
## http://boscoh.com/protein/calculating-the-solvent-accessible-surface-area-asa.html
## http://boscoh.com/protein/asapy.html
#def points_on_sphere(n):
#    pts = np.empty((int(n), 3))
#    n = float(n)
#    
#    inc = math.pi * (3 - math.sqrt(5))
#    off = 2 / n
#    
#    for k in range(0, int(n)):
#        
#        y = k * off - 1 + (off / 2)
#        r = math.sqrt(1 - y*y)
#        phi = k * inc
#        pts[k][0] = math.cos(phi) * r
#        pts[k][1] = y
#        pts[k][2] = math.sin(phi) * r
#        
#    return pts

# CONTACT FUNCTIONS

def is_hbond(donor, acceptor):
    '''
    Feed me BioPython atoms.
    '''
    
    for hydrogen_coord in donor.h_coords:
        
        h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)
        
        if h_dist <= VDW_RADII['H'] + acceptor.vdw_radius + VDW_COMP_FACTOR:
            
            if get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= CONTACT_TYPES['hbond']['angle rad']:
                
                return 1
    
    return 0

def is_weak_hbond(donor, acceptor):
    '''
    Feed me BioPython atoms.
    '''
    
    for hydrogen_coord in donor.h_coords:
        
        h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)
        
        if h_dist <= VDW_RADII['H'] + acceptor.vdw_radius + VDW_COMP_FACTOR:
            
            if get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= CONTACT_TYPES['weak hbond']['angle rad']:
                
                return 1
    
    return 0

def is_halogen_weak_hbond(donor, halogen, ob_mol):
    '''
    Feed me BioPython atoms and the OpenBabel molecule.
    '''
    
    # `nbr` WILL BE A BIOPYTHON ATOM
    # ... HOPEFULLY
    nbr = ob_to_bio[get_single_bond_neighbour(ob_mol.GetAtomById(bio_to_ob[halogen])).GetId()]
    
    for hydrogen_coord in donor.h_coords:
        
        h_dist = np.linalg.norm(halogen.coord - hydrogen_coord)
        
        if h_dist <= VDW_RADII['H'] + halogen.vdw_radius + VDW_COMP_FACTOR:
            
            if CONTACT_TYPES['weak hbond']['cx angle min rad'] <= get_angle(nbr.coord, halogen.coord, hydrogen_coord) <= CONTACT_TYPES['weak hbond']['cx angle max rad']:
                
                return 1
    
    return 0

def is_xbond(donor, acceptor, ob_mol):
    '''
    Feed me BioPython atoms and the OpenBabel molecule.
    '''
    
    # `nbr` WILL BE A BIOPYTHON ATOM
    # ... HOPEFULLY
    nbr = ob_to_bio[get_single_bond_neighbour(ob_mol.GetAtomById(bio_to_ob[donor])).GetId()]
    theta = get_angle(nbr.coord, donor.coord, acceptor.coord)
    
    if (theta >= CONTACT_TYPES['xbond']['angle theta 1 rad']):
        return 1
    
    return 0

def update_atom_fsift(atom, addition, contact_type='INTER'):
    '''
    '''
    
    atom.actual_fsift = [x | y for x, y in zip(atom.actual_fsift, addition)]
    
    if contact_type == 'INTER':
        atom.actual_fsift_inter_only = [x | y for x, y in zip(atom.actual_fsift_inter_only, addition)]
    
    if 'INTRA' in contact_type:
        atom.actual_fsift_intra_only = [x | y for x, y in zip(atom.actual_fsift_intra_only, addition)]
        
    if 'WATER' in contact_type:
        atom.actual_fsift_water_only = [x | y for x, y in zip(atom.actual_fsift_water_only, addition)]

def sift_xnor(sift1, sift2):
    '''
    '''
    
    out = []
    
    for x, y in zip(sift1, sift2):
        
        if x and not y: # TF
            out.append(0)
        
        elif not x and not y: # FF
            out.append(1)
            
        elif x and y: # TT
            out.append(1)
            
        elif not x and y: # FT
            out.append(0)
        
        else:
            raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

    return out

def sift_match_base3(sift1, sift2):
    '''
    0 = UNMATCHED
    1 = MATCHED
    2 = MATCH NOT POSSIBLE
    
    Assuming that sift1 is the potential SIFt, and sift2 is the actual SIFt.
    '''
    
    out = []
    
    for x, y in zip(sift1, sift2):
        
        if x and not y: # TF
            out.append(0) # UNMATCHED
        
        elif not x and not y: # FF
            out.append(2) # MATCH NOT POSSIBLE
            
        elif x and y: # TT
            out.append(1) # MATCHED
            
        elif not x and y: # FT
            out.append(2) # SHOULD NEVER HAPPEN
        
        else:
            raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

    return out

def human_sift_match(sift_match, feature_sift=FEATURE_SIFT):
    '''
    Takes a base-3 SIFt indicating contact matched-ness and converts it to a human readable form.
    '''
    
    terms = []
    
    for e, fp in enumerate(sift_match):
        
        if fp == 2: # MATCH NOT POSSIBLE
            continue
        
        elif fp == 1:
            terms.append('Matched {}'.format(feature_sift[e]))
        
        elif fp == 0:
            terms.append('Unmatched {}'.format(feature_sift[e]))
        
        else:
            raise ValueError
    
    terms.sort()
    return ':'.join(terms)

########
# MAIN #
########

if __name__ == '__main__':
    
    # ARGUMENT PARSING
    parser = argparse.ArgumentParser(description='''

#############
# OPENCREDO #
#############

A program for calculating CREDO interactions,
using only Open Source dependencies.

Dependencies:
- Python (v2.7)
- Numpy
- BioPython (>= v1.60)
- OpenBabel (with Python bindings)

''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('pdb', type=str, help='Path to the PDB file to be analysed.')
    
    selection_group = parser.add_mutually_exclusive_group(required=False)
    selection_group.add_argument('-s', '--selection', type=str, nargs='+', help='Select the "ligand" for interactions, using selection syntax: /<chain_id>/<res_num>[<ins_code>]/<atom_name> . Fields can be omitted.')
    selection_group.add_argument('-sf', '--selection-file', type=str, help='Selections as above, but listed in a file.')

    
    parser.add_argument('-wh', '--write-hydrogenated', action='store_true', help='Write a PDB file including the added hydrogen coordinates.')
    parser.add_argument('-mh', '--minimise-hydrogens', action='store_true', help='Energy minimise OpenBabel added hydrogens.')
    parser.add_argument('-ms', '--minimisation-steps', type=int, default=50, help='Number of hydrogen minimisation steps to perform.')
    parser.add_argument('-mf', '--minimisation-forcefield', type=str, choices=('MMFF94', 'UFF', 'Ghemical'), default='MMFF94', help='Choose the forcefield to minimise hydrogens with. Ghemical is not recommended.')
    parser.add_argument('-mm', '--minimisation-method', type=str, choices=('DistanceGeometry', 'SteepestDescent', 'ConjugateGradients'), default='ConjugateGradients', help='Choose the method to minimise hydrogens with. ConjugateGradients is recommended.')
    parser.add_argument('-co', '--vdw-comp', type=float, default=0.1, help='Compensation factor for VdW radii dependent interaction types.')    
    parser.add_argument('-i', '--interacting', type=float, default=5.0, help='Distance cutoff for grid points to be \'interacting\' with the entity.')
    #parser.add_argument('-sr', '--solvent_radius', type=float, default=1.4, help='Radius of solvent probe for accessibility calculations.')
    #parser.add_argument('-ssp', '--solvent-sphere-points', type=int, default=960, help='Number of points to use for solvent shell spheres for accessibility calculations.')
    #parser.add_argument('-st', '--sasa-threshold', type=float, default=1.0, help='Floating point solvent accessible surface area threshold (squared Angstroms) for considering an atom as \'accessible\' or not.')
    #parser.add_argument('-ca', '--consider-all', action='store_true', help='Consider all entity/selection atoms, not just solvent accessible ones. If this is set, SASAs won\'t be calculated.')
    #parser.add_argument('-spdb', '--sasa-pdb', action='store_true', help='Store a PDB with atom b-factors set based on boolean solvent accessibility.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Be chatty.')
 
    args = parser.parse_args()
    
    pdb_filename = args.pdb
    
    VDW_COMP_FACTOR = args.vdw_comp
    INTERACTING_THRESHOLD = args.interacting
    #SOLVENT_RADIUS = args.solvent_radius
    #NUM_SOLVENT_SPHERE_POINTS = args.solvent_sphere_points
    #SASA_THRESHOLD = args.sasa_threshold
    
    # LOGGING
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format='%(levelname)s//%(asctime)s//%(message)s', datefmt='%H:%M:%S')
    else:
        logging.basicConfig(level=logging.WARN, format='%(levelname)s//%(asctime)s//%(message)s', datefmt='%H:%M:%S')
    
    logging.info('Program begin.')
    
    # LOAD STRUCTURE (BIOPYTHON)
    pdb_parser = PDBParser()
    s = pdb_parser.get_structure('structure', pdb_filename)
    s_atoms = list(s.get_atoms())

    logging.info('Loaded PDB structure (BioPython)')
    
    # CHECK FOR HYDROGENS AND RAISE ERROR IF THEY
    # EXIST. HYDROGENS INTERFERE WITH CONTACT CALCULATION.
    # COULD MAKE THIS A SMOOTHER PROCESS IN THE FUTURE.
    # FOR NOW REMOVE ANY HYDROGENS FROM THE INPUT FILE.
    
    hydrogens = [x for x in s_atoms if x.element == 'H']
 
    if hydrogens:
        raise HydrogenError
    
    # LOAD STRUCTURE (OPENBABEL)
    ob_conv = ob.OBConversion()
    ob_conv.SetInFormat('pdb')
    mol = ob.OBMol()
    ob_conv.ReadFile(mol, pdb_filename)
    
    logging.info('Loaded PDB structure (OpenBabel)')
    
    # MAPPING OB ATOMS TO BIOPYTHON ATOMS AND VICE VERSA
    
    # FIRST MAP PDB SERIAL NUMBERS TO BIOPYTHON ATOMS FOR SPEED LATER
    # THIS AVOIDS LOOPING THROUGH `s_atoms` MANY TIMES
    serial_to_bio = {x.serial_number: x for x in s_atoms}
    
    # DICTIONARIES FOR CONVERSIONS
    ob_to_bio = {}
    bio_to_ob = {}
    
    for ob_atom in ob.OBMolAtomIter(mol):
        
        serial = ob_atom.GetResidue().GetSerialNum(ob_atom)
        
        # MATCH TO THE BIOPYTHON ATOM BY SERIAL NUMBER
        try:
            biopython_atom = serial_to_bio[serial]
        
        except KeyError:
            # ERRORWORTHY IF WE CAN'T MATCH AN OB ATOM TO A BIOPYTHON ONE
            raise OBBioMatchError
        
        # `Id` IS A UNIQUE AND STABLE ID IN OPENBABEL
        # CAN RECOVER THE ATOM WITH `mol.GetAtomById(id)`
        ob_to_bio[ob_atom.GetId()] = biopython_atom
        bio_to_ob[biopython_atom] = ob_atom.GetId()
    
    logging.info('Mapped OB to BioPython atoms and vice-versa.')
    
    # ADD EMPTY DATA STRUCTURES FOR TAGGED ATOM DATA
    # IN A SINGLE ITERATION
    for atom in s_atoms:
        
        # FOR ATOM TYPING VIA OPENBABEL
        atom.atom_types = set([])
        
        # LIST FOR EACH ATOM TO STORE EXPLICIT HYDROGEN COORDINATES
        atom.h_coords = []
        
        # DETECT METALS
        if atom.element.upper() in METALS:
            atom.is_metal = True
        else:
            atom.is_metal = False
        
        # DETECT HALOGENS
        if atom.element.upper() in HALOGENS:
            atom.is_halogen = True
        else:
            atom.is_halogen = False
    
    # ATOM TYPING VIA OPENBABEL
    # ITERATE OVER ATOM TYPE SMARTS DEFINITIONS
    for atom_type, smartsdict in ATOM_TYPES.items():
        
        # FOR EACH ATOM TYPE SMARTS STRING
        for smarts in smartsdict.values():
            
            # GET OPENBABEL ATOM MATCHES TO THE SMARTS PATTERN
            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(str(smarts))
            
            ob_smart.Match(mol)
            matches = [x for x in ob_smart.GetMapList()]
            
            if matches:
                
                # REDUCE TO A SINGLE LIST
                matches = set(reduce(operator.add, matches))
                
                for match in matches:
                    
                    atom = mol.GetAtom(match)
                    ob_to_bio[atom.GetId()].atom_types.add(atom_type)
                    
    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in (x for x in s_atoms if 'W' in x.get_full_id()[3][0]):
        atom.atom_types.add('hbond acceptor')
        atom.atom_types.add('hbond donor')
    
    logging.info('Typed atoms.')
    
    # INITIALISE ATOM FEATURE SIFT AND
    # DETERMINE ATOM POTENTIAL FEATURE SIFT
    
    # 5: 0: HBOND
    # 6: 1: WEAK_HBOND
    # 7: 2: HALOGEN_BOND
    # 8: 3: IONIC
    # 9: 4: METAL_COMPLEX
    #10: 5: AROMATIC
    #11: 6:  HYDROPHOBIC
    #12: 7: CARBONYL
    
    for atom in s_atoms:
        atom.potential_fsift = [0] * 8
        atom.actual_fsift = [0] * 8
        atom.actual_fsift_inter_only = [0] * 8
        atom.actual_fsift_intra_only = [0] * 8
        atom.actual_fsift_water_only = [0] * 8
        
        # 0: HBOND
        if 'hbond acceptor' in atom.atom_types or 'hbond donor' in atom.atom_types:
            atom.potential_fsift[0] = 1
        
        # 1: WEAK HBOND
        if 'weak hbond acceptor' in atom.atom_types or 'weak hbond donor' in atom.atom_types or 'hbond donor' in atom.atom_types:
            atom.potential_fsift[1] = 1
            
        # 2: HALOGEN BOND
        if 'xbond acceptor' in atom.atom_types or 'xbond donor' in atom.atom_types:
            atom.potential_fsift[2] = 1
        
        # 3: IONIC
        if 'pos ionisable' in atom.atom_types or 'neg ionisable' in atom.atom_types:
            atom.potential_fsift[3] = 1
        
        # 4: METAL COMPLEX
        if 'hbond acceptor' in atom.atom_types or atom.is_metal:
            atom.potential_fsift[4] = 1
        
        # 5: AROMATIC
        if 'aromatic' in atom.atom_types:
            atom.potential_fsift[5] = 1
        
        # 6: HYDROPHOBIC
        if 'hydrophobe' in atom.atom_types:
            atom.potential_fsift[6] = 1
        
        # 7: CARBONYL
        if 'carbonyl oxygen' in atom.atom_types or 'carbonyl carbon' in atom.atom_types:
            atom.potential_fsift[7] = 1
    
    # PERCIEVE AROMATIC RINGS
    s.rings = OrderedDict()
    
    for e, ob_ring in enumerate(mol.GetSSSR()):
        
        if not ob_ring.IsAromatic():
            continue
        
        center = ob.vector3()
        normal = ob.vector3()
        normal_opp = ob.vector3()
        
        ob_ring.findCenterAndNormal(center, normal, normal_opp)
        
        # CONVERT CENTER AND NORMALS TO NUMPY
        center = np.array([center.GetX(), center.GetY(), center.GetZ()])
        normal = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])
        normal_opp = np.array([normal_opp.GetX(), normal_opp.GetY(), normal_opp.GetZ()])
        
        # STORE RING
        s.rings[e] = {'ring_id': e,
                      'center': center,
                      'normal': normal,
                      'normal_opp': normal_opp}
    
    logging.info('Percieved and stored rings.')
    
    # ADD EXPLICIT HYDROGEN COORDS FOR H-BONDING INTERACTIONS
    # ADDING HYDROGENS DOESN'T SEEM TO INTERFERE WITH ATOM SERIALS (THEY GET ADDED AS 0)
    # SO WE CAN STILL GET BACK TO THE PERSISTENT BIOPYTHON ATOMS THIS WAY.
    mol.AddHydrogens()
    
    logging.info('Added hydrogens.')
    
    # FOR PDB OUTPUT OF OB STRUCTURES
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('pdb', 'pdb')
    
    if args.minimise_hydrogens:
        # MINIMIZE HYDROGENS ONLY
        # `https://www.mail-archive.com/openbabel-discuss@lists.sourceforge.net/msg02216.html`
        # `https://www.mail-archive.com/openbabel-discuss@lists.sourceforge.net/msg02220.html`
        # `https://github.com/dlonie/OpenBabel-BFGS/blob/master/scripts/python/examples/minSimple.py`
        
        # CONSTRAIN NON-HYDROGEN ATOMS IN THE MINIMISATION
        logging.info('Beginning hydrogen minimisation.')
        
        constraints = ob.OBFFConstraints()
        
        for ob_atom in ob.OBMolAtomIter(mol):
                
            if not ob_atom.IsHydrogen():
                constraints.AddAtomConstraint(ob_atom.GetIdx())
        
        logging.info('Constrained non-hydrogen atoms.')
        
        # INITIALISE THE FORCEFIELD
        ff = ob.OBForceField.FindForceField(args.minimisation_forcefield) # MMFF94, UFF, Ghemical
        
        logging.info('Initialised forcefield.')
        
        # TODO: MAKE LOGGING A COMMAND LINE FLAG
        ff.SetLogLevel(ob.OBFF_LOGLVL_NONE) # OBFF_LOGLVL_LOW
        ff.SetLogToStdErr()
        
        if ff.Setup(mol, constraints) == 0: #, constraints)
            logging.warn('Could not setup the hydrogen minimisation forcefield. Skipping hydrogen minimisation.')
        else:
            
            # DOESN'T WORK
            #for ob_atom in ob.OBMolAtomIter(mol):
            #    
            #    if not ob_atom.IsHydrogen():
            #        ff.SetFixAtom(ob_atom.GetIdx())
            
            # TODO: MAKE THIS COMMAND-LINE CHANGABLE
            if args.minimisation_method == 'ConjugateGradients':
                ff.ConjugateGradients(args.minimisation_steps)
            elif args.minimisation_method == 'SteepestDescent':
                ff.SteepestDescent(args.minimisation_steps)
            elif args.minimisation_method == 'DistanceGeometry':
                ff.DistanceGeometry()
            
            ff.GetCoordinates(mol)
            
            logging.info('Minimised hydrogens.')
    
    # OUTPUT OB STRUCTURE WITH HYDROGENS IF REQUESTED
    if args.write_hydrogenated:
        conv.WriteFile(mol, pdb_filename.replace('.pdb', '_hydrogenated.pdb'))
        logging.info('Wrote hydrogenated PDB file.')
    
    # ADD HYDROGENS TO BIOPYTHON ATOMS
    
    # ITERATE THROUGH THE OBMOL
    for atom in ob.OBMolAtomIter(mol):
        
        # IF THE ATOM HAS EXPLICIT HYDROGENS
        if atom.ExplicitHydrogenCount() > 0:
            
            biopython_atom = ob_to_bio[atom.GetId()]

            # GET THE BONDED ATOMS OF THE OBATOM
            for atom_neighbour in ob.OBAtomAtomIter(atom):
                if atom_neighbour.IsHydrogen():

                    # APPEND THE HYDROGEN COORDINATES TO THE BIOPYTHON ATOM 'h_coords' ATTRIBUTE
                    biopython_atom.h_coords.append(np.array([atom_neighbour.GetX(), atom_neighbour.GetY(), atom_neighbour.GetZ()]))
                
    logging.info('Added hydrogens to BioPython atoms.')
    
    # CHOOSE ENTITY FOR CATMAP CALCULATION
    # AN ENTITY IS A LIST OF BIOPYTHON ATOMS
    # USING ALL THE BIOPYTHON ATOMS AS E FOR NOW
    e = list(s_atoms)
    
    # ADD VDW RADII TO ENTITY ATOMS
    # USING OPENBABEL VDW RADII
    for atom in e:
        atom.vdw_radius = ob.etab.GetVdwRad(mol.GetAtomById(bio_to_ob[atom]).GetAtomicNum())
    
    logging.info('Added VdW radii.')
    
    # ADD COVALENT RADII TO ENTITY ATOMS
    # USING OPENBABEL VDW RADII
    for atom in e:
        atom.cov_radius = ob.etab.GetCovalentRad(mol.GetAtomById(bio_to_ob[atom]).GetAtomicNum())
    
    logging.info('Added covalent radii.')
    
    # NEIGHBORSEARCH
    ns = NeighborSearch(e)
    
    logging.info('Completed NeighborSearch.')
    
    # ADD "LIGAND" SELECTION (SUBSET OF ATOMS) FROM
    # WITHIN THE ENTITY FOR CONTACT CALCULATION
    selection = e[:]
    
    if args.selection:
        selection = selection_parser(args.selection, e)
    
    elif args.selection_file:
        with open(args.selection_file, 'rb') as fo:
            selection = selection_parser([line for line in fo], e)
    
    if len(selection) == 0:
        
        logging.error('Selection was empty.')
        sys.exit()
    
    logging.info('Made selection.')
    
    selection_set = set(selection)
    
    # EXPAND THE SELECTION TO INCLUDE THE BINDING SITE
    selection_plus = set(selection)
    
    if args.selection:
        
        for atom_bgn, atom_end in ns.search_all(6.0):
            
            if atom_bgn in selection_set or atom_end in selection_set:
                selection_plus.add(atom_bgn)
                selection_plus.add(atom_end)
        
        selection_plus = list(selection_plus)
        
        logging.info('Expanded to binding site.')
        
        # NEW NEIGHBOURSEARCH
        ns = NeighborSearch(selection_plus)
        
        logging.info('Completed new NeighbourSearch.')
    
    #if not args.consider_all:
    #    # CALCULATE PER-ATOM SOLVENT ACCESSIBILITY
    #    #
    #    # Using the Shrake-Rupley Algorithm.
    #    # - http://en.wikipedia.org/wiki/Accessible_surface_area
    #    # - http://boscoh.com/protein/calculating-the-solvent-accessible-surface-area-asa.html
    #    # - http://www.xsi-blog.com/archives/115 (Golden Section Spiral)
    #    
    #    # GENERATE SPHERES FOR ATOM SOLVATION
    #    # USING THE GOLDEN SECTION SPIRAL ALGORITHM
    #    gss = points_on_sphere(NUM_SOLVENT_SPHERE_POINTS)
    #    
    #    # CONSTANT FOR CONVERTING FROM SPHERE POINTS TO AREA
    #    const = 4.0 * np.pi / len(gss)
    #    
    #    # LOOP THROUGH ATOMS OF ENTITY OF INTEREST
    #    for atom in e:
    #        atom.num_accessible_points = 0
    #        
    #        # GET NEARBY ATOMS
    #        atom.neighbours = ns.search(atom.coord, 5.0) # NEARBY ATOM DISTANCE HARD-CODED FOR NOW, SHOULDN'T NEED TO BE CHANGED
    #        
    #        # ATOM RADII FOR SASA PURPOSES
    #        radius = atom.vdw_radius + SOLVENT_RADIUS
    #        neighbour_radii = np.array([x.vdw_radius + SOLVENT_RADIUS for x in atom.neighbours])
    #        
    #        # TRANSLATE THE SOLVENT SPHERE TO SURROUND THE ATOM
    #        atom.solvent_sphere = gss[:] * radius + atom.coord
    #        
    #        # CALCULATE A DISTANCE MATRIX OF SPHERE POINTS TO NEIGHBOUR ATOMS
    #        # ESSENTIALLY A 2D ARRAY: NEIGHBOURS ALONG X-AXIS, POINTS ALONG Y-AXIS
    #        # THUS EACH ELEMENT IS AN ARRAY FOR A POINT, IN WHICH EACH ELEMENT IS A DISTANCE TO A NEIGHBOUR
    #        n_dm = cdist(atom.solvent_sphere, np.array([x.coord for x in atom.neighbours]))
    #        
    #        # SUBTRACT THE NEIGHBOUR RADII FROM THEIR DISTANCE TO EACH POINT
    #        # NEGATIVE NUMBERS IN A SUB-ARRAY INDICATE INACCESSIBILITY
    #        # THEREFORE IF ANY OF THE 'ROW' VALUES ARE NEGATIVE, THE POINT IS NOT ACCESSIBLE
    #        n_dm_sub = n_dm - neighbour_radii
    #        
    #        # USE NUMPY ANY TO REDUCE EACH POINT ROW TO ONE TRUE OR FALSE VALUE
    #        # THEN SUM THESE TO GET THE NUMBER OF INACCESSIBLE POINTS
    #        # THEN TAKE THAT SUM AWAY FROM THE TOTAL NUMBER OF POINTS TO GET THE NUMBER OF ACCESSIBLE POINTS
    #        atom.num_accessible_points = NUM_SOLVENT_SPHERE_POINTS - sum(np.any(n_dm_sub < 0.0, axis=1))
    #        
    #        # CONVERT THE NUMBER OF ACCESSIBLE POINTS TO ACCESSIBILITY IN A^2
    #        atom.area = const * atom.num_accessible_points * radius * radius
    #        
    #        # BOOLEANISE THE ACCESSBILITY
    #        atom.is_accessible = atom.area > SASA_THRESHOLD
    #    
    #    logging.info('Calculated per-atom SASA.')
    
    # CALCULATE PAIRWISE CONTACTS
    with open(pdb_filename.replace('.pdb', '.contacts'), 'wb') as fo:
        for atom_bgn, atom_end in ns.search_all(INTERACTING_THRESHOLD):
           
            # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/contacts.py?at=default`
            # SIFT:
            # 0: CLASH
            # 1: COVALENT
            # 2: VDW_CLASH
            # 3: VDW
            # 4: PROXIMAL
            # 5: HBOND
            # 6: WEAK_HBOND
            # 7: HALOGEN_BOND
            # 8: IONIC
            # 9: METAL_COMPLEX
            #10: AROMATIC
            #11: HYDROPHOBIC
            #12: CARBONYL
            
            # IGNORE CONTACTS THAT EITHER:
            # - DON'T INVOLVE THE SELECTION
            # - DON'T INVOLVE WATER
            #if not (atom_bgn in selection_set or atom_end in selection_set):
            #    if not ('W' in atom_bgn.get_full_id()[3][0] or 'W' in atom_end.get_full_id()[3][0]):
            #        continue
            
            # DETERMINE CONTACT TYPE
            contact_type = ''
            
            if not atom_bgn in selection_set and not atom_end in selection_set:
                contact_type = 'INTRA_NON_SELECTION'
                
            if atom_bgn in selection_set and atom_end in selection_set:
                contact_type = 'INTRA_SELECTION'
                
            if (atom_bgn in selection_set and not atom_end in selection_set) or (atom_end in selection_set and not atom_bgn in selection_set):
                contact_type = 'INTER'
                
            if (atom_bgn in selection_set and 'W' in atom_end.get_full_id()[3][0]) or (atom_end in selection_set and 'W' in atom_bgn.get_full_id()[3][0]):
                contact_type = 'SELECTION_WATER'
                
            if (atom_bgn not in selection_set and 'W' in atom_end.get_full_id()[3][0]) or (atom_end not in selection_set and 'W' in atom_bgn.get_full_id()[3][0]):
                contact_type = 'NON_SELECTION_WATER'
                
            if 'W' in atom_bgn.get_full_id()[3][0] and 'W' in atom_end.get_full_id()[3][0]:
                contact_type = 'WATER_WATER'
                
            if not contact_type:
                logging.error('Could not assign a contact type for {}:{}'.format(atom_bgn, atom_end))
                raise
            
            sum_cov_radii = atom_bgn.cov_radius + atom_end.cov_radius
            sum_vdw_radii = atom_bgn.vdw_radius + atom_end.vdw_radius
            
            ob_atom_bgn = mol.GetAtomById(bio_to_ob[atom_bgn])
            ob_atom_end = mol.GetAtomById(bio_to_ob[atom_end])
            
            SIFt = [0] * 13
            
            # IGNORE INTRA-RESIDUE CONTACTS
            if atom_bgn.get_parent() == atom_end.get_parent():
                continue
            
            #print contact_type
            
            distance = np.linalg.norm(atom_bgn.coord - atom_end.coord)
            
            # COVALENT
            is_covalent = False
            
            for cov_bonded in ob.OBAtomAtomIter(ob_atom_bgn):
                
                if cov_bonded.GetId() == ob_atom_end.GetId():
                    is_covalent = True
                    break
            
            if is_covalent:
                SIFt[1] = 1
            
            # CLASHES
            elif distance < sum_cov_radii:
                SIFt[0] = 1
                
            # VDW CLASH
            elif distance < sum_vdw_radii:
                SIFt[2] = 1
            
            # VDW
            elif distance <= sum_vdw_radii + VDW_COMP_FACTOR:
                SIFt[3] = 1
                
            # PROXIMAL
            else:
                SIFt[4] = 1
            
            # FEATURE CONTACTS
            if not any(SIFt[:2]) and distance < CONTACT_TYPES_DIST_MAX:
                
                # HBOND
                
                # NO NEED TO USE HYDROGENS FOR WATERS
                if 'W' in atom_bgn.get_full_id()[3][0] and distance <= (sum_vdw_radii + VDW_COMP_FACTOR):
                    if 'hbond acceptor' in atom_end.atom_types or 'hbond donor' in atom_end.atom_types:
                        SIFt[5] = 1
                
                elif 'W' in atom_end.get_full_id()[3][0] and distance <= (sum_vdw_radii + VDW_COMP_FACTOR):
                    if 'hbond acceptor' in atom_bgn.atom_types or 'hbond donor' in atom_bgn.atom_types:
                        SIFt[5] = 1
                
                else:
                    
                    # ATOM_BGN IS DONOR
                    if 'hbond donor' in atom_bgn.atom_types and 'hbond acceptor' in atom_end.atom_types:
                        
                        SIFt[5] = is_hbond(atom_bgn, atom_end)
                    
                    # ATOM_END IS DONOR
                    elif 'hbond donor' in atom_end.atom_types and 'hbond acceptor' in atom_bgn.atom_types:
                        
                        SIFt[5] = is_hbond(atom_bgn, atom_end)
                    
                # WEAK HBOND
                
                # ATOM_BGN IS ACCEPTOR, ATOM_END IS CARBON
                if 'hbond acceptor' in atom_bgn.atom_types and 'weak hbond donor' in atom_end.atom_types:
                    SIFt[6] = is_weak_hbond(atom_end, atom_bgn)
                
                # ATOM_BGN IS CARBON, ATOM_END IS ACCEPTOR
                elif 'weak hbond donor' in atom_bgn.atom_types and 'hbond acceptor' in atom_end.atom_types:
                    SIFt[6] = is_weak_hbond(atom_bgn, atom_end)
                    
                # ATOM_BGN IS HALOGEN WEAK ACCEPTOR
                elif 'weak hbond acceptor' in atom_bgn.atom_types and atom_bgn.is_halogen and ('hbond donor' in atom_end.atom_types or 'weak hbond donor' in atom_end.atom_types):
                    SIFt[6] = is_halogen_weak_hbond(atom_end, atom_bgn, mol)
                
                # ATOM END IS HALOGEN WEAK ACCEPTOR
                elif 'weak hbond acceptor' in atom_end.atom_types and atom_end.is_halogen and ('hbond donor' in atom_bgn.atom_types or 'weak hbond donor' in atom_bgn.atom_types):
                    SIFt[6] = is_halogen_weak_hbond(atom_bgn, atom_end, mol)
                    
                # XBOND
                if distance <= sum_vdw_radii + VDW_COMP_FACTOR:
                    
                    if 'xbond donor' in atom_bgn.atom_types and 'xbond acceptor' in atom_end.atom_types:
                        SIFt[7] = is_xbond(atom_bgn, atom_end, mol)
                    
                    elif 'xbond donor' in atom_end.atom_types and 'xbond acceptor' in atom_bgn.atom_types:
                        SIFt[7] = is_xbond(atom_end, atom_bgn, mol)
                
                # IONIC
                if distance <= CONTACT_TYPES['ionic']['distance']:
                    
                    if 'pos ionisable' in atom_bgn.atom_types and 'neg ionisable' in atom_end.atom_types:
                        SIFt[8] = 1
                    
                    elif 'neg ionisable' in atom_bgn.atom_types and 'pos ionisable' in atom_end.atom_types:
                        SIFt[8] = 1
                    
                # METAL COMPLEX
                if distance <= CONTACT_TYPES['metal']['distance']:
                    
                    if 'hbond acceptor' in atom_bgn.atom_types and atom_end.is_metal:
                        SIFt[9] = 1
                        
                    elif 'hbond acceptor' in atom_end.atom_types and atom_bgn.is_metal:
                        SIFt[9] = 1
                
                # CARBONYL
                if distance <= CONTACT_TYPES['carbonyl']['distance']:
                    
                    if 'carbonyl oxygen' in atom_bgn.atom_types and 'carbonyl carbon' in atom_end.atom_types:
                        SIFt[12] = 1
                        
                    elif 'carbonyl oxygen' in atom_end.atom_types and 'carbonyl carbon' in atom_bgn.atom_types:
                        SIFt[12] = 1
                        
                # AROMATIC
                if 'aromatic' in atom_bgn.atom_types and 'aromatic' in atom_end.atom_types and distance <= CONTACT_TYPES['aromatic']['distance']:
                    SIFt[10] = 1
                
                # HYDROPHOBIC
                if 'hydrophobe' in atom_bgn.atom_types and 'hydrophobe' in atom_end.atom_types and distance <= CONTACT_TYPES['hydrophobic']['distance']:
                    SIFt[11] = 1
            
            # UPDATE ATOM FEATURE SIFts
            fsift = SIFt[5:]
            update_atom_fsift(atom_bgn, fsift, contact_type)
            update_atom_fsift(atom_end, fsift, contact_type)
            
            # WRITE OUT CONTACT SIFT TO FILE
            fo.write('{}\n'.format('\t'.join([str(x) for x in [make_pymol_string(atom_bgn), make_pymol_string(atom_end)] + SIFt])))
        
        logging.info('Calculated pairwise contacts.')
    
    # WRITE OUT SIFT MATCHING
    # LIGAND AND BINDING SITE (`selection_plus`)
    with open(pdb_filename.replace('.pdb', '.siftmatch'), 'wb') as fo, open(pdb_filename.replace('.pdb', '.specific.siftmatch'), 'wb') as specific_fo:
        for atom in selection_plus:
            
            sift_match = sift_match_base3(atom.potential_fsift, atom.actual_fsift) # WHICH SIFT TO USE?
            
            sift_match_inter = sift_match_base3(atom.potential_fsift, atom.actual_fsift_inter_only)
            sift_match_intra = sift_match_base3(atom.potential_fsift, atom.actual_fsift_intra_only)
            sift_match_water = sift_match_base3(atom.potential_fsift, atom.actual_fsift_water_only)
            
            human_readable = human_sift_match(sift_match)
            
            # SUBJECT TO CHANGE
            fo.write('{}\n'.format('\t'.join([str(x) for x in [make_pymol_string(atom)] + sift_match + [int3(sift_match)] + [human_readable]])))
            
            specific_fo.write('{}\n'.format('\t'.join([str(x) for x in [make_pymol_string(atom)] +
                                                       sift_match_inter +
                                                       [human_sift_match(sift_match_inter)] +
                                                       sift_match_intra +
                                                       [human_sift_match(sift_match_intra)] +
                                                       sift_match_water +
                                                       [human_sift_match(sift_match_water)]])))
    
    # RING-RING INTERACTIONS
    # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`
    # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/sql/populate.sql?at=default`
    # `http://marid.bioc.cam.ac.uk/credo/about`
    with open(pdb_filename.replace('.pdb', '.ri'), 'wb') as fo:
        for ring in s.rings:
            
            ring_key = ring
            ring = s.rings[ring]
            
            for ring2 in s.rings:
                
                ring_key2 = ring2
                ring2 = s.rings[ring2]
                
                # NO SELFIES
                if ring_key == ring_key2:
                    continue
                
                distance = np.linalg.norm(ring['center'] - ring2['center'])
                
                if distance > CONTACT_TYPES['aromatic']['centroid_distance']:
                    continue
                
                theta_point = ring['center'] - ring2['center']
                #iota_point = ring2['center'] - ring['center']
                
                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                dihedral = abs(ring_ring_angle(ring, ring2, True, True))
                theta = abs(ring_angle(ring, theta_point, True, True))
                
                #logging.info('Dihedral = {}     Theta = {}'.format(dihedral, theta))
                
                int_type = ''
                
                if dihedral <= 30.0 and theta <= 30.0 :
                    int_type = 'FF'
                if dihedral <= 30.0 and theta <= 60.0 :
                    int_type = 'OF'
                if dihedral <= 30.0 and theta <= 90.0 :
                    int_type = 'EE'

                if dihedral > 30.0 and dihedral <= 60.0 and theta <= 30.0:
                    int_type = 'FT'
                if dihedral > 30.0 and dihedral <= 60.0 and theta <= 60.0 :
                    int_type = 'OT'
                if dihedral > 30.0 and dihedral <= 60.0 and theta <= 90.0 :
                    int_type = 'ET'

                if dihedral > 60.0 and dihedral <= 90.0 and theta <= 30.0 :
                    int_type = 'FE'
                if dihedral > 60.0 and dihedral <= 90.0 and theta <= 60.0 :
                    int_type = 'OE'
                if dihedral > 60.0 and dihedral <= 90.0 and theta <= 90.0 :
                    int_type = 'EF'
                    
                # WRITE RING INTERACTION TO FILE
                output = [
                    ring['ring_id'],
                    list(ring['center']),
                    ring2['ring_id'],
                    list(ring2['center']),
                    int_type
                ]
                
                fo.write('{}\n'.format('\t'.join([str(x) for x in output])))
        
    # RINGS AND ATOM-RING INTERACTIONS
    with open(pdb_filename.replace('.pdb', '.ari'), 'wb') as fo, open(pdb_filename.replace('.pdb', '.rings'), 'wb') as ring_fo:
        for ring in s.rings:
            
            ring_key = ring
            ring = s.rings[ring]
            
            for atom in ns.search(ring['center'], CONTACT_TYPES['aromatic']['atom_aromatic_distance']):
                
                distance = np.linalg.norm(atom.coord - ring['center'])
                
                if distance > CONTACT_TYPES['aromatic']['atom_aromatic_distance'] or 'aromatic' in atom.atom_types:
                    continue
                
                potential_interactions = set([])
                
                if atom.element == 'C' and 'weak hbond donor' in atom.atom_types:
                    potential_interactions.add('CARBONPI')
                
                if 'pos ionisable' in atom.atom_types:
                    potential_interactions.add('CATIONPI')
                    
                if 'hbond donor' in atom.atom_types:
                    potential_interactions.add('DONORPI')
                
                if 'xbond donor' in atom.atom_types:
                    potential_interactions.add('HALOGENPI')
                
                if not potential_interactions:
                    continue
                
                # N.B.: NOT SURE WHY ADRIAN WAS USING SIGNED, BUT IT SEEMS
                #       THAT TO FIT THE CRITERIA FOR EACH TYPE OF INTERACTION
                #       BELOW, SHOULD BE UNSIGNED, I.E. `abs()`
                theta = abs(ring_angle(ring, atom.coord, True, True))
                
                if theta <= 30.0:
                    
                    #logging.info('Atom: <{}>     Theta = {}'.format(atom.get_full_id(), theta))
                    
                    # WRITE ATOM-RING INTERACTION TO FILE
                    output = [
                        make_pymol_string(atom),
                        ring['ring_id'],
                        list(ring['center']),
                        sorted(list(potential_interactions))
                    ]
                    
                    fo.write('{}\n'.format('\t'.join([str(x) for x in output])))
            
            # WRITE RING OUT TO RING FILE
            output = [
                ring['ring_id'],
                list(ring['center'])
            ]
            
            ring_fo.write('{}\n'.format('\t'.join([str(x) for x in output])))
    
    logging.info('Program End. Maximum memory usage was {}.'.format(max_mem_usage()))
