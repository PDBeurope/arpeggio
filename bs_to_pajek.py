#!/usr/bin/python

###########
# IMPORTS #
###########

import argparse
import collections
import logging
import operator
import os
import sys

from collections import OrderedDict

import igraph

#############
# FUNCTIONS #
#############

def sift_to_bool(x):
    
    return bool(int(x))

def add_vertex_without_duplication(g, name):
    '''
    By default, igraph allows you to add nodes with a name that already exists
    in the graph.
    
    This function adds a node with a name to a graph, but ensures that the name
    isn't already given to an existing node.
    '''
    
    try:
        g.vs.find(name=name)
    except (ValueError, KeyError):
        g.add_vertex(name=name)

########
# MAIN #
########

if __name__ == '__main__':
    
    # ARGUMENT PARSING
    parser = argparse.ArgumentParser(description='''

#############################
# ARPEGGIO TO GRAPH NETWORK #
#############################

A program for taking binding site contact output from
Arpeggio, and converting to Pajek graph network format.

Dependencies:
- Python (v2.7)
- python-igraph

''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('pdb', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-ai', '--allow-intra', action='store_true', help='Allow intra contacts (selection and non selection) to be output.')
    parser.add_argument('-ip', '--include-proximals', action='store_true', help='Include proximal contacts.')
    parser.add_argument('-icl', '--include-covalent-clashes', action='store_true', help='Include covalent clashes.')
    parser.add_argument('-ico', '--include-covalents', action='store_true', help='Include covalent contacts.')
    parser.add_argument('-io', '--include-others', action='store_true', help='Include contacts which aren\'t polar or apolar, e.g. atom-atom aromatic.')
    parser.add_argument('-op', type=str, default='', help='The output postfix you used with Arpeggio (if any).')
    parser.add_argument('-op2', type=str, default='', help='An additional output postfix for your .net files, if you want.')
    
    args = parser.parse_args()
    
    pdb_filename = args.pdb
    
    # ALLOWABLE TYPES
    allowable_types = set(['INTER'])
    
    if args.allow_intra:
        allowable_types.add('INTRA_NON_SELECTION')
        allowable_types.add('INTRA_SELECTION')
    
    #allowable_types.add('SELECTION_WATER')
    #allowable_types.add('NON_SELECTION_WATER')
    #allowable_types.add('WATER_WATER')
    
    # CHANGE TO AN ABSOLUTE PATH
    pdb_filename = os.path.abspath(pdb_filename)
    
    output_postfix = args.op
    output_postfix2 = args.op2
    
    contacts_extension = '.bs_contacts'
    
    #if args.use_all_binding_site:
    #    contacts_extension = '.bs_contacts'
    
    contacts_filename = pdb_filename.replace('.pdb', output_postfix + contacts_extension)
    
    # MAKE DATA STRUCTURE FOR PAJEK FORMAT
    contacts = igraph.Graph()
    polar_contacts = igraph.Graph()
    apolar_contacts = igraph.Graph()
    
    # GET CONTACTS
    with open(contacts_filename, 'rb') as fo:
        
        for line in fo:
            
            line = line.strip().split('\t')
            label = []
            
            contact_type = line[17]
            
            if contact_type not in allowable_types:
                continue
            
            node_bgn, node_end = line[0], line[1]
            
            # GET CONTACT SIFT
            clash = sift_to_bool(line[2])
            covalent = sift_to_bool(line[3])
            vdw_clash = sift_to_bool(line[4])
            vdw = sift_to_bool(line[5])
            proximal = sift_to_bool(line[6])
            hbond = sift_to_bool(line[7])
            weak_hbond = sift_to_bool(line[8])
            halogen_bond = sift_to_bool(line[9])
            ionic = sift_to_bool(line[10])
            metal_complex = sift_to_bool(line[11])
            aromatic = sift_to_bool(line[12])
            hydrophobic = sift_to_bool(line[13])
            carbonyl = sift_to_bool(line[14])
            polar = sift_to_bool(line[15])
            weak_polar = sift_to_bool(line[16])
            
            # AVOID CLASHES, COVALENTS AND PROXIMALS UNLESS THEY'RE ASKED FOR
            if not args.include_covalents:
                if covalent:
                    continue
                
            if not args.include_covalent_clashes:
                if clash:
                    continue
            
            if not args.include_proximals:
                if proximal:
                    continue
            
            # POLAR
            # USING HBOND, WEAK HBOND, HALOGEN BOND, IONIC, CARBONYL
            # METAL COMPLEX, POLAR, WEAK POLAR
            polars = [hbond, weak_hbond, halogen_bond, ionic,
                      metal_complex, carbonyl, polar, weak_polar]
            
            if any(polars):
                add_vertex_without_duplication(polar_contacts, name=node_bgn)
                add_vertex_without_duplication(polar_contacts, name=node_end)
                polar_contacts.add_edge(node_bgn, node_end) #, label='polar')
                
                label.append('polar')
            
            # HYDROPHOBIC
            if hydrophobic:
                add_vertex_without_duplication(apolar_contacts, name=node_bgn)
                add_vertex_without_duplication(apolar_contacts, name=node_end)
                apolar_contacts.add_edge(node_bgn, node_end) #, label='apolar')
                
                label.append('apolar')
            
            # DON'T OUTPUT NON-POLAR OR HYDROPHOBIC CONTACTS
            if not args.include_others:
                if not label:
                    continue
            
            # ALL WITH LABELS
            label_string = '_'.join(label)
            
            add_vertex_without_duplication(contacts, name=node_bgn)
            add_vertex_without_duplication(contacts, name=node_end)
            contacts.add_edge(node_bgn, node_end, label=label_string)
            
    contacts.write(contacts_filename.replace(contacts_extension, '_all{}.net'.format(output_postfix2)), format='pajek')
    contacts.write(contacts_filename.replace(contacts_extension, '_all{}.gml'.format(output_postfix2)), format='gml')
    
    polar_contacts.write(contacts_filename.replace(contacts_extension, '_polar{}.net'.format(output_postfix2)), format='pajek')
    polar_contacts.write(contacts_filename.replace(contacts_extension, '_polar{}.gml'.format(output_postfix2)), format='gml')
    
    apolar_contacts.write(contacts_filename.replace(contacts_extension, '_apolar{}.net'.format(output_postfix2)), format='pajek')
    apolar_contacts.write(contacts_filename.replace(contacts_extension, '_apolar{}.gml'.format(output_postfix2)), format='gml')
    
    