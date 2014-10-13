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

###########
# CLASSES #
###########

class Pajek(object):
    '''
    '''
    
    def __init__(self):
        
        self._nodes = OrderedDict()
        self._edges = OrderedDict()
    
    def add_node(self, node_name):
        
        if node_name not in self._nodes:
            
            # ADD NODE NAME AS KEY WITH NODE NUMBER AS VALUE
            self._nodes[node_name] = len(self._nodes) + 1
    
    def add_edge(self, node_bgn_name, node_end_name, weight=""):
        
        # ADDS NODES IF THEY DON'T EXIST
        self.add_node(node_bgn_name)
        self.add_node(node_end_name)
        
        if (node_bgn_name, node_end_name) not in self._edges:
            
            # KEYED WITH NODES, VALUE IS WEIGHT
            self._edges[(node_bgn_name, node_end_name)] = weight

    def write_pajek_string(self):
        
        output_list = []
        
        # NODES
        output_list.append('*Vertices {}'.format(len(self._nodes)))
        
        for node_name, node_num in self._nodes.iteritems():
            output_list.append('{} "{}"'.format(node_num, node_name))
        
        # EDGES
        output_list.append('*Edges')
        
        for edge in self._edges:
            
            node_bgn_num = self._nodes[edge[0]]
            node_end_num = self._nodes[edge[1]]
            
            weight = ''
            
            if self._edges[edge] != '':
                weight = ' {}'.format(weight)
            
            
            output_list.append('{} {}{}'.format(node_bgn_num, node_end_num, weight))
        
        return '\n'.join(output_list)
    
    def write_pajek_file(self, filename):
        
        with open(filename, 'wb') as fo:
            fo.write(self.write_pajek_string())

#############
# FUNCTIONS #
#############

def sift_to_bool(x):
    
    return bool(int(x))

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

''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('pdb', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-ai', '--allow-intra', action='store_true', help='Allow intra contacts (selection and non selection) to be output.')
    parser.add_argument('-ip', '--include-proximals', action='store_true', help='Include proximal contacts.')
    parser.add_argument('-icl', '--include-covalent-clashes', action='store_true', help='Include covalent clashes.')
    parser.add_argument('-ico', '--include-covalents', action='store_true', help='Include covalent contacts.')
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
    all_contacts = Pajek()
    polar_contacts = Pajek()
    hydrophobic_contacts = Pajek()
    
    # GET CONTACTS
    with open(contacts_filename, 'rb') as fo:
        
        for line in fo:
            
            line = line.strip().split('\t')
            
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
            
            all_contacts.add_edge(node_bgn, node_end)
            
            # POLAR
            # USING HBOND, WEAK HBOND, HALOGEN BOND, IONIC, CARBONYL
            # METAL COMPLEX, POLAR, WEAK POLAR
            polars = [hbond, weak_hbond, halogen_bond, ionic,
                      metal_complex, carbonyl, polar, weak_polar]
            
            if any(polars):
                polar_contacts.add_edge(node_bgn, node_end)
            
            # HYDROPHOBIC
            if hydrophobic:
                hydrophobic_contacts.add_edge(node_bgn, node_end)
            
    all_contacts.write_pajek_file(contacts_filename.replace(contacts_extension, '_all{}.net'.format(output_postfix2)))
    polar_contacts.write_pajek_file(contacts_filename.replace(contacts_extension, '_polar{}.net'.format(output_postfix2)))
    hydrophobic_contacts.write_pajek_file(contacts_filename.replace(contacts_extension, '_hydrophobic{}.net'.format(output_postfix2)))
    
    