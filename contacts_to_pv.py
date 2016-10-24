import argparse
import os
import re
import sys

from show_contacts import pymol_config

# do('set dash_radius, 0.15, *PI')
# do('set dash_gap, 0.1, *PI')

# do('color white, CARBONPI')
# do('color blue, DONORPI')
# do('color green, HALOGENPI')
# do('color red, CATIONPI')
# do('color yellow, METSULPHURPI')

# do('set dash_radius, 0.25, amide-ring')
# do('set dash_gap, 0.2, amide-ring')
# do('set dash_length, 0.5, amide-ring')
# do('color white, amide-ring')

# do('set dash_radius, 0.25, amide-amide')
# do('set dash_gap, 0.2, amide-amide')
# do('set dash_length, 0.5, amide-amide')
# do('color blue, amide-amide')

# do('set dash_radius, 0.25, ring-*')
# do('set dash_gap, 0.1, ring-*')
# do('color white, ring-*')

group_pymol_config = {

    'dashcolor': {
        'CARBONPI': 'white',
        'DONORPI': 'blue',
        'HALOGENPI': 'green',
        'CATIONPI': 'red',
        'METSULPHURPI': 'yellow',
        'amide-ring': 'white',
        'amide-amide': 'blue',
        'ring-FF': 'white',
        'ring-OF': 'white',
        'ring-EE': 'white',
        'ring-FT': 'white',
        'ring-OT': 'white',
        'ring-ET': 'white',
        'ring-FE': 'white',
        'ring-OE': 'white',
        'ring-EF': 'white',
    },

    'dashradius': {
        'CARBONPI': 0.15,
        'DONORPI': 0.15,
        'HALOGENPI': 0.15,
        'CATIONPI': 0.15,
        'METSULPHURPI': 0.15,
        'amide-ring': 0.25,
        'amide-amide': 0.25,
        'ring-FF': 0.25,
        'ring-OF': 0.25,
        'ring-EE': 0.25,
        'ring-FT': 0.25,
        'ring-OT': 0.25,
        'ring-ET': 0.25,
        'ring-FE': 0.25,
        'ring-OE': 0.25,
        'ring-EF': 0.25,
    },

    'dashgap': {
        'CARBONPI': 0.1,
        'DONORPI': 0.1,
        'HALOGENPI': 0.1,
        'CATIONPI': 0.1,
        'METSULPHURPI': 0.1,
        'amide-ring': 0.2,
        'amide-amide': 0.2,
        'ring-FF': 0.1,
        'ring-OF': 0.1,
        'ring-EE': 0.1,
        'ring-FT': 0.1,
        'ring-OT': 0.1,
        'ring-ET': 0.1,
        'ring-FE': 0.1,
        'ring-OE': 0.1,
        'ring-EF': 0.1,
    },

    'dashlength': {
        'CARBONPI': 0.5,
        'DONORPI': 0.5,
        'HALOGENPI': 0.5,
        'CATIONPI': 0.5,
        'METSULPHURPI': 0.5,
        'amide-ring': 0.5,
        'amide-amide': 0.5,
        'ring-FF': 0.5,
        'ring-OF': 0.5,
        'ring-EE': 0.5,
        'ring-FT': 0.5,
        'ring-OT': 0.5,
        'ring-ET': 0.5,
        'ring-FE': 0.5,
        'ring-OE': 0.5,
        'ring-EF': 0.5,
    },

}

# MAIN
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''
    
    ############
    # ARPEGGIO #
    ############
    
    Interaction Viewer output to PV Viewer Format.
    
    A program for calculating interatomic interactions,
    using only Open Source dependencies.
    
    Dependencies:
    - Python (v2.7)
    - PV Viewer ()
    
    **You must have already run your structure with Arpeggio to generate the required output files,
      and also run show_contacts.py and made a .pml file.**.
    Be careful about absolute and relative paths.
    
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('pdb', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-op', type=str, default='', help='The output postfix you used with Arpeggio (if any).')

    args = parser.parse_args()
    
    pdb_filename = args.pdb
    
    # CHANGE TO AN ABSOLUTE PATH
    pdb_filename = os.path.abspath(pdb_filename)
    short_pdb_filename = os.path.split(pdb_filename)[1]
    
    output_postfix = args.op

    script_filename = pdb_filename.replace('.pdb', output_postfix + '.pml')
    js_filename = script_filename.replace('.pml', '.js')

    # SANITISE NAME FOR JS FUNCTION NAME
    sanitised_name = os.path.split(pdb_filename)[1].replace('.pdb', output_postfix)
    sanitised_name = re.sub('[^a-zA-Z0-9]', '', sanitised_name)
    function_name = 'load' + sanitised_name

    # STORE PSUEDOATOM POINTS
    pseudoatoms = {}

    js = []

    # DEFINE LOAD FUNCTION
    js.append('function {}() {{'.format(function_name))
    js.append("pv.io.fetchPdb('{}', function(structure) {{".format(short_pdb_filename))

    js.append('''

        viewer.sline('protein', structure, {
            color: color.byChain()
        });
        viewer.ballsAndSticks('ligands', structure.select('ligand'));
        viewer.centerOn(structure);

''')

    with open(script_filename, 'rb') as fo:

        for line in fo:

            if line.startswith('pseudoatom'):

                line = line.strip()

                #pseudoatom ring_centers, pos=[0.31440000000000001, 7.6326000000000009, -10.023800000000001]
                pseudoatom_name = line.split()[1].rstrip(',')
                pseudoatom_pos = line.split('pos=')[1]

                pseudoatoms[pseudoatom_name] = pseudoatom_pos

                continue

            if line.startswith('distance'):

                # PSEUDOATOM-PSEUDOATOM CONTACTS
                m = re.search(r'^distance ([a-zA-Z_-]+), (pt1), (pt2)', line)

                if m:
                    dist_name = m.group(1)
                    dist_pos1 = pseudoatoms['pt1']
                    dist_pos2 = pseudoatoms['pt2']

                    js.append('''
drawContact(viewer, '{}',
{},
{},
'{}',
{},
{});
'''.format(dist_name,
           dist_pos1,
           dist_pos2,
           group_pymol_config['dashcolor'][dist_name],
           group_pymol_config['dashradius'][dist_name],
           group_pymol_config['dashgap'][dist_name] / 3.0,
           group_pymol_config['dashlength'][dist_name]))

                    continue

                # PSEUDOATOM-ATOM CONTACTS
                m = re.search(r'^distance ([a-zA-Z_-]+), (pt1), ([A-Za-z0-9]{1})/([0-9]+)([A-Za-z]*)/([A-Za-z0-9]+)', line)

                if m:
                    dist_name = m.group(1)
                    dist_pos = pseudoatoms['pt1']
                    chain = m.group(3)
                    res = m.group(4)
                    ins = m.group(5)
                    atm = m.group(6)

                    js.append('''
drawAtomPositionContactFromPredicate(structure, viewer, '{}',
{{chain: '{}', rnum: {}, aname: '{}'}},
{},
'{}',
{},
{});
'''.format(dist_name,
           chain,
           res,
           atm,
           dist_pos,
           group_pymol_config['dashcolor'][dist_name],
           group_pymol_config['dashradius'][dist_name],
           group_pymol_config['dashgap'][dist_name] / 3.0,
           group_pymol_config['dashlength'][dist_name]))

                    continue

                # NORMAL CONTACTS
                m = re.search(r'^distance ([a-zA-Z_-]+), ([A-Za-z0-9]{1})/([0-9]+)([A-Za-z]*)/([A-Za-z0-9]+), ([A-Za-z0-9]{1})/([0-9]+)([A-Za-z]*)/([A-Za-z0-9]+)', line)

                if not m:
                    continue

                dist_name = m.group(1)
                dist_feat, dist_dist = dist_name.split('-')
                
                if dist_feat == 'undefined' and dist_dist == 'proximal':
                    continue

                chain_bgn = m.group(2)
                res_bgn = m.group(3)
                ins_bgn = m.group(4)
                atm_bgn = m.group(5)

                chain_end = m.group(6)
                res_end = m.group(7)
                ins_end = m.group(8)
                atm_end = m.group(9)

                # SHOW STICKS FOR BINDING RESIDUES
                js.append("viewer.ballsAndSticks('binding_site', structure.select({{'cname': '{}', 'rnum': {}}}));".format(chain_bgn, res_bgn))
                js.append("viewer.ballsAndSticks('binding_site', structure.select({{'cname': '{}', 'rnum': {}}}));".format(chain_end, res_end))

                # DRAW CONTACT
                js.append('''
drawAtomAtomContactFromPredicates(structure, viewer, '{}',
{{chain: '{}', rnum: {}, aname: '{}'}},
{{cname: '{}', rnum: {}, aname: '{}'}},
'{}',
{},
{});
'''.format(dist_name,
           chain_bgn,
           res_bgn,
           atm_bgn,
           chain_end,
           res_end,
           atm_end,
           pymol_config['dashcolor'][dist_feat][dist_dist],
           pymol_config['dashradius'][dist_feat][dist_dist],
           pymol_config['dashgap'][dist_feat][dist_dist] / 3.0,
           pymol_config['dashlength'][dist_feat][dist_dist]))

                # TODO: ADD FOR RING INTERACTIONS (APPEND TO DICTIONARY
                #       AT THE TOP OF THIS SCRIPT)
                # except:
                #     pass


    # CLOSE LOAD FUNCTION
    js.append('});')
    js.append('}')

    # LOAD STRUCTURE WHEN DOM READY
    js.append('''document.addEventListener('DOMContentLoaded', {});'''.format(function_name))

    with open(js_filename, 'wb') as fo:
        fo.write('\n'.join(js))