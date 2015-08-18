import argparse
import os
import re
import sys

from show_contacts import pymol_config



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

            if line.startswith('distance'):

                # try:
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