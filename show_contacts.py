import argparse
import os
import sys
import xmlrpclib

host='localhost'
port=9123

# WITH HELP FROM
# `https://bitbucket.org/blundell/credoscript/src/0975daeea0cbd97a298daf471550794941b5aed0/credoscript/pymol.py?at=default`
# `https://bitbucket.org/blundell/credoscript/src/0975daeea0cbd97a298daf471550794941b5aed0/credoscript/support/pymolviewer.py?at=default`
# `https://bitbucket.org/blundell/credoscript/src/0975daeea0cbd97a298daf471550794941b5aed0/credoscript/config-default.json?at=default`

pymol_config = {
    
        "dashcolor":
        {
            "normal": "cyan",
            "disordered_region": "deepsalmon",

            "xbond":
            {
                "vdwclash": "blue",
                "vdw": "tv_blue",
                "proximal": "slate"
            },

            "hbond":
            {
                "vdwclash": "red",
                "vdw": "deepsalmon",
                "proximal": "salmon"
            },
            
            "polar":
            {
                "vdwclash": "red",
                "vdw": "deepsalmon",
                "proximal": "salmon"
            },

            "weakhbond":
            {
                "vdwclash": "orange",
                "vdw": "amber",
                "proximal": "peach"
            },
            
            "weakpolar":
            {
                "vdwclash": "orange",
                "vdw": "amber",
                "proximal": "peach"
            },

            "metalcomplex":
            {
                "covalent": "purple",
                "vdwclash": "purple",
                "vdw": "medium_purple"
            },

            "ionic":
            {
                "vdwclash": "yellow",
                "vdw": "maize",
                "proximal": "cream"
            },

            "aromatic":
            {
                "vdwclash": "cyan",
                "vdw": "electric_blue",
                "proximal": "baby_blue"
            },

            "hydrophobic":
            {
                "vdwclash": "green",
                "vdw": "emerald",
                "proximal": "moss_green"
            },

            "carbonyl":
            {
                "vdwclash": "rose",
                "vdw": "brilliant_rose",
                "proximal": "thulian_pink"
            },

            "undefined":
            {
                "covalent": "white",
                "vdwclash": "grey90",
                "vdw": "grey60",
                "proximal": "grey30"
            }
        },

        "dashradius":
        {
            "normal": 0.02,

            "xbond":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "hbond":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "weakhbond":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },
            
            "polar":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "weakpolar":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "metalcomplex":
            {
                "covalent": 0.16,
                "vdwclash": 0.12,
                "vdw": 0.08
            },

            "ionic":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "aromatic":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "hydrophobic":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "carbonyl":
            {
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            },

            "undefined":
            {
                "covalent": 0.16,
                "vdwclash": 0.12,
                "vdw": 0.08,
                "proximal": 0.04
            }
        },

        "dashgap":
        {
            "normal": 0.0,

            "xbond":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "hbond":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "weakhbond":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },
            
            "polar":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "weakpolar":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "metalcomplex":
            {
                "covalent": 0.20,
                "vdwclash": 0.25,
                "vdw": 0.35
            },

            "ionic":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "aromatic":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "hydrophobic":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "carbonyl":
            {
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            },

            "undefined":
            {
                "covalent": 0.0,
                "vdwclash": 0.25,
                "vdw": 0.35,
                "proximal": 0.45
            }
        },

        "dashlength":
        {
            "normal": 1.0,

            "xbond":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "hbond":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "weakhbond":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },
            
            "polar":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "weakpolar":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "metalcomplex":
            {
                "covalent": 0.10,
                "vdwclash": 0.08,
                "vdw": 0.06
            },

            "ionic":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "aromatic":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "hydrophobic":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "carbonyl":
            {
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            },

            "undefined":
            {
                "covalent": 1.0,
                "vdwclash": 0.08,
                "vdw": 0.06,
                "proximal": 0.04
            }
        },

        "colors":
        {
            "carbon": "grey",
            "non-std-res": "green",
            "variations": "magenta"
        }
}

# MAIN
if __name__ == '__main__':
    
    # ARGUMENT PARSING
    parser = argparse.ArgumentParser(description='''
    
    ############
    # ARPEGGIO #
    ############
    
    Interaction Viewer
    
    A program for calculating interatomic interactions,
    using only Open Source dependencies.
    
    Dependencies:
    - Python (v2.7)
    - PyMOL, run as `pymol -R` to enable the XML-RPC server.
    
    **You must have already run your structure with Arpeggio to generate the required output files**.
    Be careful about absolute and relative paths, this script might not mind but PyMOL probably will.
    
    ''', formatter_class=argparse.RawTextHelpFormatter)
        
    parser.add_argument('pdb', type=str, help='Path to the PDB file to be analysed.')
    parser.add_argument('-op', type=str, default='', help='The output postfix you used with Arpeggio (if any).')
    parser.add_argument('-xml', '--xml-rpc', action='store_true', help='Interact with PyMOL by XML-RPC server (`pymol -R`).')
    parser.add_argument('-s', '--script', action='store_true', help='Output a PyMOL script with the relevant commands.')
    parser.add_argument('-bs', '--use-all-binding-site', action='store_true', help='Use all the binding site contacts instead of just the INTER ones. Only works if you used a selection when running Arpeggio.')
    
    args = parser.parse_args()
    
    pdb_filename = args.pdb
    
    # CHANGE TO AN ABSOLUTE PATH
    pdb_filename = os.path.abspath(pdb_filename)
    
    output_postfix = args.op
    
    contacts_extension = '.contacts'
    
    if args.use_all_binding_site:
        contacts_extension = '.bs_contacts'
    
    contacts_filename = pdb_filename.replace('.pdb', output_postfix + contacts_extension)
    rings_filename = pdb_filename.replace('.pdb', output_postfix + '.rings')
    ari_filename = pdb_filename.replace('.pdb', output_postfix + '.ari') # ATOM-RING INTERACTIONS
    ri_filename = pdb_filename.replace('.pdb', output_postfix + '.ri') # RING-RING INTERACTIONS
    
    script_filename = pdb_filename.replace('.pdb', output_postfix + '.pml')
    
    if args.xml_rpc:
        srv = xmlrpclib.Server('http://{}:{}'.format(host, port))
    
    do = None
    
    if args.script:
        
        final_output = ''
        
        def do(command):
            
            global final_output # DON'T FORGET I EXIST!
            final_output = final_output + command + '\n'
    
    elif args.xml_rpc:
        do = srv.do
    
    else:
        print 'Please select XML-RPC (-xml) or script (-s) output.'
        sys.exit(1)
    
    # PYMOL SETUP
    do('reinitialize')

    # SET PYMOL VARIABLES FOR PRETTIER MOLECULES
    do('set valence, 1')
    do('set stick_rad, 0.15')
    do('set line_width, 1')
    do('set mesh_width, 0.3')

    # DECREASE FONT SIZE
    do('set label_size, 10')

    # MAKE SPHERES LOOK PRETTIER (WATER, IONS...')
    do('set sphere_scale, 0.15')

    # DASH PROPERTIES
    do('set dash_round_ends, 0')
    do('set dash_gap, 0.25')
    do('set dash_length, 0.05')

    # COLORS
    do('set_color maize, (251, 236, 93)') # IONIC VDW
    do('set_color cream, (255, 253, 208)') # IONIC PROXIMAL

    do('set_color electric_purple, (191, 0, 255)') # METAL COMPLEX VDW CLASH
    do('set_color medium_purple, (147, 112, 219)') # METAL COMPLEX VDW

    do('set_color amber, (255, 191, 0)')
    do('set_color peach, (255, 229, 180)')

    do('set_color electric_blue, (125, 249, 255)') #
    do('set_color baby_blue, (224, 255, 255)')

    do('set_color emerald, (80, 200, 120)')
    do('set_color moss_green, (173, 200, 173)')

    do('set_color rose, (255, 0, 127)')
    do('set_color brilliant_rose, (246, 83, 166)')
    do('set_color thulian_pink, (222, 111, 161)')

    do('set_color neon_red, (255, 0, 102)')

    # DNA
    do('set cartoon_ladder_mode, 1')
    do('set cartoon_ring_finder, 1')
    do('set cartoon_ring_mode, 3')
    do('set cartoon_nucleic_acid_mode, 4')
    do('set cartoon_ring_transparency, 0.5')
    
    # QUALITY
    do('set line_smooth, 1')
    do('set antialias, 2')
    do('set cartoon_fancy_helices, 1')
    do('set depth_cue, 1')
    do('set specular,1')
    do('set surface_quality, 1')
    do('set stick_quality, 15')
    do('set sphere_quality, 2')
    do('set cartoon_sampling, 14')
    do('set ribbon_sampling, 10')
    do('set transparency_mode, 2')
    do('set stick_ball, 1')
    do('set stick_ball_ratio, 1.5')
    do("rebuild")
    
    # LOAD STRUCTURE
    do('load {}'.format(pdb_filename))
    
    do('select binding_site, None')
    
    # DEFER UPDATE
    do('set defer_update, 1')
    
    # GET CONTACTS
    with open(contacts_filename, 'rb') as fo:
        
        used_interaction_types = set([])
        
        for n_lines, line in enumerate(fo):
            
            line = line.strip().split('\t')
            
            atom_bgn = line[0]
            atom_end = line[1]
            SIFt = [int(x) for x in line[2:17]]
            
            interactions = []
            dist_flag = ''
            
            # `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/contacts.py?at=default`
            # SIFT:
            # 0: CLASH
            # 1: COVALENT
            # 2: VDW_CLASH
            # 3: VDW
            # 4: PROXIMAL
            # 5: HBOND # 0
            # 6: WEAK_HBOND # 1
            # 7: HALOGEN_BOND # 2
            # 8: IONIC # 3
            # 9: METAL_COMPLEX # 4
            #10: AROMATIC # 5
            #11: HYDROPHOBIC # 6
            #12: CARBONYL # 7
            
            # MUTALLY EXCLUSIVE TYPES
            
            # CLASH
            if SIFt[0]:
                dist_flag = 'clash'
                
                # SKIPPING CLASHES FOR NOW
                continue
                
            # COVALENT
            elif SIFt[1]:
                
                dist_flag = 'covalent'
                
                # SKIPPING COVALENTS FOR NOW
                continue
            
            # VDW_CLASH
            elif SIFt[2]:
                dist_flag = 'vdwclash'
                
            # VDW
            elif SIFt[3]:
                dist_flag = 'vdw'
                
            # PROXIMAL
            elif SIFt[4]:
                dist_flag = 'proximal'
                
            else:
                dist_flag = 'undefined'
            
            # FEATURE CONTACTS
            
            interaction_type = ''
            
            for e, interaction in enumerate(SIFt[5:]):
                
                if interaction:
                    
                    if e == 0:
                        interaction_type = 'hbond'
                    
                    elif e == 1:
                        interaction_type = 'weakhbond'
                    
                    elif e == 2:
                        interaction_type = 'xbond'
                        
                    elif e == 3:
                        interaction_type = 'ionic'
                        
                    elif e == 4:
                        interaction_type = 'metalcomplex'
                    
                    elif e == 5:
                        interaction_type = 'aromatic'
                    
                    elif e == 6:
                        interaction_type = 'hydrophobic'
                    
                    elif e == 7:
                        interaction_type = 'carbonyl'
                        
                    elif e == 8:
                        interaction_type = 'polar'
                        
                    elif e == 9:
                        interaction_type = 'weakpolar'
                    
                    interactions.append((interaction_type, dist_flag))
                
            if not interaction_type:
                interactions.append(('undefined', dist_flag))
            
            for interaction_type, flag in interactions:
                
                do('distance {}-{}, {}, {}'.format(interaction_type, flag, atom_bgn, atom_end))
                
                do('show sticks, byres {}'.format(atom_bgn))
                do('show sticks, byres {}'.format(atom_end))
                
                do('select binding_site, binding_site + byres {}'.format(atom_bgn))
                do('select binding_site, binding_site + byres {}'.format(atom_end))
                
                used_interaction_types.add((interaction_type, flag))
            
            if not n_lines % 500:
                print 'Drew {} contacts.'.format(n_lines)
        
        print 'Drew {} total contacts.'.format(n_lines)
        
        for interaction_type, flag in used_interaction_types:
            
            label = '{}-{}'.format(interaction_type, flag)
            
            do('color {}, {}'.format(pymol_config['dashcolor'][interaction_type][flag], label))
            do('set dash_radius, {}, {}'.format(pymol_config['dashradius'][interaction_type][flag], label))
            do('set dash_gap, {}, {}'.format(pymol_config['dashgap'][interaction_type][flag], label))
            do('set dash_length, {}, {}'.format(pymol_config['dashlength'][interaction_type][flag], label))
    
    # RING INTERACTIONS
    
    # RING CENTROIDS
    with open(rings_filename, 'rb') as fo:
        for line in fo:
            line = line.strip().split('\t')
            
            do('pseudoatom ring_centers, pos={}'.format(line[2]))
    
    do('hide everything, ring_centers')
    #do('show spheres, ring_centers')
    
    # RING-RING INTERACTIONS
    with open(ri_filename, 'rb') as fo:
        for line in fo:
            line = line.strip().split('\t')
            
            if line[8] == 'INTRA_NON_SELECTION':
                continue
            
            do('pseudoatom pt1, pos={}'.format(line[2]))
            do('pseudoatom pt2, pos={}'.format(line[5]))
            
            do('distance ring-{}, pt1, pt2'.format(line[6]))
            
            do('delete pt1')
            do('delete pt2')
    
    do('set dash_radius, 0.25, ring-*')
    do('set dash_gap, 0.1, ring-*')
    do('color white, ring-*')
    
    # ATOM-RING INTERACTIONS
    with open(ari_filename, 'rb') as fo:
        for line in fo:
            line = line.strip().split('\t')
            
            do('pseudoatom pt1, pos={}'.format(line[3]))
            
            if line[6] == 'INTRA_NON_SELECTION':
                continue
            
            if type(eval(line[4])) is list:
            
                for itype in eval(line[4]):
                    do('distance {}, pt1, {}'.format(itype, line[0]))
            
            do('delete pt1')
    
    do('set dash_radius, 0.15, *PI')
    do('set dash_gap, 0.1, *PI')
    
    do('color white, CARBONPI')
    do('color blue, DONORPI')
    do('color green, HALOGENPI')
    do('color red, CATIONPI')
    do('color yellow, METSULPHURPI')
    
    # PRETTINESS
    do('hide labels')
    do('util.cbaw')
    do('bg_color white')
    do('show cartoon')
    do('set cartoon_side_chain_helper, 1')
    do('hide lines')
    do('hide everything, het')
    do('show sticks, het')
    do('show spheres, het')
    do('disable undefined-proximal')
    
    # UPDATE PYMOL NOW
    do('set defer_update, 0')
    
    if args.script:
        with open(script_filename, 'wb') as fo:
            fo.write(final_output)
        