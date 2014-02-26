from collections import OrderedDict

DEFAULT_SIFT = OrderedDict()

DEFAULT_SIFT['proximal'] = 0
DEFAULT_SIFT['pos_ionic'] = 0
DEFAULT_SIFT['neg_ionic'] = 0
DEFAULT_SIFT['carbonyl_o'] = 0
DEFAULT_SIFT['carbonyl_c'] = 0
DEFAULT_SIFT['aromatic'] = 0
DEFAULT_SIFT['hydrophobic'] = 0
DEFAULT_SIFT['xbond_donor'] = 0
DEFAULT_SIFT['xbond_acceptor'] = 0
DEFAULT_SIFT['weak_hbond_donor'] = 0
DEFAULT_SIFT['weak_hbond_acceptor'] = 0
DEFAULT_SIFT['hbond_donor'] = 0
DEFAULT_SIFT['hbond_acceptor'] = 0
DEFAULT_SIFT['metal_metal'] = 0
DEFAULT_SIFT['metal_acceptor'] = 0
DEFAULT_SIFT['ring_face'] = 0
DEFAULT_SIFT['ring_offset'] = 0
DEFAULT_SIFT['ring_edge'] = 0
DEFAULT_SIFT['ring_donor'] = 0 # COVERS CARBONPI, CATIONPI, DONORPI, HALOGENPI

#http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
VDW_RADII = {
             'H': 1.2
            }

METALS = set(['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
              'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
              'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
              'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
              'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf'])

HALOGENS = set(['F', 'CL', 'BR', 'I', 'AT'])

ATOM_TYPES = {
                
        "hbond acceptor":
        {
            "acceptor": "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]"
        },

        "hbond donor":
        {
            "donor": "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]"
        },

        "xbond acceptor":
        {
            "acceptor": "[#7&!$([nX3,#7v5]),#8,#16&!$([#16v4,#16v6]);!$([*+1,*+2,*+3])]"
        },

        "xbond donor":
        {
            "donor": "[F,Cl,Br,I;X1;$([F,Cl,Br,I]-[#6,#8]);!$([F,Cl,Br,I]C[F,Cl,Br,I])]"
        },

        "weak hbond acceptor":
        {
            "c-x halogens": "[F,Cl,Br,I;X1;$([F,Cl,Br,I]-[#6,#8])]"
        },

        "weak hbond donor":
        {
            "donor": "[#6!H0]"
        },

        "pos ionisable":
        {
            "basic group": "[NH0+0$(*(-[C!$(*=*)])(-[C!$(*=*)])-[C!$(*=*)])!$(*-[a]),NH+0$(*(-[C!$(*=*)])-[C!$(*=*)])!$(*-[a]),NH2+0$(*-[C!$(*=*)])!$(*-[a])]",
            "imidazole": "n1cncc1",
            "guanidine": "N=C(-N)-N",
            "posn": "[#7+]"
        },

        "neg ionisable":
        {
            "acidic group": "[OH,OH0-]-[C,S]=[O,P,S]",
            "hydroxylic acid": "[OH$(*-[A]=[!#6])]"
        },

        "hydrophobe":
        {
            "hydrophobe": "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]"
        },

        "carbonyl oxygen":
        {
            "oxygen": "[$([OH0]=[CX3,c])]"
        },

        "carbonyl carbon":
        {
            "carbon": "[$([CX3,c](=[OH0]))]"
        },

        "aromatic":
        {
            "arom_4": "[a;r4,!R1&r3]1:[a;r4,!R1&r3]:[a;r4,!R1&r3]:[a;r4,!R1&r3]:1",
            "arom_5": "[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1",
            "arom_6": "[a;r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:1",
            "arom_7": "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
            "arom_8": "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1"
        }
}

CONTACT_TYPES_DIST_MAX = 4.5

CONTACT_TYPES = {
        "hbond":
        {
            "distance": 3.9,
            "angle rad": 1.57,
            "angle degree": 90.0
        },

        "weak hbond":
        {
            "distance": 3.6,
            "angle rad": 2.27,
            "angle degree": 130.0,
            "cx angle min rad": 0.52,
            "cx angle min degree": 30.0,
            "cx angle max rad": 2.62,
            "cx angle max degree": 150.0
        },

        "aromatic":
        {
            "distance": 4.0,
            "centroid_distance": 6.0,
            "atom_aromatic_distance": 4.5
        },

        "xbond":
        {
            "catmap distance": 1.85, # SAME AS BROMINE 
            "angle theta 1 rad": 2.09,
            "angle theta 1 degree": 120.0,
            "angle theta 2 min rad": 1.22,
            "angle theta 2 max rad": 2.97,
            "angle theta 2 min degree": 70.0,
            "angle theta 2 max degree": 170.0
        },

        "ionic":
        {
            "distance": 4.0
        },

        "hydrophobic":
        {
            "distance": 4.5
        },

        "carbonyl":
        {
            "distance": 3.6
        },

        "metal":
        {
            "distance": 2.8
        }
}