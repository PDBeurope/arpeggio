from collections import OrderedDict

DEFAULT_SIFT = OrderedDict()

# FEATURE SIFTS
    # 5: 0: HBOND
    # 6: 1: WEAK_HBOND
    # 7: 2: HALOGEN_BOND
    # 8: 3: IONIC
    # 9: 4: METAL_COMPLEX
    #10: 5: AROMATIC
    #11: 6:  HYDROPHOBIC
    #12: 7: CARBONYL
    
    #8: POLAR - HBOND WITHOUT ANGLES
    #9: WEAK POLAR - WEAK HBOND WITHOUT ANGLES
    
FEATURE_SIFT = ('H-Bond', 'Weak H-Bond', 'Halogen Bond', 'Ionic', 'Metal Complex', 'Aromatic', 'Hydrophobic', 'Carbonyl', 'Polar', 'Weak Polar')

#http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
VDW_RADII = {
             'H': 1.2
            }

METALS = set(['LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI',
              'CU', 'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
              'IN', 'SN', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO',
              'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI',
              'PO', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF'])

HALOGENS = set(['F', 'CL', 'BR', 'I', 'AT'])

MAINCHAIN_ATOMS = set(['N', 'C', 'CA', 'O', 'OXT'])

# `https://github.com/openbabel/openbabel/blob/master/src/atom.cpp`
# THE NUMBER OF VALENCE ELECTRONS IN A FREE ATOM
VALENCE = [0,1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,
11,12,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,1,2,
4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,
8,1,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,5,6,7,8,9,10,11,12]

# THE NUMBER OF ELECTRONS REQUIRED TO MAKE UP A FULL VALENCE SHELL (I.E TO MAKE IT HAPPY)
SHELL = [0,2,2,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,18,18,18,18,18,18,
18,18,18,18,8,8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,8,
8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,8,8,18,18,18,18,
18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]

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
        
        # SEE RDKIT `BaseFeatures.fdef`
        "pos ionisable":
        {
            "rdkit basic group": "[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
            #"basic group": "[NH0+0$(*(-[C!$(*=*)])(-[C!$(*=*)])-[C!$(*=*)])!$(*-[a]),NH+0$(*(-[C!$(*=*)])-[C!$(*=*)])!$(*-[a]),NH2+0$(*-[C!$(*=*)])!$(*-[a])]",
            "imidazole": "n1cncc1",
            "guanidine": "N=C(-N)-N",
            #"posn": "[#7+]"
            "rdkit posn": "[#7;+;!$([N+]-[O-])]"
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
            "polar distance": 3.5,
            "angle rad": 1.57,
            "angle degree": 90.0
        },

        "weak hbond":
        {
            "distance": 3.6,
            "weak polar distance": 3.5,
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
            "atom_aromatic_distance": 4.5,
            "met_sulphur_aromatic_distance": 6.0
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

THETA_REQUIRED = set(['CARBONPI', 'CATIONPI', 'DONORPI', 'HALOGENPI'])