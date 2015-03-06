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

AMIDE_SMARTS = '[NX3][CX3](=[OX1])[#6]' # DEFINITION FROM `http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html`

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
            "acceptor"          : "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
            "enol"              : "[$([nH]:@c(=O))]",
            "tautomeric nH"     : "[$([n;H1;v3;!$([nH]cccc)])]",
            # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
            "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]"
        },

        "hbond donor":
        {
            "donor"             : "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]",
            "oxygen acid"       : "[$([O;H0;$(O=C([OH])-*)])]",
            "tautomer nH"       : "[$(n:a:[nH])]",
            # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES 
            "oxygen amide term" : "[$([O;H0;$(O=C-[NH2])])]"
                   
        },

        "xbond acceptor":
        {
            # SAME AS HBA
            "acceptor"          : "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
            "enol"              : "[$([nH]:@c(=O))]",
            "tautomeric nH"     : "[$([n;H1;v3;!$([nH]cccc)])]",
            # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
            "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]"
        },

        "xbond donor":
        {
            "donor"             : "[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]"
        },

        "weak hbond acceptor":
        {
            # SAME AS HBA
            "acceptor"          : "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
            "enol"              : "[$([nH]:@c(=O))]",
            "tautomeric nH"     : "[$([n;H1;v3;!$([nH]cccc)])]",
            # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
            "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]",
            "c-x halogens"      : "[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]"
        },

        "weak hbond donor":
        {
            "donor"             : "[#6!H0]"
        },
        
        # SEE RDKIT `BaseFeatures.fdef`
        "pos ionisable":
        {
            "rdkit basic group" : "[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
            "imidazole"         : "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1",
            "guanidine amidine" : "NC(=N)",
            "rdkit posn"        : "[#7;+;!$([N+]-[O-])]",
            "cations"           : "[$([*+1,*+2,*+3]);!$([N+]-[O-])]",
            "metals"            : "[Li,Be,Na,Mg,Al,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf]"
        },

        "neg ionisable":
        {
            "O acidic group"    : "[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]",
            "anions"            : "[*-1,*-2]"
        },

        "hydrophobe":
        {
            "hydrophobe"        : "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]"
        },

        "carbonyl oxygen":
        {
            "oxygen"            : "[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]"
        },

        "carbonyl carbon":
        {
            "carbon"            : "[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]"
        },

        "aromatic":
        {
            "arom_4"            : "[a;r4,!R1&r3]1:[a;r4,!R1&r3]:[a;r4,!R1&r3]:[a;r4,!R1&r3]:1",
            "arom_5"            : "[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1",
            "arom_6"            : "[a;r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:1",
            "arom_7"            : "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
            "arom_8"            : "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1"
        }
}

STD_RES = set(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO',
               'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'])

PROT_ATOM_TYPES = {
        "hbond acceptor":       [
                                "ALAO",         # all the carbonyl Oxygens in the main chain
                                "ARGO",
                                "ASNO",
                                "ASPO",
                                "CYSO",
                                "GLNO",
                                "GLUO",
                                "GLYO",
                                "HISO",
                                "ILEO",
                                "LEUO",
                                "LYSO",
                                "METO",
                                "PHEO",
                                "PROO",
                                "SERO",
                                "THRO",
                                "TRPO",
                                "TYRO",
                                "VALO",
                                "ALAOXT",         # all the carbonyl Oxygens terminals
                                "ARGOXT",
                                "ASNOXT",
                                "ASPOXT",
                                "CYSOXT",
                                "GLNOXT",
                                "GLUOXT",
                                "GLYOXT",
                                "HISOXT",
                                "ILEOXT",
                                "LEUOXT",
                                "LYSOXT",
                                "METOXT",
                                "PHEOXT",
                                "PROOXT",
                                "SEROXT",
                                "THROXT",
                                "TRPOXT",
                                "TYROXT",
                                "VALOXT",
                                "ASNOD1",  
                                "ASNND2",       #for the ambiguity of the position of the N and O
                                "ASPOD1",
                                "ASPOD2",
                                "GLNOE1"
                                "GLNNE2",       #for the ambiguity of the position of the N and O
                                "GLUOE1",
                                "GLUOE2",
                                "HISND1",       #for the ambiguity of the position of the N/C
                                "HISCE1",       #for the ambiguity of the position of the N/C
                                "HISNE2",       #for the ambiguity of the position of the N/C         
                                "HISCD2",       #for the ambiguity of the position of the N/C
                                "METSD",        #http://pubs.acs.org/doi/abs/10.1021/jz300207k and pubid 19089987
                                "CYSSG",        #pubid 19089987, also when they from di-sulfide (Cys-Cys, fig 8 paper)
                                "SEROG",        #isostar plots
                                "THROG1",       #isostar plots
                                "TYROH"         #isostar plots
                                ]
        ,
        "hbond donor":          [
                                "ALAN",          # all the amide nitrogens in the main chain except proline
                                "ARGN",
                                "ASNN",
                                "ASPN",
                                "CYSN",
                                "GLNN",
                                "GLUN",
                                "GLYN",
                                "HISN",
                                "ILEN",
                                "LEUN",
                                "LYSN",
                                "METN",
                                "PHEN",
                                "SERN",
                                "THRN",
                                "TRPN",
                                "TYRN",
                                "VALN",
                                "ARGNE",
                                "ARGNH1",
                                "ARGNH2",
                                "ASNND2",
                                "ASNOD1",       #for the ambiguity of the position of the N and O
                                "CYSSG",        #http://www.ncbi.nlm.nih.gov/pubmed/19089987
                                "GLNNE2",
                                "GLNOE1",       #for the ambiguity of the position of N/O
                                "HISND1",       #for the ambiguity of the position of the N/C
                                "HISCE1",       #for the ambiguity of the position of the N/C
                                "HISNE2",       #for the ambiguity of the position of the N/C         
                                "HISCD2",       #for the ambiguity of the position of the N/C
                                "LYSNZ",
                                "SEROG",        
                                "THROG1",       
                                "TRPNE1",
                                "TYROH"         
                                ]
                
        ,
        "xbond acceptor": [
                                "ALAO",         # all the carbonyl Oxygens in the main chain
                                "ARGO",
                                "ASNO",
                                "ASPO",
                                "CYSO",
                                "GLNO",
                                "GLUO",
                                "GLYO",
                                "HISO",
                                "ILEO",
                                "LEUO",
                                "LYSO",
                                "METO",
                                "PHEO",
                                "PROO",
                                "SERO",
                                "THRO",
                                "TRPO",
                                "TYRO",
                                "VALO",
                                "ALAOXT",         # all the carbonyl Oxygens terminals
                                "ARGOXT",
                                "ASNOXT",
                                "ASPOXT",
                                "CYSOXT",
                                "GLNOXT",
                                "GLUOXT",
                                "GLYOXT",
                                "HISOXT",
                                "ILEOXT",
                                "LEUOXT",
                                "LYSOXT",
                                "METOXT",
                                "PHEOXT",
                                "PROOXT",
                                "SEROXT",
                                "THROXT",
                                "TRPOXT",
                                "TYROXT",
                                "VALOXT",
                                "ASNOD1",  
                                "ASNND2",       #for the ambiguity of the position of the N and O
                                "ASPOD1",
                                "ASPOD2",
                                "GLNOE1"
                                "GLNNE2",       #for the ambiguity of the position of the N and O
                                "GLUOE1",
                                "GLUOE2",
                                "HISND1",       #for the ambiguity of the position of the N/C
                                "HISCE1",       #for the ambiguity of the position of the N/C
                                "HISNE2",       #for the ambiguity of the position of the N/C         
                                "HISCD2",       #for the ambiguity of the position of the N/C
                                "METSD",        #http://pubs.acs.org/doi/abs/10.1021/jz300207k and pubid 19089987
                                "CYSSG",        #pubid 19089987, also when they from di-sulfide (Cys-Cys, fig 8 paper)
                                "SEROG",        #isostar plots
                                "THROG1",       #isostar plots
                                "TYROH"         #isostar plots
                                ]
        ,
        "weak hbond acceptor": [
                                "ALAO",         # all the carbonyl Oxygens in the main chain
                                "ARGO",
                                "ASNO",
                                "ASPO",
                                "CYSO",
                                "GLNO",
                                "GLUO",
                                "GLYO",
                                "HISO",
                                "ILEO",
                                "LEUO",
                                "LYSO",
                                "METO",
                                "PHEO",
                                "PROO",
                                "SERO",
                                "THRO",
                                "TRPO",
                                "TYRO",
                                "VALO",
                                "ALAOXT",         # all the carbonyl Oxygens terminals
                                "ARGOXT",
                                "ASNOXT",
                                "ASPOXT",
                                "CYSOXT",
                                "GLNOXT",
                                "GLUOXT",
                                "GLYOXT",
                                "HISOXT",
                                "ILEOXT",
                                "LEUOXT",
                                "LYSOXT",
                                "METOXT",
                                "PHEOXT",
                                "PROOXT",
                                "SEROXT",
                                "THROXT",
                                "TRPOXT",
                                "TYROXT",
                                "VALOXT",
                                "ASNOD1",  
                                "ASNND2",       #for the ambiguity of the position of the N and O
                                "ASPOD1",
                                "ASPOD2",
                                "GLNOE1"
                                "GLNNE2",       #for the ambiguity of the position of the N and O
                                "GLUOE1",
                                "GLUOE2",
                                "HISND1",       #for the ambiguity of the position of the N/C
                                "HISCE1",       #for the ambiguity of the position of the N/C
                                "HISNE2",       #for the ambiguity of the position of the N/C         
                                "HISCD2",       #for the ambiguity of the position of the N/C
                                "METSD",        #http://pubs.acs.org/doi/abs/10.1021/jz300207k and pubid 19089987
                                "CYSSG",        #pubid 19089987, also when they from di-sulfide (Cys-Cys, fig 8 paper)
                                "SEROG",        #isostar plots
                                "THROG1",       #isostar plots
                                "TYROH"         #isostar plots
                                ]
        ,
        "weak hbond donor": [
                                "ALACA",         # all the c-alphas 
                                "ARGCA",
                                "ASNCA",
                                "ASPCA",
                                "CYSCA",
                                "GLNCA",
                                "GLUCA",
                                "GLYCA",
                                "HISCA",
                                "ILECA",
                                "LEUCA",
                                "LYSCA",
                                "METCA",
                                "PHECA",
                                "PROCA",
                                "SERCA",
                                "THRCA",
                                "TRPCA",
                                "TYRCA",
                                "VALCA",
                                "ALACB",         #cb and further down
                                "ARGCB",
                                "ARGCG",
                                "ARGCD",
                                "ASNCB",
                                "ASPCB",
                                "CYSCB",
                                "GLNCB",
                                "GLNCG",
                                "GLUCB",
                                "GLUCG",
                                "GLNCB",
                                "HISCB",
                                "ILECB",
                                "ILECG1",
                                "ILECD1",
                                "ILECG2",
                                "LEUCB",
                                "LEUCG",
                                "LEUCD1",
                                "LEUCD2",
                                "LYSCB",
                                "LYSCG",
                                "LYSCD",
                                "LYSCE",
                                "METCB",
                                "METCG",
                                "METCE",
                                "PHECB",
                                "PHECG",
                                "PHECD1",
                                "PHECD2",
                                "PHECE1",
                                "PHECE2",
                                "PHECZ",
                                "PROCB",
                                "PROCG",
                                "PROCD",
                                "SERCB",
                                "THRCB",
                                "THRCG2",
                                "TRPCB",
                                "TRPCD1"
                                "TRPCE3",
                                "TRPCZ3",
                                "TRPCH2",
                                "TRPCZ2",                                
                                "TYRCB",
                                "TYRCD1",
                                "TYRCD2",
                                "TYRCE1",
                                "TYRCE2",
                                "TRYCB",
                                "VALCB",
                                "VALCG1",
                                "VALCG2"
        ]
        ,
        "pos ionisable":   [
                                "ARGNE",
                                "ARGCZ",
                                "ARGNH1",
                                "ARGNH2",
                                "HISCG",
                                "HISND1",
                                "HISCE1",
                                "HISNE2",
                                "HISCD2",
                                "LYSNZ"
        ]
        ,
        "neg ionisable": [
                                "ASPOD1",
                                "ASPOD2",
                                "GLUOE1",
                                "GLUOE2"
        ]
        ,
        "hydrophobe": [
                                "ALACB",        
                                "ARGCB",
                                "ARGCG",
                                "ASNCB",
                                "ASPCB",
                                "CYSCB",     #sulfur in Cys has an Hydrogen, it is polarised
                                "GLNCB",
                                "GLNCG",
                                "GLUCB",
                                "GLUCG",
                                "GLNCB",
                                "HISCB",
                                "ILECB",
                                "ILECG1",
                                "ILECD1",
                                "ILECG2",
                                "LEUCB",
                                "LEUCG",
                                "LEUCD1",
                                "LEUCD2",
                                "LYSCB",
                                "LYSCG",
                                "LYSCD",
                                "METCB",
                                "METCG",
                                "METSD",
                                "METCE",
                                "PHECB",
                                "PHECG",
                                "PHECD1",
                                "PHECD2",
                                "PHECE1",
                                "PHECE2",
                                "PHECZ",
                                "PROCB",
                                "PROCG",
                                "THRCG2",
                                "TRPCB",
                                "TRPCG",
                                "TRPCD2",
                                "TRPCE3",
                                "TRPCZ3",
                                "TRPCH2",
                                "TRPCZ2",                                
                                "TRYCB",
                                "TYRCG",
                                "TYRCD1",
                                "TYRCD2",
                                "TYRCE1",
                                "TYRCE2",
                                "VALCB",
                                "VALCG1",
                                "VALCG2"
        ]
        ,
        "carbonyl oxygen":  [
                                "ALAO",         # all the carbonyl Oxygens in the main chain
                                "ARGO",
                                "ASNO",
                                "ASPO",
                                "CYSO",
                                "GLNO",
                                "GLUO",
                                "GLYO",
                                "HISO",
                                "ILEO",
                                "LEUO",
                                "LYSO",
                                "METO",
                                "PHEO",
                                "PROO",
                                "SERO",
                                "THRO",
                                "TRPO",
                                "TYRO",
                                "VALO"
        ]
        ,
        "carbonyl carbon": [
                                "ALAC",         
                                "ARGC",
                                "ASNC",
                                "ASPC",
                                "CYSC",
                                "GLNC",
                                "GLUC",
                                "GLYC",
                                "HISC",
                                "ILEC",
                                "LEUC",
                                "LYSC",
                                "METC",
                                "PHEC",
                                "PROC",
                                "SERC",
                                "THRC",
                                "TRPC",
                                "TYRC",
                                "VALC"
        ]
        ,
        "aromatic": [
                                "HISCG",
                                "HISND1",
                                "HISCE1",
                                "HISNE2",
                                "HISCD2",
                                "PHECG",
                                "PHECD1",
                                "PHECD2",
                                "PHECE1",
                                "PHECE2",
                                "PHECZ",
                                "TRPCG",
                                "TRPCD1",
                                "TRPCD2",
                                "TRPNE1",
                                "TRPCE2",
                                "TRPCE3",
                                "TRPCZ2",
                                "TRPCZ3",
                                "TRPCH2",
                                "TYRCG",
                                "TYRCD1",
                                "TYRCD2",
                                "TYRCE1",
                                "TYRCE2",
                                "TYRCZ"
  
        ]
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
        
        "amide":
        {
             "centroid_distance": 6.0,
             "angle degree": 30.0,
             "angle rad": 0.52
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

COMMON_SOLVENTS = set(["MA4", "EOH", "AGC", "EOM", "3CN", "CPS", "CPT", "TAU", "TAR", "TAS", "SPK",
"1FH", "SPM", "VA3", "SPD", "XPE", "GD", "GA", "EHN", "PE8", "LAK", "OES",
"MAC", "LAT", "CP2", "MAN", "TSM", "DTN", "TSD", "TSE", "YBT", "EDO", "PCL",
"CB5", "BEQ", "F09", "DTV", "NO3", "NI1", "IUM", "ZN", "MN3", "NFC", "NFB",
"NFO", "K", "NFS", "NFR", "DVT", "HSG", "MG", "NOE", "MN", "PC4", "CBM", "SE",
"CBX", "MXE", "BE7", "GLO", "MM4", "GLC", "DOD", "SMO", "GLY", "DOX", "FE",
"12P", "543", "FCY", "MDD", "AU3", "FCL", "FCO", "MW3", "AUC", "VER", "F",
"FMT", "FMS", "JEF", "V", "SX", "SBT", "SR", "SBY", "LMT", "LMU", "T3A", "SM",
"SBE", "SB", "SBO", "MMC", "YH", "DPR", "YB", "CAT", "DPG", "DPO", "CAC", "NHE",
"PGE", "LA", "PGO", "PGL", "HG2", "LI", "LU", "PGR", "PGQ", "144", "PG6", "NH4",
"PG4", "PG5", "ZRC", "PG0", "HGI", "HGC", "CFN", "CFO", "CFM", "TFH", "NCP",
"OC5", "TFA", "OC7", "NCO", "Y1", "CFT", "TFP", "SF4", "ZNO", "SF3", "FNE",
"EGL", "7PE", "RU", "LBT", "SAL", "RE", "RB", "BF4", "URE", "MRY", "CHM", "MRD",
"IPA", "IPH", "ZN3", "ZN2", "LI1", "4OA", "RHD", "OCM", "OCL", "OCO", "OCN",
"NGN", "P2K", "HDD", "PEG", "HDN", "MG8", "BMA", "2FU", "WO5", "WO4", "YT3",
"WO2", "CEQ", "2FH", "AF3", "ENC", "GCD", "B3P", "CE1", "EU", "XCC", "MGF",
"TBR", "TBU", "WCC", "OC8", "BTC", "BTB", "OC1", "OC3", "I42", "OC4", "XL1",
"OC6", "3OF", "ATO", "B7G", "EMC", "SEA", "SEK", "MO1", "EU3", "CLN", "CLO",
"DUM", "UNX", "BEF", "CLF", "PDT", "MO6", "PDO", "TLA", "FRU", "BBX", "BBU",
"CLP", "OMO", "ZEM", "KR", "P6G", "2ME", "2MO", "CYS", "POR", "POP", "PON",
"2BM", "FEO", "CYN", "PNL", "GOL", "202", "F50", "DHE", "DHD", "PO2", "PO3",
"MLP", "OTE", "KO4", "E1H", "NRU", "3CO", "BGC", "2PA", "2PE", "NI2", "NI3",
"DZZ", "2PO", "2PN", "TCN", "MQD", "HSH", "TZZ", "OF2", "OF3", "OF1", "HSW",
"C2O", "TCH", "VXA", "GPX", "BNG", "NIK", "C2C", "F2O", "SYL", "NAW", "FDE",
"FDD", "DMF", "MO3", "GBL", "CRY", "MO7", "NFE", "MO5", "MO4", "DMS", "DMU",
"DMX", "RGI", "AEM", "PS5", "PR", "LDA", "PT", "PB", "ZO3", "PD", "MH2", "MH3",
"MHM", "BU2", "BU3", "DDQ", "BU1", "MHA", "S", "DDH", "T1A", "LCP", "MOS",
"O4M", "MOH", "MOO", "BCN", "UMQ", "HNI", "PEU", "FSO", "TMA", "ETF", "CD1",
"PEJ", "CD3", "FSX", "COH", "ETI", "COM", "RU7", "CON", "HE5", "HE6", "BF2",
"16D", "16A", "P33", "ND4", "HEV", "CO", "CN", "CM", "HES", "CA", "HEZ", "CD",
"HEA", "HEB", "HEC", "CS", "CR", "CU", "PIS", "CO5", "PSL", "PE4", "PE5", "PE7",
"PE3", "FS1", "FS2", "FS3", "FS4", "PC3", "OEC", "NMO", "NML", "DIO", "DIA",
"C8E", "MSF", "C10", "C15", "DIS", "SOG", "PPK", "IR3", "SOM", "TBA", "PPM",
"SOR", "GAI", "IR", "3NI", "IRI", "SO4", "MTL", "IN", "VX", "PP9", "MTF", "CM5",
"MTD", "ARS", "SDS", "BA", "CCN", "CCH", "BR", "EPE", "O2", "BOG", "PIN", "MYQ",
"MYR", "POL", "PIG", "N8E", "F3S", "IOH", "IOD", "CD5", "BO4", "DR6", "NA2",
"NA6", "CXE", "NA5", "OS", "REO", "HLT", "U1", "FLH", "HO", "FLO", "MP1", "EEE",
"HG", "PTN", "SCN", "SFO", "ETN", "MNR", "MNH", "OSM", "ROP", "SFN", "MNC",
"PO4", "3BR", "FMI", "OX", "TRE", "MN5", "O", "PBM", "TRS", "PT4", "TRT", "HTG",
"MW2", "MPD", "HTO", "MPO", "MTO", "6WO", "MPR", "MPS", "CNB", "ACA", "GD3",
"SUC", "IDT", "SUL", "P4C", "ACN", "2OF", "S0H", "IDO", "HFM", "ACU", "ACT",
"ACY", "2OS", "NI", "CNM", "CU1", "NO", "NA", "MW1", "TEE", "FE2", "TEP", "1BO",
"FEC", "FEA", "DET", "FEL", "CUZ", "FES", "CUA", "CUO", "CUN", "CUM", "CUL",
"PTL", "MN6", "BRP", "BRO", "BRM", "BRJ", "MES", "AZI", "VO3", "1PG", "1PE",
"CIT", "1PS", "CN1", "CL", "TL", "SGM", "HO3", "TE", "TB", "ICI", "AG", "FPO",
"AL", "CNF", "AR", "AU", "CE", "SRM", "ICT", "15P", "ZH3", "RTC", "2NO", "DXG",
"DXE", "VSO"])

STANDARD_NUCLEOTIDES = set(['A', 'C', 'G', 'I', 'U', 'DA', 'DC', 'DG', 'DI', 'DT', 'DU', 'N'])