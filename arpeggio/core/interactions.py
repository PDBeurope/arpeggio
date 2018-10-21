import collections
import csv
import logging
import operator
import os
from functools import reduce

import numpy as np
import openbabel as ob
from Bio.PDB import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Atom import DisorderedAtom

from arpeggio.core import (AtomSerialError, OBBioMatchError,
                           SiftMatchError, config, protein_reader, utils)

AtomContacts = collections.namedtuple('AtomContacts',
                                      ['start_atom', 'end_atom', 'sifts', 'contact_type', 'distance'])

Parameters = collections.namedtuple('Parameters',
                                    ['vdw_comp_factor', 'interacting_threshold', 'has_hydrogens',
                                     'ph', ])


class InteractionComplex:
    def __init__(self, filename, vdw_comp=0.1, interacting=5.0, ph=7.4):
        """Creates interaction complex to be used by Arperggio for
        calculations.

        Args:
            filename (str): Path to the structure to be processed.
            vdw_comp (float, optional): Defaults to 0.1. Compensation
                factor for VdW radii dependent interaction types.
            interacting (float, optional): Defaults to 5.0. Distance
                cutoff for grid points to be \'interacting\' with
                the entity.
            ph (float, optional): Defaults to 7.4. Ph to be used if
                structure does not containg hydrogens.
        """
        self.id = os.path.basename(filename).split('.')[0]
        # Biopython part
        self.biopython_str = self._read_in_biopython(filename)
        self.s_atoms = list()
        for atom in list(self.biopython_str.get_atoms()):
            if type(atom) == DisorderedAtom:
                for altloc, alt_atom in atom.child_dict.items():
                    self.s_atoms.append(alt_atom)
            else:
                self.s_atoms.append(atom)

        # obabel data
        self.ob_mol = self._read_openbabel(filename)

        self.ob_to_bio = {}  # openbabel to biopython atom mapping
        self.bio_to_ob = {}  # biopython to openbabel atom mapping
        self.selection = []  # user atom selection
        self.selection_plus = []  # extended user atom selection
        self.selection_plus_residues = []  # residues
        self.ns = None  # neighbourhood search
        self.contacts = []  # arperggio contacts as a list of AtomContacts

        self.params = Parameters(
            vdw_comp_factor=vdw_comp,
            interacting_threshold=interacting,
            has_hydrogens=any(x.element == 'H' for x in self.s_atoms),
            ph=ph)

        self._establish_structure_mappping()

        if self.params.has_hydrogens:
            self.input_has_hydrogens = True
            logging.debug(("Detected that the input structure contains hydrogens. "
                          "Hydrogen addition will be skipped."))
        else:
            self.ob_mol.AddHydrogens(False, True, ph)
            self.input_has_hydrogens = False
            logging.debug('Added hydrogens.')

    # region public methods

    def structure_checks(self):
        """Check if structure is properly formated

        Raises:
            AtomSerialError: in case there are duplicities in atom ids.
        """
        all_serials = [x.serial_number for x in self.s_atoms]

        if len(all_serials) > len(set(all_serials)):
            raise AtomSerialError

    def address_ambiguities(self):
        """Remove ambiguous definitions from the config
        """
        # REMOVE FROM SMARTS DEFINITIONS
        config.ATOM_TYPES['hbond acceptor'].pop('NH2 terminal amide', None)
        config.ATOM_TYPES['hbond donor'].pop('oxygen amide term', None)
        config.ATOM_TYPES['xbond acceptor'].pop('NH2 terminal amide', None)
        config.ATOM_TYPES['weak hbond acceptor'].pop('NH2 terminal amide', None)

        # REMOVE FROM PROTEIN ATOM DEFINITIONS
        config.PROT_ATOM_TYPES['hbond acceptor'] = [x for x in config.PROT_ATOM_TYPES['hbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]
        config.PROT_ATOM_TYPES['hbond donor'] = [x for x in config.PROT_ATOM_TYPES['hbond donor'] if x not in ('ASNOD1', 'GLNOE1', 'HISCE1', 'HISCD2')]
        config.PROT_ATOM_TYPES['xbond acceptor'] = [x for x in config.PROT_ATOM_TYPES['xbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]
        config.PROT_ATOM_TYPES['weak hbond acceptor'] = [x for x in config.PROT_ATOM_TYPES['weak hbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]

    def write_atom_types(self, wd):
        """Write out atom types in a csv format.

        Args:
            wd (str): Working directory where the output will be stored.
        """
        filepath = os.path.join(wd, self.id + '_atomtypes.csv')
        with open(filepath, 'w') as fo:
            writer = csv.writer(fo, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['atom', 'atom_types'])

            for atom in self.s_atoms:
                writer.writerow([utils.make_pymol_string(atom), sorted(tuple(atom.atom_types))])

        logging.debug('Typed atoms.')

    def write_contacts(self, selection, wd):
        """Write out contacts detected by arpergio for the structure
        and the selection.

        Args:
            selection (list of Atoms): Selection perceived by Arpergio
            wd (str): Working directory
        """
        if len(selection) == 0:
            contacts = os.path.join(wd, self.id + '_contacts.csv')
            self.__write_contact_file(contacts, self.contacts)
        else:
            contacts = os.path.join(wd, self.id + '_bs_contacts.csv')
            temp_contacts = filter(lambda l: l.contact_type in ('INTER', 'SELECTION_WATER', 'WATER_WATER'), self.contacts)
            self.__write_contact_file(contacts, temp_contacts)

            contacts = os.path.join(wd, self.id + '_contacts.csv')
            self.__write_contact_file(contacts, self.contacts)

    def get_contacts(self):
        """Get json information for contacts.

        Returns:
            [list of dict]: Json-like structure of contacts data.
        """
        temp_contacts = filter(lambda l: l.contact_type in ('INTER', 'SELECTION_WATER', 'WATER_WATER'), self.contacts)
        contacts = ['clash', 'covalent', 'vdw_clash', 'vdw', 'proximal', 'hbond', 'weak_hbond',
                    'xbond', 'ionic', 'metal_complex', 'aromatic', 'hydrophobic', 'carbonyl',
                    'polar', 'weak_polar']
        result_bag = []
        for contact in temp_contacts:
            result_entry = {}
            result_entry['atom_bgn'] = utils.make_pymol_json(contact.start_atom)
            result_entry['atom_end'] = utils.make_pymol_json(contact.end_atom)
            result_entry['distance'] = round(np.float64(contact.distance), 2)
            result_entry['contact'] = [k for k, v in zip(contacts, contact.sifts) if v == 1]
            result_entry['interacting_entities'] = contact.contact_type

            result_bag.append(result_entry)

        return result_bag

    def minimize_hydrogens(self, minimisation_forcefield, minimisation_method, minimisation_steps):
        """Minimize structures hydrogen.

        Args:
            minimisation_forcefield (str): One of the following: MMFF94, UFF, Ghemical.
            minimisation_method (str): One of the following: ConjugateGradients, SteepestDescent,
                DistanceGeometry.
            minimisation_steps (int): Number of minimisation steps
        """
        # MINIMIZE HYDROGENS ONLY
        # `https://www.mail-archive.com/openbabel-discuss@lists.sourceforge.net/msg02216.html`
        # `https://www.mail-archive.com/openbabel-discuss@lists.sourceforge.net/msg02220.html`
        # `https://github.com/dlonie/OpenBabel-BFGS/blob/master/scripts/python/examples/minSimple.py`

        # CONSTRAIN NON-HYDROGEN ATOMS IN THE MINIMISATION
        logging.debug('Beginning hydrogen minimisation.')

        constraints = ob.OBFFConstraints()

        for ob_atom in ob.OBMolAtomIter(self.ob_mol):

            if not ob_atom.IsHydrogen():
                constraints.AddAtomConstraint(ob_atom.GetIdx())

        logging.debug('Constrained non-hydrogen atoms.')

        # INITIALISE THE FORCEFIELD
        ff = ob.OBForceField.FindForceField(minimisation_forcefield)  # MMFF94, UFF, Ghemical

        logging.debug('Initialised forcefield.')

        # TODO: MAKE LOGGING A COMMAND LINE FLAG
        ff.SetLogLevel(ob.OBFF_LOGLVL_NONE)  # OBFF_LOGLVL_LOW
        ff.SetLogToStdErr()

        if ff.Setup(self.ob_mol, constraints) == 0:  # , constraints)
            logging.warn('Could not setup the hydrogen minimisation forcefield. Skipping hydrogen minimisation.')
        else:

            # DOESN'T WORK
            # for ob_atom in ob.OBMolAtomIter(mol):
            #
            #    if not ob_atom.IsHydrogen():
            #        ff.SetFixAtom(ob_atom.GetIdx())

            if minimisation_method == 'ConjugateGradients':
                ff.ConjugateGradients(minimisation_steps)
            elif minimisation_method == 'SteepestDescent':
                ff.SteepestDescent(minimisation_steps)
            elif minimisation_method == 'DistanceGeometry':
                ff.DistanceGeometry()

            ff.GetCoordinates(self.ob_mol)

            logging.debug('Minimised hydrogens.')

    def write_hydrogenated(self, wd, input_structure):
        """Writes out structure with added H atoms.

        Args:
            wd (str): working directory
        """
        filetype = self._setup_filetype(input_structure)
        conv = ob.OBConversion()
        conv.SetInAndOutFormats(filetype, filetype)
        conv.WriteFile(self.ob_mol, os.path.join(wd, self.id + '_hydrogenated' + '.' + filetype))

        if not self.params.has_hydrogens:
            logging.debug('Wrote hydrogenated structure file. Hydrogenation was by Arpeggio using OpenBabel defaults.')

        else:
            logging.debug('Wrote hydrogenated structure file. Hydrogens were from the input file.')

    def run_arpeggio(self, user_selections, interacting_cutoff, vdw_comp, include_sequence_adjacent):
        """Run the arpergio algorithm on the structure

        Args:
            user_selections (list of Atoms): atom selection representing
                ligand or a couple of ligands
            interacting_cutoff (float): Distance cutoff for grid points
                to be \'interacting\' with the entity.
            vdw_comp (float): Compensation factor for VdW radii dependent
                interaction types.
            include_sequence_adjacent (bool): Include non-bonding
                interactions between residues that are next to each
                other in sequence
        """
        self._extend_atom_properties()
        self._ob_atom_typing()

        self._handle_hydrogens()
        logging.debug(("Determined atom explicit and implicit valences, bond orders, "
                      "atomic numbers, formal charge and number of bound hydrogens."))

        self._initialize_atom_sift()
        self._initialize_residue_sift()
        logging.debug('Initialised SIFts.')

        self._handle_chains_residues_and_breaks()
        logging.debug('Determined polypeptide residues, chain breaks, termini')  # and amide bonds.')

        self._perceive_rings()
        logging.debug('Percieved and stored rings.')

        self._perceive_amide_groups()
        logging.debug('Perceived and stored amide groups.')

        self._add_hydrogens_to_biopython()
        logging.debug('Added hydrogens to BioPython atoms.')

        self._add_atomic_radii()
        self._assign_aromatic_rings_to_residues()
        logging.debug('Assigned rings to residues.')

        self._make_selection(user_selections)
        logging.debug('Completed new NeighbourSearch.')

        self._calculate_contacts(interacting_cutoff, vdw_comp, include_sequence_adjacent)

    def write_atom_sifts(self, wd):
        """Write out per-atom SIFTS. Files: '_sifts' and '_specific_sifts'

        Args:
            wd (str): Working directory
        """
        sifts = os.path.join(wd, self.id + '_sifts.csv')
        specific_sifts = os.path.join(wd, self.id + '_specific_sifts.csv')

        sifts_content = [[utils.make_pymol_string(atom)] + atom.sift for atom in self.selection_plus]
        specific_sifts_content = [[utils.make_pymol_string(atom)]
                                  + atom.sift_inter_only
                                  + atom.sift_intra_only
                                  + atom.sift_water_only
                                  for atom in self.selection_plus]

        self.__write_atom_sifts(sifts, sifts_content)
        self.__write_atom_sifts(specific_sifts, specific_sifts_content)

    def write_binding_site_sifts(self, wd):
        """ Write out sift matching ligand and binding site (`selection_plus`)

        Args:
            wd (str): Working directory

        Raises:
            ValueError: [description]
            OBBioMatchError: [description]
            SiftMatchError: [description]
        """
        siftmatch = os.path.join(wd, self.id + '_siftmatch.csv')
        spec_siftmatch = os.path.join(wd, self.id + '_specific_siftmatch.csv')

        with open(siftmatch, 'w') as f, open(spec_siftmatch, 'w') as spec_f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            s_writer = csv.writer(spec_f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            # TODO @Lukas add headers

            for atom in self.selection_plus:
                sift_match = self._sift_match_base3(atom.potential_fsift, atom.actual_fsift)  # WHICH SIFT TO USE?

                sift_match_inter = self._sift_match_base3(atom.potential_fsift, atom.actual_fsift_inter_only)
                sift_match_intra = self._sift_match_base3(atom.potential_fsift, atom.actual_fsift_intra_only)
                sift_match_water = self._sift_match_base3(atom.potential_fsift, atom.actual_fsift_water_only)

                human_readable = self._human_sift_match(sift_match)

                writer.writerow([utils.make_pymol_string(atom)] + sift_match + [utils.int3(sift_match)] + [human_readable])
                s_writer.writerow([utils.make_pymol_string(atom)]
                                  + sift_match_inter
                                  + [self._human_sift_match(sift_match_inter)]
                                  + sift_match_intra
                                  + [self._human_sift_match(sift_match_intra)]
                                  + sift_match_water
                                  + [self._human_sift_match(sift_match_water)])

    def write_polar_matching(self, wd):
        """Write out potential polar contacts

        Args:
            wd (str): Working directory
        """
        polarmatch = os.path.join(wd, self.id + '_polarmatch.csv')
        s_polarmatch = os.path.join(wd, self.id + '_specific_polarmatch.csv')

        with open(polarmatch, 'w') as f, open(s_polarmatch, 'w') as pf:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            p_writer = csv.writer(pf, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # TODO add header
            for atom in self.selection_plus:
                writer.writerow([utils.make_pymol_string(atom)]
                                + [atom.potential_hbonds,
                                   atom.potential_polars,
                                   atom.actual_hbonds,
                                   atom.actual_polars])
                p_writer.writerow([utils.make_pymol_string(atom)]
                                  + [atom.potential_hbonds,
                                     atom.potential_polars,
                                     atom.actual_hbonds_inter_only,
                                     atom.actual_hbonds_intra_only,
                                     atom.actual_hbonds_water_only,
                                     atom.actual_polars_inter_only,
                                     atom.actual_polars_intra_only,
                                     atom.actual_polars_water_only])

    def write_residue_sifts(self, wd):
        self._calc_residue_sifts()
        residue_sifts = os.path.join(wd, self.id + '_residue_sifts.csv')
        with open(residue_sifts, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(config.RESIDUE_SIFT_HEADER)

            for residue in self.selection_plus_residues:
                output_list = [utils.make_pymol_string(residue), residue.is_polypeptide]

                for sift in (residue.sift, residue.sift_inter_only, residue.sift_intra_only, residue.sift_water_only,
                             residue.mc_sift, residue.mc_sift_inter_only, residue.mc_sift_intra_only, residue.mc_sift_water_only,
                             residue.sc_sift, residue.sc_sift_inter_only, residue.sc_sift_intra_only, residue.sc_sift_water_only,

                             residue.integer_sift, residue.integer_sift_inter_only, residue.integer_sift_intra_only, residue.integer_sift_water_only,
                             residue.mc_integer_sift, residue.mc_integer_sift_inter_only, residue.mc_integer_sift_intra_only, residue.mc_integer_sift_water_only,
                             residue.sc_integer_sift, residue.sc_integer_sift_inter_only, residue.sc_integer_sift_intra_only, residue.sc_integer_sift_water_only,

                             residue.ring_ring_inter_sift, residue.ring_atom_inter_sift, residue.atom_ring_inter_sift, residue.mc_atom_ring_inter_sift, residue.sc_atom_ring_inter_sift,
                             residue.ring_ring_inter_integer_sift, residue.ring_atom_inter_integer_sift, residue.atom_ring_inter_integer_sift, residue.mc_atom_ring_inter_integer_sift, residue.sc_atom_ring_inter_integer_sift,

                             residue.amide_ring_inter_sift, residue.ring_amide_inter_sift, residue.amide_amide_inter_sift,
                             residue.amide_ring_inter_integer_sift, residue.ring_amide_inter_integer_sift, residue.amide_amide_inter_integer_sift
                             ):
                    output_list = output_list + sift

                writer.writerow(output_list)

        # endregion

        # region private methods
    def _calc_residue_sifts(self):
        """Calculate integer sifts and flatten them to binary SIFTs to
        provide binary information on the contacts employed per residue.
        """
        for residue in self.biopython_str.get_residues():

            # WHOLE RESIDUE SIFTS
            for atom in residue.child_list:

                if hasattr(atom, 'integer_sift'):
                    residue.integer_sift = [x + y for x, y in zip(residue.integer_sift, atom.integer_sift)]

                if hasattr(atom, 'integer_sift_inter_only'):
                    residue.integer_sift_inter_only = [x + y for x, y in zip(residue.integer_sift_inter_only, atom.integer_sift_inter_only)]

                if hasattr(atom, 'integer_sift_intra_only'):
                    residue.integer_sift_intra_only = [x + y for x, y in zip(residue.integer_sift_intra_only, atom.integer_sift_intra_only)]

                if hasattr(atom, 'integer_sift_water_only'):
                    residue.integer_sift_water_only = [x + y for x, y in zip(residue.integer_sift_water_only, atom.integer_sift_water_only)]

            # FLATTEN TO BINARY SIFTS
            residue.sift = [1 if x else 0 for x in residue.integer_sift]
            residue.sift_inter_only = [1 if x else 0 for x in residue.integer_sift_inter_only]
            residue.sift_intra_only = [1 if x else 0 for x in residue.integer_sift_intra_only]
            residue.sift_water_only = [1 if x else 0 for x in residue.integer_sift_water_only]

            # MAINCHAIN/SIDECHAIN SIFTS FOR POLYPEPTIDE RESIDUES
            residue.mc_integer_sift = [0] * 15
            residue.sc_integer_sift = [0] * 15

            residue.mc_integer_sift_inter_only = [0] * 15
            residue.mc_integer_sift_intra_only = [0] * 15
            residue.mc_integer_sift_water_only = [0] * 15

            residue.sc_integer_sift_inter_only = [0] * 15
            residue.sc_integer_sift_intra_only = [0] * 15
            residue.sc_integer_sift_water_only = [0] * 15

            residue.mc_sift = [0] * 15
            residue.mc_sift_inter_only = [0] * 15
            residue.mc_sift_intra_only = [0] * 15
            residue.mc_sift_water_only = [0] * 15

            residue.sc_sift = [0] * 15
            residue.sc_sift_inter_only = [0] * 15
            residue.sc_sift_intra_only = [0] * 15
            residue.sc_sift_water_only = [0] * 15

            if residue in self.polypeptide_residues:

                for atom in residue.child_list:

                    if atom.name in config.MAINCHAIN_ATOMS:

                        if hasattr(atom, 'integer_sift'):
                            residue.mc_integer_sift = [x + y for x, y in zip(residue.mc_integer_sift, atom.integer_sift)]

                        if hasattr(atom, 'integer_sift_inter_only'):
                            residue.mc_integer_sift_inter_only = [x + y for x, y in zip(residue.mc_integer_sift_inter_only, atom.integer_sift_inter_only)]

                        if hasattr(atom, 'integer_sift_intra_only'):
                            residue.mc_integer_sift_intra_only = [x + y for x, y in zip(residue.mc_integer_sift_intra_only, atom.integer_sift_intra_only)]

                        if hasattr(atom, 'integer_sift_water_only'):
                            residue.mc_integer_sift_water_only = [x + y for x, y in zip(residue.mc_integer_sift_water_only, atom.integer_sift_water_only)]

                    else:

                        if hasattr(atom, 'integer_sift'):
                            residue.sc_integer_sift = [x + y for x, y in zip(residue.sc_integer_sift, atom.integer_sift)]

                        if hasattr(atom, 'integer_sift_inter_only'):
                            residue.sc_integer_sift_inter_only = [x + y for x, y in zip(residue.sc_integer_sift_inter_only, atom.integer_sift_inter_only)]

                        if hasattr(atom, 'integer_sift_intra_only'):
                            residue.sc_integer_sift_intra_only = [x + y for x, y in zip(residue.sc_integer_sift_intra_only, atom.integer_sift_intra_only)]

                        if hasattr(atom, 'integer_sift_water_only'):
                            residue.sc_integer_sift_water_only = [x + y for x, y in zip(residue.sc_integer_sift_water_only, atom.integer_sift_water_only)]

                # FLATTEN TO BINARY SIFTS
                residue.mc_sift = [1 if x else 0 for x in residue.mc_integer_sift]
                residue.mc_sift_inter_only = [1 if x else 0 for x in residue.mc_integer_sift_inter_only]
                residue.mc_sift_intra_only = [1 if x else 0 for x in residue.mc_integer_sift_intra_only]
                residue.mc_sift_water_only = [1 if x else 0 for x in residue.mc_integer_sift_water_only]

                residue.sc_sift = [1 if x else 0 for x in residue.sc_integer_sift]
                residue.sc_sift_inter_only = [1 if x else 0 for x in residue.sc_integer_sift_inter_only]
                residue.sc_sift_intra_only = [1 if x else 0 for x in residue.sc_integer_sift_intra_only]
                residue.sc_sift_water_only = [1 if x else 0 for x in residue.sc_integer_sift_water_only]

            # FLATTEN RING RELATED SIFTS
            residue.ring_ring_inter_sift = [1 if x else 0 for x in residue.ring_ring_inter_integer_sift]

            residue.ring_atom_inter_sift = [1 if x else 0 for x in residue.ring_atom_inter_integer_sift]
            residue.atom_ring_inter_sift = [1 if x else 0 for x in residue.atom_ring_inter_integer_sift]
            residue.mc_atom_ring_inter_sift = [1 if x else 0 for x in residue.mc_atom_ring_inter_integer_sift]
            residue.sc_atom_ring_inter_sift = [1 if x else 0 for x in residue.sc_atom_ring_inter_integer_sift]

            # FLATTEN AMIDE RELATED SIFTS
            residue.amide_ring_inter_sift = [1 if x else 0 for x in residue.amide_ring_inter_integer_sift]
            residue.ring_amide_inter_sift = [1 if x else 0 for x in residue.ring_amide_inter_integer_sift]
            residue.amide_amide_inter_sift = [1 if x else 0 for x in residue.amide_amide_inter_integer_sift]

    def __write_atom_sifts(self, path, content):
        """Routine to actually do the file writing.

        Args:
            path (str): File to store the data.
            content (list of any): Data to be written.
        """
        with open(path, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['atom',
                             'clash',
                             'covalent',
                             'vdw_clash',
                             'vdw',
                             'proximal',
                             'hbond',
                             'weak_hbond',
                             'xbond',
                             'ionic',
                             'metal_complex',
                             'aromatic',
                             'hydrophobic',
                             'carbonyl',
                             'polar',
                             'weak_polar',
                             'interacting_entities'
                             ])
            [writer.writerow(i) for i in content]

    def __write_contact_file(self, path, contacts):
        """Write out Arperggio contacts

        Args:
            path (str): Path where the data will be exported.
            contacts (list of AtomContacts): List of contacts to be exported.
        """
        with open(path, 'w') as fo:
            writer = csv.writer(fo, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['atom_bgn',
                             'atom_end',
                             'distance',
                             'clash',
                             'covalent',
                             'vdw_clash',
                             'vdw',
                             'proximal',
                             'hbond',
                             'weak_hbond',
                             'xbond',
                             'ionic',
                             'metal_complex',
                             'aromatic',
                             'hydrophobic',
                             'carbonyl',
                             'polar',
                             'weak_polar',
                             'interacting_entities'
                             ])
            for contact in contacts:
                writer.writerow([
                    utils.make_pymol_string(contact.start_atom),
                    utils.make_pymol_string(contact.end_atom),
                    contact.distance]
                    + contact.sifts
                    + [contact.contact_type]
                )

    def __get_contact_type(self, atom_bgn, atom_end, selection_set):
        """Based on the origin of the interacting guys inferre, whether
        or not this interaction is happening inside/outside the selection
        etc.

        Possible values are:
            * INTRA_NON_SELECTION
            * INTRA_SELECTION
            * INTER
            * SELECTION_WATER
            * NON_SELECTION_WATER
            * WATER_WATER

        Args:
            atom_bgn (Atom): First atom of the interaction
            atom_end (Atom): Second atom of the interaction
            selection_set (list of Atom): User selection

        Raises:
            ValueError: If contact type could not be assigned.

        Returns:
            str: Type of the contact for a given atom.
        """
        contact_type = ''

        if atom_bgn not in selection_set and atom_end not in selection_set:
            contact_type = 'INTRA_NON_SELECTION'

        if atom_bgn in selection_set and atom_end in selection_set:
            contact_type = 'INTRA_SELECTION'

        if (atom_bgn in selection_set and atom_end not in selection_set) or (atom_end in selection_set and atom_bgn not in selection_set):
            contact_type = 'INTER'

        if (atom_bgn in selection_set and atom_end.get_full_id()[3][0] == 'W') or (atom_end in selection_set and atom_bgn.get_full_id()[3][0] == 'W'):
            contact_type = 'SELECTION_WATER'

        if (atom_bgn not in selection_set and atom_end.get_full_id()[3][0] == 'W') or (atom_end not in selection_set and atom_bgn.get_full_id()[3][0] == 'W'):
            contact_type = 'NON_SELECTION_WATER'

        if atom_bgn.get_full_id()[3][0] == 'W' and atom_end.get_full_id()[3][0] == 'W':
            contact_type = 'WATER_WATER'

        if not contact_type:
            logging.error('Could not assign a contact type for {}:{}'.format(atom_bgn, atom_end))
            raise ValueError('Could not assign a contact type for {}:{}'.format(atom_bgn, atom_end))

        return contact_type

    def _calculate_contacts(self, interacting_cutoff, vdw_comp_factor, include_sequence_adjacent):
        """Calculate the actual molecular contacts and their type.

        Args:
            interacting_cutoff (float): Distance cutoff for grid points
                to be \'interacting\' with the entity.
            vdw_comp (float): Compensation factor for VdW radii dependent
                interaction types.
            include_sequence_adjacent (bool): Include non-bonding
                interactions between residues that are next to each
                other in sequence
        """
        for atom_bgn, atom_end in self.ns.search_all(interacting_cutoff):

            selection_set = set(self.selection)
            # IGNORE ANY HYDROGENS FOR THESE CONTACTS
            # IF HYDROGENS ARE PRESENT IN THE BIOPYTHON STRUCTURE FOR ANY REASON
            if atom_bgn.element.strip() == 'H' or atom_end.element.strip() == 'H':
                continue

            contact_type = self.__get_contact_type(atom_bgn, atom_end, selection_set)

            sum_cov_radii = atom_bgn.cov_radius + atom_end.cov_radius
            sum_vdw_radii = atom_bgn.vdw_radius + atom_end.vdw_radius

            ob_atom_bgn = self.ob_mol.GetAtomById(self.bio_to_ob[atom_bgn])
            ob_atom_end = self.ob_mol.GetAtomById(self.bio_to_ob[atom_end])

            SIFt = [0] * 15

            # IGNORE INTRA-RESIDUE CONTACTS
            res_bgn = atom_bgn.get_parent()
            res_end = atom_end.get_parent()

            if res_bgn is res_end:
                continue

            # IGNORE CONTACTS TO SEQUENCE-ADJACENT RESIDUES (BY DEFAULT)
            if not include_sequence_adjacent:
                if res_end.is_polypeptide and res_end.is_polypeptide:

                    if hasattr(res_bgn, 'prev_residue') and hasattr(res_bgn, 'next_residue') and \
                            hasattr(res_end, 'prev_residue') and hasattr(res_end, 'next_residue'):

                        if res_bgn.next_residue is res_end or res_bgn.prev_residue is res_end or \
                                res_end.next_residue is res_bgn or res_end.prev_residue is res_bgn:
                            continue

            # print contact_type

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
            elif distance <= sum_vdw_radii + vdw_comp_factor:
                SIFt[3] = 1

            # PROXIMAL
            else:
                SIFt[4] = 1

            # METAL COMPLEX
            # CAN BE COVALENT SO GO HERE
            if distance <= config.CONTACT_TYPES['metal']['distance']:

                if 'hbond acceptor' in atom_bgn.atom_types and atom_end.is_metal:
                    SIFt[9] = 1

                elif 'hbond acceptor' in atom_end.atom_types and atom_bgn.is_metal:
                    SIFt[9] = 1

            # FEATURE CONTACTS
            if not any(SIFt[:1]) and distance <= config.CONTACT_TYPES_DIST_MAX:

                # HBOND

                # NO NEED TO USE HYDROGENS FOR WATERS
                if atom_bgn.get_full_id()[3][0] == 'W' and distance <= (sum_vdw_radii + vdw_comp_factor):
                    if 'hbond acceptor' in atom_end.atom_types or 'hbond donor' in atom_end.atom_types:
                        SIFt[5] = 1
                        SIFt[13] = 1

                elif atom_end.get_full_id()[3][0] == 'W' and distance <= (sum_vdw_radii + vdw_comp_factor):
                    if 'hbond acceptor' in atom_bgn.atom_types or 'hbond donor' in atom_bgn.atom_types:
                        SIFt[5] = 1
                        SIFt[13] = 1

                else:

                    # ATOM_BGN IS DONOR
                    if 'hbond donor' in atom_bgn.atom_types and 'hbond acceptor' in atom_end.atom_types:

                        SIFt[5] = self._is_hbond(atom_bgn, atom_end, vdw_comp_factor)

                        # CHECK DISTANCE FOR POLARS
                        if distance <= config.CONTACT_TYPES["hbond"]["polar distance"]:
                            SIFt[13] = 1

                    # ATOM_END IS DONOR
                    elif 'hbond donor' in atom_end.atom_types and 'hbond acceptor' in atom_bgn.atom_types:

                        SIFt[5] = self._is_hbond(atom_end, atom_bgn, vdw_comp_factor)

                        # CHECK DISTANCE FOR POLARS
                        if distance <= config.CONTACT_TYPES["hbond"]["polar distance"]:
                            SIFt[13] = 1

                # UPDATE ATOM HBOND/POLAR COUNTS
                if SIFt[5]:
                    atom_bgn.actual_hbonds = atom_bgn.actual_hbonds + 1
                    atom_end.actual_hbonds = atom_end.actual_hbonds + 1

                    if 'INTRA' in contact_type:
                        atom_bgn.actual_hbonds_intra_only = atom_bgn.actual_hbonds_intra_only + 1
                        atom_end.actual_hbonds_intra_only = atom_end.actual_hbonds_intra_only + 1

                    elif 'INTER' in contact_type:
                        atom_bgn.actual_hbonds_inter_only = atom_bgn.actual_hbonds_inter_only + 1
                        atom_end.actual_hbonds_inter_only = atom_end.actual_hbonds_inter_only + 1

                    elif 'WATER' in contact_type:
                        atom_bgn.actual_hbonds_water_only = atom_bgn.actual_hbonds_water_only + 1
                        atom_end.actual_hbonds_water_only = atom_end.actual_hbonds_water_only + 1

                if SIFt[13]:
                    atom_bgn.actual_polars = atom_bgn.actual_polars + 1
                    atom_end.actual_polars = atom_end.actual_polars + 1

                    if 'INTRA' in contact_type:
                        atom_bgn.actual_polars_intra_only = atom_bgn.actual_polars_intra_only + 1
                        atom_end.actual_polars_intra_only = atom_end.actual_polars_intra_only + 1

                    elif 'INTER' in contact_type:
                        atom_bgn.actual_polars_inter_only = atom_bgn.actual_polars_inter_only + 1
                        atom_end.actual_polars_inter_only = atom_end.actual_polars_inter_only + 1

                    elif 'WATER' in contact_type:
                        atom_bgn.actual_polars_water_only = atom_bgn.actual_polars_water_only + 1
                        atom_end.actual_polars_water_only = atom_end.actual_polars_water_only + 1

                # WEAK HBOND

                # ATOM_BGN IS ACCEPTOR, ATOM_END IS CARBON
                if 'hbond acceptor' in atom_bgn.atom_types and 'weak hbond donor' in atom_end.atom_types:
                    SIFt[6] = self._is_weak_hbond(atom_end, atom_bgn, vdw_comp_factor)

                    # CHECK DISTANCE FOR WEAK POLARS
                    if distance <= config.CONTACT_TYPES["weak hbond"]["weak polar distance"]:
                        SIFt[14] = 1

                # ATOM_BGN IS CARBON, ATOM_END IS ACCEPTOR
                if 'weak hbond donor' in atom_bgn.atom_types and 'hbond acceptor' in atom_end.atom_types:
                    SIFt[6] = self._is_weak_hbond(atom_bgn, atom_end, vdw_comp_factor)

                    # CHECK DISTANCE FOR WEAK POLARS
                    if distance <= config.CONTACT_TYPES["weak hbond"]["weak polar distance"]:
                        SIFt[14] = 1

                # ATOM_BGN IS HALOGEN WEAK ACCEPTOR
                if 'weak hbond acceptor' in atom_bgn.atom_types and atom_bgn.is_halogen and ('hbond donor' in atom_end.atom_types or 'weak hbond donor' in atom_end.atom_types):
                    SIFt[6] = self._is_halogen_weak_hbond(atom_end, atom_bgn, self.ob_mol, vdw_comp_factor)

                    # CHECK DISTANCE FOR WEAK POLARS
                    if distance <= config.CONTACT_TYPES["weak hbond"]["weak polar distance"]:
                        SIFt[14] = 1

                # ATOM END IS HALOGEN WEAK ACCEPTOR
                if 'weak hbond acceptor' in atom_end.atom_types and atom_end.is_halogen and ('hbond donor' in atom_bgn.atom_types or 'weak hbond donor' in atom_bgn.atom_types):
                    SIFt[6] = self._is_halogen_weak_hbond(atom_bgn, atom_end, self.ob_mol, vdw_comp_factor)

                    # CHECK DISTANCE FOR WEAK POLARS
                    if distance <= config.CONTACT_TYPES["weak hbond"]["weak polar distance"]:
                        SIFt[14] = 1

                # XBOND
                if distance <= sum_vdw_radii + vdw_comp_factor:

                    if 'xbond donor' in atom_bgn.atom_types and 'xbond acceptor' in atom_end.atom_types:
                        SIFt[7] = self._is_xbond(atom_bgn, atom_end, self.ob_mol)

                    elif 'xbond donor' in atom_end.atom_types and 'xbond acceptor' in atom_bgn.atom_types:
                        SIFt[7] = self._is_xbond(atom_end, atom_bgn, self.ob_mol)

                # IONIC
                if distance <= config.CONTACT_TYPES['ionic']['distance']:

                    if 'pos ionisable' in atom_bgn.atom_types and 'neg ionisable' in atom_end.atom_types:
                        SIFt[8] = 1

                    elif 'neg ionisable' in atom_bgn.atom_types and 'pos ionisable' in atom_end.atom_types:
                        SIFt[8] = 1

                # CARBONYL
                if distance <= config.CONTACT_TYPES['carbonyl']['distance']:

                    if 'carbonyl oxygen' in atom_bgn.atom_types and 'carbonyl carbon' in atom_end.atom_types:
                        SIFt[12] = 1

                    elif 'carbonyl oxygen' in atom_end.atom_types and 'carbonyl carbon' in atom_bgn.atom_types:
                        SIFt[12] = 1

                # AROMATIC
                if 'aromatic' in atom_bgn.atom_types and 'aromatic' in atom_end.atom_types and distance <= config.CONTACT_TYPES['aromatic']['distance']:
                    SIFt[10] = 1

                # HYDROPHOBIC
                if 'hydrophobe' in atom_bgn.atom_types and 'hydrophobe' in atom_end.atom_types and distance <= config.CONTACT_TYPES['hydrophobic']['distance']:
                    SIFt[11] = 1

            # UPDATE ATOM INTEGER SIFTS
            self._update_atom_integer_sift(atom_bgn, SIFt, contact_type)
            self._update_atom_integer_sift(atom_end, SIFt, contact_type)

            # UPDATE ATOM SIFTS
            self._update_atom_sift(atom_bgn, SIFt, contact_type)
            self._update_atom_sift(atom_end, SIFt, contact_type)

            # UPDATE ATOM FEATURE SIFts
            fsift = SIFt[5:]
            self._update_atom_fsift(atom_bgn, fsift, contact_type)
            self._update_atom_fsift(atom_end, fsift, contact_type)

            self.contacts.append(AtomContacts(start_atom=atom_bgn, end_atom=atom_end, sifts=SIFt, contact_type=contact_type, distance=distance))

    def _make_selection(self, selections):
        """Select residues to calculate interactions

        Args:
            selections (list of str): User structure selection in the
                following formars:
                /<chain_id>/<res_num>[<ins_code>]/<atom_name>
                or RESNAME:<het_id>
        """
        entity = list(self.s_atoms)[:]
        selection = entity if len(selections) == 0 else utils.selection_parser(selections, entity)
        selection_ring_ids = list(self.biopython_str.rings)
        selection_amide_ids = list(self.biopython_str .amides)

        if len(selection) == 0:
            logging.error('Selection was empty.')
            raise AttributeError('Selection must not be empty.')

        logging.debug('Made selection.')
        selection_set = set(selection)

        # EXPAND THE SELECTION TO INCLUDE THE BINDING SITE
        selection_plus = set(selection)
        selection_plus_residues = set([x.get_parent() for x in selection_plus])
        selection_plus_ring_ids = set(selection_ring_ids)
        selection_plus_amide_ids = set(selection_amide_ids)

        # TODO @Lukas FIX if args.selection:

        # GET LIST OF RESIDUES IN THE SELECTION PLUS BINDING SITE
        selection_residues = set([x.get_parent() for x in selection])

        # MAKE A SET OF ALL RING IDS ASSOCIATED WITH THE SELECTION AND BINDING SITE
        selection_ring_ids = set([x for x in self.biopython_str.rings if self.biopython_str.rings[x]['residue'] in selection_residues])
        selection_amide_ids = set([x for x in self.biopython_str.amides if self.biopython_str.amides[x]['residue'] in selection_residues])

        # EXPAND THE SELECTION TO THE BINDING SITE
        for atom_bgn, atom_end in self.ns.search_all(6.0):

            if atom_bgn in selection_set or atom_end in selection_set:
                selection_plus.add(atom_bgn)
                selection_plus.add(atom_end)

        selection_plus = list(selection_plus)

        logging.debug('Expanded to binding site.')

        # GET LIST OF RESIDUES IN THE SELECTION PLUS BINDING SITE
        selection_plus_residues = set([x.get_parent() for x in selection_plus])

        # MAKE A SET OF ALL RING IDS ASSOCIATED WITH THE SELECTION AND BINDING SITE
        selection_plus_ring_ids = set([x for x in self.biopython_str.rings if self.biopython_str.rings[x]['residue'] in selection_plus_residues])

        # MAKE A SET OF ALL AMIDE IDS ASSOCIATED WITH THE SELECTION AND BINDING SITE
        selection_plus_amide_ids = set([x for x in self.biopython_str.amides if self.biopython_str.amides[x]['residue'] in selection_plus_residues])

        logging.debug('Flagged selection rings.')

        # NEW NEIGHBOURSEARCH
        self.ns = NeighborSearch(selection_plus)

        self.selection = selection
        self.selection_plus = selection_plus
        self.selection_plus_residues = selection_plus_residues

    def _assign_aromatic_rings_to_residues(self):
        # NEIGHBORSEARCH
        self.ns = NeighborSearch(self.s_atoms)

        logging.debug('Completed NeighborSearch.')

        # ASSIGN AROMATIC RINGS TO RESIDUES
        for ring_id in self.biopython_str.rings:

            ring_centroid = self.biopython_str.rings[ring_id]['center']
            atoms_near_ring = self.ns.search(ring_centroid, 3.0)  # MORE THAN REASONABLE FOR PICKING UP THE RESIDUE THE RING IS IN
            # UNLESS THERE ARE BIG PROBLEMS IN THE PDB, IN WHICH CASE THE RING
            # RESIDUE ASSIGNMENT IS THE LEAST OF CONCERNS ;)

            closest_atom = (None, None)  # (ATOM, DISTANCE)

            for nearby_atom in atoms_near_ring:

                distance = np.linalg.norm(nearby_atom.coord - ring_centroid)

                if closest_atom[1] is None or distance < closest_atom[1]:
                    closest_atom = (nearby_atom, distance)

            if closest_atom[0] is None:

                logging.warn('Residue assignment was not possible for ring {}.'.format(ring_id))
                self.biopython_str.rings[ring_id]['residue'] = None

            else:

                # ASSIGN RING TO RESIDUE, STORE RESIDUE IN THE RING
                closest_residue = closest_atom[0].get_parent()
                self.biopython_str.rings[ring_id]['residue'] = closest_residue
                self.biopython_str.rings[ring_id]['residue_shortest_distance'] = closest_atom[1]

                # ADD RING ID TO THE RESIDUE AS WELL
                if not hasattr(closest_residue, 'rings'):
                    closest_residue.rings = []

                closest_residue.rings.append(ring_id)

    def _add_atomic_radii(self):
        """Add atomic radii to the structure atoms.
        """
        # ADD VDW RADII TO ENTITY ATOMS
        # USING OPENBABEL VDW RADII
        for atom in self.s_atoms:
            atom.vdw_radius = ob.etab.GetVdwRad(self.ob_mol.GetAtomById(self.bio_to_ob[atom]).GetAtomicNum())

        logging.debug('Added VdW radii.')

        # ADD COVALENT RADII TO ENTITY ATOMS
        # USING OPENBABEL VDW RADII
        for atom in self.s_atoms:
            atom.cov_radius = ob.etab.GetCovalentRad(self.ob_mol.GetAtomById(self.bio_to_ob[atom]).GetAtomicNum())

        logging.debug('Added covalent radii.')

    def _add_hydrogens_to_biopython(self):
        """Add hydrogens atom to the biopython structure
        Note: TODO @Lukas Not sure why it is here.
        """
        for atom in ob.OBMolAtomIter(self.ob_mol):
            # IF THE ATOM HAS EXPLICIT HYDROGENS
            if atom.ExplicitHydrogenCount() > 0:

                biopython_atom = self.ob_to_bio[atom.GetId()]

                # GET THE BONDED ATOMS OF THE OBATOM
                for atom_neighbour in ob.OBAtomAtomIter(atom):
                    if atom_neighbour.IsHydrogen():

                        # APPEND THE HYDROGEN COORDINATES TO THE BIOPYTHON ATOM 'h_coords' ATTRIBUTE
                        biopython_atom.h_coords.append(np.array([atom_neighbour.GetX(), atom_neighbour.GetY(), atom_neighbour.GetZ()]))

    def _perceive_amide_groups(self):
        # DETECT AMIDE GROUPS
        # AMIDES FOR AMIDE-RELATED NON-BONDING INTERACTIONS
        self.biopython_str.amides = collections.OrderedDict()

        # GET OPENBABEL ATOM MATCHES TO THE SMARTS PATTERN
        ob_smart = ob.OBSmartsPattern()
        ob_smart.Init(config.AMIDE_SMARTS)
        ob_smart.Match(self.ob_mol)

        matches = [x for x in ob_smart.GetMapList()]

        for e, match in enumerate(matches):

            ob_match = [self.ob_mol.GetAtom(x) for x in match]
            bio_match = [self.ob_to_bio[x.GetId()] for x in ob_match]

            # CHECK FOR EXPECTED BEHAVIOUR
            assert len(bio_match) == 4
            assert bio_match[0].element == 'N'
            assert bio_match[1].element == 'C'  # SHOULD BE BACKBONE C WHEN IN PROTEIN MAINCHAIN
            assert bio_match[2].element == 'O'
            assert bio_match[3].element == 'C'  # SHOULD BE C-ALPHA WHEN IN PROTEIN MAINCHAIN

            # ASSIGN GROUP TO A RESIDUE
            bio_match_residues = [x.get_parent() for x in bio_match]

            # USE THE RESIDUE OF THE MAJORITY OF ATOMS
            # `http://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list`
            group_residue = max(bio_match_residues, key=bio_match_residues.count)

            # GET AMIDE BOND CENTROID
            # DETERMINED AS CENTRE OF MASS OF C-N (OR C-O-N?)
            con = np.array([bio_match[1].coord, bio_match[2].coord, bio_match[0].coord])  # C-O-N
            cn = np.array([bio_match[1].coord, bio_match[0].coord])  # C-N
            amide_centroid = con.sum(0) / float(len(con))
            bond_centroid = cn.sum(0) / float(len(cn))

            # GET AMIDE PLANE WITH SVD
            # `http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054621.html`

            cog = con - amide_centroid
            u, s_, vh = np.linalg.svd(cog)
            v = vh.conj().transpose()
            a, b, c = v[:, -1]
            #d = 0  # :S

            normal = np.array([a, b, c])
            normal_opp = -normal

            # STORE AMIDE GROUPS
            self.biopython_str.amides[e] = {
                'amide_id': e,
                'center': bond_centroid,  # amide_centroid,
                'normal': normal,
                'normal_opp': normal_opp,
                'atoms': bio_match,
                'residue': group_residue
            }

    def _handle_chains_residues_and_breaks(self):
        """Detection of polypeptides, residues, chains breaks and
        C/N termin. Create underlying data structures and such.
        """

        # DETECT POLYPEPTIDES, RESIDUES, CHAIN BREAKS AND TERMINI
        ppb = PPBuilder()
        polypeptides = ppb.build_peptides(self.biopython_str, aa_only=False)

        # CHAIN BREAKS AND TERMINI

        # MAKE DATA STRUCTURES FOR CHAIN POLYPEPTIDES
        chain_ids = set([x.id for x in self.biopython_str.get_chains()])
        chain_pieces = collections.OrderedDict()
        chain_polypeptides = collections.OrderedDict()
        chain_break_residues = collections.OrderedDict()
        chain_termini = collections.OrderedDict()
        # chain_sequences = OrderedDict()

        for chain_id in chain_ids:
            chain_pieces[chain_id] = 0
            chain_break_residues[chain_id] = []
            chain_polypeptides[chain_id] = []

        # GET THE CHAIN_ID(S) ASSOCIATED WITH EACH POLYPEPTIDE
        polypeptide_chain_id_sets = [set([k.get_parent().id for k in x]) for x in polypeptides]

        for e, polypeptide_chain_id_set in enumerate(polypeptide_chain_id_sets):

            # WARN IF NOT JUST ONE CHAIN ID ASSOCIATED WITH THE POLYPEPTIDE
            if len(polypeptide_chain_id_set) != 1:
                logging.warn('A polypeptide had {} chains associated with it: {}'.format(len(polypeptide_chain_id_set),
                                                                                         polypeptide_chain_id_set))

            for polypeptide_chain_id in polypeptide_chain_id_set:
                chain_pieces[polypeptide_chain_id] = chain_pieces[polypeptide_chain_id] + 1

                # ADD FIRST AND LAST RESIDUE TO THE CHAIN BREAK RESIDUES (POLYPEPTIDE TERMINAL RESIDUES)
                chain_break_residues[polypeptide_chain_id] = chain_break_residues[polypeptide_chain_id] + [polypeptides[e][0], polypeptides[e][-1]]
                chain_polypeptides[polypeptide_chain_id] = chain_polypeptides[polypeptide_chain_id] + [polypeptides[e]]

        # CHAIN BREAKS AND TERMINI
        for chain_id in chain_break_residues:

            try:
                # GET FIRST AND LAST ("GENUINE") TERMINI
                chain_termini[chain_id] = [chain_break_residues[chain_id][0], chain_break_residues[chain_id][-1]]
            except IndexError:
                logging.warn('Chain termini could not be determined for chain {}. It may not be a polypeptide chain.'.format(chain_id))

            try:
                # POP OUT THE FIRST AND LAST RESIDUES FROM THE CHAIN BREAK RESIDUES
                # TO REMOVE THE GENUINE TERMINI
                chain_break_residues[chain_id] = chain_break_residues[chain_id][1:-1]
            except IndexError:
                logging.warn('Chain termini could not be extracted from breaks for chain {}. It may not be a polypeptide chain.'.format(chain_id))

        all_chain_break_residues = []
        all_terminal_residues = []

        try:
            all_chain_break_residues = reduce(operator.add, list(chain_break_residues.values()))
        except TypeError:
            pass

        try:
            all_terminal_residues = reduce(operator.add, list(chain_termini.values()))
        except TypeError:
            pass

        # POLYPEPTIDE RESIDUES
        self.polypeptide_residues = set([])

        for pp in polypeptides:

            last_residue = None

            for residue in pp:

                # FLAG AS POLYPEPTIDE
                self.polypeptide_residues.add(residue)
                residue.is_polypeptide = True

                # FLAG IF CHAIN BREAK OR TERMINAL
                residue.is_chain_break = False
                residue.is_terminal = False
                residue.is_terminal_or_break = False

                if residue in all_chain_break_residues:
                    residue.is_chain_break = True

                if residue in all_terminal_residues:
                    residue.is_terminal = True

                residue.is_terminal_or_break = residue.is_terminal or residue.is_chain_break

                # DETERMINE PRECEEDING AND NEXT RESIDUES IN THE SEQUENCE
                residue.prev_residue = None
                residue.next_residue = None

                residue.prev_residue = last_residue

                if last_residue:
                    last_residue.next_residue = residue

                last_residue = residue

    def _perceive_rings(self):
        """Perceive aromatic rings
        """

        self.biopython_str.rings = collections.OrderedDict()

        for e, ob_ring in enumerate(self.ob_mol.GetSSSR()):

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
            self.biopython_str.rings[e] = {'ring_id': e,
                                           'center': center,
                                           'normal': normal,
                                           'normal_opp': normal_opp,
                                           'atoms': [],
                                           'ob_atom_ids': []}

            # GET RING ATOMS AND STORE
            for ob_atom in ob.OBMolAtomIter(self.ob_mol):

                if ob_ring.IsMember(ob_atom):

                    self.biopython_str.rings[e]['atoms'].append(self.ob_to_bio[ob_atom.GetId()])
                    self.biopython_str.rings[e]['ob_atom_ids'].append(ob_atom.GetId())

    def _initialize_atom_sift(self):
        """Initialise the following fields:
            * atom full sift
            * atom feature sift
            * determine atom potential feature sift.
            * atom number of potential hbonds/polar
            * no. actual hbonds/polar

        Details:
            5: 0: HBOND
            6: 1: WEAK_HBOND
            7: 2: HALOGEN_BOND
            8: 3: IONIC
            9: 4: METAL_COMPLEX
            10: 5: AROMATIC
            11: 6: HYDROPHOBIC
            12: 7: CARBONYL

            # 8: POLAR - H-BONDS WITHOUT ANGLES
            # 9: WEAK POLAR - WEAK H-BONDS WITHOUT ANGLES
        """
        for atom in self.s_atoms:

            # INTEGER SIFTS
            atom.integer_sift = [0] * 15

            atom.integer_sift_inter_only = [0] * 15
            atom.integer_sift_intra_only = [0] * 15
            atom.integer_sift_water_only = [0] * 15

            # BINARY SIFTS
            atom.sift = [0] * 15

            atom.sift_inter_only = [0] * 15
            atom.sift_intra_only = [0] * 15
            atom.sift_water_only = [0] * 15

            atom.potential_fsift = [0] * 10

            atom.actual_fsift = [0] * 10
            atom.actual_fsift_inter_only = [0] * 10

            atom.actual_fsift_intra_only = [0] * 10
            atom.actual_fsift_water_only = [0] * 10

            atom.potential_hbonds = 0
            atom.potential_polars = 0

            atom.actual_hbonds = 0
            atom.actual_polars = 0

            atom.actual_hbonds_inter_only = 0
            atom.actual_polars_inter_only = 0

            atom.actual_hbonds_intra_only = 0
            atom.actual_polars_intra_only = 0

            atom.actual_hbonds_water_only = 0
            atom.actual_polars_water_only = 0

            # ATOM POTENTIAL FEATURE SIFT
            # 0: HBOND
            if 'hbond acceptor' in atom.atom_types or 'hbond donor' in atom.atom_types:
                atom.potential_fsift[0] = 1

                if 'hbond acceptor' in atom.atom_types:

                    # NUMBER OF LONE PAIRS
                    try:
                        lone_pairs = config.VALENCE[atom.atomic_number] - atom.bond_order - atom.formal_charge
                    except AttributeError:
                        print(atom)

                    if lone_pairs != 0:
                        lone_pairs = lone_pairs / 2

                    atom.potential_hbonds = atom.potential_hbonds + lone_pairs
                    atom.potential_polars = atom.potential_polars + lone_pairs

                if 'hbond donor' in atom.atom_types:
                    atom.potential_hbonds = atom.potential_hbonds + atom.num_hydrogens
                    atom.potential_polars = atom.potential_polars + atom.num_hydrogens

            # 1: WEAK HBOND
            if 'weak hbond acceptor' in atom.atom_types or 'weak hbond donor' in atom.atom_types or 'hbond donor' in atom.atom_types or 'hbond acceptor' in atom.atom_types or atom.is_halogen:
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

            # 8: POLAR
            if 'hbond acceptor' in atom.atom_types or 'hbond donor' in atom.atom_types:
                atom.potential_fsift[8] = 1

            # 9: WEAK POLAR
            if 'weak hbond acceptor' in atom.atom_types or 'weak hbond donor' in atom.atom_types or 'hbond donor' in atom.atom_types or 'hbond acceptor' in atom.atom_types or atom.is_halogen:
                atom.potential_fsift[9] = 1

    def _initialize_residue_sift(self):
        """Initialize residue SIFT.
        """
        for residue in self.biopython_str.get_residues():

            # INITIALISE POLYPEPTIDE FLAG
            residue.is_polypeptide = False

            # INTEGER SIFTS
            residue.integer_sift = [0] * 15

            residue.integer_sift_inter_only = [0] * 15
            residue.integer_sift_intra_only = [0] * 15
            residue.integer_sift_water_only = [0] * 15

            # RING-RING SIFTS
            residue.ring_ring_inter_integer_sift = [0] * 9

            # ATOM-RING SIFTS
            residue.ring_atom_inter_integer_sift = [0] * 5
            residue.atom_ring_inter_integer_sift = [0] * 5
            residue.mc_atom_ring_inter_integer_sift = [0] * 5
            residue.sc_atom_ring_inter_integer_sift = [0] * 5

            # AMIDE-RING
            residue.amide_ring_inter_integer_sift = [0]
            residue.ring_amide_inter_integer_sift = [0]

            # AMIDE-AMIDE
            residue.amide_amide_inter_integer_sift = [0]

    def _handle_hydrogens(self):
        """Determine atom valences and explicit hydrogen counts.
        """
        for ob_atom in ob.OBMolAtomIter(self.ob_mol):
            if not self.input_has_hydrogens:
                if ob_atom.IsHydrogen():
                    continue

            # `http://openbabel.org/api/2.3/classOpenBabel_1_1OBAtom.shtml`
            # CURRENT NUMBER OF EXPLICIT CONNECTIONS
            valence = ob_atom.GetValence()

            # MAXIMUM NUMBER OF CONNECTIONS EXPECTED
            implicit_valence = ob_atom.GetImplicitValence()

            # BOND ORDER
            bond_order = ob_atom.BOSum()

            # NUMBER OF BOUND HYDROGENS
            num_hydrogens = ob_atom.ExplicitHydrogenCount()

            # ELEMENT NUMBER
            atomic_number = ob_atom.GetAtomicNum()

            # FORMAL CHARGE
            formal_charge = ob_atom.GetFormalCharge()

            bio_atom = self.ob_to_bio[ob_atom.GetId()]

            bio_atom.valence = valence
            bio_atom.implicit_valence = implicit_valence
            bio_atom.num_hydrogens = num_hydrogens
            bio_atom.bond_order = bond_order
            bio_atom.atomic_number = atomic_number
            bio_atom.formal_charge = formal_charge

    def _ob_atom_typing(self):
        """Atom type structure using openbabel.
        """
        for atom_type, smartsdict in list(config.ATOM_TYPES.items()):

            # logging.debug('Typing: {}'.format(atom_type))

            # FOR EACH ATOM TYPE SMARTS STRING
            for smarts in list(smartsdict.values()):

                # logging.debug('Smarts: {}'.format(smarts))

                # GET OPENBABEL ATOM MATCHES TO THE SMARTS PATTERN
                ob_smart = ob.OBSmartsPattern()
                ob_smart.Init(str(smarts))

                # logging.debug('Initialised for: {}'.format(smarts))

                ob_smart.Match(self.ob_mol)

                # logging.debug('Matched for: {}'.format(smarts))

                matches = [x for x in ob_smart.GetMapList()]

                # logging.debug('List comp matches: {}'.format(smarts))

                if matches:

                    # REDUCE TO A SINGLE LIST
                    matches = set(reduce(operator.add, matches))

                    # logging.debug('Set reduce matches: {}'.format(smarts))

                    for match in matches:

                        atom = self.ob_mol.GetAtom(match)
                        self.ob_to_bio[atom.GetId()].atom_types.add(atom_type)

        # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
        for atom in (x for x in self.s_atoms if x.get_full_id()[3][0] == 'W'):
            atom.atom_types.add('hbond acceptor')
            atom.atom_types.add('hbond donor')

        # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
        for residue in self.biopython_str.get_residues():

            if residue.resname in config.STD_RES:

                for atom in residue.child_list:

                    # REMOVE TYPES IF ALREADY ASSIGNED FROM SMARTS
                    for atom_type in list(config.PROT_ATOM_TYPES.keys()):
                        atom.atom_types.discard(atom_type)

                    # ADD ATOM TYPES FROM DICTIONARY
                    for atom_type, atom_ids in config.PROT_ATOM_TYPES.items():

                        atom_id = residue.resname.strip() + atom.name.strip()

                        if atom_id in atom_ids:
                            atom.atom_types.add(atom_type)

    def _extend_atom_properties(self):
        for atom in self.s_atoms:
            atom.atom_types = set([])
            atom.h_coords = []

            atom.is_metal = atom.element.upper() in config.METALS
            atom.is_halogen = atom.element.upper() in config.HALOGENS

    def _establish_structure_mappping(self):
        """Maps biopython atoms to openbabel ones and vice-versa.

        Raises:
            OBBioMatchError: If we cant match an OB atom to a biopython
        """
        # FIRST MAP PDB SERIAL NUMBERS TO BIOPYTHON ATOMS FOR SPEED LATER
        # THIS AVOIDS LOOPING THROUGH `s_atoms` MANY TIMES
        serial_to_bio = {x.serial_number: x for x in self.s_atoms}
        for ob_atom in ob.OBMolAtomIter(self.ob_mol):
            serial = ob_atom.GetResidue().GetSerialNum(ob_atom)

            # MATCH TO THE BIOPYTHON ATOM BY SERIAL NUMBER
            try:
                biopython_atom = serial_to_bio[serial]

            except KeyError:
                raise OBBioMatchError(serial)

            # `Id` IS A UNIQUE AND STABLE ID IN OPENBABEL
            # CAN RECOVER THE ATOM WITH `mol.GetAtomById(id)`
            self.ob_to_bio[ob_atom.GetId()] = biopython_atom
            self.bio_to_ob[biopython_atom] = ob_atom.GetId()

        logging.debug('Mapped OB to BioPython atoms and vice-versa.')

    def _read_in_biopython(self, path):
        """Reads in molecule using Biopython

        Args:
            path (str): path to the input molecule

        Returns:
            BioPython.Structure: Parsed biopython protein structure
        """
        extension = self._setup_filetype(path)

        s = (PDBParser().get_structure('structure', path)
             if extension == 'pdb'
             else protein_reader.read_mmcif_to_biopython(path))

        logging.debug('Loaded PDB structure (BioPython)')

        return s

    def _read_openbabel(self, path):
        """Reads in molecule using openbabel.

        Args:
            path (str): Path to the structure

        Returns:
            openbabel.OBMol: Parsed openbabel protein structure
        """
        filetype = self._setup_filetype(path)
        if filetype == 'pdb':
            ob_conv = ob.OBConversion()
            ob_conv.SetInFormat(filetype)
            mol = ob.OBMol()
            ob_conv.ReadFile(mol, path)

            logging.debug('Loaded PDB structure (OpenBabel)')

            return mol
        else:
            logging.debug('Loaded MMCIF structure (OpenBabel)')
            mol = protein_reader.read_mmcif_to_openbabel(path)

            return mol

    def _setup_filetype(self, filepath):
        """Sets up argument type to be consumedd by openbabel

        Args:
            filepath (str): Path to the source file

        Returns:
            str: file type
        """
        filetype = 'pdb'

        if filepath.endswith('.pdb') or filepath.endswith('.ent'):
            filetype = 'pdb'
        elif filepath.endswith('.mmcif') or filepath.endswith('.cif'):
            filetype = 'mmcif'

        return filetype

    def _get_angle(self, point_a, point_b, point_c):
        """Get the angle between three points in 3D space.
        Points should be supplied in Numpy array format.

        http://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space

        Args:
            point_a ([type]): [description]
            point_b ([type]): [description]
            point_c ([type]): [description]

        Returns:
            [type]: [description]
        """
        # In pseudo-code, the vector BA (call it v1) is:
        # v1 = {A.x - B.x, A.y - B.y, A.z - B.z}
        v1 = point_a - point_b

        # Similarly the vector BC (call it v2) is:
        # v2 = {C.x - B.x, C.y - B.y, C.z - B.z}
        v2 = point_c - point_b

        # The dot product of v1 and v2 is a function of the cosine of the angle between them
        # (it's scaled by the product of their magnitudes). So first normalize v1 and v2:

        # v1mag = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
        v1_mag = np.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])

        # v1norm = {v1.x / v1mag, v1.y / v1mag, v1.z / v1mag}
        v1_norm = np.array([v1[0] / v1_mag, v1[1] / v1_mag, v1[2] / v1_mag])

        # v2mag = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
        v2_mag = np.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])

        # v2norm = {v2.x / v2mag, v2.y / v2mag, v2.z / v2mag}
        v2_norm = np.array([v2[0] / v2_mag, v2[1] / v2_mag, v2[2] / v2_mag])

        # Then calculate the dot product:

        # res = v1norm.x * v2norm.x + v1norm.y * v2norm.y + v1norm.z * v2norm.z
        res = v1_norm[0] * v2_norm[0] + v1_norm[1] * v2_norm[1] + v1_norm[2] * v2_norm[2]

        # And finally, recover the angle:
        angle = np.arccos(res)

        if np.isnan(angle):
            logging.warn('Angle for <{}, {}, {}> was NaN, setting to Pi.'.format(point_a, point_b, point_c))
            angle = np.pi

        return angle

    def _group_angle(self, group, point_coords, degrees=False, signed=False):
        '''
        Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`

        `group` should be a dict with Numpy array 'center' and 'normal' attributes.
        `point_coords` should be a Numpy array.
        '''

        cosangle = np.dot(group['normal'], point_coords) / (np.linalg.norm(group['normal']) * np.linalg.norm(point_coords))

        # GET THE ANGLE AS RADIANS
        rad = np.arccos(cosangle)

        if not degrees:
            return rad

        # CONVERT RADIANS INTO DEGREES
        else:

            # CONVERT INTO A SIGNED ANGLE
            if signed:
                rad = rad - np.pi if rad > np.pi / 2 else rad

            # RETURN DEGREES
            return rad * 180 / np.pi

    def _group_group_angle(self, group, group2, degrees=False, signed=False):
        '''
        Adapted from CREDO: `https://bitbucket.org/blundell/credovi/src/bc337b9191518e10009002e3e6cb44819149980a/credovi/structbio/aromaticring.py?at=default`

        `group` and `group2` should be a dict with Numpy array 'center' and 'normal' attributes.
        '''

        cosangle = np.dot(group['normal'], group2['normal']) / (np.linalg.norm(group['normal']) * np.linalg.norm(group2['normal']))

        # GET THE ANGLE AS RADIANS
        rad = np.arccos(cosangle)

        if not degrees:
            return rad

        # CONVERT RADIANS INTO DEGREES
        else:

            # CONVERT INTO A SIGNED ANGLE
            if signed:
                rad = rad - np.pi if rad > np.pi / 2 else rad

            # RETURN DEGREES
            return rad * 180 / np.pi

    def _is_hbond(self, donor, acceptor, vdw_comp_factor):
        """Determines if atoms share Hbond.

        Args:
            donor (OBAtom): First atom
            acceptor (OBAtom): Second atom
            vdw_comp_factor (float): Compensation factor for VdW radii
                dependent interaction types.

        Returns:
            int: 1/0 binary information whether or not the atoms share
                Hbond.
        """

        for hydrogen_coord in donor.h_coords:

            h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)

            if h_dist <= config.VDW_RADII['H'] + acceptor.vdw_radius + vdw_comp_factor:

                if self._get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= config.CONTACT_TYPES['hbond']['angle rad']:

                    return 1

        return 0

    def _is_weak_hbond(self, donor, acceptor, vdw_comp_factor):
        """Determines if atoms share a weak Hbond.

        Args:
            donor (OBAtom): First atom
            acceptor (OBAtom): Second atom
            vdw_comp_factor (float): Compensation factor for VdW radii
                dependent interaction types.

        Returns:
            int: 1/0 binary information whether or not the atoms share
                weak Hbond.
        """

        for hydrogen_coord in donor.h_coords:

            h_dist = np.linalg.norm(hydrogen_coord - acceptor.coord)

            if h_dist <= config.VDW_RADII['H'] + acceptor.vdw_radius + vdw_comp_factor:

                if self._get_angle(donor.coord, hydrogen_coord, acceptor.coord) >= config.CONTACT_TYPES['weak hbond']['angle rad']:

                    return 1

        return 0

    def _is_halogen_weak_hbond(self, donor, halogen, ob_mol, vdw_comp_factor):
        """Determines if atoms share a weak halogen bond.

        Args:
            donor (OBAtom): Donor atom
            halogen (OBAtom): Halogen atom
            ob_mol (OBMol): [description]
            vdw_comp_factor (float): [description]

        Returns:
            int: 1/0 binary information whether or not the atoms share
                weak weak Hbond.
        """

        # `nbr` WILL BE A BIOPYTHON ATOM
        # ... HOPEFULLY
        nbr = self.ob_to_bio[utils.get_single_bond_neighbour(ob_mol.GetAtomById(self.bio_to_ob[halogen])).GetId()]

        for hydrogen_coord in donor.h_coords:

            h_dist = np.linalg.norm(halogen.coord - hydrogen_coord)

            if h_dist <= config.VDW_RADII['H'] + halogen.vdw_radius + vdw_comp_factor:

                if config.CONTACT_TYPES['weak hbond']['cx angle min rad'] <= self._get_angle(nbr.coord, halogen.coord, hydrogen_coord) <= config.CONTACT_TYPES['weak hbond']['cx angle max rad']:

                    return 1

        return 0

    def _is_xbond(self, donor, acceptor, ob_mol):
        """Determines if atoms share a weak halogen bond.

        Args:
            donor (OBAtom): Donor atom
            acceptor (OBAtom): Halogen atom
            ob_mol (OBMol): [description]

        Returns:
            int: 1/0 binary information whether or not the atoms share
                xbond.
        """

        # `nbr` WILL BE A BIOPYTHON ATOM
        # ... HOPEFULLY
        nbr = self.ob_to_bio[utils.get_single_bond_neighbour(ob_mol.GetAtomById(self.bio_to_ob[donor])).GetId()]
        theta = self._get_angle(nbr.coord, donor.coord, acceptor.coord)

        if (theta >= config.CONTACT_TYPES['xbond']['angle theta 1 rad']):
            return 1

        return 0

    def _update_atom_sift(self, atom, addition, contact_type='INTER'):
        '''
        '''

        atom.sift = [x or y for x, y in zip(atom.sift, addition)]

        if contact_type == 'INTER':
            atom.sift_inter_only = [x or y for x, y in zip(atom.sift_inter_only, addition)]

        if 'INTRA' in contact_type:
            atom.sift_intra_only = [x or y for x, y in zip(atom.sift_intra_only, addition)]

        if 'WATER' in contact_type:
            atom.sift_water_only = [x or y for x, y in zip(atom.sift_water_only, addition)]

    def _update_atom_fsift(self, atom, addition, contact_type='INTER'):
        """[summary]

        Args:
            atom ([type]): [description]
            addition ([type]): [description]
            contact_type (str, optional): Defaults to 'INTER'. [description]
        """'''
        '''

        atom.actual_fsift = [x or y for x, y in zip(atom.actual_fsift, addition)]

        if contact_type == 'INTER':
            atom.actual_fsift_inter_only = [x or y for x, y in zip(atom.actual_fsift_inter_only, addition)]

        if 'INTRA' in contact_type:
            atom.actual_fsift_intra_only = [x or y for x, y in zip(atom.actual_fsift_intra_only, addition)]

        if 'WATER' in contact_type:
            atom.actual_fsift_water_only = [x or y for x, y in zip(atom.actual_fsift_water_only, addition)]

    def _update_atom_integer_sift(self, atom, addition, contact_type='INTER'):
        """[summary]

        Args:
            atom ([type]): [description]
            addition ([type]): [description]
            contact_type (str, optional): Defaults to 'INTER'. [description]
        """'''
        '''

        atom.integer_sift = [x + y for x, y in zip(atom.sift, addition)]

        if contact_type == 'INTER':
            atom.integer_sift_inter_only = [x + y for x, y in zip(atom.sift_inter_only, addition)]

        if 'INTRA' in contact_type:
            atom.integer_sift_intra_only = [x + y for x, y in zip(atom.sift_intra_only, addition)]

        if 'WATER' in contact_type:
            atom.integer_sift_water_only = [x + y for x, y in zip(atom.sift_water_only, addition)]

    def _sift_xnor(self, sift1, sift2):
        """[summary]

        Args:
            sift1 ([type]): [description]
            sift2 ([type]): [description]

        Raises:
            ValueError: [description]

        Returns:
            [type]: [description]
        """'''
        '''

        out = []

        for x, y in zip(sift1, sift2):

            if x and not y:  # TF
                out.append(0)

            elif not x and not y:  # FF
                out.append(1)

            elif x and y:  # TT
                out.append(1)

            elif not x and y:  # FT
                out.append(0)

            else:
                raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

        return out

    def _sift_match_base3(self, sift1, sift2):
        """[summary]
        0 = UNMATCHED
        1 = MATCHED
        2 = MATCH NOT POSSIBLE

        Assuming that sift1 is the potential SIFt, and sift2 is the actual SIFt.
        Args:
            sift1 ([type]): [description]
            sift2 ([type]): [description]

        Raises:
            SiftMatchError: [description]
            ValueError: [description]

        Returns:
            [type]: [description]
        """

        out = []

        for x, y in zip(sift1, sift2):

            if x and not y:  # TF
                out.append(0)  # UNMATCHED

            elif not x and not y:  # FF
                out.append(2)  # MATCH NOT POSSIBLE

            elif x and y:  # TT
                out.append(1)  # MATCHED

            elif not x and y:  # FT
                raise SiftMatchError

            else:
                raise ValueError('Invalid SIFts for matching: {} and {}'.format(sift1, sift2))

        return out

    def _human_sift_match(self, sift_match, feature_sift=config.FEATURE_SIFT):
        '''
        Takes a base-3 SIFt indicating contact matched-ness and converts it to a human readable form.
        '''

        terms = []

        for e, fp in enumerate(sift_match):

            if fp == 2:  # MATCH NOT POSSIBLE
                continue

            elif fp == 1:
                terms.append('Matched {}'.format(feature_sift[e]))

            elif fp == 0:
                terms.append('Unmatched {}'.format(feature_sift[e]))

            else:
                raise ValueError

        terms.sort()
        return ':'.join(terms)

    # endregion
