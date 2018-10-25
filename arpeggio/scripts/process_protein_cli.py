import argparse
import logging
import os

from arpeggio.core import InteractionComplex
from arpeggio.core.utils import max_mem_usage

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def create_parser():
    """Create command line parser

    Returns:
        argsparser: args parser
    """

    parser = argparse.ArgumentParser(description='''
        ############
        # ARPEGGIO #
        ############

        A program for calculating interactions,
        using only Open Source dependencies.

        Dependencies:
        - Python (v3.6)
        - Numpy
        - BioPython (>= v1.60)
        - OpenBabel (with Python bindings)

        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('filename', type=str, help='Path to the file to be analysed.')

    selection_group = parser.add_mutually_exclusive_group(required=False)
    selection_group.add_argument('-s', '--selection', type=str, nargs='+', help='Select the "ligand" for interactions, using selection syntax: /<chain_id>/<res_num>[<ins_code>]/<atom_name> or RESNAME:<het_id>. Fields can be omitted.')
    selection_group.add_argument('-sf', '--selection-file', type=str, help='Selections as above, but listed in a file.')

    parser.add_argument('-wh', '--write-hydrogenated', action='store_true', help='Write an output file including the added hydrogen coordinates.')
    parser.add_argument('-mh', '--minimise-hydrogens', action='store_true', help='Energy minimise OpenBabel added hydrogens.')
    parser.add_argument('-ms', '--minimisation-steps', type=int, default=50, help='Number of hydrogen minimisation steps to perform.')
    parser.add_argument('-mf', '--minimisation-forcefield', type=str, choices=('MMFF94', 'UFF', 'Ghemical'), default='MMFF94', help='Choose the forcefield to minimise hydrogens with. Ghemical is not recommended.')
    parser.add_argument('-mm', '--minimisation-method', type=str, choices=('DistanceGeometry', 'SteepestDescent', 'ConjugateGradients'), default='ConjugateGradients', help='Choose the method to minimise hydrogens with. ConjugateGradients is recommended.')
    parser.add_argument('-co', '--vdw-comp', type=float, default=0.1, help='Compensation factor for VdW radii dependent interaction types.')
    parser.add_argument('-i', '--interacting', type=float, default=5.0, help='Distance cutoff for grid points to be \'interacting\' with the entity.')
    parser.add_argument('-ph', type=float, default=7.4, help='pH for hydrogen addition.')
    parser.add_argument('-sa', '--include-sequence-adjacent', action='store_true', help='For intra-polypeptide interactions, include non-bonding interactions between residues that are next to each other in sequence; this is not done by default.')
    parser.add_argument('-a', '--use-ambiguities', action='store_true', help='Turn on abiguous definitions for ambiguous contacts.')
    parser.add_argument('-o', '--output', default=None, help='Define output directory.')
    parser.add_argument('-m', '--mute', action='store_true', help='Silent mode without debug info.')

    return parser


def _setup_logging(args):
    """Set up logging and working directory

    Args:
        args (ArgumentParser): parsed arguments from the command line
    """

    logging_level = logging.WARNING if args.mute else logging.DEBUG
    logging.basicConfig(level=logging_level,
                        format='%(levelname)s//%(asctime)s.%(msecs).03d//%(message)s',
                        datefmt='%H:%M:%S')

    logger.info('Program begin.')


def main():
    """Run Arpeggio algorithm
    """
    parser = create_parser()
    args = parser.parse_args()
    _setup_logging(args)
    run_arpeggio(args)


def run_arpeggio(args):

    args.output = os.getcwd() if args.output is None else args.output
    os.makedirs(args.output, exist_ok=True)

    i_complex = InteractionComplex(args.filename, args.vdw_comp, args.interacting, args.ph)
    i_complex.structure_checks()

    selections = _parse_selection(args)

    if args.use_ambiguities:
        i_complex.address_ambiguities()

    if args.minimise_hydrogens:
        i_complex.minimize_hydrogens(args.minimisation_forcefield,
                                     args.minimisation_method,
                                     args.minimisation_steps)

    if args.write_hydrogenated:
        i_complex.write_hydrogenated(args.output, args.filename)

    i_complex.run_arpeggio(selections, args.interacting, args.vdw_comp, args.include_sequence_adjacent)

    # write out files
    i_complex.write_atom_types(args.output)  # _atomtypes
    i_complex.write_contacts(selections, args.output)  # _contacts; _bs_contacts
    i_complex.write_atom_sifts(args.output)  # _sift; _specific_sift
    i_complex.write_binding_site_sifts(args.output)  # _siftmatch; _specific_siftmatch
    i_complex.write_polar_matching(args.output)  # _polarmatch; _specific_polarmatch
    i_complex.write_residue_sifts(args.output)  # residue_sifts

    logger.info('Program End. Maximum memory usage was {}.'.format(max_mem_usage()))


def _parse_selection(args):
    """Parse user defined selection

    Args:
        args (ArgumentParser): Application arguments

    Returns:
        list of str: Selections in the textual form.
    """
    selection = []

    if args.selection:
        selection = args.selection
    else:
        with open(args.selection_file, 'r') as f:
            selection = [line for line in f]

    logger.info('Selection perceived: {}'.format(selection))
    return selection
