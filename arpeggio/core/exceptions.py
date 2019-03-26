"""Exception classes to distinguis internal errors
"""
import logging


class HydrogenError(Exception):

    def __init__(self):
        logging.error('Please remove all hydrogens from the structure then re-run.')


class OBBioMatchError(Exception):

    def __init__(self, serial=''):

        if not serial:
            logging.error('An OpenBabel atom could not be matched to a BioPython counterpart.')

        else:
            logging.error(f'OpenBabel OBAtom with serial number {serial} could not be matched to a BioPython counterpart.')


class AtomSerialError(Exception):

    def __init__(self):
        logging.error('One or more atom serial numbers are duplicated.')


class SiftMatchError(Exception):

    def __init__(self):
        logging.error('Seeing is not believing.')


class SelectionError(Exception):

    def __init__(self, selection):
        logging.error(f'Invalid selector: {selection}')
