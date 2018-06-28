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
            logging.error('OpenBabel OBAtom with serial number {} could not be matched to a BioPython counterpart.'.format(serial))


class AtomSerialError(Exception):

    def __init__(self):
        logging.error('One or more atom serial numbers are duplicated.')


class SiftMatchError(Exception):

    def __init__(self):
        logging.error('Seeing is not believing.')


class SelectionError(Exception):

    def __init__(self, selection):
        logging.error('Invalid selector: {}'.format(selection))
