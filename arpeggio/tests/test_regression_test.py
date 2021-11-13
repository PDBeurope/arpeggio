import json
import os
import shutil
import subprocess
from collections import namedtuple

import pytest

import arpeggio
from arpeggio.tests.helpers.utils import process_arpeggio_pair

"""Regression test to check if results provided by original implementation
of Arpeggio are exactly the same as with refactored version with the use
of mmcif files.
"""

AtomTypes = namedtuple('AtomTypes', ['name', 'types'])


class Config:
    """Test configuration. Python 2 environment with dependencies installed
    is hardcoded. So if required needs to be changed.
    """

    def __init__(self):
        root = os.path.join(os.path.dirname(arpeggio.__file__))

        self.structures_path = os.path.join(root, 'tests', 'test_data', 'structures')
        self.arpeggio = os.path.join(root, 'tests', 'test_data', 'old_arpeggio', 'arpeggio.py')
        self.py2 = '/Users/lpravda/anaconda3/envs/py2/bin/python'
        self.params = []
        with open(os.path.join(root, 'tests', 'test_data', 'config.json')) as f:
            temp = json.load(f)
            for k, v in temp.items():
                for selection in v:
                    self.params.append((k, selection))


class TestRegression:
    c = Config()

    @pytest.fixture(scope='class')
    def config(self):
        return Config()

    @staticmethod
    @pytest.mark.parametrize("pdb_id,selection", c.params)
    @pytest.fixture(scope='module')
    def run_arpeggio(tmpdir, config, pdb_id, selection):
        pdb = tmpdir.mkdir('py2')
        cif = tmpdir.mkdir('py3')

        pdb_file = pdb.join(pdb_id + '.pdb')
        cif_file = cif.join(pdb_id + '.cif')

        shutil.copyfile(os.path.join(config.structures_path, f'{pdb_id}.pdb'), pdb_file)
        shutil.copyfile(os.path.join(config.structures_path, f'{pdb_id}.cif'), cif_file)

        subprocess.run(['python', 'arpeggio',
                        str(cif_file), '-s', selection, '-o', os.path.dirname(str(cif_file))])
        subprocess.run([config.py2, config.arpeggio, str(pdb_file), '-s', selection])

        entries = process_arpeggio_pair(pdb, cif)

        return (str(pdb), str(cif))

    def compare_atom_types(run_arpeggio):
        pdb = run_arpeggio[0]
        cif = run_arpeggio[1]
