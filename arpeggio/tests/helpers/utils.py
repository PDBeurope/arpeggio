import os


def process_arpeggio_pair(pdb, cif):

    files = filter(lambda l: l.split('.')[1] in ('contacts', 'bs_contacts', 'atom_types'), os.listdir(pdb))
    problems = []

    for i in files:
        pdb_file = read_in_file(os.path.join(pdb, i))
        cif_file = read_in_file(os.path.join(cif, i))

        diff = pdb_file.difference(cif_file)
        if len(diff) != 0:
            problems.append(i)

    return problems


def read_in_file(path):
    with open(path) as f:
        lines = f.readlines()
        fixed = map(lambda l: l.replace('.0', ''), lines)  # for some reason some files contain doubles instead of ints (e.g. 2.0 vs 2).
        return set(fixed)
