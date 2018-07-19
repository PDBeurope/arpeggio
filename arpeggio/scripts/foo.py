import argparse

parser = argparse.ArgumentParser(description='add', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('filename', type=str, help='Path to the file to be analysed.')
parser.add_argument('-s', '--selection', type=str, nargs='+', help='Select the "ligand" for interactions, using selection syntax: /<chain_id>/<res_num>[<ins_code>]/<atom_name> or RESNAME:<het_id>. Fields can be omitted.')
args = parser.parse_args()
print(args.selection)
