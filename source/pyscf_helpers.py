import pyscf


def make_reasonable_mol_name(mol_str):
    mols = mol_str.split(';')
    symbols = [m.split()[0] for m in mols]
    # Make a dict, will remove duplicates
    mol_dict = dict.fromkeys(symbols)

    # count them
    mol_str = ''
    for key in mol_dict:
        N = symbols.count(key)
        mol_str += key + (str(N) if N > 1 else '')

    return mol_str
