import sys
import numpy as np
import scipy as sp

from pprint import pprint
from timeit import default_timer as timer

from QiskitVQEWrapper.vqe_wrapper.VQE_wrapper import VQEWrapper
from pyscf_helpers import make_reasonable_mol_name


vqe_wrapper = VQEWrapper()
vqe_wrapper.spin = 0
vqe_wrapper.charge = 0
vqe_wrapper.freeze_core = True


parameter_transfer = False
coords = np.linspace(0.1, 2.0, num=100)
mol_str = "H 0.0 0.0 0.0; H {} 0.0 0.0"
filename = "PLEASE_DELETE_ME.txt"
#filename = f'results/{make_reasonable_mol_name(mol_str)}_{parameter_transfer}.csv'
output_file = open(filename, 'w')

print(f'Results will be written to {filename}')

for coord in coords:
    vqe_wrapper.molecule_string = mol_str.format(coord)

    vqe_wrapper.initiate()
    result = vqe_wrapper.run_vqe()

    cost_function_evals = result['algorithm_result']['cost_function_evals']
    if parameter_transfer: 
        vqe_wrapper.initial_point = vqe_wrapper.vqe_algo.optimal_params

    energy = result['nuclear_repulsion_energy'] + result['computed_electronic_energy']

    print(f'{coord:.4}, {vqe_wrapper.vqe_time}, {cost_function_evals}, {energy}')
    output_file.write(f'{coord}, {vqe_wrapper.vqe_time}, {cost_function_evals}, {energy}\n')

output_file.close()
