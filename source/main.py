import sys
import numpy as np
import scipy as sp

from pprint import pprint

from timeit import default_timer as timer

from qiskit import Aer
from qiskit.aqua import QuantumInstance
from qiskit.aqua.operators import Z2Symmetries
from qiskit.aqua.algorithms.minimum_eigen_solvers import VQE
from qiskit.aqua.algorithms import ExactEigensolver
from qiskit.aqua.components.optimizers import SLSQP, L_BFGS_B
from qiskit.chemistry.core import Hamiltonian, TransformationType, QubitMappingType
from qiskit.chemistry.drivers import PySCFDriver, UnitsType
from qiskit.chemistry.components.variational_forms import UCCSD 
from qiskit.chemistry.components.initial_states import HartreeFock
from qiskit.chemistry.drivers import HFMethodType


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


params = {'molecule': None, 
          'spin': None,
          'charge': None,
          'coords': None,
          'param_trans': None}

# Parse cmd args
def parse_cmd_args(argv, params):
    for key in params: 
        if ('--' + key) == argv[0]:
            params[key] = argv[1]
            if len(argv) > 2: 
                parse_cmd_args(argv[2:], params)
                return
            # If there is no more entries in the list we are done
            else:
                return
        else:
            pass
    raise ValueError(f'ERROR: Unknown argument {argv[0]}')
            

def coords_arr(coords):
    if coords:
        start, num, stop = coords.split(':')
        return np.linspace(float(start), stop=float(stop), num=int(num))
    else:
        return None

parse_cmd_args(sys.argv[1:], params)


molecule = params['molecule']
charge = int(params['charge'])
spin = int(params['spin'])
coords = coords_arr(params['coords'])
param_trans = int(params['param_trans'])


print("Setting up system:")
print(f"  Molecule: {molecule}")
print(f"  Charge: {charge}")
print(f"  Spin: {spin}")

 
#######################################################

filename = f'results/{make_reasonable_mol_name(molecule)}_{bool(param_trans)}.csv'
output_file = open(filename, 'w')
print(f'Results will be written to {filename}')
initial_point = None
for coord in coords:
    driver = PySCFDriver(atom=molecule.format(coord), 
                         unit=UnitsType.ANGSTROM, 
                         charge=charge,
                         spin=spin,
                         hf_method=HFMethodType.UHF,
                         basis='sto3g')
    #                     basis='sto6g')
    #                     basis='631g')

    qmolecule = driver.run()

    core = Hamiltonian(transformation=TransformationType.FULL, 
                       qubit_mapping=QubitMappingType.JORDAN_WIGNER, 
                       two_qubit_reduction=False, 
                       freeze_core=True, 
                       orbital_reduction=[])

    qubit_op, _ = core.run(qmolecule)

    # optimizers
    optimizer = SLSQP(maxiter=50000)
    # optimizer = L_BFGS_B(maxiter=1000)

    # initial state
    init_state = HartreeFock(#num_qubits=qubit_op.num_qubits,
                             num_orbitals=core._molecule_info['num_orbitals'], 
                             qubit_mapping=core._qubit_mapping,
                             two_qubit_reduction=core._two_qubit_reduction, 
                             num_particles=core._molecule_info['num_particles'], 
    #                         sq_list=the_tapered_op.z2_symmetries.sq_list)
                             )

    # UCCSD Ansatz
    var_form = UCCSD(#num_qubits=qubit_op.num_qubits, 
                     #depth=1, 
                     num_orbitals=core._molecule_info['num_orbitals'], 
                     num_particles=core._molecule_info['num_particles'], 
                     active_occupied=None, 
                     active_unoccupied=None, 
                     initial_state=init_state, 
                     qubit_mapping=core._qubit_mapping, 
                     two_qubit_reduction=core._two_qubit_reduction, 
                     num_time_slices=1, 
    #                 z2_symmetries=the_tapered_op.z2_symmetries, 
                     shallow_circuit_concat=False)
    #                 force_no_tap_excitation=True, 
    #                 method_doubles='succ', 
    #                 excitation_type='d', 
    #                 same_spin_doubles=False)

    # set up VQE
    algo = VQE(qubit_op, var_form, optimizer, initial_point=initial_point)

    # Choose the backend (use Aer instead of BasicAer) 
    backend = Aer.get_backend('statevector_simulator') 
    quantum_instance = QuantumInstance(backend=backend)


    # run the algorithm
    #print("\nRunning VQE:")
    vqe_start = timer()
    algo_result = algo.run(quantum_instance)
    vqe_end = timer()

    # get the results
    result = core.process_algorithm_result(algo_result) 
    #pprint(result)
    cost_function_evals = result['algorithm_result']['cost_function_evals']
    if param_trans: 
        initial_point = algo.optimal_params

    #pprint(result)
    #statevector = result['algorithm_retvals']['eigenstate']
    #nuclear_repulsion_energy = result['nuclear_repulsion_energy']
    #hartree_fock_energy = result['hf_energy']

    #print('\n\n')
    #print('\n\n ### TIMES ### \n')
    print(f'{coord:.4}, {vqe_end - vqe_start}, {cost_function_evals}')
    output_file.write(f'{coord}, {vqe_end - vqe_start}, {cost_function_evals}\n')

