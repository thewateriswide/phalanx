'''
Phalanx, a Python3 package of a simulator for quantum annealer.

A quantum annealer is a qubit system subject to some Hamiltonian.  This 
system will stay in the ground state(s) of the Hamiltonian after an annealing
process. 

As a computer,  its Hamiltonian can be viewed as a description of a math 
problem,  and the final state is the solution.

Phalanx uses a special Ising model to simulate the quantum annealer.  To be 
specific,  it calculates the ground state of a Hamiltonian from the user.

In reality, many problems can be mapped to this type of calculation.  Such as 
optimization problems.  For more, you may like to search for

    'Adiabatic Quantum Computation'

Also the D-Wave company has many good resources on their web site for their 
powerful quantum annealing devices.  It is very much worth a read.

'''

__version__ = '1.0.0'
__author__ = 'Yunheng Ma'

from .annealer import Annealer
