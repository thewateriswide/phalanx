# phalanx
A Python3 package for simulating a quantum annealing computer/quantum annealer.

Phalanx uses a special Ising model to simulate the quantum annealer. To be specific, it calculates the ground state of a Hamiltonian given by the user.  It provides interfaces for users to specify the problem and acquire the solution.

It's encouraged to hack Phalanx by yourself,  to know the principle behind is the most important.  And you may like to replace the physical model in Phalanx to increase its power,  just try it !

## install
Place the 'phalanx' folder under your Python working directory,  then you can import it as a normal package.

## package requirement
numpy  
networkx  
matplotlib

## _what is a quantum annealer?_
A quantum annealer is a qubit system subjects to some Hamiltonian. This system will stay in the ground state(s) of the Hamiltonian after an annealing process. 

As a computer, its Hamiltonian can be viewed as a description of a math problem, and the final state is the solution.  In reality, many problems can be mapped to this type of calculation. Such as optimization problems. 
