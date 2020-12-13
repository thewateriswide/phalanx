# phalanx

A Python3 package for simulating a quantum annealing computer/quantum annealer.

Phalanx uses a Hamiltonian like conventional Ising model to describe the behavior of the quantum annealer.  It calculates the eigenvalues and eigenstates of a Hamiltonian specified by the user, and returns pseudo outcome according to the Gibbs distribution.  It provides interfaces for users to specify the problem and  acquire the outcome.

It's encouraged to hack Phalanx by yourself,  to know the principle  behind is the most important.  And you may like to replace the physical  model in Phalanx to increase its power,  just try it !

## install

Place the 'phalanx' folder under your Python working directory,  then you can import it as a normal package.

## package requirement

numpy
		networkx

## *what is a quantum annealer?*

A quantum annealer is an open system with many qubits in it.  After an annealing process, this system will stay in the so-called Gibbs state.  

As a computer, its Hamiltonian can be viewed as a description of a  math problem, and the final state is the solution.  In reality, many  problems can be mapped to this type of calculation. Such as optimization problems and sampling problems.