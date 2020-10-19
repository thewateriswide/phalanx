'''
Class for quantum annealer.
'''

__all__ = ['Annealer']

import numpy as np
import networkx as nx

# Two quantum gates in Ising model.
i_gate = np.array([[1,  0], [0,  1]])
z_gate = np.array([[1,  0], [0, -1]])

class Annealer(object):
    '''Class for simulating a quantum annealer.  
    
    A quantum annealer is a system stays at the ground state of its Hamiltonian. 
    The topology of the annealer is a simple undirected graph.  Its structure 
    is determined by its Hamiltonian.  The vertexes are qubits.  The 
    Hamiltonian governs it is of normal Ising model type,

        A[i,j]*Z_i*Z_j + B[i]*Z_i,

    where A and B are coupling and external strength.  Note that any two qubits
    could have a coupling.  The coupling and external are put together as 
    connection.

    Example:

        To create an instance:
            anl = Annealer(number_of_qubits)
        To set a connection:
            connect = [
                [1,2,3], 
                [0,2,4],
                [0,0,3]]
            anl.set_connection(connect)
        To check the connection graph:
            anl.draw_topology()
        To measure the final state of the annealing process:
            anl.measure()
        To have a graphic view for result of measurement:
            anl.draw_result()
    '''

    def __init__(self, number_of_qubits):
        '''Create an annealer object. 

        -IN:
            number_of_qubits --- Number of qubits in the annealer.
                type: integer

        -INFLUENCED:
            self.number_of_qubits --- Number of qubits in the annealer.
                type: integer
            self.connection --- Value of externals and couplings.
                type: Upper triangular matrix of float.
            self.graph --- Structure graph of the annealer.
                type: networkx.Graph
            self.hamiltonian --- Hamiltonian matrix.
                type: numpy.array of float.
            self.groundstate --- Ground state(s) of the Hamiltonian
                type: numpy.array of complex.
            self.result --- Qubits configuration after a measurement.
                type: List of integers. 0 for up, 1 for down.
        '''
        self.number_of_qubits = number_of_qubits
        self.graph = nx.Graph()
        self.connection = None
        self.hamiltonian = None
        self.groundstate = None
        self.result = None
        return None
        
    def set_connection(self, connection):
        '''Set the values of coupling and external strength.

        The connection is an upper triangular matrix.  Its diagonal elements
        are external strength values for qubits.  For example, 

            connection[2, 2]

        is the external strength on qubit 2.  Other upper right elements 
        are couplings.  For example, 

            connection[1, 3]

        is the coupling strength between qubit 1 and 3.

        -IN:
            connection --- Value of externals and couplings.  In a matrix form.
                type: matrix of float.

        -INFLUENCED:
            self.connection --- Value of externals and couplings.
                type: Upper triangular matrix of float.
            self.graph --- Structure graph of the annealer.
                type: networkx.Graph
            self.hamiltonian --- Hamiltonian matrix.
                type: numpy.array of float.
            self.groundstate --- Ground state(s) of the Hamiltonian
                type: numpy.array of complex.
        '''
        # Modify the size of the input connection matrix to fit the annealer.
        tmp = np.array(connection)
        self.connection = np.zeros((self.number_of_qubits, self.number_of_qubits))
        for i in range(min(tmp.shape[0], self.number_of_qubits)):
            for j in range(min(tmp.shape[1], self.number_of_qubits)):
                self.connection[i, j] = tmp[i, j]
        self.connection = np.triu(self.connection)

        # Modify the graph.
        self.graph = nx.Graph()
        for i in range(self.number_of_qubits):
            self.graph.add_node(i, external=self.connection[i, i])
            for j in range(i+1, self.number_of_qubits):
                if self.connection[i, j] != 0.0:
                    self.graph.add_edge(i, j, coupling=self.connection[i, j])

        # Modify the hamiltonian.
        self.hamiltonian = np.zeros((2**self.number_of_qubits, 2**self.number_of_qubits))
        
        # External field part of hamiltonian.
        for i in range(self.number_of_qubits):
            if self.connection[i, i] != 0.0:
                term_ext = 1.0
                for j in range(self.number_of_qubits):
                    if j == i:
                        term_ext = np.kron(term_ext, z_gate)
                    else:
                        term_ext = np.kron(term_ext, i_gate)
                term_ext = self.connection[i, i] * term_ext
                self.hamiltonian += term_ext

        # Coupling part of hamiltonian.
        for i in range(self.number_of_qubits):
            for j in range(i+1, self.number_of_qubits):
                if self.connection[i, j] != 0.0:
                    term_cpl = 1.0
                    for k in range(self.number_of_qubits):
                        if k == i or k == j:
                            term_cpl = np.kron(term_cpl, z_gate)
                        else:
                            term_cpl = np.kron(term_cpl, i_gate)
                    term_cpl = self.connection[i, j] * term_cpl
                    self.hamiltonian += term_cpl

        # Search for ground state(s) of the Hamiltonian.
        self.groundstate = []
        eigenvalue, eigenstate = np.linalg.eig(self.hamiltonian)
        min_val = np.min(eigenvalue)
        for i in range(len(eigenvalue)):
            if eigenvalue[i] == min_val:
                self.groundstate.append(eigenstate[:, i])
        self.groundstate = np.array(self.groundstate)
        return None

    def draw_topology(self):
        '''Draw the graph of topology.
        '''
        from matplotlib.pyplot import show
        nx.draw(self.graph, with_labels=True, node_size=800, font_weight='bold')
        show()
        return None
           
    def measure(self):
        '''Give out a measurement of the qubits on the annealer. 

        This returns a basis accorrding to the ground state of the annealer.

        -INFLUENCED:
            self.result --- Qubits configuration after a measurement.
                type: List of integers. 0 for up, 1 for down.

        -RETURN:
            --- A configuration list of qubits.
                type: List of integers. 0 for up, 1 for down.
        '''
        # Since there could have many ground states.  A final state of the annealer
        # could be a superposition of them, which we have no further information to
        # distinguish.  So we just randomly pick one as the final state.
        groundstate = self.groundstate[np.random.choice(range(self.groundstate.shape[0]))]

        # Probabilites for groundstate collapsing to an observation basis.
        prob = np.zeros(2**self.number_of_qubits)
        for i in range(2**self.number_of_qubits):
            prob[i] = np.abs(groundstate[i])**2

        # Pick a basis as the measurement, and return it.
        idx = np.random.choice(np.arange(2**self.number_of_qubits), p = prob)
        form = '0%db' % self.number_of_qubits
        result = list(format(idx, form))
        for i in range(self.number_of_qubits):
            result[i] = int(result[i])
        self.result = result
        return result

    def draw_result(self):
        '''Draw the result of the previous measurement.
        '''
        from matplotlib.pyplot import show
        color_map = []
        for i in self.graph.nodes:
            if self.result[i] == 1:
                color_map.append('#B0C4DE')  # Color for down.
            else:
                color_map.append('#4682B4')  # Color for up.
        nx.draw(self.graph, node_color=color_map, with_labels=True, node_size=800, font_weight='bold')
        show()
        return None