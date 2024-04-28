"""Provides functions for Ising Hamiltonian/montecarlo sim"""
import numpy as np

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        return "".join([str(bit) for bit in self.config])

    def __eq__(self, other):        
        return self.int() == other.int() 
    
    def __len__(self):
        return np.size(self.config)

    def on(self):
        return self.config.count(1)

    def off(self):
        return self.config.count(0)
    
    def flip_site(self,i):
        self.config[i] = int(not(self.config[i]))
    
    def int(self):
        return sum([(bit) * 2**i for i, bit in enumerate(reversed(self.config))])
 
    def set_config(self, s:list[int]):
        self.config = s
        
    def set_int_config(self, dec:int):
        bit_string = []
        for _ in range(self.N):
            bit_string.insert(0, dec % 2)
            dec //= 2
        self.set_config(bit_string)

class IsingHamiltonian:
    def __init__(self, J=[[()]], mu=np.zeros(1)):
        self.J = J
        self.mu = mu
        self.nodes = []
        self.js = []
        self.N = len(self.J)   

        for i in range(self.N):
            self.nodes.append(np.zeros(len(self.J[i]),dtype=int))
            self.js.append(np.zeros(len(self.J[i])))
            for jidx, j in enumerate(self.J[i]):
                self.nodes[i][jidx] = j[0]
                self.js[i][jidx] = j[1]

        self.mu = np.array([i for i in self.mu])
    
    def energy(self, bs: BitString):
        """Compute energy of configuration, `bs`

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        bs   : Bitstring
            input configuration
        -------
        energy  : float
            Energy of the input configuration
        """
        if (bs.N != self.N):
            error(f"Wrong Dimension, bitstring dim {bs.N}, J dim {self.N}")

        e = 0
        for i in range(bs.N):
            for j in self.J[i]:
                if j[0] < i:
                    continue
                if (bs.config[i] == bs.config[j[0]]):
                    e += j[1] 
                else:
                    e -= j[1]
        spins = [(-1)**(n+1) for n in bs.config]
        e += np.dot(self.mu, spins)
        return e

    def compute_average_values(self, T: float):
        """
        Compute the average value of Energy, Magnetization, 
        Heat Capacity, and Magnetic Susceptibility 

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        T    : float
            temperature of the system
        Returns
        -------
        energy  : float
        magnetization  : float
        heat capacity  : float
        magnetic susceptibility  : float
        """
        k = 1.38064852 * 10e-23 # Boltzmann constant

        my_bs = BitString(self.N)

        def magnetization(bs: BitString):
            mag = 0
            for n in bs.config:
                if n == 0:
                    mag -= 1
                else:
                    mag += 1
            return mag

        # Add code here to find the lowest energy configuration
        E = 0.0
        M = 0.0
        Z = 0.0
        EE = 0.0
        MM = 0.0
        for config in range(2**self.N):
            my_bs.set_int_config(config)
            e = self.energy(my_bs)
            m = magnetization(my_bs)
            Zi = np.exp(-1*e/T)
            E += e*Zi
            EE += e * e * Zi
            M += m * Zi
            MM += m * m * Zi
            Z += Zi

        E /= Z
        M /= Z
        EE /= Z
        MM /= Z

        HC = (EE-E*E)/(T*T)
        MS = (MM-M*M)/T

        return E, M, HC, MS

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pass
