"""Provides functions for Ising Hamiltonian/montecarlo sim"""

class IsingHamiltonian:
    def __init__(self, J=[[()]], mu=np.zeroes(1)):
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
                if (bs.config[i] == bs.config[j]):
                    e += J[1] 
                else:
                    e -= J[0]
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
        for config in range(2**N):
            my_bs.set_int_config(config)
            e = energy(my_bs, G)
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
