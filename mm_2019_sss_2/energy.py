import numpy as np


class EnergyFunctions:

    @staticmethod
    def LJ(r: float):
        return 4 * ((1.0 / r) ** 12
                    - (1.0 / r) ** 6)

    @property
    def tail_correction(self):
        r3 = np.power(1.0 / self.cutoff, 3)
        r9 = np.power(r3, 3)
        v = self.volume
        n2 = np.power(self.n_particles, 2)

        return ((8.0 * np.pi * n2) / (3.0 * v)) * (((1.0 / 3.0) * r9) - r3)

    # TODO: Commented out until I can check this.
    # @staticmethod
    # def Buckingham(a: float, c: float, rho: float, r: float):
    #     return a * np.exp(-r / rho) - c / r ** 6

    @staticmethod
    def minimum_image_distance(r_i: np.array, r_j: np.array,
                               box_length: float):
        """Calculates the shortest distance between a particle and another
        instance in a periodic boundary condition image.
        Parameters
        ----------
        r_i: np.array([n,3])
            The x, y, z coordinates for a particle, i.
        r_j: np.array([n,3])
            The x, y, z coordinates for a particle, j.
        box_length: float, int
            The length of a side of the side box for the periodic boundary.
        Returns
        -------
        rij2: np.array([n,3])
            The minimum image distance between the two particles, r_i and r_j.
        """
        rij = r_i - r_j
        rij -= box_length * np.round(rij / box_length)
        rij2 = np.dot(rij, rij)
        distance = np.sqrt(rij2)

        return distance

    @property
    def initial_energy(self):
        """Iterates over a set of coordinates to calculate total system energy
        This function computes the sum of all pairwise VDW energy between each
        pair of particles in the system. This is the first instance of the
        energy calculation. Subsequent uses call calculate_pair_energy.
        Parameters
        ----------
        coordinates : np.array([n,3])
            An array of atomic coordinates. Size should be [n, 3] where n is the
            number of particles.
        box_length : float
            A float indicating the size of the simulation box. Can be either
            hard-coded or calculated using num_particles and reduced_density.
        cutoff: float
            The square of the simulation_cutoff, which is the cutoff distance
            between two interacting particles.
        i_particle: int
            Intitial particle for pairwise count
        Returns
        -------
        e_total : float
            The sum of all pairwise VDW energy between each pair of particles in
            the system.
        """
        e_total = 0.0

        for i_particle in range(self.n_particles):
            for j_particle in range(i_particle):
                r_i = self.coordinates[i_particle]
                r_j = self.coordinates[j_particle]
                distance = self.minimum_image_distance(r_i, r_j,
                                                       self.box_length)
                if distance < self.cutoff:
                    e_pair = self.LJ(distance)
                    e_total += e_pair

        return e_total

    def calculate_pair_energy(self, i_particle, coordinates):
        """This function computes the sum of all pairwise VDW energy between each
            pair of particles in the system.
        Parameters
        ----------
        coordinates : np.array
            An array of atomic coordinates. Size should be (n, 3) where n is the
            number of particles.
        box_length : float
            A float indicating the size of the simulation box. Can be either
            hard-coded or calculated using num_particles and reduced_density.
        cutoff: float
            The square of the simulation_cutoff, which is the cutoff distance
            between two interacting particles.
        i_particle: integer
            Intitial particle for pairwise count
        Returns
        -------
        e_total : float
            The sum of all pairwise VDW energy between each pair of particles in
            the system.
        """

        e_total = 0.0
        i_position = self.coordinates[i_particle]

        for j_particle in range(self.n_particles):
            if i_particle != j_particle:
                j_position = coordinates[j_particle]
                distance = self.minimum_image_distance(i_position, j_position,
                                                       self.box_length)
                if distance < self.cutoff:
                    e_pair = self.LJ(distance)
                    e_total += e_pair

        return e_total
