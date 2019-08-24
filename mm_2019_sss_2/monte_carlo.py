import numpy as np
from energy import EnergyFunctions


class MonteCarlo(EnergyFunctions):

    def __init__(self, setup: object = None,
                 tune_displacement: bool = True,
                 max_displacement: (int, float) = 0.1,
                 output_freq: (int, float) = 1000):

        self.coordinates = setup.coordinates
        self.n_particles = setup.n_particles
        self.box_length = setup.box_length
        self.cutoff = setup.cutoff
        self.volume = setup.volume
        self.beta = setup.beta
        self._max_displacement = max_displacement
        self._tune_displacement = tune_displacement
        self._output_freq = output_freq

    @property
    def max_displacement(self):
        return self._max_displacement

    @property
    def tune_displacement(self):
        return self._tune_displacement

    @property
    def output_freq(self):
        return self._output_freq

    @staticmethod
    def accept_or_reject(delta_e: float, beta: float):
        """Accept or reject a move based on the energy difference and system
        temperature.
        This function uses a random numbers to adjust the acceptance criteria.
        Parameters
        ----------
        delta_e : float
            The difference between the proposed and current energies.
        beta : float
            The inverse value of the reduced temperature.
        Returns
        -------
        accept : boolean
            Either a "True" or "False" to determine whether to reject the trial.
        """
        if delta_e < 0.0:
            accept = True
        else:
            random_number = np.random.rand(1)
            p_acc = np.exp(-beta * delta_e)
            if random_number < p_acc:
                accept = True
            else:
                accept = False

        return accept

    def adjust_displacement(self, n_trials: int, n_accept: int):
        """Change the acceptance criteria to get the desired rate.
        When the acceptance rate is too high, the maximum displacement is
        adjusted to be higher. When the acceptance rate is too low, the
        maximum displacement is adjusted lower.
        Parameters
        ----------
        n_trials : integer
            The number of trials that have been performed when the function is
            initiated.
        n_accept : integer
            The current number of accepted trials when the function is
            initiated.
        max_displacement : float
            The specified maximum value for the displacement of the trial.
        Returns
        -------
        max_displacement : float
            The adjusted displacement based on the acceptance rate.
        n_trials : integer, 0
            The new number of trials.
        n_accept : integer, 0
            The new number of trials.
        """

        acc_rate = float(n_accept) / float(n_trials)

        if acc_rate < 0.38:
            self.max_displacement *= 0.8

        elif acc_rate > 0.42:
            self.max_displacement *= 1.2

        n_trials = 0
        n_accept = 0

        return n_trials, n_accept

    def run_simulation(self, num_steps: int = 6000):

        total_pair_energy = self.initial_energy
        tail_correction = self.tail_correction
        
        print(f'Total initial pair energy: {self.initial_energy}')

        energy_array = np.zeros(num_steps)

        # start the Monte Carlo iterations
        n_trials = 0
        n_accept = 0

        for i_step in range(num_steps):
            n_trials += 1

            i_particle = np.random.randint(self.n_particles)

            current_energy = self.calculate_pair_energy(i_particle,
                                                        self.coordinates)

            random_displacement = (2.0 * np.random.rand(3) - 1.0) \
                                  * self.max_displacement
            proposed_coordinates = self.coordinates.copy()
            proposed_coordinates[i_particle] += random_displacement
            proposed_coordinates -= self.box_length \
                                    * np.round(proposed_coordinates
                                               / self.box_length)
            proposed_energy = self.calculate_pair_energy(i_particle,
                                                         proposed_coordinates)

            delta_e = proposed_energy - current_energy

            if self.accept_or_reject(delta_e, self.beta):
                n_accept += 1

                total_pair_energy += delta_e

                self.coordinates[i_particle] += random_displacement

                self.coordinates -= self.box_length \
                                    * np.round(self.coordinates
                                               / self.box_length)

            total_energy = (total_pair_energy + tail_correction) \
                            / self.n_particles

            energy_array[i_step] = total_energy

            if np.mod(i_step + 1, self.output_freq) == 0:
                print(f'Step: {i_step} Energy: {total_energy}')

            # if self.tune_displacement:
            #     max_displacement, n_trials, n_accept \
            #         = self.adjust_displacement(n_trials, n_accept)

        return
