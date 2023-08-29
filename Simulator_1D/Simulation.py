import numpy as np
import math


class Simulation:

    def __init__(self, grid, potential=lambda x=0: 0):
        """
        Class that carries out the actual simulation.

        :param grid: The grid on which the simulation is carried out.
        :param potential: The potential function used for the simulation. (default 0)
        """

        self.grid = grid
        self.potential = potential

    def set_init_function(self, func):
        """
        Sets the initial function of the system.

        :param func: The function to which the initial function should be set to. (wrapper function)
        """
        self.grid.set_init_function(func)

    def plot_energy_evolution(self, save=False, log=False):
        """
        Plots the energy evolution of the system. (wrapper function)

        :param save: Set true if the plot should be saved instead of shown.
        :param log: Set true if the log should be taken for the y-axis.
        """

        self.grid.plot_energy_evolution(save, log)

    def plot_2d(self, time_step):
        """
        Makes a 2D plot of the system at a given timestep. (wrapper function)

        :param time_step: The timestep of the plot.
        """
        self.grid.plot_2d(time_step, self.potential)

    def plot_3d(self, square=True):
        """
        Plots a 3D plot of the whole system. (wrapper function)

        :param square: Set to true if the system should be squared.
        """

        self.grid.plot_3d(square)

    def gif(self):
        """
        Makes a gif of the whole system. Consisting of a 2D plot of every timestep. (wrapper function)
        """

        self.grid.gif(self.potential)

    def heatmap(self, square=True, save=False):
        """
        Makes a heatmap of the whole system. (wrapper function)

        :param square: Set to true if the system should be squared.
        :param save: Set true if the plot should be saved instead of shown.
        """

        self.grid.heatmap(square, save)

    def start_split_operator(self):
        """
        Starts the 1D simulation of the Schrödinger Poison equation using the split operator method.
        """

        self.grid.method = "Split-Time"

        dt = self.grid.grid_parameters.time_step_size

        # define operators (which are actually vectors)
        # ---------------------------------------------
        v = self.potential(self.grid.x_axis)
        opr_r = np.power(math.e, -1 * 1j * v * dt)

        # Info: From k = (2*pi)/L | The 2 vanishes since the grid is divided in half.
        space_step = self.grid.grid_parameters.space_steps
        k = (np.fft.fftfreq(space_step) * (space_step / 2) * math.pi / self.grid.x_axis[-1])
        opr_k = np.power(math.e, -0.5 * 1j * np.power(k, 2) * dt)

        k[0] = 0.0001

        # set the first energy
        # ------------------------
        self.grid.energy[0] = np.sum(np.square(np.abs(self.grid.grid[0])))

        # iterate
        # ------------------------
        for i in range(0, self.grid.grid_parameters.time_steps - 1):

            # FFT
            # ------------------------
            tmp = np.fft.fft(self.grid.grid[i])
            v = np.fft.fft(np.abs(np.power(self.grid.grid[i], 2)))

            # Poisson
            # ------------------------
            v = (-1 / (np.power(np.abs(k), 2))) * v  # (np.power(np.abs(k)

            # Momentum
            # ------------------------
            tmp = tmp * opr_k

            # IFFT
            # ------------------------
            tmp = np.fft.ifft(tmp)
            v = np.fft.ifft(v)

            # Position
            # ------------------------
            v = v * 5 + self.potential(self.grid.x_axis)
            opr_r = np.power(math.e, -1 * 1j * v * dt)
            tmp = tmp * opr_r

            self.grid.grid[i + 1] = tmp

            # save the current energy of the system
            # ------------------------
            self.grid.energy[i + 1] = np.sum(np.square(np.abs(self.grid.grid[i + 1])))
