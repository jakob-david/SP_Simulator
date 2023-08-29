import numpy as np
import math


class Simulation:

    def __init__(self, grid, gravity, potential=lambda x=0: 0):
        """
        Class that carries out the actual simulation.

        :param grid: The grid on which the simulation is carried out.
        :param potential: The potential function used for the simulation. (default 0)
        :param gravity: The "gravity" of the system. (How much the waves attract one another.)
        """

        self.grid = grid
        self.potential = potential
        self.gravity = gravity

    # Some wrapper functions
    ##################################
    def set_init_function(self, func):
        """
        Sets the initial function of the system.

        :param func: The function to which the initial function should be set to. (wrapper function)
        """

        self.grid.set_init_function(func)

    def plot_3d(self, time_step=0, square=True):
        """
        Funktion prints a 3D wireframe of the system at one specified timestep. (wrapper function)

        :param time_step: The timestep of the plot.
        :param square: Set to true if the function in the graph should be squared.
        """

        self.grid.plot_3d(time_step, square)

    def heatmap(self, time_step=0, square=True, save=False):
        """
        Makes a heatmap of the system at a given timestep. (wrapper function)

        :param time_step: The timestep which should be plotted.
        :param square: Set to true if the system should be squared.
        :param save: Set true if the plot should be saved instead of shown.
        """

        self.grid.heatmap(time_step, square, save)

    def plot_energy_evolution(self, save=False, log=False):
        """
        Plots the energy evolution of the system. (wrapper function)

        :param save: Set true if the plot should be saved instead of shown.
        :param log: Set true if the log should be taken for the y-axis.
        """

        self.grid.plot_energy_evolution(save, log)

    def gif(self):
        """
        Makes a gif of the whole system. Consisting of a 3D plot of every timestep. (wrapper function)
        """

        self.grid.gif()

    def plot_3d_potential(self, save=False):
        """
        Plots the potential as a 3D graph. (wrapper function)

        :param save: Set true if the plot should be saved instead of shown.
        """

        self.grid.plot_3d_potential(save)

    def start_split_operator(self):
        """
        Starts the 2D simulation of the Schrödinger Poison equation using the split operator method.
        """

        dt = self.grid.grid_parameters.time_step_size

        # define operators (which are actually vectors)
        # -------------------------------------------------------------
        xv, yv = np.meshgrid(self.grid.x_axis, self.grid.y_axis)
        v = self.potential(xv, yv)
        opr_r = np.power(math.e, -1 * 1j * v * dt)

        # Info: From k = (2*pi)/L | The 2 vanishes since the grid is divided in half.
        space_steps_x = self.grid.grid_parameters.space_steps_X
        kx = np.fft.fftfreq(space_steps_x) * (space_steps_x / 2) * math.pi / self.grid.x_axis[-1]

        # Info: From k = (2*pi)/L | The 2 vanishes since the grid is divided in half.
        space_steps_y = self.grid.grid_parameters.space_steps_Y
        ky = np.fft.fftfreq(space_steps_y) * (space_steps_y / 2) * math.pi / self.grid.y_axis[-1]

        k = np.meshgrid(kx, ky)
        k = np.array(k)
        k = np.sqrt(np.power(k[0], 2) + np.power(k[1], 2))
        k[0][0] = 0.0001

        opr_k = np.power(math.e, -0.5 * 1j * np.power(k, 2) * dt)

        # meshgrid for 2D-Potential
        # ------------------------
        x, y = np.meshgrid(self.grid.x_axis.real, self.grid.y_axis.real)

        # set the first energy
        # ------------------------
        self.grid.energy[0] = np.sum(np.square(np.abs(self.grid.grid[0])))

        # iterate
        # ------------------------
        for i in range(0, self.grid.grid_parameters.time_steps - 1):

            # FFT
            # ------------------------
            tmp = np.fft.fft2(self.grid.grid[i])
            v = np.fft.fft2(np.abs(np.power(self.grid.grid[i], 2)))

            # Poisson
            # ------------------------
            v = (-1 / (np.power(np.abs(k), 2))) * v  # (np.power(np.abs(k)

            # Momentum
            # ------------------------
            tmp = tmp * opr_k

            # IFFT
            # ------------------------
            tmp = np.fft.ifft2(tmp)
            v = np.fft.ifft2(v)

            # Position
            # ------------------------
            v = v * self.gravity + self.potential(x, y)
            opr_r = np.power(math.e, -1 * 1j * v * dt)
            tmp = tmp * opr_r

            self.grid.grid[i + 1] = tmp

            self.grid.energy[i + 1] = np.sum(np.square(np.abs(self.grid.grid[i + 1])))
