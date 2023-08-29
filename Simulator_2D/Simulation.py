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
    def initFunction(self, func):
        self.grid.set_init_function(func)

    def plot3D(self, time_step=0, square=True, use_for_gif=False):
        self.grid.plot_3d(time_step, square, use_for_gif)

    def heatmap(self, time_step=0, square=True, save=False):
        self.grid.heatmap(time_step, square, save)

    def plotEnergyEvolution(self, save=False, log=False):
        self.grid.plot_energy_evolution(save, log)

    def gif(self):
        self.grid.gif()

    def plot3DPotential(self, download=False):
        self.grid.plot_3d_potential(download)

    # Split-Time Poisson
    ##################################
    def SO_Poisson_2D(self, adjust_dt=False):
        self.grid.method = "Split-Time Poisson"

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
