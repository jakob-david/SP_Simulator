import numpy as np
import math

class Simulation:

    # Constructor
    ##################################
    def __init__(self, Grid, gravity, Potential=lambda x=0: 0):
        self.grid = Grid
        self.potential = Potential
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
        ###############################################################
        xv, yv = np.meshgrid(self.grid.x_axis, self.grid.y_axis)
        V = self.potential(xv, yv)
        opr_R = np.power(math.e, -1 * 1j * V * dt)

        kx = np.fft.fftfreq(self.grid.grid_parameters.space_steps_X) * (
                    self.grid.grid_parameters.space_steps_X / 2) * math.pi / self.grid.x_axis[
                 -1]  # From k = (2*pi)/L | The 2 vanishes since the grid is devided in half.
        ky = np.fft.fftfreq(self.grid.grid_parameters.space_steps_Y) * (
                    self.grid.grid_parameters.space_steps_Y / 2) * math.pi / self.grid.y_axis[
                 -1]  # From k = (2*pi)/L | The 2 vanishes since the grid is devided in half.
        k = np.meshgrid(kx, ky)
        k = np.array(k)
        k = np.sqrt(np.power(k[0], 2) + np.power(k[1], 2))
        k[0][0] = 0.0001

        opr_K = np.power(math.e, -0.5 * 1j * np.power(k, 2) * dt)

        # meshgrid for 2D-Potential
        #################################
        X, Y = np.meshgrid(self.grid.x_axis.real, self.grid.y_axis.real)

        # set the first enery
        #################################
        self.grid.energy[0] = np.sum(np.square(np.abs(self.grid.grid[0])))

        # iterate
        #####################
        for i in range(0, self.grid.grid_parameters.time_steps - 1):
            # FFT
            ############################
            tmp = np.fft.fft2(self.grid.grid[i])
            V = np.fft.fft2(np.abs(np.power(self.grid.grid[i], 2)))

            # Poisson
            ############################
            V = (-1 / (np.power(np.abs(k), 2))) * V  # (np.power(np.abs(k)

            # Momentum
            ############################
            tmp = tmp * opr_K

            # IFFT
            ############################
            tmp = np.fft.ifft2(tmp)
            V = np.fft.ifft2(V)

            # Position
            ############################
            V = V * self.gravity + self.potential(X, Y)
            opr_R = np.power(math.e, -1 * 1j * V * dt)
            tmp = tmp * opr_R

            self.grid.grid[i + 1] = tmp

            self.grid.energy[i + 1] = np.sum(np.square(np.abs(self.grid.grid[i + 1])))
