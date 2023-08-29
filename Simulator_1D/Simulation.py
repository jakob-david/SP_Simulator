import numpy as np
import math

class Simulation:

    # Constructor
    ##################################
    def __init__(self, Grid, Potential=lambda x=0: 0):
        self.grid = Grid
        self.potential = Potential

    # Some wrapper functions
    ##################################
    def initFunction(self, func):
        self.grid.set_init_function(func)

    def plotEnergyEvolution(self, save=False, log=False):
        self.grid.plot_energy_evolution(save, log)

    def plot2D(self, time_step):
        self.grid.plot_2d(time_step, self.potential)

    def plot3D(self, square=True):
        self.grid.plot_3d(square)

    def gif(self):
        self.grid.gif(self.potential)

    def heatmap(self, square=True, save=False):
        self.grid.heatmap(square, save)

    # Split-Operator method for Schrödinger Poisson
    ########################################################
    def SO_Poisson(self, adjust_dt=False):
        self.grid.method = "Split-Time"

        dt = self.grid.grid_parameters.time_step_size

        # define operators (which are actually vectors)
        ###############################################################
        V = self.potential(self.grid.x_axis)
        opr_R = np.power(math.e, -1 * 1j * V * dt)

        k = np.fft.fftfreq(self.grid.grid_parameters.space_steps) * (
                    self.grid.grid_parameters.space_steps / 2) * math.pi / self.grid.x_axis[
                -1]  # From k = (2*pi)/L | The 2 vanishes since the grid is devided in half.
        opr_K = np.power(math.e, -0.5 * 1j * np.power(k, 2) * dt)

        k[0] = 0.0001

        # set the first enery
        #################################
        self.grid.energy[0] = np.sum(np.square(np.abs(self.grid.grid[0])))

        # iterate
        #####################
        for i in range(0, self.grid.grid_parameters.time_steps - 1):
            # FFT
            ############################
            tmp = np.fft.fft(self.grid.grid[i])
            V = np.fft.fft(np.abs(np.power(self.grid.grid[i], 2)))

            # Poisson
            ############################
            V = (-1 / (np.power(np.abs(k), 2))) * V  # (np.power(np.abs(k)

            # Momentum
            ############################
            tmp = tmp * opr_K

            # IFFT
            ############################
            tmp = np.fft.ifft(tmp)
            V = np.fft.ifft(V)

            # Position
            ############################
            V = V * 5 + self.potential(self.grid.x_axis)
            opr_R = np.power(math.e, -1 * 1j * V * dt)
            tmp = tmp * opr_R

            self.grid.grid[i + 1] = tmp

            # save the current energy of the system
            ##########################################
            self.grid.energy[i + 1] = np.sum(np.square(np.abs(self.grid.grid[i + 1])))
