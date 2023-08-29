import numpy as np
import matplotlib.pyplot as plt
import warnings
import imageio.v2 as imageio
import os


class Grid:

    def __init__(self, grid_parameters):
        """
        The grid on which the simulation is performed.

        :param grid_parameters: The parameters for the grid. (dataclass)
        """

        self.grid_parameters = grid_parameters
        self.grid = np.zeros((self.grid_parameters.time_steps, self.grid_parameters.space_steps), dtype=np.complex_)
        self.energy = np.zeros(self.grid_parameters.time_steps)
        self.method = ""

        # Precalculate the x-axis since it is used quite often
        # ---------------------------------------------------------------
        self.x_axis = np.zeros(self.grid_parameters.space_steps, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps / 2) * self.grid_parameters.space_step_size
        if (self.grid_parameters.space_steps % 2) == 0:
            offset = offset + (self.grid_parameters.space_step_size / 2)
        for i in range(0, self.grid_parameters.space_steps):
            self.x_axis[i] = offset + i * self.grid_parameters.space_step_size
        # ---------------------------------------------------------------

    def set_init_function(self, func, normalize=True):
        """
        Sets the initial function.

        :param func: The initial function.
        :param normalize: Set to true if the function should be normalized.
        """

        if normalize:
            self.grid[0] = func(self.x_axis) / np.linalg.norm(func(self.x_axis))
        else:
            self.grid[0] = func(self.x_axis)

    def print_grid(self):
        """
        Prints the grid.
        """

        print(np.flipud(self.grid.round(2)))

    # Functions for plotting the system.
    # ---------------------------------------------------------------

    def plot_energy_evolution(self, save=False, log=False):
        """
        Plots the energy evolution of the system.

        :param save: Set true if the plot should be saved instead of shown.
        :param log: Set true if the log should be taken for the y-axis.
        """

        time_axis = np.zeros(self.grid_parameters.time_steps)

        for i in range(0, self.grid_parameters.time_steps):
            time_axis[i] = i * self.grid_parameters.time_step_size

        if log:
            plt.plot(time_axis, self.energy)
            plt.yscale("log")
        else:
            plt.ylim([-3, 3])  # kind of
            plt.plot(time_axis, self.energy)

        plt.xlabel("time")
        plt.ylabel("energy")
        plt.title("Energy Evolution")

        if save:

            if not os.path.exists("./pictures"):
                os.makedirs("./pictures")

            plt.savefig("./pictures/energy.png", dpi=300)
            plt.close()
        else:
            plt.show()

    def plot_2d(self, time_step, pot=lambda a, u: 0, use_for_gif=False):
        """
        Makes a 2D plot of the system at a given timestep.

        :param time_step: The timestep of the plot.
        :param pot: The potential function of the system.
        :param use_for_gif: Only for internal use. (Is set to true when generating a gif.)
        """

        # Deactivate Warnings while generating gif
        warnings.filterwarnings('ignore')

        y_pot = np.zeros(len(self.x_axis))

        for i in range(0, len(y_pot)):
            y_pot[i] = pot(self.x_axis[i], self.grid[time_step][i])

        x_plot = np.square(abs(self.grid[time_step]))
        x_plot = x_plot / np.linalg.norm(x_plot)

        x_plot_sec = (self.grid[time_step])

        fig, ax1 = plt.subplots(figsize=(10, 7), dpi=160)
        ax2 = ax1.twinx()

        ax1.plot(self.x_axis, x_plot, label='wave_squared', color="teal")
        ax1.plot(self.x_axis, x_plot_sec, label='wave', color="lightcoral")
        ax2.plot(self.x_axis, y_pot, label='potential', color="darkorange")

        handles, labels = [(a + b) for a, b in zip(ax1.get_legend_handles_labels(), ax2.get_legend_handles_labels())]
        plt.legend(handles, labels, loc='upper right')

        plt.title(self.method + "\n[dx=" + str(self.grid_parameters.time_step_size) + ", dt=" + str(
            self.grid_parameters.space_step_size) + "]   Time: " + str(
            round(time_step * self.grid_parameters.time_step_size, 3)) + "    Time-Step: " + str(time_step))
        ax1.set_ylabel("wave function")
        ax2.set_ylabel("potential")
        ax1.set_xlabel("space")

        if use_for_gif:
            plt.savefig("./gif/picture_for_giv_" + str(time_step) + ".png")
            plt.close()

        warnings.filterwarnings('default')

    def plot_3d(self, square=True):
        """
        Plots a 3D plot of the whole system.

        :param square: Set to true if the system should be squared.
        """

        # Deactivate Warnings while generating gif
        warnings.filterwarnings('ignore')

        y_axis = np.zeros(self.grid_parameters.time_steps)
        for i in range(0, self.grid_parameters.time_steps):
            y_axis[i] = i * self.grid_parameters.time_step_size

        x, y = np.meshgrid(self.x_axis.real, y_axis)

        if square:
            z = np.square(abs(self.grid))
        else:
            z = self.grid.real

        ax = plt.axes(projection='3d')
        ax.plot_wireframe(x, y, z, color='green')
        ax.set_xlabel('space')
        ax.set_ylabel('time')
        ax.set_zlabel('wave function')

        warnings.filterwarnings('default')

    def heatmap(self, square=True, save=False):
        """
        Makes a heatmap of the whole system.

        :param square: Set to true if the system should be squared.
        :param save: Set true if the plot should be saved instead of shown.
        """

        if square:
            z = np.square(abs(self.grid))
        else:
            z = self.grid.real

        plt.figure(dpi=150)
        plt.imshow(z, cmap='viridis', interpolation='nearest', origin='lower')
        plt.xlabel("grid points on the x-axis")
        plt.ylabel("time-step")
        plt.title("1D Split-Operator Method\n" + "dt=" + str(self.grid_parameters.time_step_size) + " dx=" + str(
            self.grid_parameters.space_step_size))  #
        cbar = plt.colorbar()
        cbar.set_label('squared wave function value')

        if save:

            if not os.path.exists("./pictures"):
                os.makedirs("./pictures")

            plt.savefig("./pictures/2D-heatmap.png", dpi=300)
            plt.close()
        else:
            plt.show()
            plt.close()

    def gif(self, pot=lambda a: 0):
        """
        Makes a gif of the whole system. Consisting of a 2D plot of every timestep.

        :param pot: The potential function of the system.
        """

        if not os.path.exists("./gif"):
            os.makedirs("./gif")

        for i in range(0, self.grid_parameters.time_steps):
            print("Generating picture: (" + str(i + 1) + "/" + str(self.grid_parameters.time_steps) + ")")
            self.plot_2d(i, pot, True)

        images = list()
        for i in range(0, self.grid_parameters.time_steps):
            images.append(imageio.imread("gif/picture_for_giv_" + str(i) + ".png"))
        imageio.mimsave("simulation.gif", images)
