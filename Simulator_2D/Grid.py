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
        self.grid = np.zeros(
            (self.grid_parameters.time_steps, self.grid_parameters.space_steps_X, self.grid_parameters.space_steps_Y),
            dtype=np.complex_)
        self.energy = np.zeros(self.grid_parameters.time_steps)
        self.method = ""

        # Precalculate the x-axis since it is used quite often
        # ---------------------------------------------------------------
        self.x_axis = np.zeros(self.grid_parameters.space_steps_X, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps_X / 2) * self.grid_parameters.space_step_size_X
        if (self.grid_parameters.space_steps_X % 2) == 0:
            offset = offset + (self.grid_parameters.space_step_size_X / 2)
        for i in range(0, self.grid_parameters.space_steps_X):
            self.x_axis[i] = offset + i * self.grid_parameters.space_step_size_X
        # ---------------------------------------------------------------

        # Precalculate the y-axis since it is used quite often
        # ---------------------------------------------------------------
        self.y_axis = np.zeros(self.grid_parameters.space_steps_Y, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps_Y / 2) * self.grid_parameters.space_step_size_Y
        if (self.grid_parameters.space_step_size_Y % 2) == 0:
            offset = offset + (self.grid_parameters.space_step_size_Y / 2)
        for i in range(0, self.grid_parameters.space_steps_Y):
            self.y_axis[i] = offset + i * self.grid_parameters.space_step_size_Y
        # ---------------------------------------------------------------

    def set_init_function(self, func, normalize=True):
        """
        Sets the initial function.

        :param func: The initial function.
        :param normalize: Set to true if the function should be normalized.
        """

        xv, yv = np.meshgrid(self.x_axis, self.y_axis)

        if normalize:
            self.grid[0] = func(xv, yv) / np.linalg.norm(func(xv, yv))
        else:
            self.grid[0] = func(xv, yv)

    def print_grid(self):
        """
        Prints the grid.
        """

        print(np.flipud(self.grid.round(2)))

    # Functions for plotting the system.
    # ---------------------------------------------------------------

    def plot_3d(self, time_step=0, square=True, use_for_gif=False):
        """
        Funktion prints a 3D wireframe of the system at one specified timestep.

        :param time_step: The timestep which should be plotted.
        :param square: Set to true if the function in the graph should be squared.
        :param use_for_gif: Only for internal use. (Is set to true when generating a gif.)
        """

        # Deactivate Warnings while generating gif
        warnings.filterwarnings('ignore')

        x, y = np.meshgrid(self.x_axis.real, self.y_axis.real)

        if square:
            z = np.square(abs(self.grid[time_step]))
        else:
            z = self.grid[time_step].real

        # fig = plt.figure(figsize=(10, 7), dpi=90)
        ax = plt.axes(projection='3d')
        # ax.contour3D(X, Y, Z, 50, cmap='binary')
        ax.plot_wireframe(x, y, z, color='teal')
        # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        ax.set_xlabel('x-axis')
        ax.set_ylabel('y-axis')
        ax.set_zlabel('wave function squared')
        ax.set_title(self.method + "\ntime-step: " + str(time_step) + " dt=" + str(
            self.grid_parameters.time_step_size) + " dx=" + str(self.grid_parameters.space_step_size_X) + " dy=" + str(
            self.grid_parameters.space_step_size_Y))

        if use_for_gif:
            plt.savefig("./gif/picture_for_giv_" + str(time_step) + ".png", dpi=300)
            plt.close()

        warnings.filterwarnings('default')

    def gif(self):
        """
        Makes a gif of the whole system. Consisting of a 3D plot of every timestep.
        """

        if not os.path.exists("./gif"):
            os.makedirs("./gif")

        for i in range(0, self.grid_parameters.time_steps):
            print("Generating picture: (" + str(i + 1) + "/" + str(self.grid_parameters.time_steps) + ")")
            self.plot_3d(i, True, True)

        images = list()
        for i in range(0, self.grid_parameters.time_steps):
            images.append(imageio.imread("./gif/picture_for_giv_" + str(i) + ".png"))
        imageio.mimsave("simulation.gif", images)

    # Funktion to plot one time-step as a 2D heatmap
    ##############################################################
    def heatmap(self, time_step=0, square=True, save=False):
        """
        Makes a heatmap of the system at a given timestep.

        :param time_step: The timestep which should be plotted.
        :param square: Set to true if the system should be squared.
        :param save: Set true if the plot should be saved instead of shown.
        """

        if square:
            z = np.square(abs(self.grid[time_step]))
        else:
            z = self.grid[time_step].real

        # plt.figure(figsize=(10, 7), dpi=160)
        plt.imshow(z, cmap='viridis', interpolation='nearest')
        plt.xlabel("grid points on the x-axis")
        plt.ylabel("grid points on the y-axis")
        plt.title("time-step: " + str(time_step) + "\ndt=" + str(self.grid_parameters.time_step_size) + " dx=" + str(
            self.grid_parameters.space_step_size_X) + " dy=" + str(self.grid_parameters.space_step_size_Y))
        cbar = plt.colorbar()
        cbar.set_label('squared wave function value')

        if save:

            if not os.path.exists("./pictures"):
                os.makedirs("./pictures")

            plt.savefig("./pictures/2D-heatmap-" + str(time_step) + ".png", dpi=300)
            plt.close()
        else:
            plt.show()
            plt.close()

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

    def plot_3d_potential(self, potential_function, save=False):
        """
        Plots the potential as a 3D graph.

        :param potential_function: The potential function which should be plotted.
        :param save: Set true if the plot should be saved instead of shown.
        """

        x, y = np.meshgrid(self.x_axis.real, self.y_axis.real)
        z = potential_function(x, y)

        # fig = plt.figure(figsize=(10, 7), dpi=90)
        ax = plt.axes(projection='3d')
        # ax.contour3D(X, Y, Z, 50, cmap='binary')
        ax.plot_wireframe(x, y, z, color='coral')
        # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        ax.set_xlabel('x-axis')
        ax.set_ylabel('y-axis')
        ax.set_zlabel('potential')
        ax.set_title('Potential')

        if save:

            if not os.path.exists("./pictures"):
                os.makedirs("./pictures")

            plt.savefig("./pictures/potential.png", dpi=300)
            plt.close()
        else:
            plt.show()
