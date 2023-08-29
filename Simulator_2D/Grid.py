import numpy as np
import matplotlib.pyplot as plt
import warnings
import imageio.v2 as imageio
import os


class Grid:
    def __init__(self, GridParameters):
        self.grid_parameters = GridParameters
        self.grid = np.zeros(
            (self.grid_parameters.time_steps, self.grid_parameters.space_steps_X, self.grid_parameters.space_steps_Y),
            dtype=np.complex_)
        self.energy = np.zeros(self.grid_parameters.time_steps)
        self.method = ""

        # Precalculate the x-axis since it is used quite often
        ########################################################
        self.x_axis = np.zeros(self.grid_parameters.space_steps_X, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps_X / 2) * self.grid_parameters.space_step_size_X
        if ((self.grid_parameters.space_steps_X % 2) == 0):
            offset = offset + (self.grid_parameters.space_step_size_X / 2)
        for i in range(0, self.grid_parameters.space_steps_X):
            self.x_axis[i] = offset + i * self.grid_parameters.space_step_size_X

        # Precalculate the y-axis since it is used quite often
        ########################################################
        self.y_axis = np.zeros(self.grid_parameters.space_steps_Y, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps_Y / 2) * self.grid_parameters.space_step_size_Y
        if ((self.grid_parameters.space_step_size_Y % 2) == 0):
            offset = offset + (self.grid_parameters.space_step_size_Y / 2)
        for i in range(0, self.grid_parameters.space_steps_Y):
            self.y_axis[i] = offset + i * self.grid_parameters.space_step_size_Y

    def initFunction(self, func, normalize=True):

        xv, yv = np.meshgrid(self.x_axis, self.y_axis)

        if (normalize):
            self.grid[0] = func(xv, yv) / np.linalg.norm(func(xv, yv))
        else:
            self.grid[0] = func(xv, yv)

    def printGrid(self):
        print(np.flipud(self.grid.round(2)))

    ##############################################################
    ##############################################################
    #
    # Functions for plotting the System
    #
    ##############################################################
    ##############################################################

    # Funktion to one time-step as 3D wireframe
    ##############################################################
    def plot3D(self, time_step=0, square=True, use_for_gif=False):

        # Deactivate Warnings while generating gif
        warnings.filterwarnings('ignore')

        X, Y = np.meshgrid(self.x_axis.real, self.y_axis.real)

        if (square):
            Z = np.square(abs(self.grid[time_step]))
        else:
            Z = self.grid[time_step].real

        fig = plt.figure(figsize=(10, 7), dpi=90)
        ax = plt.axes(projection='3d')
        # ax.contour3D(X, Y, Z, 50, cmap='binary')
        ax.plot_wireframe(X, Y, Z, color='teal')
        # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        ax.set_xlabel('x-axis')
        ax.set_ylabel('y-axis')
        ax.set_zlabel('wave function squared')
        ax.set_title(self.method + "\ntime-step: " + str(time_step) + " dt=" + str(
            self.grid_parameters.time_step_size) + " dx=" + str(self.grid_parameters.space_step_size_X) + " dy=" + str(
            self.grid_parameters.space_step_size_Y))

        if (use_for_gif):
            plt.savefig("./gif/picture_for_giv_" + str(time_step) + ".png", dpi=300)
            plt.close()

        warnings.filterwarnings('default')

    # Funktion to generate a gif out of the 3D wireframes.
    ##############################################################
    def gif(self):

        if not os.path.exists("./gif"):
            os.makedirs("./gif")

        for i in range(0, self.grid_parameters.time_steps):
            print("Generating picture: (" + str(i + 1) + "/" + str(self.grid_parameters.time_steps) + ")")
            self.plot3D(i, True, True)

        images = list()
        for i in range(0, self.grid_parameters.time_steps):
            images.append(imageio.imread("./gif/picture_for_giv_" + str(i) + ".png"))
        imageio.mimsave("simulation.gif", images)
        #files.download('simulation.gif')   TODO: Make download.

    # Funktion to plot one time-step as a 2D heatmap
    ##############################################################
    def heatmap(self, time_step=0, square=True, save=False):

        if (square):
            Z = np.square(abs(self.grid[time_step]))
        else:
            Z = self.grid[time_step].real

        # plt.figure(figsize=(10, 7), dpi=160)
        plt.imshow(Z, cmap='viridis', interpolation='nearest')
        plt.xlabel("grid points on the x-axis")
        plt.ylabel("grid points on the y-axis")
        plt.title("time-step: " + str(time_step) + "\ndt=" + str(self.grid_parameters.time_step_size) + " dx=" + str(
            self.grid_parameters.space_step_size_X) + " dy=" + str(self.grid_parameters.space_step_size_Y))
        cbar = plt.colorbar()
        cbar.set_label('squared wave function value')

        if (save):
            plt.savefig("2D-heatmap-" + str(time_step) + ".png", dpi=300)
            #files.download("2D-heatmap-" + str(time_step) + ".png") TODO: Make download.
            plt.close()
        else:
            plt.show()
            plt.close()

    # Funktion to plot the energy evolution of the system
    ##############################################################
    def plotEnergyEvolution(self, save=False, log=False):

        time_axis = np.zeros(self.grid_parameters.time_steps)
        for i in range(0, self.grid_parameters.time_steps):
            time_axis[i] = i * self.grid_parameters.time_step_size

        if (log):
            plt.plot(time_axis, self.energy)
            plt.yscale("log")
        else:
            plt.ylim([-3, 3])  # kind of
            plt.plot(time_axis, self.energy)

        plt.xlabel("time")
        plt.ylabel("energy")
        plt.title("Energy Evolution")

        if (save):
            plt.savefig("energy.png", dpi=300)
            #files.download("energy.png") TODO: Make download.
            plt.close()
        else:
            plt.show()

    # Funktion to plot the potential as a 3D wireframe
    ##############################################################
    def plot3DPotential(self, potential_function, download=False):

        X, Y = np.meshgrid(self.x_axis.real, self.y_axis.real)
        Z = potential_function(X, Y)

        fig = plt.figure(figsize=(10, 7), dpi=90)
        ax = plt.axes(projection='3d')
        # ax.contour3D(X, Y, Z, 50, cmap='binary')
        ax.plot_wireframe(X, Y, Z, color='coral')
        # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        ax.set_xlabel('x-axis')
        ax.set_ylabel('y-axis')
        ax.set_zlabel('potential')
        ax.set_title('Potential')

        if (download):
            plt.savefig("potential.png", dpi=300)
            #files.download("potential.png") TODO: Make download.
            plt.close()
        else:
            plt.show()
