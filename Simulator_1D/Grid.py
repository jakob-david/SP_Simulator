import numpy as np
import matplotlib.pyplot as plt
import warnings
import imageio.v2 as imageio


class Grid:
    def __init__(self, GridParameters):
        self.grid_parameters = GridParameters
        self.grid = np.zeros((self.grid_parameters.time_steps, self.grid_parameters.space_steps), dtype=np.complex_)
        self.energy = np.zeros(self.grid_parameters.time_steps)
        self.method = ""

        # Precalculate the x-axis since it is used quite often
        ########################################################
        self.x_axis = np.zeros(self.grid_parameters.space_steps, dtype=np.complex_)
        offset = -1 * int(self.grid_parameters.space_steps / 2) * self.grid_parameters.space_step_size
        if (self.grid_parameters.space_steps % 2) == 0:
            offset = offset + (self.grid_parameters.space_step_size / 2)
        for i in range(0, self.grid_parameters.space_steps):
            self.x_axis[i] = offset + i * self.grid_parameters.space_step_size

    def initFunction(self, func, normalize=True):
        if normalize:
            self.grid[0] = func(self.x_axis) / np.linalg.norm(func(self.x_axis))
        else:
            self.grid[0] = func(self.x_axis)

    def printGrid(self):
        print(np.flipud(self.grid.round(2)))

    ##############################################################
    ##############################################################
    #
    # Functions for plotting the System
    #
    ##############################################################
    ##############################################################

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
            # files.download("energy.png") TODO: file saving
            plt.close()
        else:
            plt.show()

    # Funktion to plot a 2D representation of one time-step
    ##############################################################
    def plot2D(self, time_step, pot=lambda a, u: 0, use_for_gif=False):

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

        if (use_for_gif):
            plt.savefig("picture_for_giv_" + str(time_step) + ".png")
            plt.close()

        warnings.filterwarnings('default')

    # Funktion to plot a 3D wireframe of the whole simulation.
    ##############################################################
    def plot3D(self, square=True):

        # Deactivate Warnings while generating gif
        warnings.filterwarnings('ignore')

        y_axis = np.zeros(self.grid_parameters.time_steps)
        for i in range(0, self.grid_parameters.time_steps):
            y_axis[i] = i * self.grid_parameters.time_step_size

        X, Y = np.meshgrid(self.x_axis.real, y_axis)

        if (square):
            Z = np.square(abs(self.grid))
        else:
            Z = self.grid.real

        ax = plt.axes(projection='3d')
        ax.plot_wireframe(X, Y, Z, color='green')
        ax.set_xlabel('space')
        ax.set_ylabel('time')
        ax.set_zlabel('wave function')

        warnings.filterwarnings('default')

    # Funktion to plot a 2D heatmap of the whole simulation.
    ##############################################################
    def heatmap(self, square=True, save=False):

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
            plt.savefig("2D-heatmap.png", dpi=300)
            # files.download("2D-heatmap.png") TODO: file saving
            plt.close()
        else:
            plt.show()
            plt.close()

    # Funktion to create a gif from the 2D-plots of each time-step
    ##############################################################
    def gif(self, pot=lambda a: 0):

        for i in range(0, self.grid_parameters.time_steps):
            print("Generating picture: (" + str(i + 1) + "/" + str(self.grid_parameters.time_steps) + ")")
            self.plot2D(i, pot, True)

        images = list()
        for i in range(0, self.grid_parameters.time_steps):
            images.append(imageio.imread("picture_for_giv_" + str(i) + ".png"))
        imageio.mimsave("simulation.gif", images)
        # files.download('simulation.gif') TODO: file saving
