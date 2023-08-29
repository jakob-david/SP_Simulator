from Simulator_1D.GridParameters import GridParameters
from Simulator_1D.Grid import Grid
from Simulator_1D.Simulation import Simulation

grid_parameters = GridParameters()
grid = Grid(grid_parameters)



my_sym = Simulation(grid, grid_parameters.potential_function)
my_sym.initFunction(grid_parameters.initial_function)

my_sym.SO_Poisson()
my_sym.heatmap(True)
#my_sym.plotEnergyEvolution(True)
#my_sym.gif()