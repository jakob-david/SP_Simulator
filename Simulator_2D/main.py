from Simulator_2D.GridParameters import GridParameters
from Simulator_2D.Grid import Grid
from Simulator_2D.Simulation import Simulation

grid_parameters = GridParameters()
grid = Grid(grid_parameters)

my_sym = Simulation(grid, grid_parameters.input_gravity, grid_parameters.potential_function)
my_sym.initFunction(grid_parameters.initial_function)

my_sym.SO_Poisson_2D()

my_sym.gif()

my_sym.heatmap(0, True, True)
my_sym.heatmap(100, True, True)
my_sym.heatmap(200, True, True)
my_sym.heatmap(300, True, True)
my_sym.heatmap(400, True, True)
my_sym.heatmap(500, True, True)
# my_sym.plotEnergyEvolution()
