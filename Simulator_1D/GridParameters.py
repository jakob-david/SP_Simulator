from dataclasses import dataclass
from scipy import stats


@dataclass
class GridParameters:
    """This class holds some parameters for the simulation."""

    time_steps: int = 400  # t=0 is considered as the first step
    time_step_size: float = 0.5
    space_steps: int = 500
    space_step_size: float = 0.05

    def initial_function(self, x):
        return stats.norm.pdf(x, -5, 1) + stats.norm.pdf(x, 5, 1)

    def potential_function(self, x, u=0):
        y = 0  # .05*x**2 # + u

        return y
