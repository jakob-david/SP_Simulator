from dataclasses import dataclass
from scipy import stats


@dataclass
class GridParameters:
    time_steps: int = 501         # t=0 is considered as the first step
    time_step_size: float = 0.1
    space_steps_X: int = 190
    space_step_size_X: float = 0.3
    space_steps_Y: int = 190
    space_step_size_Y: float = 0.3

    input_gravity = 10


    def initial_function(self, x, y):
        return stats.norm.pdf(x, -4, 1) * stats.norm.pdf(y, -4, 1) + stats.norm.pdf(x, 4, 1) * stats.norm.pdf(y, 4, 1)

    def potential_function(x, y, u=0):
        z = 0   # .05*x**2 + 0.05*y**2# + u
        return z
