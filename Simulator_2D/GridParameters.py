from dataclasses import dataclass
from scipy import stats


@dataclass
class GridParameters:
    """
    This class holds the parameters for the simulation.
    """

    time_steps: int = 501         # t=0 is considered as the first step
    time_step_size: float = 0.1
    space_steps_X: int = 190
    space_step_size_X: float = 0.3
    space_steps_Y: int = 190
    space_step_size_Y: float = 0.3

    input_gravity = 10

    @staticmethod
    def initial_function(x, y):
        """
        The initial function of the system.

        :param x: One parameter the function is dependent of.
        :param y: A second parameter the function is dependent of.
        """

        return stats.norm.pdf(x, -4, 1) * stats.norm.pdf(y, -4, 1) + stats.norm.pdf(x, 4, 1) * stats.norm.pdf(y, 4, 1)

    @staticmethod
    def potential_function(x, y, u=0):
        """
        The potential function of the system.

        :param x: One parameter the function is dependent of.
        :param y: A second parameter the function is dependent of.
        :param u: A third parameter the function is dependent of.
        """

        z = 0*x+0*y+u   # .05*x**2 + 0.05*y**2# + u
        return z
