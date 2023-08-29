from dataclasses import dataclass
from scipy import stats


@dataclass
class GridParameters:
    """
    This class holds the parameters for the simulation.
    """

    time_steps: int = 400  # t=0 is considered as the first step
    time_step_size: float = 0.5
    space_steps: int = 500
    space_step_size: float = 0.05

    @staticmethod
    def initial_function(x, d=0):
        """
        The initial function of the system.

        :param x: One parameter the function is dependent of.
        :param d: A second parameter the function is dependent of.
        """

        return stats.norm.pdf(x, -5, 1) + stats.norm.pdf(x, 5, 1) + d

    @staticmethod
    def potential_function(x=0, u=0):
        """
        The potential function of the system.

        :param x: One parameter the function is dependent of.
        :param u: A second parameter the function is dependent of.
        """

        return 0*x+u    # y = 0.05*x**2 # + u
