from pyCONTRA.Utilities import *
from pyCONTRA.LineSearch import *

class LBFGS:

    def __init__(self, M: int = 20, TERMINATION_RATIO: float = 1e-5, MAX_ITERATIONS: int = 1000, SMALL_STEP_RATIO: float = 1e-5, MAX_SMALL_STEPS: int = 3, MAX_STEP_NORM:float = 1e10):
        self.M = M
        self.TERMINATION_RATIO = TERMINATION_RATIO
        self.MAX_ITERATIONS = MAX_ITERATIONS
        self.SMALL_STEP_RATIO = SMALL_STEP_RATIO
        self.MAX_SMALL_STEPS = MAX_SMALL_STEPS
        self.MAX_STEP_NORM = MAX_STEP_NORM

    def Minimise(self,x0: list):
        n = len(x0)
        f = list()
        gamma = list()
        x = [list() for i in range(2)]
        x = [list() for i in range(2)]
        x = [list() for i in range(M)]
        x = [list() for i in range(M)]


