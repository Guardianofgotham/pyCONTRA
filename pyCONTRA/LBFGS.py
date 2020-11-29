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
        f = [0.0, 0.0]
        gamma = [0.0, 0.0]
        x = [[0.0 for j in range(n)] for i in range(2)]
        g = [[0.0 for j in range(n)] for i in range(2)]
        s = [[0.0 for j in range(n)] for i in range(self.M)]
        y = [[0.0 for j in range(n)] for i in range(self.M)]
        rho = list()
        gradient_ratio = None
        f0 = None

        x[0] = x0
        f[0] = f0 = ComputeFunction(x[0])

        if(f[0] > 1e20):
            Report
            return f[0]
        
        ComputeGradient(g[0], x[0])
        gradient_ratio = Norm(g[0])/max(1,Norm(x[0]))



