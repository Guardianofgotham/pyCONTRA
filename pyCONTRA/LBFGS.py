from pyCONTRA.Utilities import *
from pyCONTRA.LineSearch import *
import math

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
        f = [0, 0]
        gamma = [0, 0]
        x = [[0 for j in range(n)] for i in range(2)]
        g = [[0 for j in range(n)] for i in range(2)]
        s = [[0 for j in range(n)] for i in range(self.M)]
        y = [[0 for j in range(n)] for i in range(self.M)]
        rho = [0 for i in range(M)]
        gradient_ratio = None
        f0 = None

        x[0] = x0
        f[0] = f0 = ComputeFunction(x[0])

        if(f[0] > 1e20):
            # Report     can just print
            return f[0]
        
        ComputeGradient(g[0], x[0])
        gradient_ratio = Norm(g[0])/max(1,Norm(x[0]))
        if(gradient_ratio < self.TERMINATION_RATIO):
            print("Termination before optimization: gradient vector small {} < {}".format(gradient_ratio, self.TERMINATION_RATIO))
            return f[0]
        gamma[0] = 1/Norm(g[0]) #what is Real

        print (0, x[0], f[0], 0)

        progress_mode = False
        num_consecutive_small_steps = 0
        k = 0

        while(True):
            d = [0 for i in range(-g[k%2])]
            a = [0 for i in range(self.M)]
            for i in range(k-1,k-self.M,-1):
                a[(i+self.M) % self.M] = rho[(i+self.M) % self.M] * DotProduct(s[(i+self.M) % self.M], d)
                d += (a[(i+self.M) % self.M] - b) * s[(i+self.M) % self.M]
            
            step = DoLineSearch(x[k%2], f[k%2], g[k%2], d, x[(k+1)%2], f[(k+1)%2], g[(k+1)%2], 0, min(10, (self.MAX_STEP_NORM/max(1,Norm(d)))))

            print(k+1, x[(k+1)%2], f[(k+1)%2], step)

            if(k+1 >= self.MAX_ITERATIONS):
                print("Termination condition: maximum number of iterations reached")
                break

            gradient_ratio = Norm(g[(k+1)%2]) / max(1, Nor(x[k+1]%2))
            if(gradient_ratio < self.TERMINATION_RATIO):
                print("Termination condition: gradient vector small ({} < {})".format(gradient_ratio, self.TERMINATION_RATIO))
                break
            if(step == 0):
                num_consecutive_small_steps = self.MAX_SMALL_STEPS
            elif((f[k%2] - f[(k+1)%2]) / max(1, f0 - f[(k+1)%2]) < self.SMALL_STEP_RATIO):
                num_consecutive_small_steps+=1
            else:
                num_consecutive_small_steps = 0
                progress_mode = True
            if(num_consecutive_small_steps == self.MAX_SMALL_STEPS):
                if(M > 0 and progress_mode):
                    progress_mode = False
                    num_consecutive_small_steps = 0
                    print("Restart: Too many consecutive small steps")

                    for i in range(self.M):
                        s = [0 for i in range(len(s))]
                        y = [0 for i in range(len(y))]
                        rho[i] = 0
                else:
                    print("Termination: Too many consecutive small steps")
                    break
            
            s[k % self.M] = x[(k+1) % 2] - x[k % 2]
            y[k % self.M] = g[(k+1) % 2] - g[k % 2]
            rho[k % self.M] = 1 / DotProduct(y[k % self.M], s[k % self.M])


            if(!math.isfinite(rho[k%self.M]) or rho[k%M] <= 0):
                s = [0 for i in range(len(s))]
                y = [0 for i in range(len(y))]
                rho[k%self.M] = 0
            
            gamma[(k+1) % 2] = DotProduct(s[(k-1+self.M) %
                                            self.M], y[(k-1+self.M) % self.M])
            
            if (!math.isfinite(gamma[(k+1) % 2])):
                gamma[(k+1)%2] = gamma[k%2]
            
            k+=1    # recheck
        x0 = x[(k+1)%2]
        return f[(k+1)%2]

        






