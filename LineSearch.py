from pyCONTRA.Utilities import *
import math

class LineSearch:

    def __init__(self, T_INIT: int = 1,  MU: float = 0.001,  MIN_IMPROVEMENT_RATIO: float = 0.1,  MAX_EVALUATIONS: int = 10,  GAMMA1: float = 0.01,  GAMMA2: float = 0.8):
        self.T_INIT = T_INIT                                  # initial step size
        self.MU = MU                                      # phi(t) <= phi(0) + MU * phi'(0) t      (sufficient decrease)
        self.MIN_IMPROVEMENT_RATIO = MIN_IMPROVEMENT_RATIO,                   # minimum proportion of overall improvement required to keep going
        self.MAX_EVALUATIONS = MAX_EVALUATIONS,                          # maximum number of function evaluations
        self.GAMMA1 = GAMMA1                                 # maximum step length shrinkage 
        self.GAMMA2 = GAMMA2                            # minimum step length shrinkage 


    def DoLineSearch(self, x: list, f: float, g: list, d: list, new_x: list, new_f: float, new_g: list, T_MIN: float, T_MAX: float):

        assert (T_MIN <= T_MAX), "Line search called with T_MIN > T_MAX."
        dot_prod = DotProduct(d, g)
        sufficient_decrease = False

        t_best = 0 
        t_last = self.T_INIT 
        t_prev = 0
        t_last = max(T_MIN, t_last)
        t_last = min(T_MAX, t_last)
        f_best = f 
        f_last = ComputeFunction(x + t_last * d) 
        f_prev = 0

        if (f_last < f_best): 
            f_best = f_last
            t_best = f_last
        if (f_last <= f + self.MU * t_last * dot_prod):
            sufficient_decrease = True
        
        iteration = 2
        for iteration in range(2, self.MAX_EVALUATIONS+1):

            if (sufficient_decrease and iteration > 2 and (f_prev - f_last) / max(1, f - f_best) < self.MIN_IMPROVEMENT_RATIO):
                break

            if iteration == 2:
                a = ((f_last - f) / t_last - dot_prod) / t_last
                b = dot_prod
                t_new = -b / 2*a

            else:
                c = dot_prod
                d = f

                f_l = f_last - d - c * t_last
                f_p = f_prev - d - c * t_prev

                a = (f_l / (t_last * t_last) - f_p / (t_prev * t_prev)) / (t_last - t_prev)
                b = (f_l / (t_last * t_last * t_last) - f_p / (t_prev * t_prev * t_prev)) / (1 / t_last - 1 / t_prev)

                A = 3*a
                B = 2*b
                C = c

                t_new = (-B + math.sqrt(B*B - 4*A*C)) / (2*A);  # pick the left root

            if (not math.isfinite(t_new)):
                t_new = t_last

            lower_bound = T_MIN
            upper_bound = T_MAX

            if iteration > 3:
                lower_bound = max(lower_bound, * t_last)
                upper_bound = min(upper_bound, * t_last)


            if (lower_bound > upper_bound):
                break

            t_new = max(lower_bound, t_new)
            t_new = min(upper_bound, t_new)

            f_new = ComputeFunction(x + t_new * d)

            if (f_new < f_best): 
                f_best = f_new
                t_best = t_new

            if (f_new <= f + self.MU * t_new * dot_prod):
                sufficient_decrease = True

            t_prev = t_last; f_prev = f_last
            t_last = t_new; f_last = f_new


        new_f = f_best
        new_x = x + t_best * d
        ComputeGradient(new_g, new_x)
        return t_best;    

