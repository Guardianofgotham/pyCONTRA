from pyCONTRA.DistributedComputation import *

class ComputationEngine(DistributedComputationBase):
    def __init__(self):
        pass

    def DistributedComputation(self):
        pass

    # routine for performing an individual work unit
    def DoComputation(self, result, shared, nonshared):
        pass

    # methods to act on individual work units
    def CheckParsability(self, result, nonshared):
        pass

    def ComputeSolutionNormBound(self, result,   shared,   nonshared):
        pass

    def ComputeGradientNormBound(self, result,   nonshared):
        pass

    def ComputeLoss(self, result,   shared,   nonshared):
        pass

    def ComputeFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeMStepFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeGammaMLEFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeHessianVectorProduct(self, result,   shared,   nonshared):
        pass

    def Predict(self, result,   shared,   nonshared):
        pass

    def CheckZerosInData(self, result,   shared,   nonshared):
        pass

    def ComputeGammaMLEScalingFactor(self, result,   shared,   nonshared):
        pass

    def ComputeFunctionAndGradientSE(self, result,   shared,   nonshared, need_gradient):
        pass
