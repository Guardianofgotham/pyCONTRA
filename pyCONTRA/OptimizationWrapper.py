from pyCONTRA.ComputationWrapper import *


class OptimizationWrapper:
    def __init__(self):
        pass

    def OptimizationWrapper(self, computation_wrapper):
        pass

    def Train(self,  units,  w,  w0,   C):
        pass

    def TrainEM(self,  units,  w,   C,   train_max_iter):
        pass

    def TrainSGD(self,  units,  w,   C):
        pass

    def LearnHyperparameters(self, units,  values):
        pass

    def LearnHyperparametersEM(self, units,  values,   train_max_iter):
        pass

    def Indent(self,):
        pass

    def Unindent(self,):
        pass

    def PrMessage(self, s):
        pass

    # getters
    def GetOptions(self):
        pass

    def GetDescriptions(self):
        pass

    def GetInferenceEngine(self):
        pass

    def GetParameterManager(self):
        pass

    def GetComputationEngine(self):
        pass

    def GetComputationWrapper(self):
        pass
