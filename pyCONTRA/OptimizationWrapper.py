from pyCONTRA.ComputationWrapper import *


class OptimizationWrapper:
    def __init__(self, optimization_wrapper: OptimizationWrapper, initial_w: list, training: list, holdout:list):
        self.
        

    def OptimizationWrapper(self, computation_wrapper):
        raise Exception("Not implemented")

    def Train(self,  units,  w,  w0,   C):
        raise Exception("Not implemented")

    def TrainEM(self,  units,  w,   C,   train_max_iter):
        raise Exception("Not implemented")

    def TrainSGD(self,  units,  w,   C):
        raise Exception("Not implemented")

    def LearnHyperparameters(self, units,  values):
        raise Exception("Not implemented")

    def LearnHyperparametersEM(self, units,  values,   train_max_iter):
        raise Exception("Not implemented")

    def Indent(self,):
        raise Exception("Not implemented")

    def Unindent(self,):
        raise Exception("Not implemented")

    def PrMessage(self, s):
        raise Exception("Not implemented")

    # getters
    def GetOptions(self):
        raise Exception("Not implemented")

    def GetDescriptions(self):
        raise Exception("Not implemented")

    def GetInferenceEngine(self):
        raise Exception("Not implemented")

    def GetParameterManager(self):
        raise Exception("Not implemented")

    def GetComputationEngine(self):
        raise Exception("Not implemented")

    def GetComputationWrapper(self):
        raise Exception("Not implemented")
