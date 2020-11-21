import argparse
from pyCONTRA.ParameterManager import *
from pyCONTRA.InferenceEngine import *
from pyCONTRA.ComputationEngine import *
from pyCONTRA.ComputationWrapper import *

def RunGradientSanityCheck(args: argparse.Namespace, description: list):
    ## The architecture of the code is somewhat complicated here, so
    ## here's a quick explanation:
    ## 
    ##    ParameterManager: associates each parameter of the model
    ##                      with a name and manages hyperparameter
    ##                      groups
    ##                     
    ##    InferenceEngine: performs application-specific
    ##                     (loss-augmented) inference
    ##
    ##    ComputationEngine: makes all necessary calls to dynamic
    ##                       programming routines for processing
    ##                       individual sequences and interfaces with
    ##                       distributed computation module
    ##
    ##    ComputationWrapper: provides a high-level interface for
    ##                        performing computations on groups of
    ##                        sequences
    ##
    ##    OuterOptimizationWrapper # InnerOptimizationWrapper:
    ##                        interface between computation routines
    ##                        and optimization routines
    parameter_manager = ParameterManager()
    inference_engine = InferenceEngine()
    inference_engine.RegisterParameters(parameter_manager)
    computation_engine = ComputationEngine()
    computation_wrapper = ComputationWrapper()

    if(computation_engine.IsComputeNode()):
        computation_engine.RunAsComputeNode()
    
    w = list()

    initweights_filename = args.train_initweights_filename
    if(initweights_filename!=""):
        parameter_manager.ReadFromFile(initweights_filename, w)
    else:
        for i in range(0, parameter_manager.GetNumLogicalParameters()):
            w.append(0)
    computation_wrapper.SanityCheckGradient(computation_wrapper.GetAllUnits(), w)
    computation_engine.StopComputeNodes()
    print("inside SanityGradientCheck")