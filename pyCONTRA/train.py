import argparse
from pyCONTRA.ParameterManager import *
from pyCONTRA.InferenceEngine import *
from pyCONTRA.ComputationEngine import *
from pyCONTRA.ComputationWrapper import *
from pyCONTRA.OptimizationWrapper import *

def RunTrainingMode(args: argparse.Namespace, description: list):
    parameter_manager = ParameterManager()
    inference_engine = InferenceEngine(args.allow_noncomplementary, args.num_data_sources)
    inference_engine.RegisterParameters(parameter_manager)
    computation_engine = ComputationEngine(args, description, inference_engine, parameter_manager)
    computation_wrapper = ComputationWrapper(computation_engine)

    if(computation_engine.IsComputeNode()):
        computation_engine.RunAsComputeNode();
        return
    
    initweights_filename = args.train_initweights_filename
    priorweights_filename = args.train_priorweights_filename
    train_max_iter = args.train_max_iter
    hyperparam_data = args.hyperparam_data

    w = list()
    if(initweights_filename!=""):
        parameter_manager.ReadFromFile(initweights_filename, w)
    else:
        for i in range(0, parameter_manager.GetNumLogicalParameters()):
            w.append(0)
    
    w0 = list()
    if(priorweights_filename!=""):
        parameter_manager.ReadFromFile(priorweights_filename, w0)
    else:
        for i in range(0, parameter_manager.GetNumLogicalParameters()):
            w0.append(0)
    
    units = computation_wrapper.FilterNonparsable(computation_wrapper.GetAllUnits())
    optimization_wrapper = OptimizationWrapper()

    if(args.holdout_ratio <= 0):
        regularization_coefficients = list()
        if(args.training_mode=="em"):
            optimization_wrapper.TrainEM(units, w, regularization_coefficients, train_max_iter)
        elif args.training_mode=="em-sgd":
            print(parameter_manager.GetNumParameterGroups())
            if parameter_manager.GetNumParameterGroups()!=2:
                raise Exception("Using em-sgd with multiple hyperparameters is not supported")
            regularization_coefficients[1] = 0
            optimization_wrapper.TrainSGD(units,w, regularization_coefficients)
        else:
            optimization_wrapper.Train(units, w, w0, regularization_coefficients)
    else:
        if args.training_mode == "em":
            optimization_wrapper.LearnHyperparametersEM(units, w, train_max_iter)
        else:
            optimization_wrapper.LearnHyperparameters(units, w)
    parameter_manager.WriteToFile("optimize.params.final", w)
    computation_engine.StopComputeNodes()
    