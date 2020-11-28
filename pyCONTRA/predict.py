'''Corresponds to Contrafold.cpp'''
import argparse
from pyCONTRA.ParameterManager import *
from pyCONTRA.InferenceEngine import *
from pyCONTRA.ComputationEngine import *
from pyCONTRA.ComputationWrapper import *
from pyCONTRA.Defaults import *


def RunPredictionMode(args: argparse.Namespace, description: list):
    parameter_manager = ParameterManager()
    inference_engine = InferenceEngine()
    inference_engine.RegisterParameters(parameter_manager)
    computation_engine = ComputationEngine(args, description, inference_engine, parameter_manager)
    computation_wrapper = ComputationWrapper(computation_engine)

    if(computation_engine.IsComputeNode()):
        computation_engine.RunAsComputeNode()
        return

    output_parens_destination = args.output_parens_destination
    output_bpseq_destination = args.output_bpseq_destination
    output_posteriors_destination = args.output_posteriors_destination

    w = list()
    if(args.parameter_filename != ""):
        parameter_manager.ReadFromFile(args.parameter_filename, w)
    else:
        if(args.allow_noncomplementary):
            w = GetDefaultNoncomplementaryValues()
        else:
            w = GetDefaultComplementaryValues()
    

    if args.gamma < 0:
        if output_parens_destination != "":
            pass
        if output_bpseq_destination != "":
            pass
        if output_posteriors_destination != "":
            pass

        for k in range(-5, 11):
            gamma = 2.0**float(k)

            if(len(description) > 1):
                if output_parens_destination!="":
                    pass
                if output_bpseq_destination!="":
                    pass
                if output_posteriors_destination!="":
                    pass
            computation_wrapper.Predict(computation_engine, w, gamma, args.log_base)
    else:
        if len(description)>1:
            if output_parens_destination!="":
                pass
            if output_bpseq_destination!="":
                pass
            if output_posteriors_destination!="":
                pass
        computation_wrapper.Predict(computation_wrapper.GetAllUnits(), w, args.gamma, args.log_base) 
    print("Inside Prediction")
    computation_engine.StopComputeNodes()