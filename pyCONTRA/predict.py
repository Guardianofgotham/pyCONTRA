import argparse
from pyCONTRA.ParameterManager import *
from pyCONTRA.InferenceEngine import *
from pyCONTRA.ComputationEngine import *
from pyCONTRA.ComputationWrapper import *


def RunPredictionMode(args: argparse.Namespace, description: list):
    parameter_manager = ParameterManager()
    inference_engine = InferenceEngine()
    inference_engine.RegisterParameters(parameter_manager)
    computation_engine = ComputationEngine()
    computation_wrapper = ComputationWrapper()

    if(computation_engine.IsComputeNode()):
        computation_engine.RunAsComputeNode()

    output_parens_destination = args.output_parens_destination
    output_bpseq_destination = args.output_bpseq_destination
    output_posteriors_destination = args.output_posteriors_destination

    w = list()
    if(args.parameter_filename != ""):
        parameter_manager.ReadFromFile(args.parameter_filename, w)
    else:
        if PROFILE != 0:
            pass
        else:
            if(args.allow_noncomplementary):
                pass
            else:
                pass

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
        computation_wrapper.Predict(computation_wrapper, w, args.gamma, args.log_base)  
    computation_engine.StopComputeNodes()