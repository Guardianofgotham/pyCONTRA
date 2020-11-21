# Driver Test Code in here
from pyCONTRA import *

args = getArgumentsObject()
print(f"Training mode: {args.training_mode}")
print(f"Use constraints: {args.use_constraints}")
print(f"Use evidence: {args.use_evidence}")

description = list()
# MakeFileDescription

if(args.gradient_sanity_check):
    RunGradientSanityCheck(args)
elif(args.training_mode!=""):
    RunTrainingMode(args)
else:
    RunPredictionMode(args, description)