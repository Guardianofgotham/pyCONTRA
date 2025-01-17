# Driver Test Code in here
from pyCONTRA import *
# import pdb

# pdb.set_trace()
args, fileNames = getArgumentsObject()
print(f"Training mode: {args.training_mode}")
print(f"Use constraints: {args.use_constraints}")
print(f"Use evidence: {args.use_evidence}")

description = list()
MakeFileDescriptions(args, fileNames, description)

if(args.gradient_sanity_check):
    RunGradientSanityCheck(args, description)
elif(args.training_mode == "supervised" or args.training_mode == "em" or args.training_mode == "em-sgd"):
    RunTrainingMode(args, description)
else:
    RunPredictionMode(args, description)
