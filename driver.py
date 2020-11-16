# Driver Test Code in here
from pyCONTRA import getArgumentsObject

args = getArgumentsObject()
if(args.training_mode=="predict"):
    print("Basil: AI Champion")
elif(args.training_mode=="train"):
    print("Ritvik: Statistics God")