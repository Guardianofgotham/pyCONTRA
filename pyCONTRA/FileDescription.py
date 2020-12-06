from pyCONTRA.SStruct import *
import math
class FileDescription:
    def __init__(self, input_filename: str, allow_noncomplementary: bool, num_data_sources: int):
        self.sstruct = SStruct(input_filename, num_data_sources)
        self.input_filename = input_filename
        self.size = math.pow(self.sstruct.GetLength(), 3.0)
        self.weight = 1.0
        if not allow_noncomplementary:
            self.sstruct.RemoveNoncomplementaryPairings()
        if (self.sstruct.GetNumSequences() > 1):
            print("Folding multiple input sequences without --profile mode enabled.");

    