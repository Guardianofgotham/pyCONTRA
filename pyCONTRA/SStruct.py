class SStruct:
    UNPAIRED = -1
    UNKNOWN = -1
    UNKNOWN_POTENTIAL = -1.0

    def __init__(self):
        self.names = list()
        self.sequences = list()
        self.mapping = list()
        self.unpaired_potentials = list()
        self.has_struct = False
        self.has_evidence = False
        self.num_data_sources = 1
        self.which_evidence = list()

    def AnalyzeFormat(self, filename: str):
        pass
    def LoadFASTA(self, filename: str):   
        pass
    def LoadRAW(self, filename: str):   
        pass
    def LoadBPSEQ(self, filename: str):   
        pass
    def LoadBPP2SEQ(self, filename: str):   
        pass
    def LoadBPP2TSEQ(self, filename: str):   
        pass
    def FilterSequence(self, sequence: str):   
        pass
    def FilterParens(self, sequence: str):   
        pass
    def ConvertParensToMapping(self, parens: str):   
        pass
    def ConvertMappingToParens(self, mapping: list):   
        pass
    def ValidateMapping(self, mapping: list):   
        pass
    def Load(self, filename: str):   
        pass
    def ContainsPseudoknots(self):   
        pass
    def RemoveNoncomplementaryPairings(self, seq=0):   
        pass
    def WriteBPSEQ(self, outfile, seq=0):   
        pass
    def WriteParens(self, outfile):   
        pass
    def ComputePercentIdentity(self):   
        pass
    def ComputePositionBasedSequenceWeights(self, ):   
        pass
    def SetMapping(self, mapping):   
        pass
    def GetNames(self):   
        pass
    def GetSequence(self):
        pass
    def GetMapping(self):
        pass
    def GetUnpairedPotential(self, which: int):
        pass
    def GetPairedPotentials(self, which: int):
        pass
    def GetLength(self):
        pass
    def GetNumSequences(self):
        pass
    def HasStruct(self):
        pass
    def HasEvidence(self):
        pass