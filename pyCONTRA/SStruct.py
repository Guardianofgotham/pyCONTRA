from pyCONTRA.Utilities import *

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
    
    def Load(self, filename: str):
        FileFormat = AnalyzeFormat(filename)
        if(FileFormat=="FASTA"):
            LoadFasta(filename)
        elif(FileFormat=="BPSEQ"):
            LoadBPSEQ(filename)
        elif (FileFormat=="RAW"):
            LoadRAW(filename)
        else:
            raise Exception("Unable to determine file type.")




    def AnalyzeFormat(self, filename: str):
        try:
            data = open(filename).readlines()
        except:
            raise Exception("Unable to open input file: " + filename)
        s = None
        for i in data:
            if (len(i) > 0):
                s = i
                break
        FileFormat = None
        if(s[0]==">"):
            FileFormat = "FASTA"
        else:
            iss = s.split(" ")
            if(iss[0].isnumeric and len(iss[1]) == 1 and iss[2].isnumeric):
                FileFormat = "BPSEQ"
        else:
            FileFormat = "RAW"
        
        return format

    def LoadFASTA(self, filename: str):   
        self.names = list()
        self.sequences = list()
        self.mapping = list()

        try:
            data = open(filename).readlines()
        except:
            raise Exception("Unable to open input file: " + filename)

        for i in data:
            s = i.strip()
            if(len(s)==0):
                continue
            if(s[0] == ">"):
                self.names.append(s[1:])
                self.sequences.append("@")
            else:
                if(len(self.sequences)==0):
                    raise Exception("Expected header for FASTA file:", filename)
                for j in range(len(s)):
                    if(s[j] == " "):
                        continue
                    self.sequences[-1]+=s[j]   #recheck
        
        if(len(self.sequences)==0):
            raise Exception("No sequences read.")
        if(len(self.sequences[0])==1):
            raise Exception("Zero-length sequnce read.")
        for i in self.sequences:
            if(len(i)!=len(self.sequences[0])):
                raise Exception("Not all sequences have the same length.")
        
        consensus_found = False
        for i in sequences:
            is_consensus  = True
            for j in i[1:]:
                if(is_consensus==False):
                    break;
                if(j.isalpha()):
                    is_consensus = False
            if(is_consensus):
                if(consensus_found):
                    raise Exception("More than one consensus base-pairing structure found.")
                else:
                    self.mapping = ConvertParensToMapping(FilterParens(i))
                    self.sequences.pop(1)
                    self.names.pop(1)
                    consensus_found = True
                    continue
        
        if(consensus_found==False):
            self.mapping = [UNKNOWN for i in range(len(self.sequences[0]))]

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
