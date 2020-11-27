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
        for i in self.sequences:
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
        self.names = list()
        self.sequences = list()
        self.mapping = list()

        self.names.append(filename)
        self.sequences.append("@")
        self.mapping.append(UNKNOWN)

        try:
            data = open(filename).readlines()
        except:
            raise Exception("Unable to open input file: " + filename)
        
        row = 0
        for i in data:
            tokens = i.split()
            if(!tokens[0].isnumeric):
                raise Exception("Could not read row number:", filename)
            if(tokens[0]<=0):
                raise Exception("Row numbers must be positive:", filename)
            if(tokens[0] != (row+1)):
                raise Exception(
                    "Rows of BPSEQ file must occur in increasing order:", filename.)
            row = tokens[0]

            try:
                tokens[1]
            except:
                raise Exception ("Expected sequence letter after row number:", filename)
            if(len(tokens[1])!=1):
                raise Exception(
                    "Expected sequence letter after row number:", filename)
            ch = tokens[1][0]
            try:
                tokens[2]
            except:
                raise Exception(
                    "Expected mapping letter after sequence letter:", filename)
            if(!tokens[2].isnumeric):
                raise Exception("Could not read matching row number:", filename)
            if(tokens[2] <= -1):
                raise Exception(
                    "Matching row numbers must be greater than or equal to -1:", filename)
            self.sequences[-1].append(ch)
            self.mapping.append(tokens[2])


    def LoadBPP2SEQ(self, filename: str):   
        pass
    def LoadBPP2TSEQ(self, filename: str):   
        pass
    def FilterSequence(self, sequence: str):   
        if(sequence[0] != "@"):
            raise Exception("Improperly formatted sequence.")
        for i in range(1, len(sequence)):
            if(sequence[i] == "-"):
                sequence[i] = "."
                break
    def FilterParens(self, sequence: str):   
        if(sequence[0] != "@"):
            raise Exception("Improperly formatted sequence.")
        for i in range(1,len(sequence)):
            if(sequence[i]=="-"):
                sequence[i] = "."
                break
            elif(sequence[i] == "?" or sequence[i] == "." or sequence[i] == "(" or sequence[i] == ")"):
                break
            else:
                raise Exception("Unexpected character {} in parenthesized structure.", sequence[i])
        
        return sequence;

    def ConvertParensToMapping(self, parens: str):   
        mapping = [UNKNOWN for i in range(len(parens))]
        stack = list()

        assert parens[0] == "@", "Invalid parenthesized string."
        for i in range(1,len(parens)):
            if(parens[i] == "?"):
                break
            elif(parens[i] == "."):
                mapping[i] = UNPAIRED
                break
            elif(parens[i] == "("):
                stack.append(i)
                break
            elif(parens[i] == ")"):
                if(len(stack)==0):
                    raise Exception("Parenthesis mismatch.")
                mapping[i] = stack[-1]
                mapping[stack[-1]] = i
                stack.pop(-1)
                break
            else:
                raise Exception(
                    "Unexpected character {} in parenthesized structure.".format(parens[i]))
            if(len(stack)!=0):
                raise Exception("Parenthesis mismatch.")

            return mapping

        
                
    def ConvertMappingToParens(self, mapping: list):   
        assert !ContainsPseudoknots(), "Should not attempt to convert a mapping with pseudoknots."
        parens = "@"

        for i in range(1,len(mapping)):
            if (mapping[i] == UNKNOWN):
                parens += "?"
            else if (mapping[i] == UNPAIRED):
                parens += "."
            else if (mapping[i] > i):
                parens += "("
            else if (mapping[i] < i):
                parens += ")"
            else
                raise Exception("Invalid structure.")

    def ValidateMapping(self, mapping: list):   
        if(len(mapping)==0 or mapping[0] != UNKNOWN):
            raise Exception("Invalid mapping.")
        for i in range(1,len(mapping)):
            if(mapping[i] == UNPAIRED or mapping[i] == UNKNOWN):
                continue
            if(mapping[i] < 1 or mapping[i] >= len(mapping)):
                raise Exception("Position {} of sequence maps to invalid psotion".format(i))
            if(mapping[mapping[i]]!=i):
                raise Exception("Positions {} and {} of sequence do not map to each other.".format(i,mapping[i]))
            if(mapping[i] == i):
                raise Exception("Position {} of sequence maps to itself.".format(i))

    def ContainsPseudoknots(self):   
        stack = list()
        for i in range(1,len(self.mapping)):
            if(self.mapping[i] == UNPAIRED or self.mapping == UNKNOWN):
                continue
            if(self.mapping[i] > i):
                stack.append(i)
            elif(self.mapping < i):
                if(stack[-1] == self.mapping[i]):
                    stack.pop(-1)
                else:
                    return True
            else:
                raise Exception("Invalid structure: positions may not map to themselves.")
        
        if(len(stack)!=0):
            raise Exception(
                "Invalid structure: bad pairings found.")
        return False

    def RemoveNoncomplementaryPairings(self, seq=0):   
        if(seq < 0 or seq >= len(self.sequences)):
            raise Exception("Refernece to invalid sequence.")
        assert len(self.sequences[seq]) == len(self.mapping), "Inconsistent lengths."
        for i in range(1, len(self.mapping)):
            if(mapping[i] > i and !IsComplementary(self.sequences[seq][i], self.sequences[seq][mapping[i]])):
                self.mapping[self.mapping[i]] = UNPAIRED
                self.mapping[i] = UNPAIRED

    def WriteBPSEQ(self, outfile, seq=0):   
        if(seq<0 or seq>= len(self.sequences)):
            raise Exception("Refernece to invalid sequence.")
        assert len(self.sequences[seq]) == len(self.mapping), "Inconsistent lengths."
        for i in range(1,len(self.mapping)):
            outfile.write(str(i)+" "+self.sequences[seq][i] +" "+ self.mapping[i])

    def WriteParens(self, outfile):   
        if(ContainsPseudoknots()):
            raise Exception("Cannot write strucutre containing pseudoknot using parenthesized format.")

        for i in range(len(self.sequences)):
            outfile.write(">"+self.names[i])
            outfile.write(self.sequences[i][1:])
        
        outfile.write(">structure")
        outfile.write(ConvertMappingToParens(self.mapping)[1:])
    def ComputePercentIdentity(self):   
        pid = 0.0
        for i in range(len(self.sequences)):
            for j in range(i+1,len(self.sequences)):
                identities = 0
                len1=0
                len2 = 0
                s = self.sequences[i]
                t = self.sequences[j]

                for k in range(len(s)):
                    if(s[k].isalpha()):
                        len1+=1
                    if(t[k].isalpha()):
                        len2 += 1
                    if(s[k].isalpha() and s[k].upper() == t[k].upper()):
                        identities+=1
                
                den = min(len1,len2)
                if(den==0):
                    pairwise_pid = 0.0
                else:
                    pairwise_pid = double(identities) / den

                    

    def ComputePositionBasedSequenceWeights(self, ):   
        pass
    def SetMapping(self, mapping):   
        self.mapping = mapping
        ValidateMapping(mapping)
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