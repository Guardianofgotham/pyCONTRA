from __future__ import annotations
from pyCONTRA.config import *
from pyCONTRA.SStruct import *
from pyCONTRA.ParameterManager import *
import array as arr
from queue import Queue
import sys


class InferenceEngine:
    DATA_LOW_THRESH = 1e-7

    def __init__(self, allow_noncomplementary: bool, num_data_sources: int):
        self.allow_noncomplementary = allow_noncomplementary

        self.char_mapping = [0]*256
        self.is_complementary = []
        for i in range(M+1):
            self.is_complementary.append([0]*(M+1))
        self.cache_initialized = False
        self.parameter_manager = None
        self.num_data_sources = num_data_sources
        self.L = 0
        self.SIZE = 0

        if PROFILE == 0:
            N, SIZE2 = None, None
            A, weights = [], []

        self.s, self.offset = [], []

        self.allow_unpaired_position = []
        self.allow_unpaired, self.allow_paired = [], []
        self.loss_unpaired_position = []
        self.loss_unpaired, self.loss_paired = [], []

        # ## dynamic programming matrices
        self.FCt, self.F5t, self.FMt, self.FM1t = [], [], [], []
        self.FCv, self.F5v, self.FMv, self.FM1v = [], [], [], []
        self.FCi, self.F5i, self.FMi, self.FM1i = [], [], [], []
        self.FCo, self.F5o, self.FMo, self.FM1o = [], [], [], []

        self.FCi_ess, self.F5i_ess, self.FMi_ess, self.FM1i_ess = [], [], [], []
        self.FCo_ess, self.F5o_ess, self.FMo_ess, self.FM1o_ess = [], [], [], []

        self.posterior = []

        self.score_unpaired_position = []
        self.score_paired_position = []
        self.score_unpaired_position_raw = []
        self.score_paired_position_raw = []
        self.score_base_pair = list()
        for i in range(M+1):
            self.score_base_pair.append([])
            for j in range(M+1):
                self.score_base_pair[-1].append((0, 0))

        self.score_terminal_mismatch = list()
        self.score_internal_explicit = list()
        self.score_internal_1x1_nucleotides = list()
        self.score_helix_stacking = list()
        self.score_helix_closing = list()
        self.score_dangle_left = list()
        self.score_dangle_right = list()
        for i in range(M+1):
            self.score_terminal_mismatch.append([])
            self.score_internal_explicit.append([])
            self.score_internal_1x1_nucleotides.append([])
            self.score_helix_stacking.append([])
            self.score_helix_closing.append([])
            self.score_dangle_left.append([])
            self.score_dangle_right.append([])
            for j in range(M+1):
                self.score_terminal_mismatch[-1].append([])
                self.score_internal_explicit[-1].append([])
                self.score_internal_1x1_nucleotides[-1].append(None)
                self.score_helix_stacking[-1].append([])
                self.score_helix_closing[-1].append([])
                self.score_dangle_left[-1].append([])
                self.score_dangle_right[-1].append([])
                for k in range(M+1):
                    self.score_terminal_mismatch[-1][-1].append([])
                    # self.score_terminal_mismatch[-1][-1].append([])
                    self.score_helix_stacking[-1][-1].append([])
                    self.score_dangle_left[-1][-1].append([])
                    self.score_dangle_right[-1][-1].append([])
                    for l in range(M+1):
                        self.score_terminal_mismatch[-1][-1][-1].append((0, 0))
                        # self.score_terminal_mismatch[-1][-1][-1].append((0, 0))
                        self.score_helix_stacking[-1][-1][-1].append((0, 0))

        self.cache_score_hairpin_length = [0]*(D_MAX_HAIRPIN_LENGTH+1)
        self.score_hairpin_length_at_least = [0]*(D_MAX_HAIRPIN_LENGTH+1)
        self.score_bulge_length_at_least = [0]*(D_MAX_BULGE_LENGTH+1)
        self.score_internal_length_at_least = [0]*(D_MAX_BULGE_LENGTH+1)
        self.score_internal_symmetric_length_at_least = [
            0]*(D_MAX_BULGE_LENGTH+1)
        self.score_bulge_0x1_nucleotides = [0]*(D_MAX_BULGE_LENGTH+1)
        self.score_bulge_1x0_nucleotides = [0]*(D_MAX_BULGE_LENGTH+1)
        self.score_internal_asymmetry_at_least = [0]*(D_MAX_BULGE_LENGTH+1)
        self.score_multi_base = None
        self.score_multi_unpaired = None
        self.score_multi_paired = None
        self.score_external_unpaired = None
        self.score_external_paired = None
        self.log_score_evidence = None
        self.cache_score_single = None
        self.cache_score_helix_sums = []

    def FillScores(self):
        raise Exception("Not implemented")

    def FillCounts(self):
        raise Exception("Not implemented")

    def ComputeRowOffset(self, i: int, N: int) -> int:
        if(i < 0 and i > N):
            raise Exception("Index out-of-bonuds")
        return i*(N+N-i-1)//2
        # raise Exception("Not implemented")

    def IsComplementary(self, i: int, j: int) -> bool:
        return self.is_complementary[self.s[i]][self.s[j]]
        # raise Exception("Not implemented")

    def ScoreJunctionA(self, i: int,  j: int):
        raise Exception("Not implemented")

    def ScoreJunctionB(self, i: int, j: int):
        raise Exception("Not implemented")

    def ScoreBasePair(self, i: int,  j: int):
        raise Exception("Not implemented")

    def ScoreHairpin(self, i: int, j: int):
        raise Exception("Not implemented")

    def ScoreHelix(self, i: int,  j: int, m: int):
        raise Exception("Not implemented")

    def ScoreSingleNucleotides(self, i: int,  j: int, p: int, q: int):
        raise Exception("Not implemented")

    def ScoreSingle(self, i: int, j: int, p: int, q: int):
        raise Exception("Not implemented")

    def CountJunctionA(self, i: int,  j: int, value):
        raise Exception("Not implemented")

    def CountJunctionB(self, i: int, j: int, value):
        raise Exception("Not implemented")

    def CountBasePair(self, i: int, j: int,  value):
        raise Exception("Not implemented")

    def CountHairpin(self, i: int, j: int,  value):
        raise Exception("Not implemented")

    def CountHelix(self, i: int,  j: int,  m: int, value):
        raise Exception("Not implemented")

    def CountSingleNucleotides(self, i: int, j: int, p: int,  q: int, value):
        raise Exception("Not implemented")

    def CountSingle(self, i: int,  j: int, p: int, q: int,  value):
        raise Exception("Not implemented")

    def ScorePairedUnpositionEvidenceRaw(self, which_data: int, i: int):
        raise Exception("Not implemented")

    def ScoreUnpairedPositionEvidenceRaw(self, which_data: int, i: int):
        raise Exception("Not implemented")

    def ScoreUnpairedPositionEvidence(self, i: int):
        raise Exception("Not implemented")

    def ScorePairedPositionEvidence(self, i: int):
        raise Exception("Not implemented")

    def ScoreUnpairedEvidence(self, i: int, j: int):
        raise Exception("Not implemented")

    def ScoreBasePairEvidence(self, i: int,  j: int):
        raise Exception("Not implemented")

    def ScoreHelixEvidence(self, i: int, j: int,  m: int):
        raise Exception("Not implemented")

    def CountBasePairEvidence(self, i: int, j: int,  value):
        raise Exception("Not implemented")

    def CountHelixEvidence(self, i: int, j: int, m: int,  value):
        raise Exception("Not implemented")

    def ScoreHairpinEvidence(self, i: int, j: int):
        raise Exception("Not implemented")

    def ScoreSingleNucleotidesEvidence(self, i: int,  j: int,  p: int,  q: int):
        raise Exception("Not implemented")

    def ScoreSingleEvidence(self, i: int, j: int, p: int, q: int):
        raise Exception("Not implemented")

    def CountHairpinEvidence(self, i: int, j: int, value):
        raise Exception("Not implemented")

    def CountSingleNucleotidesEvidence(self, i: int,  j: int,  p: int,  q: int,  value):
        raise Exception("Not implemented")

    def CountSingleEvidence(self, i: int, j: int, p: int,  q: int, value):
        raise Exception("Not implemented")

    def ComputeProfileScore(self, profile_score, pos, dimensions, table):
        raise Exception("Not implemented")

    def ConvertProfileCount(self, profile_score,  pos,  dimensions: int, table):
        raise Exception("Not implemented")

    # register params with the parameter manager
    def RegisterParameters(self, parameter_manager: ParameterManager):
        buffer = ""
        buffer2 = ""
        self.cache_initialized = False
        self.parameter_manager = parameter_manager
        parameter_manager.ClearParameters()
        parameter_manager.AddParameterGroup("all_params")
        for i in range(M+1):
            for j in range(M+1):
                if (i == M or j == M):
                    self.score_base_pair[i][j] = (0, 0)
                else:
                    buffer = f"base_pair_{alphabet[i]}{alphabet[j]}"
                    buffer2 = f"base_pair_{alphabet[j]}{alphabet[i]}"
                    if ((buffer < buffer2)):
                        parameter_manager.AddParameterMapping(
                            buffer, self.score_base_pair[i][j])
                    else:
                        parameter_manager.AddParameterMapping(
                            buffer2, self.score_base_pair[i][j])

        # Complete if-endif statements
        for i1 in range(0, M+1):
            for j1 in range(0, M+1):
                for i2 in range(0, M+1):
                    for j2 in range(0, M+1):
                        if (i1 == M or j1 == M or i2 == M or j2 == M):
                            self.score_terminal_mismatch[i1][j1][i2][j2] = (0, 0)
                        else:
                            buffer = f"terminal_mismatch_{alphabet[i1]}{alphabet[j1]}{alphabet[i2]}{alphabet[j2]}"
                            parameter_manager.AddParameterMapping(buffer, self.score_terminal_mismatch[i1][j1][i2][j2])
        
        for i in range(0, D_MAX_HAIRPIN_LENGTH+1):
            buffer = f"hairpin_length_at_least_{i}"
            parameter_manager.AddParameterMapping(
                buffer, self.score_hairpin_length_at_least[i])

        for i in range(0, D_MAX_INTERNAL_EXPLICIT_LENGTH+1):
            for j in range(0, D_MAX_INTERNAL_EXPLICIT_LENGTH+1):
                if i == 0 or j == 0:
                    self.score_internal_explicit[i][j] = (0, 0)
                else:
                    buffer = f"internal_explicit_{min(i, j)}_{max(i, j)}"
                    parameter_manager.AddParameterMapping(
                        buffer, self.score_internal_explicit[i][j])

        for i in range(0, D_MAX_BULGE_LENGTH+1):
            if i == 0:
                self.score_bulge_length_at_least[i] = (0, 0)
            else:
                buffer = f"bulge_length_at_least_{i}"
                parameter_manager.AddParameterMapping(
                    buffer, self.score_bulge_length_at_least[i])
        
        for i in range(0, D_MAX_INTERNAL_LENGTH+1):
            if (i < 2):
                self.score_internal_length_at_least[i] = (0, 0)
            else:
                buffer = f"internal_length_at_least_{i}"
                parameter_manager.AddParameterMapping(
                    buffer, self.score_internal_length_at_least[i])
        
        for i in range(0,  D_MAX_INTERNAL_SYMMETRIC_LENGTH+1):
            if (i == 0):
                self.score_internal_symmetric_length_at_least[i] = (0, 0)
            else:
                buffer = f"internal_symmetric_length_at_least_{i}"
                parameter_manager.AddParameterMapping(
                    buffer, self.score_internal_symmetric_length_at_least[i])
        
        for i in range(0, D_MAX_INTERNAL_ASYMMETRY+1):
            if (i == 0):
                self.score_internal_asymmetry_at_least[i] = (0, 0)
            else:
                buffer = f"internal_asymmetry_at_least_{i}"
                parameter_manager.AddParameterMapping(
                    buffer, self.score_internal_asymmetry_at_least[i])
        
        for i1 in range(0,  M+1):
            if (i1 == M):
                self.score_bulge_0x1_nucleotides[i1] = (0, 0)
                self.score_bulge_1x0_nucleotides[i1] = (0, 0)
            else:
                buffer = f"bulge_0x1_nucleotides_{alphabet[i1]}"
                parameter_manager.AddParameterMapping(
                    buffer, self.score_bulge_0x1_nucleotides[i1])
                parameter_manager.AddParameterMapping(
                    buffer, self.score_bulge_1x0_nucleotides[i1])
        
        for i1 in range(0, M+1):
            for i2 in range(0, M+1):
                if (i1 == M or i2 == M):
                    self.score_internal_1x1_nucleotides[i1][i2] = (0, 0)
                else:
                    buffer = f"internal_1x1_nucleotides_{alphabet[i1]}{alphabet[i2]}"
                    buffer2 = f"internal_1x1_nucleotides_{alphabet[i2]}{alphabet[i1]}"
                if (buffer < buffer2):
                    parameter_manager.AddParameterMapping(
                        buffer, self.score_internal_1x1_nucleotides[i1][i2])
                else:
                    parameter_manager.AddParameterMapping(
                        buffer2, self.score_internal_1x1_nucleotides[i1][i2])
        
        for i1 in range(0, M+1):
            for j1 in range(0, M+1):
                for i2 in range(0, M+1):
                    for j2 in range(0, M+1):
                        if (i1 == M or j1 == M or i2 == M or j2 == M):
                            self.score_helix_stacking[i1][j1][i2][j2] = (0, 0)
                        else:
                            buffer = f"helix_stacking_{alphabet[i1]}{alphabet[j1]}{alphabet[i2]}{alphabet[j2]}"
                            buffer2 = f"helix_stacking_{alphabet[j2]}{alphabet[i2]}{alphabet[j1]}{alphabet[i1]}"
                        if (buffer < buffer2):
                            parameter_manager.AddParameterMapping(
                                buffer, self.score_helix_stacking[i1][j1][i2][j2])
                        else:
                            parameter_manager.AddParameterMapping(
                                buffer2, self.score_helix_stacking[i1][j1][i2][j2])
        
        print(f"PMANAGER: {parameter_manager.GetNumLogicalParameters()}")
        for i in range(0, M+1):
            for j in range(0, M+1):
                if (i == M or j == M):
                    self.score_helix_closing[i][j] = (0, 0)
                else:
                    buffer = f"helix_closing_{alphabet[i]}{alphabet[j]}"
                    parameter_manager.AddParameterMapping(
                        buffer, self.score_helix_closing[i][j])
        
        parameter_manager.AddParameterMapping(
            "multi_base", self.score_multi_base)
        parameter_manager.AddParameterMapping(
            "multi_unpaired", self.score_multi_unpaired)
        parameter_manager.AddParameterMapping(
            "multi_paired", self.score_multi_paired)
        
        for i1 in range(0, M+1):
            for j1 in range(0, M+1):
                for i2 in range(0, M+1):
                    if (i1 == M or j1 == M or i2 == M):
                        self.score_dangle_left[i1][j1][i2] = (0, 0)
                    else:
                        buffer = f"dangle_left_{alphabet[i1]}{alphabet[j1]}{alphabet[i2]}"
                        parameter_manager.AddParameterMapping(
                            buffer, self.score_dangle_left[i1][j1][i2])

        
        for i1 in range(0, M+1):
            for j1 in range(0, M+1):
                for j2 in range(0, M+1):
                    if (i1 == M or j1 == M or j2 == M):
                        self.score_dangle_right[i1][j1][j2] = (0, 0)
                    else:
                        buffer = f"dangle_right_{alphabet[i1]}{alphabet[j1]}{ alphabet[j2]}"
                        parameter_manager.AddParameterMapping(
                            buffer, self.score_dangle_right[i1][j1][j2])
        parameter_manager.AddParameterMapping(
            "external_unpaired", self.score_external_unpaired)
        parameter_manager.AddParameterMapping(
            "external_paired", self.score_external_paired)
        parameter_manager.AddParameterGroup("evidence_cpds");
        for num_data_sources_current in range(0, self.num_data_sources):
            for i0 in range(0, 2):
                for i1 in range(0, M):
                    for i2 in range(0, 2):
                        if (i0 == 0):
                            buffer = f"log_score_evidence{num_data_sources_current}_k_{alphabet[i1]}{i2}"
                        else:
                            buffer = f"log_score_evidence{num_data_sources_current}_theta_{alphabet[i1]}{i2}"
                        parameter_manager.AddParameterMapping(
                            buffer, self.log_score_evidence[num_data_sources_current][i0][i1][i2])
        # raise Exception("Not implemented")

    # load sequence
    def LoadSequence(self, sstruct):
        self.cache_initialized = False
        self.L = sstruct.GetLength()
        self.SIZE = (self.L+1)*(self.L+2) // 2
        print(self.SIZE)
        self.s = [0]*(self.L+1)
        self.offset = [0]*(self.L+1)
        self.allow_unpaired_position = [0]*(self.L+1)
        self.allow_unpaired = [0]*(self.SIZE)
        self.allow_paired = [0]*(self.SIZE)
        self.loss_unpaired_position = [0]*(self.L+1)
        self.loss_unpaired = [0]*(self.SIZE)
        self.loss_paired = [0]*(self.SIZE)

        self.cache_score_helix_sums.clear()
        self.cache_score_helix_sums = [(0, 0)]*((2*self.L+1)*self.L)
        sequence = sstruct.GetSequences()[0]
        self.s[0] = len(alphabet)
        for i in range(1, self.L+1):
            self.s[i] = self.char_mapping[ord(sequence[i])]

        for i in range(0, self.L):
            self.offset[i] = self.ComputeRowOffset(i, self.L+1)
            self.allow_unpaired_position[i] = 1
            self.loss_unpaired_position[i] = 0

        for i in range(0, self.SIZE):
            self.allow_unpaired[i] = 1
            self.allow_paired[i] = 1
            self.loss_unpaired[i] = 0
            self.loss_paired[i] = 0

        for i in range(0, self.L+1):
            self.allow_paired[self.offset[0]+i] = 0
            self.allow_paired[self.offset[i]+i] = 0
        if not self.allow_noncomplementary:
            for i in range(1, self.L+1):
                for j in range(1, self.L+1):
                    if not self.IsComplementary(i, j):
                        self.allow_paired[self.offset[i]+j] = 0
        for i in range(0, self.num_data_sources):
            self.score_unpaired_position_raw[i].clear()
            self.score_unpaired_position[i].clear()
            self.score_paired_position_raw[i].clear()
            self.score_paired_position[i].clear()
            self.score_unpaired_position_raw[i] = sstruct.GetUnpairedPotentials(
                i)
            self.score_paired_position_raw[i] = sstruct.GetPairedPotentials(i)

        # raise Exception("Not implemented")

    # load parameter values

    def LoadValues(self, values: list):
        print(len(values), self.parameter_manager.GetNumLogicalParameters())
        if len(values) != self.parameter_manager.GetNumLogicalParameters():
            raise Exception("Parameter Size MisMatch")
        self.cache_initialized = False
        for i in range(0, len(values)):
            physical_parameters = self.parameter_manager.GetPhysicalParameters(i)
            for j in range(0, len(physical_parameters)):
                physical_parameters[j] = (values[i],0)
        # raise Exception("Not implemented")

    def UseLoss(self, true_mapping,  example_loss):
        raise Exception("Not implemented")

    # use Constraints
    def UseConstraints(self, true_mapping):
        # print(len(true_mapping), self.L+1)
        if(not(len(true_mapping) == self.L+1)):
            raise Exception("Supplied mapping of incorrect length!")
        cache_initialized = False

        for i in range(1, self.L+1):
            self.allow_unpaired_position[i] = (
                true_mapping[i] == SStruct.UNKNOWN or true_mapping[i] == SStruct.UNPAIRED)

        for i in range(self.L+1):
            self.allow_unpaired[self.offset[i]+i] = 1
            self.allow_paired[self.offset[i]+i] = 0
            for j in range(1, self.L+1):
                self.allow_unpaired[self.offset[i]+j] = self.allow_unpaired[self.offset[i] +
                                                                            j-1] and self.allow_unpaired_position[j]
                self.allow_paired[self.offset[i]+j] = (i > 0 and (true_mapping[i] == SStruct.UNKNOWN or true_mapping[i] == j) and (
                    true_mapping[j] == SStruct.UNKNOWN or true_mapping[j] == i) and (self.allow_noncomplementary or self.IsComplementary(i, j)))

    # Viterbi inference

    def ComputeViterbi(self):
        raise Exception("Not implemented")

    def GetViterbiScore(self):
        raise Exception("Not implemented")

    def PredictPairingsViterbi(self):
        raise Exception("Not implemented")

    def ComputeViterbiFeatureCounts(self):
        raise Exception("Not implemented")

    # MEA inference

    def ComputeInside(self):
        self.InitializeCache()
        self.F5i=[NEG_INF]*(self.L+1)
        self.FCi=[NEG_INF]*(self.SIZE);
        self.FMi=[NEG_INF]*(self.SIZE);
        self.FM1i=[NEG_INF]*(self.SIZE);
       

        for i in range(self.L,-1,-1):#(int i = L; i >= 0; i--)
            for j in range(i,self.L+1):#(int j = i; j <= L; j++
            
                FM2i = NEG_INF
                if (i+2 <= j):
                    p1 = self.FM1i[self.offset[i]+i+1]
                    p2 = self.FMi[self.offset[i+1]+j]
                    for k in range(i+1,j):# (register int k = i+1; k < j; k++)
    
                        Fast_LogPlusEquals(FM2i, (p1) + (p2))
                        p1+=1
                        p2 += self.L-k;
                    
                if (0 < i and j < self.L and self.allow_paired[self.offset[i]+j+1]):
                
                    sum_i = (NEG_INF);
                    
                    # compute ScoreHairpin(i,j)
                    
                    if (self.allow_unpaired[self.offset[i]+j] and j-i >= C_MIN_HAIRPIN_LENGTH):
                        Fast_LogPlusEquals(sum_i, ScoreHairpin(i,j));

                    score_helix=0
                    #score_helix = (i+2 <= j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    if(i+2<=j):
                        score_helix=ScoreBasePair(i+1,j)+ScoreHelixStacking(i,j+1)
                    else:
                        score_helix=0
                    score_other = ScoreJunctionB(i,j);
                    
                    for p in range(i,(min(i+C_MAX_SINGLE_LENGTH,j)+1)):  #(int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    
                        if (p > i and not(self.allow_unpaired_position[p])):
                            break;

                        q_min = max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        FCptr = self.FCi[self.offset[p+1]-1]
                        for q in range(q_min,j-1,-1):#(int q = j; q >= q_min; q--):
                            if(q < j and not(self.allow_unpaired_position[q+1])):
                                break;
                            if(not(self.allow_paired[self.offset[p+1]+q])):
                                continue;
                            
                            if(p == i and q == j):
                                score=score_helix + FCptr[q]
                            # score = (p == i && q == j) ?
                            else:
                                score=score_other + self.cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q)

                            #     (score_helix + FCptr[q]) :
                            #     (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                            #     ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                                Fast_LogPlusEquals(sum_i, score)
                        
                            Fast_LogPlusEquals(sum_i, FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                            self.FCi[self.offset[i]+j] = sum_i;
                            if (0 < i and i+2 <= j and j < self.L):
                                sum_i = NEG_INF
                                
                                #compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                                
                                if (self.allow_paired[self.offset[i+1]+j]):
                                    Fast_LogPlusEquals(sum_i, self.FCi[self.offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j));
                                
                                #compute FM1[i+1,j] + b
                                
                                if (self.allow_unpaired_position[i+1]):
                                    Fast_LogPlusEquals(sum_i, self.FM1i[self.offset[i+1]+j] + ScoreMultiUnpaired(i+1));
                                
                                self.FM1i[self.offset[i]+j] = sum_i;
                            if ((0 < i) and (i+2 <= j) and (j < self.L)):
                                
                                sum_i = NEG_INF
                                
                                #compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                                
                                Fast_LogPlusEquals(sum_i, FM2i);
                                
                                #compute FM[i,j-1] + b
                                
                                if (self.allow_unpaired_position[j]):
                                    Fast_LogPlusEquals(sum_i, FMi[self.offset[i]+j-1] + ScoreMultiUnpaired(j));
                                
                                #compute FM1[i,j]
                                
                                Fast_LogPlusEquals(sum_i, self.FM1i[self.offset[i]+j]);
                                
                                self.FMi[self.offset[i]+j] = sum_i
            self.F5i[0] = 0
            for j in range(1,self.L+1):#(int j = 1; j <= L; j++)
                
                sum_i = NEG_INF;
                
                #compute F5[j-1] + ScoreExternalUnpaired()
                
                if (self.allow_unpaired_position[j]):
                    Fast_LogPlusEquals(sum_i, self.F5i[j-1] + ScoreExternalUnpaired(j));
                
                #compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
                
                for k in range(0,j):
                    if (self.allow_paired[self.offset[k+1]+j]):
                        Fast_LogPlusEquals(sum_i, self.F5i[k] + self.FCi[self.offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
                
                self.F5i[j] = sum_i;

                
                # compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
        # raise Exception("Not implemented")

    def ComputeLogPartitionCoefficient(self): 
        raise Exception("Not implemented")

    def ComputeOutside(self):
        self.F5o.clear(); self.F5o=[NEG_INF]*(self.L+1)
        self.FCo.clear(); self.FCo=[NEG_INF]*self.SIZE
        self.FMo.clear(); self.FMo=[NEG_INF]*self.SIZE
        self.FM1o.clear(); self.FM1o=[NEG_INF]*self.SIZE
        self.F5o[self.L] = 0
        for j in range(self.L,0,-1):
            if self.allow_unpaired_position[j]:
                Fast_LogPlusEquals(self.F5o[j-1], self.F5o[j] + ScoreExternalUnpaired(j))
        
            for k in range(0, j):
                if self.allow_paired[self.offset[k+1]+j]:
                    temp = self.F5o[j] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k)
                    Fast_LogPlusEquals(self.F5o[k], temp + self.FCi[self.offset[k+1]+j-1])
                    Fast_LogPlusEquals(self.FCo[self.offset[k+1]+j-1], temp + self.F5i[k])
        for i in range(0, self.L+1):
            for j in range(self.L, i-1, -1):
                FM2o = NEG_INF
                if(i>0 and i+2<=j and j<self.L):
                    Fast_LogPlusEquals(FM2o, self.FMo[self.offset[i]+j])
                
                    # // compute FM[i,j-1] + b
                    
                    if self.allow_unpaired_position[j]:
                        Fast_LogPlusEquals(self.FMo[self.offset[i]+j-1], self.FMo[self.offset[i]+j] + ScoreMultiUnpaired(j))
                    
                    # // compute FM1[i,j]
                    
                    Fast_LogPlusEquals(self.FM1o[self.offset[i]+j], self.FMo[self.offset[i]+j])

                    
                if (0 < i and i+2 <= j and j < self.L):
                # // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                    if (self.allow_paired[self.offset[i+1]+j]):
                        Fast_LogPlusEquals(self.FCo[self.offset[i+1]+j-1], self.FM1o[self.offset[i]+j] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j))
                    
                    # // compute FM1[i+1,j] + b
                    
                    if (self.allow_unpaired_position[i+1]):
                        Fast_LogPlusEquals(self.FM1o[self.offset[i+1]+j], self.FM1o[self.offset[i]+j] + ScoreMultiUnpaired(i+1))
                

                if (0 < i and j < self.L and self.allow_paired[self.offset[i]+j+1]):
                    score_helix = self.FCo[offset[i]+j] + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1)  if (i+2 <= j) else 0
                    score_other = self.FCo[self.offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for p in range(i, min(i+C_MAX_SINGLE_LENGTH,j)+1):
                        if (p > i and not allow_unpaired_position[p]):
                            break
                        q_min = max(p+2,p-i+j-C_MAX_SINGLE_LENGTH)
                        FCptr = self.FCo[self.offset[p+1]-1]
                        for q in range(j, q_min-1, -1):
                            if (q < j and not self.allow_unpaired_position[q+1]):
                                break
                            if (not self.allow_paired[self.offset[p+1]+q]):
                                continue;
                            Fast_LogPlusEquals(FCptr[q], 
                                            score_helix if (p == i and q == j) else score_other + self.cache_score_single[p-i][j-q][0] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q))
                    Fast_LogPlusEquals(FM2o, self.FCo[self.offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                

                
                if i+2 <= j:
                    p1i = self.FM1i[self.offset[i]+i+1]
                    p2i = self.FMi[self.offset[i+1]+j]
                    p1o = self.FM1o[self.offset[i]+i+1]
                    p2o = self.FMo[self.offset[i+1]+j]
                    for k in range(i+1, j):
                        Fast_LogPlusEquals(p1o, FM2o + p2i)
                        Fast_LogPlusEquals(p2o, FM2o + p1i)
                        p1i+=1
                        p1o+=1
                        p2i += L-k
                        p2o += L-k
        
        raise Exception("Not implemented")

    def ComputeFeatureCountExpectations(self):
        raise Exception("Not implemented")

    def ComputePosterior(self):
        raise Exception("Not implemented")

    def PredictPairingsPosterior(self,   gamma: int):
        if(not(gamma > 0)):
            print("Non-negative gamma expected.")
        unpaired_posterior = []
        score = []
        traceback = []

        for i in range(1, self.L+1):
            unpaired_posterior.append(1)
            for j in range(i):
                unpaired_posterior[i] -= self.posterior[self.offset[j]+i]
            for j in range(i+1, self.L+1):
                unpaired_posterior[i] -= self.posterior[self.offset[i]+j]

        for i in range(1, self.L+1):
            unpaired_posterior[i] /= 2 * gamma
        score = [-1]*self.SIZE
        traceback = [-1]*self.SIZE
        print(self.SIZE)
        # DP

        for i in range(self.L, -1, -1):
            for j in range(i, self.L+1):
                this_score = score[self.offset[i]+j]
                this_traceback = traceback[self.offset[i]+j]
                if(i == j):
                    UPDATE_MAX(this_score, this_traceback, 0, 0)
                else:
                    if(self.allow_unpaired_position[i+1]):
                        UPDATE_MAX(
                            this_score, this_traceback, unpaired_posterior[i+1] + score[self.offset[i+1]+j], 1)
                    if (self.allow_unpaired_position[j]):
                        UPDATE_MAX(
                            this_score, this_traceback, unpaired_posterior[j] + score[self.offset[i]+j-1], 2)
                    if(i+2 <= j):
                        if(self.allow_paired[self.offset[i+1]+j]):
                            UPDATE_MAX(
                                this_score, this_traceback, self.posterior[self.offset[i+1]+j] + score[self.offset[i+1]+j-1], 3)

                        p1 = score[self.offset[i]+i+1]
                        p2 = score[self.offset[i+1]+j]
                        for k in range(i+1, j):
                            UPDATE_MAX(this_score, this_traceback,
                                       (p1) + (p2), k+4)
                            p1 += 1
                            p2 += self.L-k

        solution = [SStruct.UNPAIRED]*(self.L+1)
        solution[0] = SStruct.UNKNOWN
        traceback_queue = Queue(0)
        traceback_queue.put((0, self.L))

        while(traceback_queue.qsize() != 0):
            t = traceback_queue.get()
            i = t[0]
            j = t[1]

            if(traceback[self.offset[i]+j] == -1):
                print("should not get here ")

            elif(traceback[self.offset[i]+j] == 0):
                raise Exception("Not implemented")
            elif(traceback[self.offset[i]+j] == 1):
                traceback_queue.put((i+1, j))

            elif(traceback[self.offset[i]+j] == 2):
                traceback_queue.put((i, j-1))
            elif(traceback[self.offset[i]+j] == 3):
                solution[i+1] = j
                solution[j] = i+1
                traceback_queue.put((i+1, j-1))

            else:
                k = traceback[self.offset[i]+j]-4
                traceback_queue.put((i, k))
                traceback_queue.put((k, j))

        return solution

    def PredictPairingsPosteriorCentroid(self, gamma):
        raise Exception("Not implemented")

    def GetPosterior(self,   posterior_cutoff):
        raise Exception("Not implemented")

    # EM inference

    def ComputeInsideESS(self):
        raise Exception("Not implemented")

    def ComputeOutsideESS(self):
        raise Exception("Not implemented")

    def ComputeFeatureCountExpectationsESS(self):
        raise Exception("Not implemented")

    def ComputePosteriorESS(self):
        raise Exception("Not implemented")

    def ComputeLogPartitionCoefficientESS(self):
        raise Exception("Not implemented")

    def ComputeESS(self):
        raise Exception("Not implemented")

    def ComputeGammaMLESum(self, ev_cpd_id,  ignorePairing,  usePosterior,  which_data):
        raise Exception("Not implemented")

    def ComputeGammaMLESS(self, ev_cpd_id,  ignoreZeros,  useMLE,  which_data):
        raise Exception("Not implemented")

    def ComputeGammaMLEESS(self, ev_cpd_id,  ignoreZeros,  useMLE,  which_data):
        raise Exception("Not implemented")

    def GetNumExamplesSeqPairing(self, ev_cpd_id,  ignoreZeros,  which_data):
        raise Exception("Not implemented")

    def GetNumExamplesSeq(self, ev_cpd_id,  ignoreZeros,  which_data):
        raise Exception("Not implemented")

    def AreZerosInSeqPairing(self,  id_base,  id_pairings,  which_base):
        raise Exception("Not implemented")

    def AreZerosInSeq(self,  id_base,  which_base):
        raise Exception("Not implemented")

    def LogGammaProb(self,  data,  which_data,  isUnpaired,  seq):
        raise Exception("Not implemented")

    def UpdateEvidenceStructuresK(self,  which_data):
        self.score_unpaired_position[which_data].clear()
        self.score_paired_position[which_data].clear()

        for i in range(0, self.L):
            score_unpaired = 0  if self.score_unpaired_position_raw[which_data][i] == SStruct.UNKNOWN_POTENTIAL else self.LogGammaProb(self.score_unpaired_position_raw[which_data][i],which_data,1,self.s[i+1])
            self.score_unpaired_position[which_data].append(score_unpaired)

            score_paired = 0 if self.score_paired_position_raw[which_data][i] == SStruct.UNKNOWN_POTENTIAL else self.LogGammaProb(score_paired_position_raw[which_data][i],which_data,0,s[i+1])
            self.score_paired_position[which_data].append(score_paired)
        

    def UpdateEvidenceStructures(self):
        for i in range(0, self.num_data_sources):
            UpdateEvidenceStructuresK(i)
