from __future__ import annotations
from pyCONTRA.config import *
from pyCONTRA.SStruct import *
from queue import Queue


class InferenceEngine:
    DATA_LOW_THRESH = 1e-7

    def __init__(self):
        self.allow_noncomplementary = None
        self.char_mapping = "0"*256
        self.is_complementary = []
        for i in range(M+1):
            self.is_complementary.append([0]*(M+1))
        self.cache_initialized = None
        self.parameter_manager = None
        self.num_data_sources = None
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
        FCt, F5t, FMt, FM1t = [], [], [], []
        FCv, F5v, FMv, FM1v = [], [], [], []
        FCi, F5i, FMi, FM1i = [], [], [], []
        FCo, F5o, FMo, FM1o = [], [], [], []

        FCi_ess, F5i_ess, FMi_ess, FM1i_ess = [], [], [], []
        FCo_ess, F5o_ess, FMo_ess, FM1o_ess = [], [], [], []

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

        score_terminal_mismatch = list()
        for i in range(M+1):
            score_terminal_mismatch.append([])
            for j in range(M+1):
                score_terminal_mismatch[-1].append([])
                for k in range(M+1):
                    score_terminal_mismatch[-1][-1].append([])
                    for l in range(M+1):
                        score_terminal_mismatch[-1][-1][-1].append((0, 0))

        score_hairpin_length_at_least = [0]*(D_MAX_HAIRPIN_LENGTH+1)
        cache_score_hairpin_length = [0]*(D_MAX_HAIRPIN_LENGTH+1)
        score_internal_explicit = None
        score_bulge_length_at_least = None
        score_internal_length_at_least = None
        score_internal_symmetric_length_at_least = None
        score_internal_asymmetry_at_least = None
        score_bulge_0x1_nucleotides = None
        score_bulge_1x0_nucleotides = None
        score_internal_1x1_nucleotides = None
        score_helix_stacking = None
        score_helix_closing = None
        score_multi_base = None
        score_multi_unpaired = None
        score_multi_paired = None
        score_dangle_left = None
        score_dangle_right = None
        score_external_unpaired = None
        score_external_paired = None
        log_score_evidence = None
        cache_score_single = None
        cache_score_helix_sums = None
        pass

    def FillScores(self):
        pass

    def FillCounts(self):
        pass

    def ComputeRowOffset(self, i: int, N: int) -> int:
        pass

    def IsComplementary(self, i: int, j: int) -> bool:
        pass

    def ScoreJunctionA(self, i: int,  j: int):
        pass

    def ScoreJunctionB(self, i: int, j: int):
        pass

    def ScoreBasePair(self, i: int,  j: int):
        pass

    def ScoreHairpin(self, i: int, j: int):
        pass

    def ScoreHelix(self, i: int,  j: int, m: int):
        pass

    def ScoreSingleNucleotides(self, i: int,  j: int, p: int, q: int):
        pass

    def ScoreSingle(self, i: int, j: int, p: int, q: int):
        pass

    def CountJunctionA(self, i: int,  j: int, value):
        pass

    def CountJunctionB(self, i: int, j: int, value):
        pass

    def CountBasePair(self, i: int, j: int,  value):
        pass

    def CountHairpin(self, i: int, j: int,  value):
        pass

    def CountHelix(self, i: int,  j: int,  m: int, value):
        pass

    def CountSingleNucleotides(self, i: int, j: int, p: int,  q: int, value):
        pass

    def CountSingle(self, i: int,  j: int, p: int, q: int,  value):
        pass

    def ScorePairedUnpositionEvidenceRaw(self, which_data: int, i: int):
        pass

    def ScoreUnpairedPositionEvidenceRaw(self, which_data: int, i: int):
        pass

    def ScoreUnpairedPositionEvidence(self, i: int):
        pass

    def ScorePairedPositionEvidence(self, i: int):
        pass

    def ScoreUnpairedEvidence(self, i: int, j: int):
        pass

    def ScoreBasePairEvidence(self, i: int,  j: int):
        pass

    def ScoreHelixEvidence(self, i: int, j: int,  m: int):
        pass

    def CountBasePairEvidence(self, i: int, j: int,  value):
        pass

    def CountHelixEvidence(self, i: int, j: int, m: int,  value):
        pass

    def ScoreHairpinEvidence(self, i: int, j: int):
        pass

    def ScoreSingleNucleotidesEvidence(self, i: int,  j: int,  p: int,  q: int):
        pass

    def ScoreSingleEvidence(self, i: int, j: int, p: int, q: int):
        pass

    def CountHairpinEvidence(self, i: int, j: int, value):
        pass

    def CountSingleNucleotidesEvidence(self, i: int,  j: int,  p: int,  q: int,  value):
        pass

    def CountSingleEvidence(self, i: int, j: int, p: int,  q: int, value):
        pass

    def ComputeProfileScore(self, profile_score, pos, dimensions, table):
        pass

    def ConvertProfileCount(self, profile_score,  pos,  dimensions: int, table):
        pass

    # register params with the parameter manager
    def RegisterParameters(self, parameter_manager):
        pass

    # load sequence
    def LoadSequence(self, sstruct):
        pass

    # load parameter values
    def LoadValues(self, values):
        # load loss function
        pass

    def UseLoss(self, true_mapping,  example_loss):
        pass

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
        pass

    def GetViterbiScore(self):
        pass

    def PredictPairingsViterbi(self):
        pass

    def ComputeViterbiFeatureCounts(self):
        pass

    # MEA inference

    def ComputeInside(self):
        pass

    def ComputeLogPartitionCoefficient(self):
        pass

    def ComputeOutside(self):
        pass

    def ComputeFeatureCountExpectations(self):
        pass

    def ComputePosterior(self):
        pass

    def PredictPairingsPosterior(self,   gamma: int):
        if(not(gamma>0)):
            print("Non-negative gamma expected.")
        unpaired_posterior =[]
        score=[]
        traceback=[]
        
        for i in range(1,self.L+1):
            unpaired_posterior.append(1)
            for j in range(i):
                unpaired_posterior[i] -= self.posterior[self.offset[j]+i]
            for j in range(i+1,self.L+1):
                unpaired_posterior[i] -= self.posterior[self.offset[i]+j]

        for i in range(1,self.L+1):
            unpaired_posterior[i] /= 2 * gamma
        score=[-1]*self.SIZE
        traceback=[-1]*self.SIZE

        #DP

        for i in range(self.L,-1,-1):
            for j in range(i,self.L+1):
                this_score=score[self.offset[i]+j]
                this_traceback = traceback[self.offset[i]+j];
                if(i==j):
                    UPDATE_MAX(this_score, this_traceback,0, 0)
                else:
                    if(self.allow_unpaired_position[i+1]):
                        UPDATE_MAX(this_score, this_traceback, unpaired_posterior[i+1] + score[self.offset[i+1]+j], 1)
                    if (self.allow_unpaired_position[j]):
                        UPDATE_MAX(this_score, this_traceback, unpaired_posterior[j] + score[self.offset[i]+j-1], 2)
                    if(i+2 <= j):
                        if(self.allow_paired[self.offset[i+1]+j]):
                            UPDATE_MAX(this_score, this_traceback, self.posterior[self.offset[i+1]+j] + score[self.offset[i+1]+j-1], 3)



                        p1=score[self.offset[i]+i+1]
                        p2=score[self.offset[i+1]+j]
                        for k in range(i+1,j):
                            UPDATE_MAX(this_score, this_traceback, (p1) + (p2), k+4)
                            p1+=1
                            p2 += self.L-k

        solution=[SStruct.UNPAIRED]*(self.L+1)
        solution[0] = SStruct.UNKNOWN
        traceback_queue = Queue(0)
        traceback_queue.put((0,self.L))

        while(traceback_queue.qsize()!=0):
            t=traceback_queue.get()
            i=t[0]
            j=t[1]

            if(traceback[self.offset[i]+j]==-1):
                print("should not get here ")
                
            elif(traceback[self.offset[i]+j]==0):
                pass
            elif(traceback[self.offset[i]+j]==1):
                traceback_queue.put((i+1,j))
                
            elif(traceback[self.offset[i]+j]==2):
                traceback_queue.put((i,j-1))
            elif(traceback[self.offset[i]+j]==3):
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.put((i+1,j-1))

            else:
                k=traceback[self.offset[i]+j]-4
                traceback_queue.put((i,k))
                traceback_queue.put((k,j))

        return solution

        

    def PredictPairingsPosteriorCentroid(self, gamma):
        pass

    def GetPosterior(self,   posterior_cutoff):
        pass

    # EM inference

    def ComputeInsideESS(self):
        pass

    def ComputeOutsideESS(self):
        pass

    def ComputeFeatureCountExpectationsESS(self):
        pass

    def ComputePosteriorESS(self):
        pass

    def ComputeLogPartitionCoefficientESS(self):
        pass

    def ComputeESS(self):
        pass

    def ComputeGammaMLESum(self, ev_cpd_id,  ignorePairing,  usePosterior,  which_data):
        pass

    def ComputeGammaMLESS(self, ev_cpd_id,  ignoreZeros,  useMLE,  which_data):
        pass

    def ComputeGammaMLEESS(self, ev_cpd_id,  ignoreZeros,  useMLE,  which_data):
        pass

    def GetNumExamplesSeqPairing(self, ev_cpd_id,  ignoreZeros,  which_data):
        pass

    def GetNumExamplesSeq(self, ev_cpd_id,  ignoreZeros,  which_data):
        pass

    def AreZerosInSeqPairing(self,  id_base,  id_pairings,  which_base):
        pass

    def AreZerosInSeq(self,  id_base,  which_base):
        pass

    def LogGammaProb(self,  data,  which_data,  isUnpaired,  seq):
        pass

    def UpdateEvidenceStructures(self,  which_data):
        pass

    def UpdateEvidenceStructures(self):
        pass
