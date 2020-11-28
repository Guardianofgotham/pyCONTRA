from __future__ import annotations
from pyCONTRA.config import *
from pyCONTRA.SStruct import *


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

        s, offset = [], []

        allow_unpaired_position = []
        allow_unpaired, allow_paired = [], []
        loss_unpaired_position = []
        loss_unpaired, loss_paired = [], []

        # ## dynamic programming matrices
        FCt, F5t, FMt, FM1t = [], [], [], []
        FCv, F5v, FMv, FM1v = [], [], [], []
        FCi, F5i, FMi, FM1i = [], [], [], []
        FCo, F5o, FMo, FM1o = [], [], [], []

        FCi_ess, F5i_ess, FMi_ess, FM1i_ess = [], [], [], []
        FCo_ess, F5o_ess, FMo_ess, FM1o_ess = [], [], [], []

        posterior = []

        score_unpaired_position = []
        score_paired_position = []
        score_unpaired_position_raw = []
        score_paired_position_raw = []
        score_base_pair = list()
        for i in range(M+1):
            score_base_pair.append([])
            for j in range(M+1):
                score_base_pair[-1].append((0, 0))

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

    def PredictPairingsPosterior(self,   gamma):
        pass

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
