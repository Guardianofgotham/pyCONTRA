from __future__ import annotations
from pyCONTRA.config import *


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
        num_data_sources = None
        L = None
        SIZE = None

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
    def LoadValues(self, alues):
        pass

    # load loss function
    def UseLoss(self, true_mapping,  example_loss):
        pass

    # use raints
    def Useraints(self, true_mapping):
        pass

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
