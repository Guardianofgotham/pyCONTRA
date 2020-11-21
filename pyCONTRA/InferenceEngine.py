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
        
        # // dynamic programming matrices
        FCt, F5t, FMt, FM1t = [], [], [], []
        FCv, F5v, FMv, FM1v = [], [], [], []
        FCi, F5i, FMi, FM1i = [], [], [], []
        FCo, F5o, FMo, FM1o = [], [], [], []

        FCi_ess, F5i_ess, FMi_ess, FM1i_ess = [], [], [], []
        FCo_ess, F5o_ess, FMo_ess, FM1o_ess = [], [], [], []

        posterior = []

        score_unpaired_position=[]
        score_paired_position=[]
        score_unpaired_position_raw=[]
        score_paired_position_raw=[]
        pass
