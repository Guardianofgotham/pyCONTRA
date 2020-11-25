from pyCONTRA.config import *

class SharedInfo:
    def __init__(self):
        self.command = None
        self.w = [0]*SHARED_PARAMETER_SIZE
        self.v = [0]*SHARED_PARAMETER_SIZE
        self.use_nonsmooth = False
        self.use_loss = False
        self.gamma = None
        self.log_base = None
        self.evidence_data_scale = None
        self.id_base=None
        self.id_pairing=None
        self.which_data=None
        self.areZeros=None
        self.hyperparam_data=None
    
    