from pyCONTRA.ComputationEngine import *
from pyCONTRA.SharedInfo import *
from pyCONTRA.NonSharedInfo import *
from pyCONTRA.ProcessingType import *


class ComputationWrapper(object):
    def __init__(self, computation_engine: ComputationEngine):
        self.computation_engine = computation_engine
        self.shared_info = SharedInfo()
        self.nonshared_info = list()
        self.cached_toggle_use_nonsmooth = False
        self.cached_toggle_use_loss = False
        self.cached_units = list()
        self.cached_w = list()
        self.cached_function = list()
        self.cached_function_gammamle = list()
        self.cached_gradient_gammamle = list()

    # retrieve list of work units
    def GetAllUnits(self):
        ret=list()
        for i in range(len(self.GetDescriptions())):
            ret.append(i)
        return ret

    # methods to act on vectors of work units
    def FilterNonparsable(self,  units: list):
        ret = []
        if(not self.computation_engine.IsMasterNode()):
            raise Exception("Routine should only be called by master process.")
        parsable=[]
        self.shared_info.command=ProcessingType.CHECK_PARSABILITY

        #self.nonshared_info.resize(units)
        for i in range(len(units)):
            self.nonshared_info.append(NonSharedInfo())
            self.nonshared_info[i].index=units[i]
        
        parsable = self.computation_engine.DistributeComputation(parsable,self.shared_info,self.nonshared_info)
        for i in range(len(units)):
            print(parsable, units)
            # assert not(units[i] >= 0 and units[i] < len(parsable) ), "Out-of-bounds index."
            if(parsable[units[i]]):
                ret.append(units[i])
            else:
                print("No valid parse for file: ", self.GetDescriptions()[units[i]].input_filename )
        return units

    def ComputeSolutionNormBound(self,  units,   C,  log_base):
        raise Exception("Not implemented")

    def ComputeGradientNormBound(self,  units,   C,  log_base):
        raise Exception("Not implemented")

    def Predict(self,  units: list,   w: list,  gamma: float,  log_base: float):
        if not self.computation_engine.IsMasterNode():
            raise Exception("Routine should only be called by master process.")
        if len(w) > SHARED_PARAMETER_SIZE:
            raise Exception(f"SHARED_PARAMETER_SIZE in Config.hpp too small; increase to at least {len(w)}.")
        if self.GetOptions().verbose:
            print(f"Performing predictions with gamma={gamma}...")
        ret = list()
        self.shared_info.command = ProcessingType.PREDICT
        for i in range(0, len(w)):
            self.shared_info.w[i] = w[i]
        self.shared_info.w = self.shared_info.w[:len(w)]
        self.shared_info.gamma=gamma
        self.shared_info.log_base=log_base
        self.nonshared_info.clear()
        for i in range(len(units)):
            self.nonshared_info.append(NonSharedInfo())
            self.nonshared_info[i].index = units[i]
        self.computation_engine.DistributeComputation(ret, self.shared_info, self.nonshared_info)

    def ComputeLoss(self,  units,   w,  log_base):
        raise Exception("Not implemented")

    def ComputeFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        raise Exception("Not implemented")

    def ComputeGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        raise Exception("Not implemented")

    def ComputeEMFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base):
        raise Exception("Not implemented")

    def ComputeEMGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base):
        raise Exception("Not implemented")

    def ComputeFunctionSE(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        raise Exception("Not implemented")

    def ComputeGradientSE(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_base):
        raise Exception("Not implemented")

    def ComputeGammaMLEFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  evidence_cpd_id1,  evidence_cpd_id2,  evidence_cpd_id3,  evidence_data_scale,  which_data):
        raise Exception("Not implemented")

    def ComputeGammaMLEGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  evidence_cpd_id1,  evidence_cpd_id2,  evidence_cpd_id3,  evidence_data_scale,  which_data):
        raise Exception("Not implemented")

    def FindZerosInData(self,  units,  evidence_cpd_id1,  evidence_cpd_id2,  which_data):
        raise Exception("Not implemented")

    def ComputeGammaMLEScalingFactor(self,  units,   w,  evidence_cpd_id1,  evidence_cpd_id2,  which_data):
        raise Exception("Not implemented")

    def ComputeHessianVectorProduct(self,  units,   w,   v,  toggle_use_loss,  log_base):
        raise Exception("Not implemented")

    # for debugging

    def SanityCheckGradient(self,  units,   w):
        raise Exception("Not implemented")
    # getters
    def GetOptions(self):
        return self.computation_engine.GetOptions()
    
    def GetDescriptions(self):
        return self.computation_engine.GetDescriptions()

    def GetInferenceEngine(self):
        return self.computation_engine.GetInferenceEngine()

    def GetParameterManager(self):
        return self.computation_engine.GetParameterManager()

    def GetComputationEngine(self):
        return self.computation_engine
