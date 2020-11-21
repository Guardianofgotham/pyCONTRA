class ComputationWrapper():
    def __init__(self):
        pass

    # retrieve list of work units
    def GetAllUnits(self):
        return

    # methods to act on vectors of work units
    def FilterNonparsable(self,  units):
        pass

    def ComputeSolutionNormBound(self,  units,   C,  log_base):
        pass

    def ComputeGradientNormBound(self,  units,   C,  log_base):
        pass

    def Predict(self,  units,   w,  gamma,  log_base):
        pass

    def ComputeLoss(self,  units,   w,  log_base):
        pass

    def ComputeFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        pass

    def ComputeGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        pass

    def ComputeEMFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base):
        pass

    def ComputeEMGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base):
        pass

    def ComputeFunctionSE(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_data):
        pass

    def ComputeGradientSE(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  hyperparam_base):
        pass

    def ComputeGammaMLEFunction(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  evidence_cpd_id1,  evidence_cpd_id2,  evidence_cpd_id3,  evidence_data_scale,  which_data):
        pass

    def ComputeGammaMLEGradient(self,  units,   w,  toggle_use_nonsmooth,  toggle_use_loss,  log_base,  evidence_cpd_id1,  evidence_cpd_id2,  evidence_cpd_id3,  evidence_data_scale,  which_data):
        pass

    def FindZerosInData(self,  units,  evidence_cpd_id1,  evidence_cpd_id2,  which_data):
        pass

    def ComputeGammaMLEScalingFactor(self,  units,   w,  evidence_cpd_id1,  evidence_cpd_id2,  which_data):
        pass

    def ComputeHessianVectorProduct(self,  units,   w,   v,  toggle_use_loss,  log_base):
        pass

    # for debugging

    def SanityCheckGradient(self,  units,   w):
        pass
    # getters
    # def  Options GetOptions()  { return computation_engine.GetOptions(): }
    # def  GetDescriptions()  { return computation_engine.GetDescriptions(): }
    # def InferenceEngine<> GetInferenceEngine() { return computation_engine.GetInferenceEngine(): }
    # def ParameterManager<> GetParameterManager() { return computation_engine.GetParameterManager(): }
    # def ComputationEngine<> GetComputationEngine() { return computation_engine: }
