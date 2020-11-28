from pyCONTRA.ComputationEngine import *
class ComputationWrapper():
    def __init__(self, computation_engine: ComputationEngine):
        self.computation_engine = computation_engine
        self.shared_info = SharedInfo()
        self.nonshared_info = NonSharedInfo()
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
        for i in range(len(GetDescriptions())):
            ret.append(i)
        return ret

    # methods to act on vectors of work units
    def FilterNonparsable(self,  units: list):
        ret = []
        if(not self.computation_engine.IsMasterNode()):
            print("Routine should only be called by master process.")
        parsable=[]
        self.shared_info.command=CHECK_PARSABILITY

        #self.nonshared_info.resize(units)
        for i in range(len(units)):
            self.nonshared_info[i].index=units[i]
        
        self.computation_engine.DistributeComputation(parsable,self.shared_info,self.nonshared_info)
        for i in range(len(units)):
            if(not(units[i] >= 0 and units[i] < int(parsable.size()))):
                print("Out-of-bounds index.")
            if(parsable[units[i]]):
                ret.push_back(units[i])
            else:
                print("No valid parse for file: ", GetDescriptions()[units[i]].input_filename )
        return ret

    def ComputeSolutionNormBound(self,  units,   C,  log_base):
        pass

    def ComputeGradientNormBound(self,  units,   C,  log_base):
        pass

    def Predict(self,  units: list,   w: list,  gamma,  log_base):
        if not self.computation_engine.IsMasterNode():
            raise Exception("Routine should only be called by master process.")
        if len(w) > SHARED_PARAMETER_SIZE:
            raise Exception(f"SHARED_PARAMETER_SIZE in Config.hpp too small; increase to at least {len(w)}.")
        if self.GetOptions().verbose_output:
            print(f"Performing predictions with gamma={gamma}...")
        ret = list()
        self.shared_info.command = ProcessingType.PREDICT
        for i in range(0, len(w)):
            self.shared_info.w[i] = w[i]
        self.shared_info.gamma=gamma
        self.shared_info.log_base=log_base
        self.nonshared_info = self.nonshared_info[:len(units)]
        for i in range(len(units)):
            self.nonshared_info[i].index = units[i]
        self.computation_engine.DistributeComputation(ret, self.shared_info, self.nonshared_info)

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
