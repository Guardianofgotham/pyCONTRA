from __future__ import annotations
import argparse
from pyCONTRA.DistributedComputation import *
from pyCONTRA.ProcessingType import *
from pyCONTRA.NonSharedInfo import *
from pyCONTRA.config import *
from pyCONTRA.SStruct import *
from pyCONTRA.SharedInfo import *
from pyCONTRA.InferenceEngine import *
from pyCONTRA.ParameterManager import *


class ComputationEngine(DistributedComputationBase):
    def __init__(self, options: argparse.Namespace, description: list, inference_engine: InferenceEngine, parameter_manager: ParameterManager):
        self.options = options
        self.description = description
        self.inference_engine = inference_engine
        self.parameter_manager = parameter_manager

    def DistributedComputation(self):
        pass

    # routine for performing an individual work unit
    def DoComputation(self, result: list, shared: SharedInfo, nonshared: NonSharedInfo):
        if shared.commmand == ProcessingType.CHECK_PARSABILITY:
            self.CheckParsability(result, nonshared)
        elif shared.command == ProcessingType.COMPUTE_SOLUTION_NORM_BOUND:
            self.ComputeSolutionNormBound(result, shared, nonshared)
        elif self.command == ProcessingType.COMPUTE_GRADIENT_NORM_BOUND:
            self.ComputeGradientNormBound(result, nonshared)
        elif self.command == ProcessingType.COMPUTE_LOSS:
            self.ComputeLoss(result, shared, nonshared)
        elif self.command == ProcessingType.COMPUTE_FUNCTION:
            self.ComputeFunctionAndGradient(result, shared, nonshared, False)
        elif self.command == ProcessingType.COMPUTE_GRADIENT:
            self.ComputeFunctionAndGradient(result, shared, nonshared, True)
        elif self.command == ProcessingType.COMPUTE_MSTEP_FUNCTION:
            self.ComputeMStepFunctionAndGradient(
                result, shared, nonshared, False)
        elif self.command == ProcessingType.COMPUTE_MSTEP_GRADIENT:
            self.ComputeMStepFunctionAndGradient(
                result, shared, nonshared, True)
        elif self.command == ProcessingType.COMPUTE_GAMMAMLE_FUNCTION:
            self.ComputeGammaMLEFunctionAndGradient(
                result, shared, nonshared, False)
        elif self.command == ProcessingType.COMPUTE_GAMMAMLE_GRADIENT:
            self.ComputeGammaMLEFunctionAndGradient(
                result, shared, nonshared, True)
        elif self.command == ProcessingType.COMPUTE_GAMMAMLE_SCALING_FACTOR:
            self.ComputeGammaMLEScalingFactor(result, shared, nonshared)
        elif self.command == ProcessingType.CHECK_ZEROS_IN_DATA:
            self.CheckZerosInData(result, shared, nonshared)
        elif self.command == ProcessingType.COMPUTE_FUNCTION_SE:
            self.ComputeFunctionAndGradientSE(result, shared, nonshared, False)
        elif self.command == ProcessingType.COMPUTE_GRADIENT_SE:
            self.ComputeFunctionAndGradientSE(result, shared, nonshared, True)
        elif self.command == ProcessingType.COMPUTE_HV:
            self.ComputeHessianVectorProduct(result, shared, nonshared)
        elif self.command == ProcessingType.PREDICT:
            Predict(result, shared, nonshared)
        else:
            assert False, "Unknown command type"
    # methods to act on individual work units

    def CheckParsability(self, result: list, nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        self.inference_engine.LoadValues(
            self.parameter_manager.GetNumLogicalParameters())
        self.inference_engine.UseConstraints(sstruct.GetMapping())
        self.inference_engine.UpdateEvidenceStructures()
        self.inference_engine.ComputeViterbi()
        conditional_score = self.inference_engine.GetViterbiScore()
        result.clear()
        result = [0]*len(self.description)
        if conditional_score < NEG_INF/2:
            result[nonshared.index] = 0
        else:
            result[nonshared.index] = 1

    def ComputeSolutionNormBound(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        max_entropy = 0
        max_loss = 0
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        w = [0]*self.parameter_manager.GetNumLogicalParameters()
        self.inference_engine.LoadValues(w)
        self.inference_engine.UpdateEvidenceStructures()
        if SMOOTH_MAX_MARGIN != 0:
            if self.options.viterbi_parsing:
                self.inference_engine.ComputeInside()
                max_entropy += self.inference_engine.ComputeLogPartitionCoefficient()
        # if defined(HAMMING_LOSS)
        # inference_engine.UseLoss(sstruct.GetMapping(), RealT(HAMMING_LOSS));
        # inference_engine.ComputeViterbi();
        # max_loss += inference_engine.GetViterbiScore();
        # endif
        result.clear()
        result = [0]*len(self.description)
        result[nonshared.index] = max_entropy/shared.log_base +max_loss
        result *= self.description[nonshared.index].weight

    def ComputeGradientNormBound(self, result: list,   nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)

        w= [1]*self.parameter_manager.GetNumLogicalParameters()
        self.inference_engine.LoadValues(w)
        self.inference_engine.UpdateEvidenceStructures()

        self.inference_engine.ComputeViterbi()
        max_L1_norm = self.inference_engine.GetViterbiScore()
        result.clear()
        result = [0]*len(self.description)
        result[nonshared.index] = max_L1_norm 

    def ComputeLoss(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        w = [shared.w + self.parameter_manager.GetNumLogicalParameters()]+shared.w
        self.inference_engine.LoadValues(w*shared.log_base)
        self.inference_engine.UpdateEvidenceStructures()
        solution=None
        if(self.options.viterbi_parsing):
            self.inference_engine.ComputeViterbi()
            solution = SStruct(sstruct)
            solution.SetMapping(self.inference_engine.PredictPairingsViterbi())
        else:
            self.inference_engine.ComputeInside()
            self.inference_engine.ComputeOutside()
            self.inference_engine.ComputePosterior()
            solution = SStruct(sstruct)
            solution.SetMapping(self.inference_engine.PredictPairingsPosterior(shared.gamma))
        if not shared.use_less:
            raise Exception("Must be using loss function in order to compute loss.")
        self.inference_engine.LoadValues([0]*len(w))
        self.inference_engine.UseConstraints(solution.GetMapping())
        self.inference_engine.UpdateEvidenceStructures()
        self.inference_engine.ComputeViterbi()

        result.clear()
        result.append(self.inference_engine.GetViterbiScore())

        result *= self.description[nonshared.index].weight
        result[-1] /= shared.log_base

    def ComputeFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeMStepFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeGammaMLEFunctionAndGradient(self, result,   shared,   nonshared, need_gradient):
        pass

    def ComputeHessianVectorProduct(self, result,   shared,   nonshared):
        pass

    def Predict(self, result,   shared,   nonshared):
        pass

    def CheckZerosInData(self, result,   shared,   nonshared):
        pass

    def ComputeGammaMLEScalingFactor(self, result,   shared,   nonshared):
        pass

    def ComputeFunctionAndGradientSE(self, result,   shared,   nonshared, need_gradient):
        pass
