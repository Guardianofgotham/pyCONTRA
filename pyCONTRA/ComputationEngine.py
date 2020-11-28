from __future__ import annotations
import argparse
import math
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
        super().__init__(options.verbose)
        self.options = options
        self.description = description
        self.inference_engine = inference_engine
        self.parameter_manager = parameter_manager

    def DistributedComputation(self):
        pass

    # routine for performing an individual work unit
    def DoComputation(self, result: list, shared: SharedInfo, nonshared: NonSharedInfo):
        print("-"*20)
        print(shared.command)
        print("-"*20)
        if shared.command == ProcessingType.CHECK_PARSABILITY:
            self.CheckParsability(result, nonshared)
        elif shared.command == ProcessingType.COMPUTE_SOLUTION_NORM_BOUND:
            self.ComputeSolutionNormBound(result, shared, nonshared)
        elif shared.command == ProcessingType.COMPUTE_GRADIENT_NORM_BOUND:
            self.ComputeGradientNormBound(result, nonshared)
        elif shared.command == ProcessingType.COMPUTE_LOSS:
            self.ComputeLoss(result, shared, nonshared)
        elif shared.command == ProcessingType.COMPUTE_FUNCTION:
            self.ComputeFunctionAndGradient(result, shared, nonshared, False)
        elif shared.command == ProcessingType.COMPUTE_GRADIENT:
            self.ComputeFunctionAndGradient(result, shared, nonshared, True)
        elif shared.command == ProcessingType.COMPUTE_MSTEP_FUNCTION:
            self.ComputeMStepFunctionAndGradient(
                result, shared, nonshared, False)
        elif shared.command == ProcessingType.COMPUTE_MSTEP_GRADIENT:
            self.ComputeMStepFunctionAndGradient(
                result, shared, nonshared, True)
        elif shared.command == ProcessingType.COMPUTE_GAMMAMLE_FUNCTION:
            self.ComputeGammaMLEFunctionAndGradient(
                result, shared, nonshared, False)
        elif shared.command == ProcessingType.COMPUTE_GAMMAMLE_GRADIENT:
            self.ComputeGammaMLEFunctionAndGradient(
                result, shared, nonshared, True)
        elif shared.command == ProcessingType.COMPUTE_GAMMAMLE_SCALING_FACTOR:
            self.ComputeGammaMLEScalingFactor(result, shared, nonshared)
        elif shared.command == ProcessingType.CHECK_ZEROS_IN_DATA:
            self.CheckZerosInData(result, shared, nonshared)
        elif shared.command == ProcessingType.COMPUTE_FUNCTION_SE:
            self.ComputeFunctionAndGradientSE(result, shared, nonshared, False)
        elif shared.command == ProcessingType.COMPUTE_GRADIENT_SE:
            self.ComputeFunctionAndGradientSE(result, shared, nonshared, True)
        elif shared.command == ProcessingType.COMPUTE_HV:
            self.ComputeHessianVectorProduct(result, shared, nonshared)
        elif shared.command == ProcessingType.PREDICT:
            self.Predict(result, shared, nonshared)
        else:
            assert False, "Unknown command type"
    # methods to act on individual work units

    def CheckParsability(self, result: list, nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        self.inference_engine.LoadValues(
            self.parameter_manager.GetNumLogicalParameters())
        # print(sstruct.GetMapping())
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
        result[nonshared.index] = max_entropy/shared.log_base + max_loss
        result *= self.description[nonshared.index].weight

    def ComputeGradientNormBound(self, result: list,   nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)

        w = [1]*self.parameter_manager.GetNumLogicalParameters()
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
        w = [shared.w + self.parameter_manager.GetNumLogicalParameters()] + \
            shared.w
        self.inference_engine.LoadValues(w*shared.log_base)
        self.inference_engine.UpdateEvidenceStructures()
        solution = None
        if(self.options.viterbi_parsing):
            self.inference_engine.ComputeViterbi()
            solution = SStruct(sstruct)
            solution.SetMapping(self.inference_engine.PredictPairingsViterbi())
        else:
            self.inference_engine.ComputeInside()
            self.inference_engine.ComputeOutside()
            self.inference_engine.ComputePosterior()
            solution = SStruct(sstruct)
            solution.SetMapping(
                self.inference_engine.PredictPairingsPosterior(shared.gamma))
        if not shared.use_less:
            raise Exception(
                "Must be using loss function in order to compute loss.")
        self.inference_engine.LoadValues([0]*len(w))
        self.inference_engine.UseConstraints(solution.GetMapping())
        self.inference_engine.UpdateEvidenceStructures()
        self.inference_engine.ComputeViterbi()

        result.clear()
        result.append(self.inference_engine.GetViterbiScore())

        result *= self.description[nonshared.index].weight
        result[-1] /= shared.log_base

    def ComputeFunctionAndGradient(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo, need_gradient: bool):
        sstruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)

        w = [shared.w + self.parameter_manager]*shared.w
        self.inference_engine.LoadValues(w*shared.log_base)
        # if defined(HAMMING_LOSS)
        # if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
        # endif
        unconditional_score = None
        unconditional_counts = list()

        if shared.use_nonsmooth:
            self.inference_engine.ComputeViterbi()
            unconditional_score = self.inference_engine.GetViterbiScore()
            if need_gradient:
                unconditional_counts = self.inference_engine.ComputeViterbiFeatureCounts()
        else:
            self.inference_engine.ComputeInside()
            unconditional_score = self.inference_engine.ComputeLogPartitionCoefficient()
            if need_gradient:
                self.inference_engine.ComputeOutside()
                unconditional_counts = self.inference_engine.ComputeFeatureCountExpectations()

        conditional_score = None
        conditional_counts = list()

        self.inference_engine.UseConstraints(sstruct.GetMapping())
        if shared.use_nonsmooth:
            self.inference_engine.ComputeViterbi()
            conditional_score = self.inference_engine.GetViterbiScore()
            if need_gradient:
                conditional_counts = self.inference_engine.ComputeViterbiFeatureCounts()
        else:
            self.inference_engine.ComputeInside()
            conditional_score = self.inference_engine.ComputeLogPartitionCoefficient()
            if need_gradient:
                self.inference_engine.ComputeOutside()
                conditional_counts = self.inference_engine.ComputeFeatureCountExpectations()

        result.clear()
        if need_gradient:
            result = unconditional_counts-conditional_counts
        if(conditional_score > unconditional_score):
            raise Exception(
                "Conditional score cannot exceed unconditional score.")
        result.append(unconditional_score-conditional_score)

        if conditional_score < NEG_INF/2:
            print(
                f"Unexpected bad parse for file: {self.descriptions[nonshared.index].input_filename}")
            result = [0]*len(result)
        if NONCONVEX_MULTIPLIER != 0:
            pass

        if result[-1] < 0:
            if result[-1] < 1e-6:
                print(
                    f"Encountered negative function value for {self.descriptions[nonshared.index].input_filename}: {result[-1]}")
                self.parameter_manager.WriteToFile(None, w)
        result *= self.description[nonshared.index].weight
        result[-1] /= shared.log_base

    def ComputeMStepFunctionAndGradient(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo, need_gradient: bool):
        sstruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        w = [shared.w+self.parameter_manager]*shared.w
        self.inference_engine.LoadValues(w*shared.log_base)
        self.inference_engine.UpdateEvidenceStructures()

        # if defined(HAMMING_LOSS)
        # Error("HAMMING_LOSS not implemented within EM training");
        # if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
        # endif
        unconditional_score = None
        unconditional_counts = list()
        if shared.use_nonsmooth:
            raise Exception(
                "Viterbi training not supported within EM training")
            self.inference_engine.ComputeViterbi()
            unconditional_score = self.inference_engine.GetViterbiScore()
            if need_gradient:
                unconditional_counts = self.inference_engine.ComputeViterbiFeatureCounts()
        else:
            self.inference_engine.ComputeInside()
            unconditional_score = self.inference_engine.ComputeLogPartitionCoefficient()
            if need_gradient:
                self.inference_engine.ComputeOutside()
                unconditional_counts = self.inference_engine.ComputeFeatureCountExpectations()
        conditional_score = None
        conditional_counts = list()
        if shared.use_nonsmooth:
            raise Exception(
                "Viterbi training not supported within EM training")
            self.inference_engine.ComputeViterbi()
            conditional_score = self.inference_engine.GetViterbiScore()
            if need_gradient:
                conditional_counts = self.inference_engine.ComputeViterbiFeatureCounts()
        else:
            if not sstruct.HasStruct():
                self.inference_engine.ComputeInsideESS()
                self.inference_engine.ComputeOutsideESS()
                conditional_counts = self.inference_engine.ComputeFeatureCountExpectationsESS()
                conditional_score = DotProduct(w, conditional_counts)
            else:
                self.inference_engine.UseConstraints(sstruct.GetMapping())
                self.inference_engine.UpdateEvidenceStructures()
                self.inference_engine.ComputeInside()
                conditional_score = self.inference_engine.ComputeLogPartitionCoefficient()
                if need_gradient:
                    self.inference_engine.ComputeOutside()
                    conditional_score = self.inference_engine.ComputeFeatureCountExpectations()
        result.clear()
        if need_gradient:
            result = unconditional_counts-conditional_counts

        if(unconditional_score < conditional_score):
            raise Exception(
                "Conditional score cannot exceed unconditional score.")
        result.append(unconditional_score-conditional_score)

        if conditional_score < NEG_INF/2:
            print(
                f"Unexpected bad parse for file: {self.descriptions[nonshared.index].input_filename}")
            result = [0]*len(result)

        if NONCONVEX_MULTIPLIER != 0:
            raise Exception(
                "Nonconvex training not supported within EM training")

        if result[-1] < 0:
            if result[-1] < 1e-6:
                print(
                    f"Encountered negative function value for {self.descriptions[nonshared.index].input_filename}: {result[-1]}")
                self.parameter_manager.WriteToFile(None, w)
                exit(-255)
            result = [0]*len(result)

        result *= self.description[nonshared.index].weight
        result[-1] /= shared.log_base

    def ComputeGammaMLEFunctionAndGradient(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo, need_gradient: bool):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        which_data = shared.which_data
        if not sstruct.HasEvidence(which_data):
            result.clear()
            if need_gradient:
                result.extend([0, 0, 0])
            result.append(0)
            return
        self.inference_engine.LoadSequence(sstruct)
        self.inference_engine.UseConstraints(sstruct.GetMapping())
        w = [shared.w + self.parameter_manager.GetNumLogicalParameters()] * \
            shared.w
        self.inference_engine.LoadValues(w*shared.log_base)

        self.inference_engine.UpdateEvidenceStructures()
        update_gammamle_sssum = None
        update_gammamle_sssumlog = None
        update_gammamle_num_examples = None
        update_gammamle_ssq = None
        ll = None
        stats = [0]*3
        j = shared.id_base
        k = shared.id_pairing
        areZeros = shared.areZeros
        scale = shared.evidence_data_scale
        evidence_cpd_id = list()
        evidence_cpd_id.extend([j, k, areZeros])
        current_k = math.exp(w[self.parameter_manager.GetLogicalIndex(
            self.inference_engine.GetLogScoreEvidence(0, j, k, which_data))])
        current_theta = math.exp(w[self.parameter_manager.GetLogicalIndex(
            self.inference_engine.GetLogScoreEvidence(1, j, k, which_data))])
        if shared.use_nonsmooth:
            raise Exception(
                "Viterbi training not supported within EM training")
        else:
            if not sstruct.HasStruct():
                self.inference_engine.ComputeInsideESS()
                self.inference_engine.ComputeOutsideESS()
                self.inference_engine.ComputePosteriorESS()
                stats = self.inference_engine.ComputeGammaMLEESS(
                    evidence_cpd_id, not areZeros, not areZeros, which_data)
                update_gammamle_num_examples = self.inference_engine.GetNumExamplesSeq(
                    evidence_cpd_id, False, which_data)
                update_gammamle_sssum = stats[0]
                update_gammamle_ssq = stats[2]
                update_gammamle_sssumlog = stats[1]
                if not areZeros:
                    ll = (current_k-1)*update_gammamle_sssumlog - update_gammamle_sssum / current_theta - update_gammamle_num_examples * \
                        current_k * \
                        math.log(current_theta) - \
                        update_gammamle_num_examples * lgamma(current_k)
                else:
                    stats = self.inference_engine.ComputeGammaMLEESS(
                        evidence_cpd_id, True, True, which_data)
                    update_gammamle_sssumlog_nozeros = stats[1]
                    update_gammamle_num_examples_nozeros = self.inference_engine.GetNumExamplesSeq(
                        evidence_cpd_id, True, which_data)
                    ll = (current_k - 1) * update_gammamle_sssumlog_nozeros - update_gammamle_sssum / current_theta - \
                        update_gammamle_num_examples_nozeros * current_k * \
                        math.log(
                            current_theta) - update_gammamle_num_examples_nozeros * lgamma(current_k)
                update_gammamle_sssum = update_gammamle_sssum / scale
                update_gammamle_sssumlog = update_gammamle_sssumlog - \
                    update_gammamle_ssq * math.log(scale)
            else:
                stats = self.inference_engine.ComputeGammaMLESS(
                    evidence_cpd_id, not areZeros, not areZeros, which_data)
                update_gammamle_num_examples = inference_engine.GetNumExamplesSeqPairing(
                    evidence_cpd_id, False, which_data)

                update_gammamle_sssum = stats[0]
                update_gammamle_sssumlog = stats[1]

                if not areZeros:
                    ll = (current_k - 1) * update_gammamle_sssumlog - update_gammamle_sssum / current_theta - \
                        update_gammamle_num_examples * current_k * \
                        math.log(current_theta) - \
                        update_gammamle_num_examples * lgamma(current_k)
                else:
                    stats = self.inference_engine.ComputeGammaMLESS(
                        evidence_cpd_id, True, True, which_data)
                    update_gammamle_sssumlog_nozeros = stats[1]
                    update_gammamle_num_examples_nozeros = self.inference_engine.GetNumExamplesSeqPairing(
                        evidence_cpd_id, True, which_data)
                    ll = (current_k - 1) * update_gammamle_sssumlog_nozeros - update_gammamle_sssum / current_theta - \
                        update_gammamle_num_examples_nozeros * current_k * \
                        math.log(
                            current_theta) - update_gammamle_num_examples_nozeros * lgamma(current_k)
                update_gammamle_sssum = update_gammamle_sssum / scale
                update_gammamle_sssumlog = update_gammamle_sssumlog - \
                    update_gammamle_num_examples * math.log(scale)

        result.clear()
        if need_gradient:
            result.extend(
                [update_gammamle_sssum, update_gammamle_sssumlog, update_gammamle_num_examples])
        result.append(ll)
        result *= self.description[nonshared.index].weight

    def ComputeHessianVectorProduct(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        w = [shared.w+self.parameter_manager.GetNumLogicalParameters()] * \
            shared.w
        v = [shared.v+self.parameter_manager.GetNumLogicalParameters()] * \
            shared.v
        if self.options.viterbi_parsing:
            raise Exception(
                "Should not use Hessian-vector products with Viterbi parsing.")
        EPSILON = 1e-8
        shared_temp = SharedInfo(shared)
        result2 = list()

        for i in range(0, self.parameter_manager.GetNumLogicalParameters()):
            shared_temp.w[i] = shared.w[i] + EPSILON * v[i]
        ComputeFunctionAndGradient(result, shared_temp, nonshared, True)

        for i in range(0, self.parameter_manager.GetNumLogicalParameters()):
            shared_temp.w[i] = shared.w[i] - EPSILON * v[i]
        ComputeFunctionAndGradient(result2, shared_temp, nonshared, True)
        result = (result-result2)/(2*EPSILON)

    def Predict(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        result.clear()
        sstruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)
        if (self.options.use_constraints):
            self.inference_engine.UseConstraints(sstruct.GetMapping())
        # // load parameters
        w = shared.w[:]
        # print(w)
        self.inference_engine.LoadValues(w)# * shared.log_base)
        self.inference_engine.UpdateEvidenceStructures()
        # // perform inference
        solution = None
        if (self.options.viterbi_parsing):
            if (self.options.use_evidence):
                raise Exception(
                    "Viterbi parsing is not supported with evidence yet")
            self.inference_engine.ComputeViterbi()
            if self.options.partition_function_only:
                print(
                    f"Viterbi score for {self.descriptions[nonshared.index].input_filename} : {inference_engine.GetViterbiScore()}")
                return
            solution = SStruct(sstruct)
            solution.SetMapping(self.inference_engine.PredictPairingsViterbi())
        else:
            if self.options.use_evidence:
                self.inference_engine.ComputeInsideESS()
                if self.options.partition_function_only:
                    print(
                        f"Log partition coefficient for {self.descriptions[nonshared.index].input_filename} : {inference_engine.ComputeLogPartitionCoefficientESS()}")
                    return
                self.inference_engine.ComputeOutsideESS()
                self.inference_engine.ComputePosteriorESS()
            else:
                self.inference_engine.ComputeInside()
                if self.options.partition_function_only:
                    print(
                        f"Log partition coefficient for {self.description[nonshared.index].input_filename} :")
                    return
                self.inference_engine.ComputeOutside()
                self.inference_engine.ComputePosterior()
            solution = SStruct(sstruct, None)
            if self.options.centroid_estimator:
                print(f"Predicting using centroid estimator.")
                solution.SetMapping(
                    inference_engine.PredictPairingsPosteriorCentroid(shared.gamma))
            else:
                print(f"Predicting using MEA estimator.")
                solution.SetMapping(
                    self.inference_engine.PredictPairingsPosterior(shared.gamma))
        if self.options.output_parens_destination == "":
            filename = MakeOutputFilename(self.description[nonshared.index].input_filename,
                                          self.options.output_parens_destination, self.options.gamma < 0, shared.gamma)
            pass
        if self.options.output_bpseq_destination != "":
            filename = MakeOutputFilename(self.descriptions[nonshared.index].input_filename,
                                          self.options.output_bpseq_destination, self.options.gamma < 0, shared.gamma)
            pass
        if self.options.output_posteriors_destination != "":
            pass
        return

    def CheckZerosInData(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        which_data = shared.which_data
        sstruct: SStruct = self.description[nonshared.index].sstruct
        self.inference_engine.LoadSequence(sstruct)

        self.inference_engine.LoadValues(
            [0]*self.parameter_manager.GetNumLogicalParameters())
        self.inference_engine.UseConstraints(sstruct.GetMapping())
        result.clear()
        result = [0]*len(self.description)
        if not sstruct.HasEvidence(which_data):
            result[nonshared.index] = 0
            return
        self.inference_engine.UpdateEvidenceStructures(which_data)
        areZeroes = None
        if not sstruct.HasEvidence():
            areZeroes = self.inference_engine.AreZerosInSeq(
                shared.id_base, shared.which_data)
        else:
            areZeroes = self.inference_engine.AreZerosInSeqPairing(
                shared.id_base, shared.id_pairing, shared.which_data)
        result[nonshared.index] = areZeroes

    def ComputeGammaMLEScalingFactor(self, result: list,   shared: SharedInfo,   nonshared: NonSharedInfo):
        sstruct: SStruct = self.description[nonshared.index].sstruct
        which_data = shared.which_data

        if not sstruct.HasEvidence(which_data):
            result.clear()
            result.extend([0,0])
        self.inference_engine.LoadSequence(sstruct)
        self.inference_engine.UseConstraints(sstruct.GetMapping())
        w = [shared.w+self.parameter_manager]*shared.w
        self.inference_engine.LoadValues(w*shared.log_base)
        self.inference_engine.UpdateEvidenceStructures()
        j = shared.id_base
        k = shared.id_pairing
        evidence_cpd_id = [j, k]
        if not sstruct.HasStruct():
            self.inference_engine.ComputeInsideESS()
            self.inference_engine.ComputeOutsideESS()
            self.inference_engine.ComputePosteriorESS()
            update_gammamle_sssum = self.inference_engine.ComputeGammaMLESum(evidence_cpd_id, True, True, which_data)  # ignore pairing since don't know it, use posterior
            update_gammamle_num_examples = self.inference_engine.GetNumExamplesSeq(evidence_cpd_id,False, which_data)  # don't ignore zeros
        else:
            update_gammamle_sssum = self.inference_engine.ComputeGammaMLESum(evidence_cpd_id, False, False, which_data)  #don't ignore pairing since we know it, don't use posterior
            update_gammamle_num_examples = self.inference_engine.GetNumExamplesSeqPairing(evidence_cpd_id,False, which_data)  # don't ignore zeros

        result.clear()
        result.extend([update_gammamle_sssum, update_gammamle_num_examples])


    def ComputeFunctionAndGradientSE(self, result,   shared,   nonshared, need_gradient):
        pass
    
    def GetOptions(self):
        return self.options

    def GetDescriptions(self):
        return self.description

    def GetInferenceEngine(self):
        return self.inference_engine
    
    def GetParameterManager(self):
        return self.parameter_manager