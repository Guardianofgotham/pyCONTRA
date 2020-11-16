import argparse
from .helper import *
import os
GAMMA_DEFAULT = 1.0
REGULARIZATION_DEFAULT = 1.0
TRAIN_MAX_ITER_DEFAULT = 10
HYPERPARAM_DATA_DEFAULT = 1


def verfiyFileNames(args: argparse.Namespace):
    if(args.parameter_filename != ""):
        checkFileExists(args.parameter_filename)
    if(args.train_examplefile != ""):
        checkFileExists(args.train_examplefile)
    if(args.train_initweights_filename != ""):
        checkFileExists(args.train_initweights_filename)
    if(args.train_priorweights_filename != ""):
        checkFileExists(args.train_priorweights_filename)
    return;


def getArgumentsObject():
    parseObject = argparse.ArgumentParser(
        description="Arguments passed to pyContra")
    parseObject.add_argument("--training_mode", type=str, default="",
                             required=True, help="predict | supervised | em | em-sgd")

    parseObject.add_argument("--verbose", type=bool, default=False,
                             required=False, help="show detailed console output: bool")

    parseObject.add_argument("--log_base", type=float, default=1.0,
                             required=False, help="set base of log-sum-exp: float")

    parseObject.add_argument("--viterbi_parsing", type=bool, default=False,
                             required=False, help="use Viterbi instead of posterior decoding for prediction, or max-margin instead of log-likelihood for training: bool")

    parseObject.add_argument("--allow_noncomplementary", type=bool, default=False,
                             required=False, help="allow non-{AU,CG,GU} pairs: bool")

    parseObject.add_argument("--parameter_filename", type=str, default="",
                             required=False, help="use particular model parameters: str")

    parseObject.add_argument("--use_constraints", type=bool, default=False,
                             required=False, help="use existing constraints (requires BPSEQ or FASTA format input): bool")

    parseObject.add_argument("--centroid_estimator", type=bool, default=False,
                             required=False, help="use centroid estimator (as opposed to MEA estimator): bool")

    parseObject.add_argument("--use_evidence", type=bool, default=False,
                             required=False, help="use experimental evidence (requires BPSEQ format input): bool")

    parseObject.add_argument("--gamma", type=float, default=GAMMA_DEFAULT,
                             required=False, help="set sensivity/specificity tradeoff parameter (default: GAMMA=6)\
                             if GAMMA > 1, emphasize sensitivity\
                             if 0 <= GAMMA <= 1, emphasize specificity\
                             if GAMMA < 0, try tradeoff parameters of 2^k for k = -5,...,10: float")

    parseObject.add_argument("--output_parens_destination", type=str, default="",
                             required=False, help="rite parenthesized output to file or directory: str")

    parseObject.add_argument("--output_bpseq_destination", type=str, default="",
                             required=False, help="write BPSEQ output to file or directory: str")

    parseObject.add_argument("--output_posteriors_cutoff", type=float, default=0,
                             required=False, help="posterior pairing probabilities Cutoffs: float")

    parseObject.add_argument("--output_posteriors_destination", type=str, default="",
                             required=False, help="write posterior pairing probabilities to file or directory: str")

    parseObject.add_argument("--gradient_sanity_check", type=bool, default=False,
                             required=False, help="perform gradient sanity check: bool")

    parseObject.add_argument("--partition_function_only", type=bool, default=False,
                             required=False, help="partition_function_only: bool")

    parseObject.add_argument("--holdout_ratio", type=float, default=0,
                             required=False, help="use fraction F of training data for holdout cross-validation: float")

    parseObject.add_argument("--regularization_coefficient", type=float, default=REGULARIZATION_DEFAULT,
                             required=False, help="perform BFGS training, using a single regularization coefficient C: float")

    parseObject.add_argument("--train_examplefile", type=str, default="",
                             required=False, help="read list of input files from provided text file (instead of as arguments): str")

    parseObject.add_argument("--train_max_iter", type=int, default=TRAIN_MAX_ITER_DEFAULT,
                             required=False, help="for single regularization coefficient the max number of iterations: int")

    parseObject.add_argument("--train_initweights_filename", type=str, default="",
                             required=False, help="train_initweights_filename: str")

    parseObject.add_argument("--train_priorweights_filename", type=str, default="",
                             required=False, help="train_priorweights_filename: str")

    parseObject.add_argument("--num_data_sources", type=int, default=0,
                             required=False, help="num_data_sources: int")

    parseObject.add_argument("--batch_size", type=int, default=1,
                             required=False, help="mini-batch size for stochastic gradient training: int")

    parseObject.add_argument("--s0", type=float, default=0.0001,
                             required=False, help="stepsize for SGD is s0/(1+iter)^s1: float")

    parseObject.add_argument("--s1", type=float, default=0,
                             required=False, help="stepsize for SGD is s0/(1+iter)^s1: float")

    parseObject.add_argument("--hyperparam_data", type=float, default=HYPERPARAM_DATA_DEFAULT,
                             required=False, help="weight on data-only examples: float")

    args = parseObject.parse_args()
    verfiyFileNames(args)
    return args
