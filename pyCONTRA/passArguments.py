import argparse
GAMMA_DEFAULT = 1.0
REGULARIZATION_DEFAULT = 1.0
TRAIN_MAX_ITER_DEFAULT = 10
HYPERPARAM_DATA_DEFAULT = 1

def getArgumentsObject():
    parseObject = argparse.ArgumentParser(
        description="Arguments passed to pyContra")
    parseObject.add_argument("--training_mode", type=str, default=-1,
                             required=True, help="predict or train")

    parseObject.add_argument("--verbose", type=bool, default=False,
                             required=False, help="True, For More Detailed Output.")

    parseObject.add_argument("--log_base", type=float, default=1.0,
                             required=False, help="Log Base: Real Value")

    parseObject.add_argument("--viterbi_parsing", type=bool, default=False,
                             required=False, help="viterbi_parsing: bool")

    parseObject.add_argument("--allow_noncomplementary", type=bool, default=False,
                             required=False, help="allow_noncomplementary: bool")

    parseObject.add_argument("--parameter_filename", type=str, default="",
                             required=False, help="parameter_filename: str")

    parseObject.add_argument("--use_constraints", type=bool, default=False,
                             required=False, help="use_constraints: bool")

    parseObject.add_argument("--centroid_estimator", type=bool, default=False,
                             required=False, help="centroid_estimator: bool")

    parseObject.add_argument("--use_evidence", type=bool, default=False,
                             required=False, help="use_evidence: bool")

    parseObject.add_argument("--gamma", type=float, default=GAMMA_DEFAULT,
                             required=False, help="centroid_estimator: float")

    parseObject.add_argument("--output_parens_destination", type=str, default="",
                             required=False, help="output_parens_destination: str")

    parseObject.add_argument("--output_bpseq_destination", type=str, default="",
                             required=False, help="output_bpseq_destination: str")

    parseObject.add_argument("--output_posteriors_cutoff", type=float, default=0,
                             required=False, help="output_posteriors_cutoff: float")

    parseObject.add_argument("--output_posteriors_destination", type=str, default="",
                             required=False, help="output_posteriors_destination: str")

    parseObject.add_argument("--gradient_sanity_check", type=bool, default=False,
                             required=False, help="gradient_sanity_check: bool")

    parseObject.add_argument("--partition_function_only", type=bool, default=False,
                             required=False, help="partition_function_only: bool")

    parseObject.add_argument("--holdout_ratio", type=float, default=0,
                             required=False, help="holdout_ratio: float")

    parseObject.add_argument("--regularization_coefficient", type=float, default=REGULARIZATION_DEFAULT,
                             required=False, help="regularization_coefficient: float")

    parseObject.add_argument("--train_examplefile", type=str, default="",
                             required=False, help="train_examplefile: str")

    parseObject.add_argument("--train_max_iter", type=int, default=TRAIN_MAX_ITER_DEFAULT,
                             required=False, help="train_max_iter: int")

    parseObject.add_argument("--train_initweights_filename", type=str, default="",
                             required=False, help="train_initweights_filename: str")

    parseObject.add_argument("--train_priorweights_filename", type=str, default="",
                             required=False, help="train_priorweights_filename: str")

    parseObject.add_argument("--num_data_sources", type=int, default=0,
                             required=False, help="num_data_sources: int")

    parseObject.add_argument("--batch_size", type=int, default=1,
                             required=False, help="batch_size: int")

    parseObject.add_argument("--s0", type=float, default=0.0001,
                             required=False, help="s0: float")

    parseObject.add_argument("--s1", type=float, default=0,
                             required=False, help="s0: float")

    parseObject.add_argument("--hyperparam_data", type=float, default=HYPERPARAM_DATA_DEFAULT,
                             required=False, help="hyperparam_data: float")

    args = parseObject.parse_args()
    return args
