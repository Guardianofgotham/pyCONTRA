# ## Global configuration file.
COMMENT = 0

# ######################################################################
# ## Miscellaneous options
# ######################################################################

# ## upper bound on the number of logical parameters in the model; the
# ## program will fail to operate properly if this value is set too low
SHARED_PARAMETER_SIZE = 5000

# ## Options related to general inference
# ## showing timings for inference routines
SHOW_TIMINGS = 0

# ## use straightforward calculation for FM2 matrix
SIMPLE_FM2 = 0

# ## use candidate list optimization for Viterbi parsing
CANDIDATE_LIST = 1

# ## use unrolled computation for single branch loops
FAST_SINGLE_BRANCH_LOOPS = 1

# ## use caching algorithm for fast helix length scores
FAST_HELIX_LENGTHS = 1

# ## Options related to training mode configuration
STOCHASTIC_GRADIENT = 0

# ## (A) Options related to max-margin training

# ## the maximum loss DELTA(y,y') allocated to each training example; if
# ## this symbol is undefined, then the DELTA(y,y') loss function is not
# ## included
# ##    -- for a straight CRF, this value should be undefined
# ##    -- for a max-margin model, this value should be set to 1
# HAMMING_LOSS=1

# ## multiplier used in the iterative convex-concave procedure (CCCP) for
# ## improving the solution of a max-margin model
# ##    -- for a regular max-margin model, this value should be set to 0
# ##    -- for a standard nonconvex model, this value should be set to 1
NONCONVEX_MULTIPLIER = 0.0

# ## number of steps of CCCP; there is no need to change this to 1 in the
# ## case that NONCONVEX_MULTIPLIER == 0, as the code will detect this
# ## and abort after the first CCCP iteration by default
NUM_CCCP_STEPS = 5

# ## use smooth approximation of max-margin algorithm for inference
# ## during training
SMOOTH_MAX_MARGIN = 0

# ######################################################################
# ## (B) Regularization type
# ######################################################################

SINGLE_HYPERPARAMETER = 1
MULTIPLE_HYPERPARAMETERS = 0
ARD_HYPERPARAMETERS = 0

######################################################################
# (C) Options related to regularization hyperparameter estimation
######################################################################

# Three possible modes:
# -- holdout cross-validation via grid search
# -- holdout cross-validation via gradient-based optimization
##    -- majorization-minimization

HYPERPARAMETER_GRID_SEARCH = 1
HYPERPARAMETER_GRADIENT_OPTIMIZATION = 0
HYPERPARAMETER_MAJORIZATION_MINIMIZATION = 0

######################################################################
# (C1) Grid-search options
######################################################################

# use logloss instead of regular holdout loss for holdout cross-validation
CROSS_VALIDATE_USING_LOGLOSS = 1

######################################################################
# (C2) Gradient-based optimization options
######################################################################

# starting regularization parameter
INITIAL_LOG_C = 5.0

######################################################################
# (C3) majorization-minimization-options
######################################################################

# number of iterative relinearization steps if using
# majorization-minmimization algorithm
NUM_ITERATIVE_RELINEARIZATION_STEPS = 5

# smoothing used for majorization-minimization algorithm
MM_SMOOTHING = 1.0

######################################################################
# (D) Input type
######################################################################

PROFILE = 0

######################################################################
# (E) Used parameter groups
######################################################################

# *
# define PARAMS_BASE_PAIR                           1
# define PARAMS_BASE_PAIR_DIST                      1
# define PARAMS_TERMINAL_MISMATCH                   1
# define PARAMS_HAIRPIN_LENGTH                      1
# define PARAMS_HAIRPIN_3_NUCLEOTIDES               1
# define PARAMS_HAIRPIN_4_NUCLEOTIDES               1
# define PARAMS_HELIX_LENGTH                        1
# define PARAMS_ISOLATED_BASE_PAIR                  1
# define PARAMS_INTERNAL_EXPLICIT                   1
# define PARAMS_BULGE_LENGTH                        1
# define PARAMS_INTERNAL_LENGTH                     1
# define PARAMS_INTERNAL_SYMMETRY                   1
# define PARAMS_INTERNAL_ASYMMETRY                  1
# define PARAMS_BULGE_0x1_NUCLEOTIDES               1
# define PARAMS_BULGE_0x2_NUCLEOTIDES               1
# define PARAMS_BULGE_0x3_NUCLEOTIDES               1
# define PARAMS_INTERNAL_1x1_NUCLEOTIDES            1
# define PARAMS_INTERNAL_1x2_NUCLEOTIDES            1
# define PARAMS_INTERNAL_2x2_NUCLEOTIDES            1
# define PARAMS_HELIX_STACKING                      1
# define PARAMS_HELIX_CLOSING                       1
# define PARAMS_MULTI_LENGTH                        1
# define PARAMS_DANGLE                              1
# define PARAMS_EXTERNAL_LENGTH                     1
# *

PARAMS_BASE_PAIR = 1
PARAMS_BASE_PAIR_DIST = 0
PARAMS_TERMINAL_MISMATCH = 1
PARAMS_HAIRPIN_LENGTH = 1
PARAMS_HAIRPIN_3_NUCLEOTIDES = 0
PARAMS_HAIRPIN_4_NUCLEOTIDES = 0
PARAMS_HELIX_LENGTH = 0
PARAMS_ISOLATED_BASE_PAIR = 0
PARAMS_INTERNAL_EXPLICIT = 1
PARAMS_BULGE_LENGTH = 1
PARAMS_INTERNAL_LENGTH = 1
PARAMS_INTERNAL_SYMMETRY = 1
PARAMS_INTERNAL_ASYMMETRY = 1
PARAMS_BULGE_0x1_NUCLEOTIDES = 1
PARAMS_BULGE_0x2_NUCLEOTIDES = 0
PARAMS_BULGE_0x3_NUCLEOTIDES = 0
PARAMS_INTERNAL_1x1_NUCLEOTIDES = 1
PARAMS_INTERNAL_1x2_NUCLEOTIDES = 0
PARAMS_INTERNAL_2x2_NUCLEOTIDES = 0
PARAMS_HELIX_STACKING = 1
PARAMS_HELIX_CLOSING = 1
PARAMS_MULTI_LENGTH = 1
PARAMS_DANGLE = 1
PARAMS_EXTERNAL_LENGTH = 1
PARAMS_EVIDENCE = 1

# *
# define PARAMS_BASE_PAIR                           1
# define PARAMS_BASE_PAIR_DIST                      0
# define PARAMS_TERMINAL_MISMATCH                   0
# define PARAMS_HAIRPIN_LENGTH                      0
# define PARAMS_HAIRPIN_3_NUCLEOTIDES               0
# define PARAMS_HAIRPIN_4_NUCLEOTIDES               0
# define PARAMS_HELIX_LENGTH                        0
# define PARAMS_ISOLATED_BASE_PAIR                  0
# define PARAMS_INTERNAL_EXPLICIT                   0
# define PARAMS_BULGE_LENGTH                        0
# define PARAMS_INTERNAL_LENGTH                     0
# define PARAMS_INTERNAL_SYMMETRY                   0
# define PARAMS_INTERNAL_ASYMMETRY                  0
# define PARAMS_BULGE_0x1_NUCLEOTIDES               0
# define PARAMS_BULGE_0x2_NUCLEOTIDES               0
# define PARAMS_BULGE_0x3_NUCLEOTIDES               0
# define PARAMS_INTERNAL_1x1_NUCLEOTIDES            0
# define PARAMS_INTERNAL_1x2_NUCLEOTIDES            0
# define PARAMS_INTERNAL_2x2_NUCLEOTIDES            0
# define PARAMS_HELIX_STACKING                      0
# define PARAMS_HELIX_CLOSING                       0
# define PARAMS_MULTI_LENGTH                        0
# define PARAMS_DANGLE                              0
# define PARAMS_EXTERNAL_LENGTH                     0
# *

######################################################################
# (F) Miscellaneous model constants
######################################################################

C_MIN_HAIRPIN_LENGTH = 0
C_MAX_SINGLE_LENGTH = 30
D_MAX_HAIRPIN_LENGTH = 30
D_MAX_BP_DIST_THRESHOLDS = 10
D_MAX_BULGE_LENGTH = 30
D_MAX_INTERNAL_LENGTH = 30
D_MAX_INTERNAL_SYMMETRIC_LENGTH = 15
D_MAX_INTERNAL_ASYMMETRY = 28
D_MAX_INTERNAL_EXPLICIT_LENGTH = 4
D_MAX_HELIX_LENGTH = 30
BP_DIST_LAST_THRESHOLD = 132
BP_DIST_THRESHOLDS = [3, 9, 12, 16, 21, 26, 34, 47, 71, BP_DIST_LAST_THRESHOLD]
alphabet = "ACGU"  # allowed symbols -- all other letters ignored
M = 4  # number of alphabet symbols
MAX_DIMENSIONS = 4
# debugging statements in local functions (does not currently appear)
VERBOSE_LOCAL = False

#####################################################################
# (G) BMRM stuff
######################################################################

# define BMRM_AVAILABLE                              0
# define DAIFLETCHER

# endif
