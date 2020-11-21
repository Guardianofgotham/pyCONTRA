######################################################################
# DistributedComputation.hpp
##
# This is a class for performing distributed computation via MPI.  In
# general, suppose you have an indexed family of functions f_i(x) for
# i belonging to some set S.  This class will allow you to run f_i(x)
# for each i in S efficiently by distributing the work among a
# collection of processors.  The result returned will be
##
# sum_i f_i(x)
##
# One particular assumption made by the DistributedComputation class
# is that the results returned will *always* be a vector of doubles.
# This assumption is made to ensure that each of the "reductions" is
# efficient.
##
# Examples:
##
# (1) Computation of the gradient for a batch mode machine learning
# algorithm.  In this case, the set S can represent the set of
# training examples, each i corresponds to a single training
# example, and x is a parameter set, and f_i(x) computes the
# gradient of some error function for the ith training example at
# the parameters x.
##
# (2) Processing a collection of files in some way.  Here, each i in
# S is some file descriptor, and x contains any shared data that
# must be provided to do the processing.  Note that it may be the
# case that the return value of each f_i is irrelevant, and that
# the main purpose of the distributed computation is simply to
# "run" each f_i(x) to ensure that all files are processed.
##
# In order to use this class, you must
##
# (1) Create a class for storing shared data which is to be
# "broadcast" to all processors.  This class should not contain
# any pointers, as the data in this class will be distributed to
# other processors by a direct memcpy().  In the following, we
# will assume that the name of this class is SharedData.
##
# (2) Create a class for storing nonshared data which is to be passed
# to each processor individually.  This class should not contain
# any pointers, as the data in this class will be distributed to
# other processors by a direct memcpy().  In the following, we
# will assume that the name of this class is NonSharedData.
##
# (3) Create a subclass of DistributedComputation which implements
# the method DoComputation().  Essentially, the DoComputation()
# method implements the functionality for a single f_i(x).
# Observe that the DoComputation() method takes essentially three
# arguments,
##
# std::vector<double> &result
# const SharedData &shared_data
# const NonSharedData &nonshared_data
##
# The latter two arguments correspond to the "x" and "i"
# arguments of each functon f_i(x).  The first argument is where
# the results of the computation (if any) should be stored.
##
# (5) Call the RunAsComputeNode() from all compute nodes and
# the DistributeComputation() method from the master node to
# perform a distributed computation.  The arguments to the
# DistributeComputation() method are:
##
# std::vector<double> &result
# const Shared Data &shared_data
# const std::vector<const NonSharedData> &nonshared_data>
##
# The DistributedComputation class will take care of the details
# to ensure that the data is shuffled to the compute nodes in an
# efficient manner.  Note that the work units associated with
# each of the nonshared data entries is allocated in the order
# supplied in the nonshared_data[] vector.  In general, for
# efficiency, it makes sense to sort the entries of
# nonshared_data[] in order of decreasing expected time to
# completion.
##
# (5) Call the StopComputeNodes() routine from the master node to
# ensure that the compute nodes stop running.
##
# That's it!
######################################################################

class DistributedComputationBase:
    def __init__(self):
        pass

    def DoComputation(self, result, shared_data, nonshared_data):
        pass
    
    ## start and stop compute nodes
    def RunAsComputeNode(self):
        pass

    def StopComputeNodes(self):
        pass


    ## perform distributed computation (to be called by master node)
    def DistributeComputation(self, result, shared_data, nonshared_data):
        pass

    ## some simple routines for dealing with node IDs
    def IsComputeNode(self):
        pass
    def IsMasterNode(self):
        pass
    def GetNumNodes(self):
        pass
    def GetNodeID(self):
        pass

    ## query statistics regarding the efficiency of distributed computation
    def GetEfficiency(self):
        pass
    def ResetEfficiency(self):
        pass