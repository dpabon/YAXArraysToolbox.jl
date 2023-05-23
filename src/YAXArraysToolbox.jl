module YAXArraysToolbox

import Dates
import StatsBase, Distributions, Test, LinearAlgebra, GLM, NamedArrays, Combinatorics, YAXArrays, YAXArrayBase
# aggregate time
include("aggregate_time.jl")

export aggregate_time

# Space4time

include("space4time.jl")

export space4time_proc

# Masking
include("masking_altitude.jl") # TO IMPROVE function names!!
include("masking_general.jl") # Maybe TO MERGE with masking_altitude.

export altitude_mask_results_proc, altitude_masking_proc
export masking_proc
# Creating time dimension on input cube based on another cube time series TO IMPROVE!!

include("filling_time.jl")



end
