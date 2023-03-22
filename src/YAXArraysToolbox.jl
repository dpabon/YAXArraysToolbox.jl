module YAXArraysToolbox

import Dates
import StatsBase, Distributions, Test, LinearAlgebra, GLM, NamedArrays, Combinatorics, YAXArrays, YAXArrayBase
# Statistics by index

include("statistics_by_index.jl") # to add max, min, mean, and sd.
include("dates_builder.jl")

export median_by_index
export dates_builder
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
