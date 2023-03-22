module YAXArraysToolbox

import Dates
import StatsBase, Distributions, Test, LinearAlgebra, GLM, NamedArrays, Combinatorics, YAXArrays, YAXArrayBase
# Statistics by index

include("statistics_by_index.jl") # to add max, min, mean, and sd.
include("dates_builder.jl")

# Space4time

include("space4time.jl")


# Masking
include("masking_altitude.jl") # TO IMPROVE function names!!
include("masking_general.jl") # Maybe TO MERGE with masking_altitude.

# Creating time dimension on input cube based on another cube time series TO IMPROVE!!

include("filling_time.jl")

end
