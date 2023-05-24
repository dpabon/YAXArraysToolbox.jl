module YAXArraysToolbox

import Dates
import StatsBase,
    Distributions,
    Test,
    LinearAlgebra,
    GLM,
    NamedArrays,
    Combinatorics,
    YAXArrays,
    YAXArrayBase,
    TimeSeries


###### Basic functions ######

# aggregate time
include("aggregate_time.jl")

export aggregate_time

# Masking
include("masking_altitude.jl") # TO IMPROVE function names!!
include("masking_general.jl") # Maybe TO MERGE with masking_altitude.

export altitude_mask_results_proc, altitude_masking_proc
export masking_proc
# filling time

include("filling_time.jl")

# plot time
include("plot_time.jl")
export plot_time

# Space4time

include("space4time.jl")

export space4time_proc




end
