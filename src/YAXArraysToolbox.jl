module YAXArraysToolbox

import StatsBase,
    Distributions,
    Test,
    LinearAlgebra,
    GLM,
    NamedArrays,
    Combinatorics,
    YAXArrays,
    YAXArrayBase,
    TimeSeries,
    Dates,
    DataFrames,
    CSV,
    Random,
    MLUtils


include("aggregate_time.jl")

export aggregate_time

include("filling_time.jl")

export filling_time

include("plot_time.jl")
export plot_time

include("plot_space.jl")
export plot_space

include("masking_time.jl")
export masking_time

include("masking_space.jl")
export masking_space

include("spacetime_folds.jl")
export spacetime_folds

# Masking 4 space4time
include("masking_altitude.jl") # TO IMPROVE function names!!
include("masking_general.jl") # Maybe TO MERGE with masking_altitude.

export altitude_mask_results_proc, altitude_masking_proc
export masking_proc

# Space4time

include("space4time.jl")

export space4time_proc

export space4time_proc_space

end
