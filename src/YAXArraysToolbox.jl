module YAXArraysToolbox

import StatsBase,
    LinearAlgebra,
    GLM,
    NamedArrays,
    Combinatorics,
    YAXArrays,
    TimeSeries,
    Dates,
    CairoMakie,
    GeoMakie,
    DimensionalData

using YAXArrays: YAXArray


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


# Masking 4 space4time
include("masking_altitude.jl") # TO IMPROVE function names!!
include("masking_general.jl") # Maybe TO MERGE with masking_altitude.

export altitude_mask_results_proc, altitude_masking_proc
export masking_proc

export altitude_mask_results_proc2
export altitude_mask_results2

# Space4time

include("space4time.jl")

export space4time_proc

export space4time_proc_space

export space4time_proc_old

end
