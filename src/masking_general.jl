using YAXArrays, StatsBase, DimensionalData


function masking_stats(cube_out, cube_in_to_mask, cube_summary_stats; rsquared_thr)
    cube_out .= NaN
    for i in eachindex(cube_in_to_mask)

        if cube_summary_stats[i] < rsquared_thr

            cube_out[i] = NaN

        else
            cube_out[i] = cube_in_to_mask[i]

        end

    end
end

function masking_co_occurrence(
    cube_out,
    cube_in_to_mask,
    cube_co_occurence;
    co_occurence_thr,
)

    cube_out .= NaN

    for i in eachindex(cube_in_to_mask)

        if cube_co_occurence[i] < co_occurence_thr

            cube_out[i] = NaN

        else

            cube_out[i] = cube_in_to_mask[i]

        end

    end

end


function masking_delta(cube_out, cube_in_to_mask, cube_delta; minmax)

    cube_out .= cube_in_to_mask

    if typeof(minmax) == Tuple{Float64,Float64} || typeof(minmax) == Tuple{Int64,Int64}
        if any(!isnan, cube_in_to_mask)

            for i in eachindex(cube_in_to_mask)

                if cube_delta[i] < minmax[1] || cube_delta[i] > minmax[2]

                    cube_out[i] = NaN
                end
            end
        end

    elseif typeof(minmax) == Tuple{Float64,Nothing} ||
           typeof(minmax) == Tuple{Int64,Nothing}
        if any(!isnan, cube_in_to_mask)

            for i in eachindex(cube_in_to_mask)

                if cube_delta[i] < minmax[1]

                    cube_out[i] = NaN
                end
            end
        end

    elseif typeof(minmax) == Tuple{Nothing,Float64} ||
           typeof(minmax) == Tuple{Nothing,Int64}
        if any(!isnan, cube_in_to_mask)

            for i in eachindex(cube_in_to_mask)

                if cube_delta[i] > minmax[2]

                    cube_out[i] = NaN
                end
            end
        end
    end
end




@doc raw"""

# Masking processor

## Arguments:

- ```cube_in_to_mask```: YAXArray cube to be masked.

- ```cube_rsquare```: Nothing, or YAXArray cube with the ``R^{2}`` variable. If set to ```nothing``` no mask is applied
- ```rsquare_thr```: Float64. ``R^{2}`` threshold. All values lower than ```rsquare_thr``` are set to ```NaN```

- ```cube_co_occurrence```: Nothing, or YAXArray cube with the co-occurrence variable. If set to ```nothing``` no mask is applied.
- ```co_occurence_thr```: Float64. Co-occurence threshold. All values lower than ```co_occurence_thr``` are set to ```NaN```

- ```cube_delta```: Nothing, or YAXArray cube with delta variable. If set to ```nothing``` no mask is applied.
- ```minmax_delta```: Tuple. Minimum and maximum thresholds of delta variable. Values lower and higher than the thresholds are set to ```NaN```. It is also possible to set any of the thresholds as ```nothing``` e.g. ```(-1, nothing)``` or ```(nothing, 1)``` in these cases only one threshold is applied.

- ``` time_dim```: Nothing, or String. Name of the time dimension. This dimensions needs to be present in all the cubes. If set to ```nothing``` no time dimension considered (It can result in slower computation time!). By default ```time_dim = time```

- ```showprog```: Boolean. Show progress bar. By default ```showprog = true```

## Output:

- YAXArray cube masked.

"""
function masking_proc(
    cube_in_to_mask::YAXArray;
    cube_rsquared = nothing,
    rsquared_thr = nothing,
    cube_co_occurrence = nothing,
    co_occurence_thr = nothing,
    cube_delta = nothing,
    minmax_delta = nothing,
    time_dim = :Ti,
    showprog = true,
)

    if !isnothing(time_dim)

        indims = (InDims(time_dim), InDims(time_dim))
        outdims = OutDims(Dim{time_dim}(lookup(cube_in_to_mask, time_dim).data))

    else
        indims = (InDims(), InDims())
        outdims = OutDims()
    end

    if !isnothing(cube_rsquared)
        results_int = mapCube(
            masking_stats,
            (cube_in_to_mask, cube_rsquared),
            indims = indims,
            outdims = outdims,
            showprog = showprog;
            rsquared_thr = rsquared_thr,
        )
    end

    if !isnothing(cube_co_occurrence)
        if @isdefined results
            results_int = mapCube(
                masking_co_occurrence,
                (results_int, cube_co_occurrence),
                indims = indims,
                outdims = outdims,
                showprog = showprog;
                co_occurence_thr = co_occurence_thr,
            )
        else
            results_int = mapCube(
                masking_co_occurrence,
                (cube_in_to_mask, cube_co_occurrence),
                indims = indims,
                outdims = outdims,
                showprog = showprog;
                co_occurence_thr = co_occurence_thr,
            )
        end
    end

    if !isnothing(cube_delta)
        if @isdefined results
            results_int = mapCube(
                masking_delta,
                (results_int, cube_delta),
                indims = indims,
                outdims = outdims,
                showprog = showprog;
                minmax = minmax_delta,
            )
        else
            results_int = mapCube(
                masking_delta,
                (cube_in_to_mask, cube_delta),
                indims = indims,
                outdims = outdims,
                showprog = showprog;
                minmax = minmax_delta,
            )
        end
    end


    return results_int
end
