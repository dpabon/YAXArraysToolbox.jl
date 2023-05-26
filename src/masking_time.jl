using YAXArrays, Statistics, Zarr, NetCDF

#=
axlist = [
    RangeAxis("time", range(1, 20, length = 20)),
    RangeAxis("x", range(1, 10, length = 10)),
    RangeAxis("y", range(1, 5, length = 15)),
    CategoricalAxis("Variable", ["var1", "var2"]),
]

data = rand(20, 10, 15, 2)


props = Dict(
    "time" => "days",
    "x" => "lon",
    "y" => "lat",
    "var1" => "one of your variables",
    "var2" => "your second variable",
)

ds = YAXArray(axlist, data, props)

cube_in = ds
=#

function masking_time_comp_1(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(>=(val), cube_mask)) > 0
        if length(findall(>=(val), cube_mask)) > 1
            cube_out[findall(>=(val), cube_mask)] .= NaN
        else
            cube_out[findall(>=(val), cube_mask)] = NaN
        end

    end
end


function masking_time_comp_2(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(>(val), cube_mask)) > 0
        if length(findall(>(val), cube_mask)) > 1
            cube_out[findall(>(val), cube_mask)] .= NaN
        else
            cube_out[findall(>(val), cube_mask)] = NaN
        end


    end
end


function masking_time_comp_3(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(<=(val), cube_mask)) > 0
        if length(findall(<=(val), cube_mask)) > 1
            cube_out[findall(<=(val), cube_mask)] .= NaN
        else
            cube_out[findall(<=(val), cube_mask)] = NaN
        end

    end
end

function masking_time_comp_4(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(<(val), cube_mask)) > 0
        if length(findall(<(val), cube_mask)) > 1
            cube_out[findall(<(val), cube_mask)] .= NaN
        else
            cube_out[findall(<(val), cube_mask)] = NaN
        end
    end
end


function masking_time_comp_5(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(>=(val), cube_mask)) > 0
        if length(findall(>=(val), cube_mask)) > 1
            cube_out[findall(>=(val), cube_mask)] .= NaN
        else
            cube_out[findall(>=(val), cube_mask)] = NaN
        end
    end
end

function masking_time_comp_6(cube_out, cube_in, cube_mask; val = val)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if length(findall(!=(val), cube_mask)) > 0
        if length(findall(!=(val), cube_mask)) > 1
            cube_out[findall(!=(val), cube_mask)] .= NaN
        else
            cube_out[findall(!=(val), cube_mask)] = NaN
        end
    end
end

function masking_time_comp_1_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(>=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(>=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(>=(quantile(filter(!isnan, cube_out), p)), cube_out)] .=
                    NaN
            else
                cube_out[findall(>=(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end

function masking_time_comp_2_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(>(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(>(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(>(quantile(filter(!isnan, cube_out), p)), cube_out)] .= NaN
            else
                cube_out[findall(>(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end


function masking_time_comp_3_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(<=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(<=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(<=(quantile(filter(!isnan, cube_out), p)), cube_out)] .=
                    NaN
            else
                cube_out[findall(<=(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end

function masking_time_comp_4_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(<(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(<(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(<(quantile(filter(!isnan, cube_out), p)), cube_out)] .= NaN
            else
                cube_out[findall(<(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end

function masking_time_comp_5_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(==(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(==(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(==(quantile(filter(!isnan, cube_out), p)), cube_out)] .=
                    NaN
            else
                cube_out[findall(==(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end

function masking_time_comp_6_quant(cube_out, cube_in; p = p)
    cube_out .= cube_in
    if any(ismissing, cube_out)
        replace!(cube_out, missing => NaN)
    end
    if !all(isnan, cube_out)
        if length(findall(!=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 0
            if length(findall(!=(quantile(filter(!isnan, cube_out), p)), cube_out)) > 1

                cube_out[findall(!=(quantile(filter(!isnan, cube_out), p)), cube_out)] .=
                    NaN
            else
                cube_out[findall(!=(quantile(filter(!isnan, cube_out), p)), cube_out)] = NaN
            end
        end
    end
end

"""

# Masking using time dimension.

The function implements two methods:
1. Masking based on a threshold value for one of the variables presented in the cube. e.g., masking the values of all the variables presented in the cube where radiation is lower than X.
2. Masking based on the quantile threshold, where the quantile is estimated using the time series for each one of the variables presented in the cube.

The masked vales are set as ```NaN```.

## Arguments:

- ```cube_in``` YAXArray Cube.
- ```time_axis```: String. Name of the time axis.
- ```var_axis```: String. Name of the axis containing the variables.
- ```var_mask```: String or nothing. Name of the variable to be used to mask the other variables. If String ```val``` must be an Int64 or Float64 number. If nothing, ```val``` must be nothing and ```p``` must be a Float64 in the interval [0,1].
- ```val```: Float64 or nothing. The value of the threshold in ```var_mask``` to be used to mask all the variables in the cube. If ```var_mask = nothing``` then, ```val=nothing```
- ```p```: Float64 or nothing. Quantile value used as a threshold to mask the variables.
- ```comp```: String. Standard comparison operation between the threshold value and each one of the elements. ```comp``` Must be one of the following: "==", "!=" "<", "<=", ">", ">=".
- ```showprog```: Boolean. Progress Bar.
- ```max_cache```: String. Maximum cache to read the data. It must be in MB e.g. "100MB" or in GB "10GB".

## Examples

"""
function masking_time(
    cube_in;
    time_axis = "time",
    var_axis = "Variable",
    var_mask = "var1",
    val = 0.1,
    p = nothing,
    comp = ">",
    showprog = true,
    max_cache = "100MB",
)

    if typeof(cube_in) == Dataset
        error(
            "$cube_in is a Dataset, use Cube() function to transform a Dataset into a cube.",
        )
    end

    if last(max_cache, 2) == "MB"
        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-6)

    elseif last(max_cache, 2) == "GB"

        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-9)

    else
        error("only MB or GB values are accepted for max_cache")
    end

    if typeof(var_mask) == String && typeof(p) == Nothing && typeof(val) != Nothing

        # In this case one variable is used to mask the others. e.g. Mask the values of all the variables presented in the cube where radiation is lower than X.

        kwarg = (; Symbol(var_axis) => var_mask)

        cube_mask = getindex(cube_in; kwarg...)

        indims = (InDims(time_axis), InDims(time_axis))

        outdims = OutDims()

        if comp == ">="
            return mapCube(
                masking_time_comp_1,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == ">"
            return mapCube(
                masking_time_comp_2,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "<="
            return mapCube(
                masking_time_comp_3,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "<"
            return mapCube(
                masking_time_comp_4,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "=="
            return mapCube(
                masking_time_comp_5,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )
        elseif comp == "!="
            return mapCube(
                masking_time_comp_6,
                (cube_in, cube_mask),
                indims = indims,
                outdims = outdims;
                val = val,
                showprog = showprog,
                max_cache = max_cache,
            )
        else
            error("incorrect 'comp' value. 'comp' must be '>=', '>', '<=', '<', '=='")
        end

    elseif typeof(var_mask) == Nothing && typeof(p) != Nothing && typeof(val) == Nothing

        # in this case the mask is applied based on the quantile value (p) of the time series for each one of the variables in the cube.

        indims = InDims(time_axis)

        outdims = OutDims()

        if comp == ">="
            return mapCube(
                masking_time_comp_1_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == ">"
            return mapCube(
                masking_time_comp_2_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "<="
            return mapCube(
                masking_time_comp_3_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "<"
            return mapCube(
                masking_time_comp_4_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "=="
            return mapCube(
                masking_time_comp_5_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif comp == "!="
            return mapCube(
                masking_time_comp_6_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        else
            error("incorrect 'comp' value. 'comp' can be '>=', '>', '<=', '<', '=='")
        end

    else

    end
end
