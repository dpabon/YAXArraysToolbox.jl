
function dates_builder_month(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2]))
    end

    return out
end

function dates_builder_day(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2], x[i][3]))
    end

    return out
end


### median ###

function median_by_index_1(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
               if !all(isnan, view(xin, index_list[i]))
                xout[i] = median(skipmissing(view(xin, index_list[i])))
               end
            end
        end
    end
end

function median_by_index_2(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                    xout[i] = median(skipmissing(filter(!isnan, view(xin, index_list[i]))))
                   end
                end
            end 
        end
    end
end

function median_by_index_3(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = median(filter(!isnan, view(xin, index_list[i])))
                   end
                end
            end
        end
    end
end

function median_by_index_4(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)

        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = median(view(xin, index_list[i]))
                   end
                    
                end
            end
        end
    end
end

###################################################################

#### mean ####

function mean_by_index_1(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
                if !all(isnan, view(xin, index_list[i]))
                    xout[i] = mean(skipmissing(view(xin, index_list[i])))
                end
            end
        end
    end
end

function mean_by_index_2(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = mean(skipmissing(filter(!isnan, view(xin, index_list[i]))))
                   end
                end
            end
        end
    end
end

function mean_by_index_3(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)

        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = mean(filter(!isnan, view(xin, index_list[i])))
                   end
                end
            end
        end
    end
end

function mean_by_index_4(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = mean(view(xin, index_list[i]))
                   end
                end
            end
        end
    end
end

###############################################################

#### std ####

function std_by_index_1(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
               if !all(isnan, view(xin, index_list[i]))
                    xout[i] = std(skipmissing(view(xin, index_list[i])))
               end
            end
        end
    end
end

function std_by_index_2(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = std(skipmissing(filter(!isnan, view(xin, index_list[i]))))
                   end
                end
            end
        end
    end
end

function std_by_index_3(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = std(filter(!isnan, view(xin, index_list[i])))
                   end
                end
            end
        end
    end
end

function std_by_index_4(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = std(view(xin, index_list[i]))
                   end
                end
            end
        end
    end
end

###################################################################################################

#### variance ####

function var_by_index_1(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
               if !all(isnan, view(xin, index_list[i]))
                    xout[i] = var(skipmissing(view(xin, index_list[i])))
               end
            end
        end
    end
end

function var_by_index_2(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = var(skipmissing(filter(!isnan, view(xin, index_list[i]))))
                   end
                end
            end 
        end
    end
end

function var_by_index_3(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = var(filter(!isnan, view(xin, index_list[i])))
                   end
                end
            end
        end
    end
end

function var_by_index_4(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = var(view(xin, index_list[i]))
                   end
                end
            end
        end
    end
end

####################################################################################################

#### sum ####

function sum_by_index_1(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
               if !all(isnan, view(xin, index_list[i]))
                xout[i] = sum(skipmissing(view(xin, index_list[i])))
               end
            end
        end
    end
end

function sum_by_index_2(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = sum(skipmissing(filter(!isnan, view(xin, index_list[i]))))
                   end
                end
            end
        end
    end
end

function sum_by_index_3(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = sum(filter(!isnan, view(xin, index_list[i])))
                   end
                end
            end 
        end
    end
end

function sum_by_index_4(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = sum(view(xin, index_list[i]))
                   end
                end
            end
        end
    end
end

##########################################################################################################

#### quantile ####

function quant_by_index_1(xout, xin; index_list = time_to_index, p = p)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(ismissing, view(xin, index_list[i]))
               if !all(isnan, view(xin, index_list[i]))
                xout[i] = quantile(skipmissing(view(xin, index_list[i])), p)
               end
            end
        end
    end
end

function quant_by_index_2(xout, xin; index_list = time_to_index, p = p)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = quantile(skipmissing(filter(!isnan, view(xin, index_list[i]))), p)
                   end
                end
            end
        end
    end
end

function quant_by_index_3(xout, xin; index_list = time_to_index, p = p)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = quantile(filter(!isnan, view(xin, index_list[i])), p)
                   end
                end
            end 
        end
    end
end

function quant_by_index_4(xout, xin; index_list = time_to_index, p = p)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        if !all(isnan, xin)
            for i in eachindex(index_list)
                if !all(ismissing, view(xin, index_list[i]))
                   if !all(isnan, view(xin, index_list[i]))
                        xout[i] = quantile(view(xin, index_list[i]), p)
                   end
                end
            end
        end
    end
end

###############################################################################################################

"""

# Aggregate by time

## Arguments:

- ```cube_in``` YAXArray Cube.
- ```time_axis```: String. Name of the time axis.
- ```new_resolution```: String. New temporal resolution can be ```"day"```, ```"month"```, ```"year"```.
- ```new_time_step```: Int64. Time step to be computed in the new time series. e.g. ```new_resolution="day", new_time_step=8``` will compute the function each 8 days. The new time dimension will only contain the days corresponding to the 8th day.
- ```fun```: String. Function to be applied to aggregate the time. It can be "median", "mean", "std", "var", "sum", "quant", "min", "max".
- ```p```: Float64 in the interval [0,1]. If ```fun=quant``` p is the value of the quantile. 
- ```skipMissing```: Boolean. Skip missing values when aggregating the data. If all values are missing, NaN is returned.
- ```skipnan```: Boolean. Skip NaN values when aggregating the data. If all values are NaN, NaN is returned.
- ```showprog```: Boolean. Progress Bar.
- ```max_cache```: String. Maximum cache to read the data. It must be in MB e.g. "100MB" or in GB "10GB".

## Examples

```julia

using YAXArrays, Zarr, DimensionalData, YAXArraysToolbox

esds = open_dataset("https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr")
esdc = Cube(esds)

# Estimating the monthly LAI

lai_month = aggregate_time(esdc[Variable = At("leaf_area_index")]; time_axis = :Ti, new_resolution = "month", new_time_step=1, fun="mean", p=nothing, skipMissing=true, skipnan=true, showprog=true, max_cache="1GB")


```

"""
function aggregate_time(
    cube_in;
    time_axis = :Ti,
    new_resolution = "month",
    new_time_step = 1,
    fun = "median",
    p = Nothing,
    skipMissing = true,
    skipnan = true,
    showprog = true,
    max_cache = "100MB",
)

    if last(max_cache, 2) == "MB"
        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-6)

    elseif last(max_cache, 2) == "GB"

        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-9)

    else
        error("only MB or GB values are accepted for max_cache")
    end


    if typeof(cube_in) == Dataset
        error(
            "$cube_in is a Dataset, use Cube() function to transform a Dataset into a cube.",
        )
    end

    time_org = collect(lookup(cube_in, time_axis))

    if new_resolution == "year"
        time_index = year.(time_org)
        new_dates = unique(time_index)
        new_dates_axis = DateTime.(new_dates)

    elseif new_resolution == "month"

        time_index = yearmonth.(time_org)
        new_dates = unique(time_index)
        new_dates_axis = dates_builder_month(new_dates)

    elseif new_resolution == "day"

        time_index = yearmonthday.(time_org)
        new_dates = unique(time_index)
        new_dates_axis = dates_builder_day(new_dates)

    end

    if length(new_dates) % new_time_step !== 0
        temp = length(new_dates)
        error(
            "$new_time_step is not a multiple of $temp, try defining a different new_time_step.",
        )
    end

    if new_time_step > 1

        # OutDims definition

        outdims =
            OutDims(Dim{Symbol("Ti")}(new_dates_axis[new_time_step:new_time_step:end]))

        # indicies in the cube to be aggregated
        index_in_cube_temp = [findall(==(i), time_index) for i in new_dates]
        count = 1
        index_in_cube = Vector{Int64}[]
        for i in 1:(Int(length(index_in_cube_temp) / new_time_step))
            #println(count)
            push!(
                index_in_cube,
                reduce(vcat, index_in_cube_temp[count:(count+new_time_step-1)]),
            )
            count += new_time_step
        end
    else
        # outdims definition

        outdims =
            OutDims(Dim{Symbol("Ti")}(new_dates_axis))

        # indices in the cube to be aggregated
        index_in_cube = [findall(==(i), time_index) for i in new_dates]
    end


    # indims definition
    indims = InDims(time_axis)

    if fun == "median"

        if skipMissing == true && skipnan == false
            return mapCube(
                median_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                median_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                median_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                median_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    elseif fun == "mean"

        if skipMissing == true && skipnan == false
            return mapCube(
                mean_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                mean_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                mean_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                mean_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    elseif fun == "std"

        if skipMissing == true && skipnan == false
            return mapCube(
                std_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                std_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                std_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                std_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end


    elseif fun == "var"

        if skipMissing == true && skipnan == false
            return mapCube(
                var_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                var_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                var_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                var_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    elseif fun == "sum"

        if skipMissing == true && skipnan == false

            return mapCube(
                sum_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                sum_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                sum_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                sum_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    elseif fun == "quant"

        if skipMissing == true && skipnan == false
            return mapCube(
                quant_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                quant_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true
            return mapCube(
                quant_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                quant_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    elseif fun == "max"

        if skipMissing == true && skipnan == false

            return mapCube(
                max_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                max_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true

            return mapCube(
                max_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                max_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end


    elseif fun == "min"

        if skipMissing == true && skipnan == false

            return mapCube(
                min_by_index_1,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == true && skipnan == true

            return mapCube(
                min_by_index_2,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == true

            return mapCube(
                min_by_index_3,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif skipMissing == false && skipnan == false

            return mapCube(
                min_by_index_4,
                cube_in,
                indims = indims,
                outdims = outdims;
                index_list = index_in_cube,
                showprog = showprog,
                max_cache = max_cache,
            )

        end

    else
        error(
            "$fun is not a valid argument for 'fun' try with: 'median', 'mean', 'std', 'var', 'sum', 'quant', 'min', or 'max'",
        )

    end

end
