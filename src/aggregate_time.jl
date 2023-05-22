using Dates, YAXArrays, Zarr, Statistics

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



function median_by_index(xout, xin; index_list=time_to_index, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = median(skipmissing(xin[index_list[i]]))

                else
                    xout[i] = median(xin[index_list[i]])
                end
                if skipnan == true
                    xout[i] = median(filter(!isnan, xin[index_list[i]]))
                else
                    xout[i] = median(xin[index_list[i]])

                end
            end
        end
    end
end

function mean_by_index(xout, xin; index_list=time_to_index, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin) || !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = mean(skipmissing(xin[index_list[i]]))
                else
                    xout[i] = mean(xin[index_list[i]])
                end
                if skipnan == true
                    xout[i] = mean(filter(!isnan, xin[index_list[i]]))
                else
                    xout[i] = mean(xin[index_list[i]])
                end
            end
        end
    end
end

function std_by_index(xout, xin; index_list=time_to_index, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin) || !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = std(skipmissing(xin[index_list[i]]))

                else
                    xout[i] = std(xin[index_list[i]])
                end
                if skipnan == true
                    xout[i] = std(filter(!isnan, xin[index_list[i]]))
                else
                    xout[i] = std(xin[index_list[i]])
                end
            end
        end
    end
end

function var_by_index(xout, xin; index_list=time_to_index, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin) || !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = var(skipmissing(xin[index_list[i]]))

                else
                    xout[i] = var(xin[index_list[i]])
                end
                if skipnan == true
                    xout[i] = var(filter(!isnan, xin[index_list[i]]))
                else
                    xout[i] = var(xin[index_list[i]])
                end
            end
        end
    end
end

function sum_by_index(xout, xin; index_list=time_to_index, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin) || !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = sum(skipmissing(xin[index_list[i]]))

                else
                    xout[i] = sum(xin[index_list[i]])
                end
                if skipnan == true
                    xout[i] = sum(filter(!isnan, xin[index_list[i]]))
                else
                    xout[i] = sum(xin[index_list[i]])
                end
            end
        end
    end
end


function quant_by_index(xout, xin; index_list=time_to_index, p=p, skipMissing=skipMissing, skipnan=skipnan)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin) || !all(ismissing, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                if skipMissing == true
                    xout[i] = quantile(skipmissing(xin[index_list[i]]), p=p)

                else
                    xout[i] = quantile(xin[index_list[i]], p=p)
                end
                if skipnan == true
                    xout[i] = quantile(filter(!isnan, xin[index_list[i]]), p=p)
                else
                    xout[i] = quantile(xin[index_list[i]], p=p)
                end
            end
        end
    end
end


@doc raw"""

# Aggregate by time

## Arguments:

- ```cube_in``` YAXArray Cube.
- ```time_axis_name```: String. Name of the time axis
- ```new_resolution```: String. New temporal resolution can be ```"day"```, ```"month"```, ```"year"```.
- ```new_time_step```: Int64. Time step to be computed in the new time series. e.g. ```new_resolution="day", new_time_step=8``` will compute the function each 8 days. The new time dimenssion will only contain the days corresponding to the 8th day.
- ```fun```: String. Function to be applied to aggregate the time. It can be "median", "mean", "std", "var", "sum", "quant".
- ```p```: Float64 on the interval [0,1]. If ```fun=quant``` p is the value of the quantile. 
- ```skipMissing```: Boolean. Skip missing values when aggregating the data.
- ```skipnan```: Boolean. Skip NaN values when aggregating the data
- ```showprog```: Boolean. Progress Bar.
- ```max_cache```: String. Maximum cache to read the data. It needs to be expressed in MB e.g. "100MB" or in GB "10GB"

## Examples

"""
function aggregate_time_proc(cube_in; time_axis="time", new_resolution="month", new_time_step=1, fun="median", p=Nothing, skipMissing=true, skipnan=true, showprog=true, max_cache="100MB")

    if last(max_cache, 2) == "MB"
        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-6)

    elseif last(max_cache, 2) == "GB"

        max_cache = parse(Float64, max_cache[begin:end-2]) / (10^-9)

    else
        error("only MB or GB values are accepted for max_cache")
    end


    if typeof(cube_in) == Dataset
        error("$cube_in is a Dataset, use Cube() function to transform a Dataset into a cube.")
    end

    time_org = getAxis(time_axis, cube_in).values

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
        error("$new_time_step is not a multiple of $temp, try defining a different new_time_step.")
    end

    if new_time_step > 1

        # OutDims definition

        outdims = OutDims(RangeAxis(time_axis, new_dates_axis[new_time_step:new_time_step:end]))

        # indicies in the cube to be aggregated
        index_in_cube_temp = [findall(==(i), time_index) for i in new_dates]
        count = 1
        index_in_cube = Vector{Int64}[]
        for i in 1:(Int(length(index_in_cube_temp) / new_time_step))
            #println(count)
            push!(index_in_cube, reduce(vcat, index_in_cube_temp[count:(count+new_time_step-1)]))
            count += new_time_step
        end
    else
        # outdims definition

        outdims = OutDims(RangeAxis(time_axis, new_dates_axis))

        # indices in the cube to be aggregated
        index_in_cube = [findall(==(i), time_index) for i in new_dates]
    end


    # indims definition
    indims = InDims(time_axis)

    if fun == "median"

        return mapCube(median_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)

    elseif fun == "mean"
        return mapCube(mean_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)

    elseif fun == "std"
        return mapCube(std_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)

    elseif fun == "var"

        return mapCube(var_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)

    elseif fun == "sum"

        return mapCube(sum_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)

    elseif fun == "quant"

        return mapCube(quant_by_index, cube_in, indims=indims, outdims=outdims; index_list=index_in_cube, p=p, skipMissing=skipMissing, skipnan=skipnan, showprog=showprog, max_cache=max_cache)
    end


end

