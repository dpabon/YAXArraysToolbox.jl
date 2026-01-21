

function collapse_space_mean(cube_out, cube_in)

    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= mean(filter(!isnan, skipmissing(cube_in)))
    end

end


function collapse_space_median(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= median(filter(!isnan, skipmissing(cube_in)))
    end

end

function collapse_space_std(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= std(filter(!isnan, skipmissing(cube_in)))
    end
end

function collapse_space_var(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= var(filter(!isnan, skipmissing(cube_in)))
    end
end

function collapse_space_sum(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= sum(filter(!isnan, skipmissing(cube_in)))
    end
end

function collapse_space_quant(cube_out, cube_in; p = p)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= quantile(filter(!isnan, skipmissing(cube_in)), p)
    end
end

function collapse_space_min(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= minimum(filter(!isnan, skipmissing(cube_in)))
    end
end

function collapse_space_max(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in) || !all(ismissing, cube_in)
        cube_out .= maximum(filter(!isnan, skipmissing(cube_in)))
    end
end


"""

# Plot time

The function allow to plot the time series of a given variables in a cube or all the variables present in a cube. As is expected that cubes contain spatial dimensions the spatial dimensions are collapsed using a function e.g., estimating the mean of the variable using the pixels of a certain area for each time step.
    
    ## Arguments:
    
    - ```cube_in``` YAXArray Cube.
    - ```time_axis```: String. Name of the time axis.
    - ```var_axis```: String. Name of the axis containing the variables.
    - ```var```: String or nothing. Name of the variable to be plotted. If nothing all the variables presented in the cube are plotted.
    - ```lat_axis```: String. Name of the latitude axis.
    - ```lon_axis```: String. Name of the longitute axis.
    - ```fun```: String. Name of the function used to collapse the spatial dimensions. It must be "median", "mean", "std", "var", "sum", "quant", "min", or "max".
    - ```plot_type```: String. Name of the plot type. By default: "lines". It can also be "scatter".
    - ```p```: Float64. in the interval [0,1]. If ```fun=quant``` p is the value of the quantile. 
    - ```resolution```: Tuple. Plot resolution. By default ```resolution = (600, 400)```. 
    - ```ncol```: Number of plots by column. By default ```ncol = 1```.
    - ```nrow```: Number of plots by row. By default ```ncol = 1```.
    - ```showprog```: Boolean. Progress Bar.
    - ```max_cache```: String. Maximum cache to read the data. It must be in MB e.g. "100MB" or in GB "10GB".
    
    
    ## Examples
    
    ```julia
    using YAXArrays, Zarr, CairoMakie, GeoMakie, Statistics, DimensionalData

    metric = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]


    cube_in = open_dataset(
        "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr",
    )

    cube_in = Cube(cube_in)
    cube_in.Variable

    cube_in = cube_in[
        lon = (-9.0 .. 0.0),
        lat = (35 .. 40),
        Ti = (Date(2010) .. Date(2014)),
        Variable = At(["leaf_area_index", "sensible_heat"]),
    ]



    plot_time(
        cube_in;
        time_axis = :Ti,
        var_axis = :Variable,
        lon_axis = :lon,
        lat_axis = :lat,
        var = "sensible_heat",
        fun = "median",
        p = 0.2,
        showprog = true,
        max_cache = "100MB",
    )

    plot_time(
        cube_in;
        time_axis = :Ti,
        var_axis = :Variable,
        lon_axis = :lon,
        lat_axis = :lat,
        var = nothing,
        fun = "median",
        resolution = (900, 600),
        p = 0.2,
        showprog = true,
        max_cache = "100MB",
        ncol = 2,
    )

    for i in eachindex(metric)
        println(metric[i])
        plot_time(
            cube_in;
            time_axis = :Ti,
            var_axis = :Variable,
            lon_axis = :lon,
            lat_axis = :lat,
            var = "sensible_heat",
            fun = metric[i],
            p = 0.2,
            showprog = true,
            max_cache = "100MB",
        )
    end

    ```
    """
    function plot_time(
        cube_in;
        time_axis = :Ti,
        var_axis = :Variable,
        var = nothing,
        lat_axis = :lat,
        lon_axis = :lon,
        fun = "mean",
        plot_type = "lines",
        p = nothing,
        resolution = (600, 400),
        ncol = 1,
        nrow = 1,
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
        
        if typeof(var) != Nothing
            kwarg = (; Symbol(var_axis) => At(var))
            
            cube_in = getindex(cube_in; kwarg...)
            
            indims = InDims(lat_axis), InDims(lon_axis)
            outdims = OutDims()
            
            if fun == "mean"
                
                temp_cube = mapCube(
                collapse_space_mean,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "median"

            temp_cube = mapCube(
                collapse_space_median,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "std"

            temp_cube = mapCube(
                collapse_space_std,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "var"

            temp_cube = mapCube(
                collapse_space_var,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "sum"

            temp_cube = mapCube(
                collapse_space_sum,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "quant"

            temp_cube = mapCube(
                collapse_space_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "min"

            temp_cube = mapCube(
                collapse_space_min,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "max"

            temp_cube = mapCube(
                collapse_space_max,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
                )
                
            end
            dates = lookup(temp_cube, time_axis).data
            vals = temp_cube.data[:]
            
            ta = TimeArray(dates, rand(length(dates)))
            
            tempo = string.(timestamp(ta))
            lentime = length(tempo)
            slice_dates = range(1, lentime, step = lentime ÷ 8)
            
            if plot_type == "lines"
                fig = Figure(resolution = resolution)
                ax = Axis(fig[1, 1], xlabel = "Date", ylabel = fun, title = var)
                line1 = lines!(ax, 1:lentime, vals; color = :black, linewidth = 0.85)
                ax.xticks = (slice_dates, tempo[slice_dates])
                ax.xticklabelrotation = π / 4
                ax.xticklabelalign = (:right, :center)
                return fig
                
            elseif plot_type == "scatter"
                
                fig = Figure(resolution = resolution)
                ax = Axis(fig[1, 1], xlabel = "Date", ylabel = fun, title = var)
                line1 = scatter!(ax, 1:lentime, vals; color = :black)
                ax.xticks = (slice_dates, tempo[slice_dates])
                ax.xticklabelrotation = π / 4
                ax.xticklabelalign = (:right, :center)
                return fig
                
            else
                error("wrong value for plot_type")
                
            end
            
            
        else
            
            indims = InDims(lat_axis, lon_axis)
            outdims = OutDims()
            
            if fun == "mean"
                
                temp_cube = mapCube(
                collapse_space_mean,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "median"

            temp_cube = mapCube(
                collapse_space_median,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "std"

            temp_cube = mapCube(
                collapse_space_std,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "var"

            temp_cube = mapCube(
                collapse_space_var,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "sum"

            temp_cube = mapCube(
                collapse_space_sum,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "quant"

            temp_cube = mapCube(
                collapse_space_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "min"

            temp_cube = mapCube(
                collapse_space_min,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "max"

            temp_cube = mapCube(
                collapse_space_max,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
                )
                
            end
            
            dates = lookup(temp_cube, time_axis).data
            
            ta = TimeArray(dates, rand(length(dates)))
            
            tempo = string.(timestamp(ta))
            lentime = length(tempo)
            slice_dates = range(1, lentime, step = lentime ÷ 8)
            variables_loc = lookup(cube_in, var_axis).data
            
            
            
            
            if length(variables_loc) > 6
                
                @warn "There are more than 6 variables to be plotted. Consider using nrow and ncol accordingly."
                
            end
            
            fig = Figure(resolution = resolution)
            
            if plot_type == "lines"
                
                if ncol == 1 && (nrow * ncol) < length(variables_loc)
                    
                    nrow = length(variables_loc)
                    
                    
                    for j = 1:nrow
                        
                        ax =
                        Axis(fig[j, 1], xlabel = "Date", ylabel = fun, title = variables_loc[j])
    
                    kwarg = (; Symbol(var_axis) => variables_loc[j])
    
                    temp_cube2 = getindex(temp_cube; kwarg...)
                    vals = temp_cube2.data
    
                    line1 = lines!(ax, 1:lentime, vals; color = :black, linewidth = 0.85)
                    ax.xticks = (slice_dates, tempo[slice_dates])
                    ax.xticklabelrotation = π / 4
                    ax.xticklabelalign = (:right, :center)
    
                end
                return fig
    
    
            elseif ncol != 1 && (nrow * ncol) < length(variables_loc)
                error(
                    "Number of rows and columns is less than the number of variables to be plotted.",
                    )
                else
                    
                    if plot_type == "lines"
                        
                        counter = 1
                        
                        for i = 1:ncol
                            for j = 1:nrow
                                
                                ax = Axis(
                                fig[j, i],
                                xlabel = "Date",
                                ylabel = fun,
                                title = variables_loc[counter],
                                )
                                
                                kwarg = (; Symbol(var_axis) => At(variables_loc[counter]))
                                
                                temp_cube2 = getindex(temp_cube; kwarg...)
                                
                                vals = temp_cube2.data
                                
                                line1 = lines!(ax, 1:lentime, vals; color = :black, linewidth = 0.85)
                                ax.xticks = (slice_dates, tempo[slice_dates])
                                ax.xticklabelrotation = π / 4
                                ax.xticklabelalign = (:right, :center)
                                
                                counter += 1
                                if counter > length(variables_loc)
                                    break
                                end
                                
                            end
                            
                        end
                        return fig
                        
                    elseif plot_type == "scatter"
                        
                        counter = 1
                        
                        for i = 1:ncol
                            for j = 1:nrow
                                
                                ax = Axis(
                                fig[j, i],
                                xlabel = "Date",
                                ylabel = fun,
                                title = variables_loc[counter],
                                )
                                
                                kwarg = (; Symbol(var_axis) => variables_loc[counter])
                                
                                temp_cube2 = getindex(temp_cube; kwarg...)
                                
                                vals = temp_cube2.data
                                
                                line1 = scatter!(ax, 1:lentime, vals; color = :black)
                                ax.xticks = (slice_dates, tempo[slice_dates])
                                ax.xticklabelrotation = π / 4
                                ax.xticklabelalign = (:right, :center)
                                
                                counter += 1
                                if counter > length(variables_loc)
                                    break
                                end
                                
                            end
                            
                        end
                        return fig
                        
                    else
                        error("wrong value for plot_type")
                        
                    end 
                    
                end
                
            end
        end
    end