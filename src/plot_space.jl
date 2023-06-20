using YAXArrays, Zarr, CairoMakie, GeoMakie, Statistics


function collapse_time_mean(cube_out, cube_in)

    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= mean(filter(!isnan, skipmissing(cube_in)))
        end
    end

end


function collapse_time_median(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= median(filter(!isnan, skipmissing(cube_in)))
        end
    end

end

function collapse_time_std(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= std(filter(!isnan, skipmissing(cube_in)))
        end
    end
end

function collapse_time_var(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= var(filter(!isnan, skipmissing(cube_in)))
        end
    end
end

function collapse_time_sum(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= sum(filter(!isnan, skipmissing(cube_in)))
        end
    end
end

function collapse_time_quant(cube_out, cube_in; p = p)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= quantile(filter(!isnan, skipmissing(cube_in)), p)
        end
    end
end

function collapse_time_min(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= minimum(filter(!isnan, skipmissing(cube_in)))
        end
    end
end

function collapse_time_max(cube_out, cube_in)
    cube_out .= NaN
    if !all(isnan, cube_in)
        if !all(ismissing, cube_in)
            cube_out .= maximum(filter(!isnan, skipmissing(cube_in)))
        end
    end
end


"""
# Plot Space/Maps


## Arguments 

- ```cube_in```: YAXArray Cube.
- ```time_axis```: String. Name of the time axis.
- ```var_axis```: String. Name of the axis containing the variables.
- ```var```: String or nothing. Name of the variable to be plotted. If nothing all the variables presented in the cube are plotted.
- ```lat_axis```: String. Name of the latitude axis.
- ```lon_axis```: String. Name of the longitute axis.
- ```fun```: String. Name of the function used to collapse the spatial dimensions. It must be "median", "mean", "std", "var", "sum", "quant", "min", or "max".
- ```p```: Float64. in the interval [0,1]. If ```fun=quant``` p is the value of the quantile.
- ```colormap```: Color Map. By default: ```colormap = Reverse(:batlow)```
- ```resolution```: Plot resolution. By default ```resolution = (800, 300)```.
- ```ncol```: Number of plots by column. By default ```ncol = 1```.
- ```nrow```: Number of plots by row. By default ```ncol = 1```.
- ```showprog```: Boolean. Progress Bar.
- ```max_cache```: String. Maximum cache to read the data. It must be in MB e.g. "100MB" or in GB "10GB".

## Examples

```julia

cube_in = open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr",
)

cube_in = Cube(cube_in)


cube_in = cube_in[
    lon = (-9.0, 0.0),
    lat = (35, 40),
    time = (Date(2010), Date(2014)),
    Variable = ["leaf_area_index", "sensible_heat"],
]

plot_space(cube_in; time_axis = "time", resolution = (900, 600), var_axis = "Variable", var =  "leaf_area_index", fun = "median")


metric = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]


for i in eachindex(metric)
    println(metric[i])
    plot_space(
        cube_in;
        time_axis = "time",
        var_axis = "Variable",
        lon_axis = "lon",
        lat_axis = "lat",
        var = "sensible_heat",
        fun = metric[i],
        p = 0.2,
        showprog = true,
        max_cache = "100MB",
    )
end



plot_space(
    cube_in;
    time_axis = "time",
    var_axis = "Variable",
    lon_axis = "lon",
    lat_axis = "lat",
    var = nothing,
    fun = "median",
    resolution = (1200, 600),
    p = 0.2,
    showprog = true,
    max_cache = "100MB",
    ncol = 2,
)

```

"""
function plot_space(
    cube_in::YAXArray;
    time_axis = "time",
    var_axis = "Variable",
    var = nothing,
    lat_axis = "lat",
    lon_axis = "lon",
    fun = "mean",
    p = nothing,
    colormap = Reverse(:batlow),
    resolution = (800, 300),
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
        kwarg = (; Symbol(var_axis) => var)

        cube_in = getindex(cube_in; kwarg...)

        indims = InDims(time_axis)
        outdims = OutDims(outtype=Float32)

        if fun == "mean"

            temp_cube = mapCube(
                collapse_time_mean,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "median"

            temp_cube = mapCube(
                collapse_time_median,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "std"

            temp_cube = mapCube(
                collapse_time_std,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "var"

            temp_cube = mapCube(
                collapse_time_var,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "sum"

            temp_cube = mapCube(
                collapse_time_sum,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "quant"

            temp_cube = mapCube(
                collapse_time_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "min"

            temp_cube = mapCube(
                collapse_time_min,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "max"

            temp_cube = mapCube(
                collapse_time_max,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        end

        lon = getAxis(lon_axis, temp_cube).values
        lat = getAxis(lat_axis, temp_cube).values

        fig = Figure(resolution = resolution)

        ga = GeoAxis(
            fig[1, 1],
            source = "+proj=longlat +datum=WGS84",
            dest = "+proj=longlat",
            coastlines = true,
            lonlims = (minimum(lon), maximum(lon)),
            latlims = (minimum(lat), maximum(lat)),
            title = var *
                    " \n " *
                    string(first(getAxis(time_axis, cube_in).values)) *
                    " / " *
                    string(last(getAxis(time_axis, cube_in).values)),
        )
        map1 = CairoMakie.heatmap!(ga, lon, lat, temp_cube[:, :], colormap = colormap)
        cbar1 = Colorbar(
            fig[1, 2],
            map1,
            label = fun,
            ticklabelsize = 18,
            labelpadding = 5,
            width = 10,
        )

        return fig

    else

        indims = InDims(time_axis)
        outdims = OutDims(outtype=Float32)

        if fun == "mean"

            return temp_cube = mapCube(
                collapse_time_mean,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "median"

            temp_cube = mapCube(
                collapse_time_median,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "std"

            temp_cube = mapCube(
                collapse_time_std,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "var"

            temp_cube = mapCube(
                collapse_time_var,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "sum"

            temp_cube = mapCube(
                collapse_time_sum,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "quant"

            temp_cube = mapCube(
                collapse_time_quant,
                cube_in,
                indims = indims,
                outdims = outdims;
                p = p,
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "min"

            temp_cube = mapCube(
                collapse_time_min,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        elseif fun == "max"

            temp_cube = mapCube(
                collapse_time_max,
                cube_in,
                indims = indims,
                outdims = outdims;
                showprog = showprog,
                max_cache = max_cache,
            )

        end


        variables_loc = getAxis(var_axis, cube_in).values


        if length(variables_loc) > 6

            @warn "There are more than 6 variables to be plotted. Consider using nrow and ncol accordingly."

        end

        fig = Figure(resolution = resolution)

        lon = getAxis(lon_axis, temp_cube).values
        lat = getAxis(lat_axis, temp_cube).values




        if ncol == 1 && (nrow * ncol) < length(variables_loc)

            nrow = length(variables_loc)

            for j = 1:nrow


                kwarg = (; Symbol(var_axis) => variables_loc[j])

                temp_cube2 = getindex(temp_cube; kwarg...)

                ga = GeoAxis(
                    fig[j, 1],
                    source = "+proj=longlat +datum=WGS84",
                    dest = "+proj=longlat",
                    coastlines = true,
                    lonlims = (minimum(lon), maximum(lon)),
                    latlims = (minimum(lat), maximum(lat)),
                    title = variables_loc[j] *
                            " \n " *
                            string(first(getAxis(time_axis, cube_in).values)) *
                            " / " *
                            string(last(getAxis(time_axis, cube_in).values)),
                )

                map1 =
                    CairoMakie.heatmap!(ga, lon, lat, temp_cube2[:, :], colormap = colormap)

                cbar1 = Colorbar(
                    fig[j, 2],
                    map1,
                    label = fun,
                    ticklabelsize = 18,
                    labelpadding = 5,
                    width = 10,
                )

            end
            return fig

        elseif ncol != 1 && (nrow * ncol) < length(variables_loc)
            error(
                "Number of rows and columns is less than the number of variables to be plotted.",
            )
        else

            init_row = 1
            init_col = 1

            for i = 1:nrow
                for j = 1:ncol
                    kwarg = (; Symbol(var_axis) => variables_loc[init_row])

                    temp_cube2 = getindex(temp_cube; kwarg...)
                    ga = GeoAxis(
                        fig[i, init_col],
                        source = "+proj=longlat +datum=WGS84",
                        dest = "+proj=longlat",
                        coastlines = true,
                        lonlims = (minimum(lon), maximum(lon)),
                        latlims = (minimum(lat), maximum(lat)),
                        title = variables_loc[init_row] *
                                " \n " *
                                string(first(getAxis(time_axis, cube_in).values)) *
                                " / " *
                                string(last(getAxis(time_axis, cube_in).values)),
                    )
                    map1 = CairoMakie.heatmap!(
                        ga,
                        lon,
                        lat,
                        temp_cube2[:, :],
                        colormap = Reverse(:batlow),
                    )
                    init_row += 1
                    init_col += 1
                    cbar1 = Colorbar(
                        fig[i, init_col],
                        map1,
                        label = fun,
                        ticklabelsize = 18,
                        labelpadding = 5,
                        width = 10,
                    )
                    init_col += 1
                    if init_row > length(variables_loc)
                        break
                    end
                end
                init_col = 1
            end
            return fig
        end

    end
end
