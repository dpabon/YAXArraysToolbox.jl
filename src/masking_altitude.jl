using StatsBase,
    Distributions,
    Test,
    LinearAlgebra,
    GLM,
    NamedArrays,
    Combinatorics,
    YAXArrays,
    YAXArrayBase,
    DimensionalData



function altitude_mask_results(cube_out, cube_in; center_coord)

    cube_out .= NaN

    if (any(!isnan, cube_in))

        # v1
        cube_out[1] = mean(filter(!isnan, vec(cube_in[2, :, :])))

        # v2
        cube_out[2] = abs(
            cube_in[1, center_coord, center_coord] -
            mean(filter(!isnan, vec(cube_in[1, :, :]))),
        )

        # v3
        cube_out[3] = abs(
            cube_in[2, center_coord, center_coord] -
            mean(filter(!isnan, vec(cube_in[2, :, :]))),
        )
    end
end


@doc raw"""
# Topographical variability processor

## Arguments:
- ```cube_in_altitude``` : Altitude YAXARRAY with two variables mean, and sd.

- ```lon_axis_name``` : String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lon"```

- ```lat_axis_name``` :  String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lat"```

- ```variable_name``` : String. Name of the Variable containing the variables "mean", and "sd".

- ```winsize```: Edge size of the moving window on pixels. By default winsize = 5. E.g. ```winsize = 5``` will produce a moving window with 5^2 pixels.

- ```showprog```: Show progress bar. By default ```showprog = true```

## Output:

The Topographical variability processor produces a YAXARRAY.cube with three Indicators:
- ```v1```: ``v_{1}=\dfrac{1}{n}\sum_{i = 1}^{n}\sigma_{h,i}``
High values of ``v_{1}`` indicate hilly terrain over the considered scale, which should be discarded from the analysis.
- ```v2```: ``v_{2}=\vert \mu_{h} - \dfrac{1}{n} \sum_{i = 1}^{n}\mu_{h,i} \vert``
``v_{2}`` indicates how different the mean elevation within the central pixel is from the average elevation in the local window.

- ```v3```: ``v_{3}= \vert \sigma_{h} - \dfrac{1}{n}\sum_{i = 1}^{n}\sigma_{h,i} \vert``
``v_{3}``: compares the central pixel's standard deviation of elevation with the standard deviation across the moving window.

## See also:

- ```altitude_masking_proc``` TO ADD LINK!!

## Bibliography:
- Duveiller, G., Hooker, J., & Cescatti, A. (2018). A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5(1), Article 1. [https://doi.org/10.1038/sdata.2018.14](https://doi.org/10.1038/sdata.2018.14)

"""
function altitude_mask_results_proc(
    cube_in_altitude;
    lon_axis_name = "lon",
    lat_axis_name = "lat",
    variable_name = "Variable",
    winsize = 5,
    showprog = true,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end

    center_coord = Int(pre_step + 1)

    indims_altitude = InDims(
        variable_name,
        MovingWindow(lon_axis_name, pre_step, after_step),
        MovingWindow(lat_axis_name, pre_step, after_step),
        window_oob_value = NaN,
    )

    outdims_altitude = OutDims(CategoricalAxis("Indicators", ["v1", "v2", "v3"]))

    result_cube = mapCube(
        altitude_mask_results,
        cube_in_altitude,
        indims = indims_altitude,
        outdims = outdims_altitude,
        showprog = showprog;
        center_coord = center_coord,
    )

    return (result_cube)

end


function altitude_masking(cube_out, cube_altitude, cube_to_mask; v1_thr, v2_thr, v3_thr)
    # println("size of cube out is ", length(cube_out))
    # println("size of cube altitude is ", length(cube_altitude))
    # println("size of cube to mask is ", length(cube_to_mask))
    cube_out .= cube_to_mask
    if any(!isnan, cube_altitude)
        if cube_altitude[1] >= v1_thr &&
           cube_altitude[2] >= v2_thr &&
           cube_altitude[3] >= v3_thr
            if length(cube_to_mask) == 1
                cube_out[] = NaN
            else
                cube_out .= NaN
            end
        end
    end
end

@doc raw"""
# Topographical masking processor

## Arguments:
- ```cube_in_to_mask```: YAXArray Cube to be masked.

- ```cube_in_altitude``` : Altitude YAXARRAY with two variables mean, and sd.

- ```lon_axis_name``` : String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lon"```

- ```lat_axis_name``` :  String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lat"```

- ```variable_name``` : String. Name of the Variable containing the variables "mean", and "sd".

- ```time_axis_name``` : String or NaN. It is strongly recommended to pass this parameter if the cube to be masked contains a time dimension, otherwise ```nothing```.

- ```winsize```: Edge size of the moving window on pixels. By default ```winsize = 5```. E.g. ```winsize = 5``` will produce a moving window with 5^2 pixels.

- ```v1_thr``` : Float. Threshold to mask values using ``v_{1}`` indicator. All values higer or equal to ```v1_thr``` are set to NaN. By default ```v1_thr = 50``` 

- ```v2_thr``` : Float. Threshold to mask values using ``v_{2}`` indicator. All values higer or equal to ```v2_thr``` are set to NaN. By detault ```v2_thr = 100```

- ```v3_thr``` : Float. Threshold to mask values using ``v_{3}`` indicator. All values higer or equal to ```v3``` are set to NaN. By default ```v3_thr = 100```

- ```showprog```: Boolean. Show progress bar. By default ```showprog = true```.

## Output:

- YAXArray Datase with two variables:
    - ```cube masked```: YAXArray Cube with same dimensions as ```cube_in_to_mask```.
    - ```masked_pixels```: YAXArray Cube with same lat, lon, dimensions as cube_masked but with a single boolean variable indicating if the pixel was masked or not.
## Topographical variability indicators

- ```v1```: ``v_{1}=\dfrac{1}{n}\sum_{i = 1}^{n}\sigma_{h,i}``
High values of ``v_{1}`` indicate hilly terrain over the considered scale, which should be discarded from the analysis.
    
- ```v2```: ``v_{2}=\vert \mu_{h} - \dfrac{1}{n} \sum_{i = 1}^{n}\mu_{h,i} \vert``
``v_{2}`` indicates how different the mean elevation within the central pixel is from the average elevation in the local window.

- ```v3```: ``v_{3}= \vert \sigma_{h} - \dfrac{1}{n}\sum_{i = 1}^{n}\sigma_{h,i} \vert``
``v_{3}``: compares the central pixel's standard deviation of elevation with the standard deviation across the moving window.
## See also

- ```altitude_mask_results_proc``` function TO ADD LINK!!

## Bibliography
- Duveiller, G., Hooker, J., & Cescatti, A. (2018). A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5(1), Article 1. [https://doi.org/10.1038/sdata.2018.14](https://doi.org/10.1038/sdata.2018.14)



"""
function altitude_masking_proc(
    cube_in_to_mask,
    cube_in_altitude;
    lon_axis_name = "lon",
    lat_axis_name = "lat",
    variable_name = "Variable",
    time_axis_name = nothing,
    winsize = 5,
    v1_thr = 50,
    v2_thr = 100,
    v3_thr = 100,
    showprog = true,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end

    center_coord = Int(pre_step + 1)

    indims_altitude = InDims(
        variable_name,
        MovingWindow(lon_axis_name, pre_step, after_step),
        MovingWindow(lat_axis_name, pre_step, after_step),
        window_oob_value = NaN,
    )

    outdims_altitude = OutDims(CategoricalAxis("Indicators", ["v1", "v2", "v3"]))

    result_cube = mapCube(
        altitude_mask_results,
        cube_in_altitude,
        indims = indims_altitude,
        outdims = outdims_altitude,
        showprog = showprog;
        center_coord = center_coord,
    )

    indims_altitude_results = InDims("Indicators")


    if isnothing(time_axis_name)
        indims_to_mask = InDims()
        outdims_to_mask = OutDims()

        result = mapCube(
            altitude_masking,
            (result_cube, cube_in_to_mask),
            indims = (indims_altitude_results, indims_to_mask),
            outdims = outdims_to_mask,
            showprog = showprog;
            v1_thr = v1_thr,
            v2_thr = v2_thr,
            v3_thr = v3_thr,
        )


    else
        indims_to_mask = InDims(time_axis_name)
        outdims_to_mask = OutDims(time_axis_name)

        result = mapCube(
            altitude_masking,
            (result_cube, cube_in_to_mask),
            indims = (indims_altitude_results, indims_to_mask),
            outdims = outdims_to_mask,
            showprog = showprog;
            v1_thr = v1_thr,
            v2_thr = v2_thr,
            v3_thr = v3_thr,
        )

    end

    indims_altitude_results = InDims("Indicators")

    masked_pixels_outdims = OutDims(CategoricalAxis("Masked_pixel", ["value"]))

    function local_fun(cube_out, cube_in; v1_thr, v2_thr, v3_thr)
        cube_out .= false
        if any(!isnan, cube_in)
            if cube_in[1] >= v1_thr && cube_in[2] >= v2_thr && cube_in[3] >= v3_thr
                cube_out[] = true
            end
        end
    end

    result2 = mapCube(
        local_fun,
        result_cube,
        indims = indims_altitude_results,
        outdims = masked_pixels_outdims,
        showprog = showprog;
        v1_thr = v1_thr,
        v2_thr = v2_thr,
        v3_thr = v3_thr,
    )

    return Dataset(; cube_masked = result, masked_pixels = result2)

end

function altitude_mask_results2(cube_out, cube_in; center_coord)

    cube_out .= NaN

    if (any(!isnan, cube_in))

        # v1
        cube_out[1] = std(filter(!isnan, vec(cube_in[:, :])))

        # v2
        cube_out[2] = abs(
            cube_in[center_coord, center_coord] -
            mean(filter(!isnan, vec(cube_in[:, :]))),
        )
    end
end


function altitude_mask_results_proc2(
    cube_in_altitude;
    lon_axis_name = "lon",
    lat_axis_name = "lat",
    variable_name = "Variable",
    winsize = 5,
    showprog = true,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end

    center_coord = Int(pre_step + 1)

    indims_altitude = InDims(
        variable_name,
        MovingWindow(lon_axis_name, pre_step, after_step),
        MovingWindow(lat_axis_name, pre_step, after_step),
        window_oob_value = NaN,
    )

    outdims_altitude = OutDims(CategoricalAxis("Indicators", ["v1", "v2"]))

    result_cube = mapCube(
        altitude_mask_results2,
        cube_in_altitude,
        indims = indims_altitude,
        outdims = outdims_altitude,
        showprog = showprog;
        center_coord = center_coord,
    )

    return (result_cube)

end