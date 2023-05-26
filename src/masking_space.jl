using YAXArrays, Zarr

#=

# Masking based on the spatial-dimensions.

One method:

1. Masking based in a static layer YAXArray without time dimension but with the same spatial-extension.

Where values that are NaN or missing in the static layer are set to NaN in the datacube

=#

function masking_space_int(cube_out, cube_in; mask)

    cube_out .= cube_in
    if length(mask) > 0
        cube_out[mask] .= NaN

    end
end



"""

# Masking using spatial dimension


The masked vales are set as ```NaN```!!.

## Arguments:

- ```cube_in``` YAXArray Cube to be masked.
-```mask``` YAXArray Cube without time dimension and with a single variable to be used as mask. All values equal to NaN or missing will be masked in cube_in. The mask will be applied to all the variables and time steps presented in ```cube_in```.
- ```lat_axis```: String. Name of the latitude axis.
- ```lon_axis```: String. Name of the longitude axis.
- ```val_mask```: NaN or missing. Value present in ```mask``` to be used as reference to mask ```cube_in```. Must be NaN or missing.
- ```showprog```: Boolean. Progress Bar.
- ```max_cache```: String. Maximum cache to read the data. It must be in MB e.g. "100MB" or in GB "10GB".

## Examples

```julia
using YAXArrays, Zarr


axlist = [
RangeAxis("time", range(1, 20, length = 20)),
RangeAxis("x", range(1, 10, length = 10)),
RangeAxis("y", range(1, 5, length = 15)),
CategoricalAxis("Variable", ["var1", "var2"]),
]


data = rand(20, 10, 15, 2)


ds = YAXArray(axlist, data, props)

axlist = [
RangeAxis("x", range(1, 10, length = 10)),
RangeAxis("y", range(1, 5, length = 15)),
CategoricalAxis("Variable", ["var1"]),
]


data = rand(10, 15, 1)

ds_mask = YAXArray(axlist, data)

masking_space(ds, ds_mask; lat_axis = "x", lon_axis = "y")
```

"""
function masking_space(
    cube_in,
    mask;
    lat_axis,
    lon_axis,
    val_mask = NaN,
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

    indims = InDims(lon_axis, lat_axis)

    outdims = OutDims(
        RangeAxis(lon_axis, getAxis(lon_axis, cube_in).values),
        RangeAxis(lat_axis, getAxis(lat_axis, cube_in).values),
    )

    temp_mask = mask.data

    if isnan(val_mask)

        to_mask = findall(isnan, temp_mask)

    elseif ismissing(val_mask)
        to_mask = findall(ismissing, temp_mask)

    else
        error("$val_mask is not a valid mask value only missing or NaN is possible.")
    end

    if length(to_mask) == 0
        @warn "The mask cube does not contain NaN or missing values. No mask is applied."
    end

    return mapCube(
        masking_space_int,
        cube_in,
        indims = indims,
        outdims = outdims;
        mask = to_mask,
        showprog = showprog,
        max_cache = max_cache,
    )
end
