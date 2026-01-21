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
using YAXArrays, Zarr, DimensionalData, Test
axlist = (
    Dim{:Ti}(range(1, 20, length = 20)),
    Dim{:x}(range(1, 10, length = 10)),
    Dim{:y}(range(1, 5, length = 15)),
    Dim{:Variable}(["var1", "var2"]),
    )
    
    
    data = rand(20, 10, 15, 2)
    
    
    ds = YAXArray(axlist, data)
    
    axlist = (
    Dim{:x}(range(1, 10, length = 10)),
    Dim{:y}(range(1, 5, length = 15)),
    Dim{:Variable}(["var1"]),
    )
    
    
    data = rand(10, 15, 1)
    
    data[3,5,1] = NaN
    
    data[1,10,1] = NaN
    
    
    data[9,5,1] = NaN
    
    ds_mask = YAXArray(axlist, data)
    
    
    
    test_cube = masking_space(ds, ds_mask; lat_axis = :x, lon_axis = :y)
```
"""
function masking_space(
    cube_in,
    mask;
    lat_axis = :lat,
    lon_axis = :lon,
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
        Dim{lon_axis}(lookup(cube_in, lon_axis).data),
        Dim{lat_axis}(lookup(cube_in, lat_axis).data),
    )

    temp_mask = mask.data

    if ismissing(val_mask)

        to_mask = findall(ismissing, temp_mask)

    elseif isnan(val_mask)
        to_mask = findall(isnan, temp_mask)

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
