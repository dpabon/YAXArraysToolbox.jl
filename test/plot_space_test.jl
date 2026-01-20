using YAXArraysToolbox
using CairoMakie
using Statistics
using GeoMakie
using YAXArrays
using DimensionalData


cube_in = open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr",
)

cube_in = Cube(cube_in)


cube_in = cube_in[
    lon=(-9.0..0.0),
    lat=(35..40),
    Ti=(Date(2010)..Date(2014)),
    Variable=At(["leaf_area_index", "sensible_heat"]),
]

plot_space(
    cube_in;
    time_axis = :Ti,
    resolution = (900, 500),
    xticklabel_pad = 25,
    yticklabel_pad = 25,
    var_axis = :Variable,
    var = "leaf_area_index",
    fun = "median",
)


metric = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]


for i in eachindex(metric)
    println(metric[i])
    plot_space(
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



plot_space(
    cube_in;
    time_axis = :Ti,
    var_axis = :Variable,
    lon_axis = :lon,
    lat_axis = :lat,
    var = nothing,
    fun = "median",
    resolution = (1200, 300),
    p = 0.2,
    showprog = true,
    max_cache = "100MB",
    ncol = 2,
)
