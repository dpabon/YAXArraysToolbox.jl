metric = ["median", "mean", "std", "var", "sum", "quant", "min", "max"]


cube_in = open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-1x720x1440-2.1.1.zarr",
)

cube_in = Cube(cube_in)
cube_in.Variable
cube_in = cube_in[
    lon = (-9.0, 0.0),
    lat = (35, 40),
    time = (Date(2010), Date(2014)),
    Variable = ["leaf_area_index", "sensible_heat"],
]

for i in eachindex(metric)
    println(metric[i])
    plot_time(
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

plot_time(
    cube_in;
    time_axis = "time",
    var_axis = "Variable",
    lon_axis = "lon",
    lat_axis = "lat",
    var = "sensible_heat",
    fun = "median",
    p = 0.2,
    showprog = true,
    max_cache = "100MB",
)

plot_time(
    cube_in;
    time_axis = "time",
    var_axis = "Variable",
    lon_axis = "lon",
    lat_axis = "lat",
    var = nothing,
    fun = "median",
    resolution = (900, 600),
    p = 0.2,
    showprog = true,
    max_cache = "100MB",
    ncol = 2
)

test_out
