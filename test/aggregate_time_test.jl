using YAXArrays, Zarr, DimensionalData

esds = open_dataset(
    "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v2.1.1/esdc-8d-0.25deg-184x90x90-2.1.1.zarr",
)
esdc = Cube(esds)

# Estimating the monthly LAI

lai_month = aggregate_time(
    esdc[Variable=At("leaf_area_index")];
    time_axis = :Ti,
    new_resolution = "month",
    new_time_step = 1,
    fun = "mean",
    p = nothing,
    skipMissing = true,
    skipnan = true,
    showprog = true,
    max_cache = "1GB",
)
