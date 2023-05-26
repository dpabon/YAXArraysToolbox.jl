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
