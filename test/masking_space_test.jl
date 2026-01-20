using YAXArrays, Zarr, DimensionalData, Test

@testset "masking_space" begin
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

    data[3, 5, 1] = NaN

    data[1, 10, 1] = NaN


    data[9, 5, 1] = NaN

    ds_mask = YAXArray(axlist, data)



    test_cube = masking_space(ds, ds_mask; lat_axis = :x, lon_axis = :y)

    @test all(isnan, test_cube[3, 5, :, :].data)
    @test all(isnan, test_cube[1, 10, :, :].data)
    @test all(isnan, test_cube[9, 5, :, :].data)


    data = convert(Array{Union{Float64,Missing},3}, data)

    data[3, 5, 1] = missing

    data[1, 10, 1] = missing


    data[9, 5, 1] = missing

    ds_mask = YAXArray(axlist, data)



    test_cube = masking_space(ds, ds_mask; lat_axis = :x, lon_axis = :y, val_mask = missing)

    @test all(isnan, test_cube[3, 5, :, :].data)
    @test all(isnan, test_cube[1, 10, :, :].data)
    @test all(isnan, test_cube[9, 5, :, :].data)

end
