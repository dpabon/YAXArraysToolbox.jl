using MLJ
using DataFrames, CSV
using MLUtils
using Random

cookfarm = CSV.read("/home/dpabon/Desktop/cookfarm.csv", DataFrame)

x = cookfarm

spacevar = "SOURCEID"
timevar = "Date"
k = 10
class = nothing
seed = 23

if !isnothing(class)


end

if !isnothing(spacevar)
    if k > length(unique(x[!, spacevar]))
        k = length(unique(x[!, spacevar]))
        @warn "k is higher than number of unique locations. k is set to $k."

    end
end

if !isnothing(timevar)
    if k > length(unique(x[!, timevar]))
        k = length(unique(x[!, timevar]))
        @warn " k is higher than number of unique points in time. k is set to $k"

    end
end
#split space into k folds
if !isnothing(spacevar)
    train_folds_space = Vector[]
    test_folds_space = Vector[]

    Random.seed!(seed)

    for (x_train, x_val) in kfolds(shuffleobs(1:length(unique(x[!, spacevar]))), k = k)
        push!(train_folds_space, unique(x[!, spacevar])[x_train])
        push!(test_folds_space, unique(x[!, spacevar])[x_val])
    end
end
#split time into k folds
if !isnothing(timevar)
    train_folds_time = Vector[]
    test_folds_time = Vector[]

    Random.seed!(seed)

    for (x_train, x_val) in kfolds(shuffleobs(1:length(unique(x[!, timevar]))), k = k)
        push!(train_folds_time, unique(x[!, timevar])[x_train])
        push!(test_folds_time, unique(x[!, timevar])[x_val])
    end
end

# combine space and time folds

cvindices_train = Vector[]
cvindices_test = Vector[]

for i = 1:k
    if !isnothing(timevar) && !isnothing(spacevar)
        push!(
            cvindices_train,
            intersect(
                reduce(
                    vcat,
                    [
                        findall(==(train_folds_space[i][i_l]), x[!, spacevar]) for
                        i_l in eachindex(train_folds_space[i])
                    ],
                ),
                reduce(
                    vcat,
                    [
                        findall(==(train_folds_time[i][i_l]), x[!, timevar]) for
                        i_l in eachindex(train_folds_time[i])
                    ],
                ),
            ),
        )

        push!(
            cvindices_test,
            setdiff(
                1:nrow(x),
                intersect(
                    reduce(
                        vcat,
                        [
                            findall(==(train_folds_space[i][i_l]), x[!, spacevar]) for
                            i_l in eachindex(train_folds_space[i])
                        ],
                    ),
                    reduce(
                        vcat,
                        [
                            findall(==(train_folds_time[i][i_l]), x[!, timevar]) for
                            i_l in eachindex(train_folds_time[i])
                        ],
                    ),
                ),
            ),
        )

    end

    if isnothing(timevar) && !isnothing(spacevar)

        push!(
            cvindices_train,
            reduce(
                vcat,
                [
                    findall(==(train_folds_space[i][i_l]), x[!, spacevar]) for
                    i_l in eachindex(train_folds_space[i])
                ],
            ),
        )

        push!(
            cvindices_test,
            setdiff(
                1:nrow(x),
                reduce(
                    vcat,
                    [
                        findall(==(train_folds_space[i][i_l]), x[!, spacevar]) for
                        i_l in eachindex(train_folds_space[i])
                    ],
                ),
            ),
        )
    end

    if !isnothing(timevar) && isnothing(spacevar)
        push!(
            cvindices_train,
            reduce(
                vcat,
                [
                    findall(==(train_folds_time[i][i_l]), x[!, timevar]) for
                    i_l in eachindex(train_folds_time[i])
                ],
            ),
        )

        push!(
            cvindices_test,
            setdiff(
                1:nrow(x),
                reduce(
                    vcat,
                    [
                        findall(==(train_folds_time[i][i_l]), x[!, timevar]) for
                        i_l in eachindex(train_folds_time[i])
                    ],
                ),
            ),
        )
    end

end
