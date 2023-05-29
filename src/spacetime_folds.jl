using DataFrames, CSV
using MLUtils
using Random

#=
Original code from:
CAST R package
Version: 0.8.0
Authors@R: c(person("Hanna", "Meyer", email = "hanna.meyer@uni-muenster.de", role = c("cre", "aut")),
             person("Carles", "Milà", role = c("aut")),
             person("Marvin", "Ludwig", role = c("aut")),
             person("Jan", "Linnenbrink", role = c("aut")),
             person("Philipp", "Otto", role = c("ctb")),
             person("Chris", "Reudenbach", role = c("ctb")),
             person("Thomas", "Nauss", role = c("ctb")),
             person("Edzer", "Pebesma", role = c("ctb")))
Author: Hanna Meyer [cre, aut], Carles Milà [aut], Marvin Ludwig [aut], Jan Linnenbrink [aut], Philipp Otto [ctb], Chris Reudenbach [ctb], Thomas Nauss [ctb], Edzer Pebesma [ctb]
Maintainer: Hanna Meyer <hanna.meyer@uni-muenster.de>
Description: Supporting functionality to run 'caret' with spatial or spatial-temporal data. 'caret' is a frequently used package for model training and prediction using machine learning. CAST includes functions to improve spatial or spatial-temporal modelling tasks using 'caret'. It includes the newly suggested 'Nearest neighbor distance matching' cross-validation to estimate the performance of spatial prediction models and allows for spatial variable selection to selects suitable predictor variables in view to their contribution to the spatial model performance. CAST further includes functionality to estimate the (spatial) area of applicability of prediction models. Methods are described in Meyer et al. (2018) <doi:10.1016/j.envsoft.2017.12.001>; Meyer et al. (2019) <doi:10.1016/j.ecolmodel.2019.108815>; Meyer and Pebesma (2021) <doi:10.1111/2041-210X.13650>; Milà et al. (2022) <doi:10.1111/2041-210X.13851>; Meyer and Pebesma (2022) <doi:10.1038/s41467-022-29838-9>.
License: GPL (>= 2)
URL: https://github.com/HannaMeyer/CAST,
    https://hannameyer.github.io/CAST/


=#

"""

# Create Space-time Folds

Create spatial, temporal or spatio-temporal Folds for cross validation based on pre-defined groups.

## Arguments:

- ```x``` DataFrame containing spatio-temporal data.
- ```spacevar```: String. which column of x identifies the spatial units (e.g. ID of weather stations).
- ```timevar```: String. which column of x identifies the temporal units (e.g. the day of the year).
- ```k```: Int64. Number of folds. If spacevar or timevar is nothing and a leave one location out or leave one time step out cv should be performed, set k to the number of unique spatial or temporal units.
- ```class```: String. which column of x identifies a class unit (e.g. land cover) NOT IMPLEMENTED YET!!.
- ```seed```: Int64 or Float64, See ?Random.seed!().

## Return

```cv_indices_train, cv_indices_test = spacetime_folds(x;spacevar="var1", timevar="var2", k=10, class=nothing, seed=23)```

## References

Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., Nauß, T. (2018): Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software 101: 1-9.
"""
function spacetime_folds(
    x;
    spacevar=nothing,
    timevar=nothing,
    k=10,
    class=nothing,
    seed=23,
)
### if classification is used, make sure that classes are equally distributed across folds
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

    return cvindices_train, cvindices_test
end

