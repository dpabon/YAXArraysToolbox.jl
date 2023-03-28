using YAXArraysToolbox
using Test
using NeutralLandscapes
using YAXArrays
using Zarr
using Random
using TiledViews
using Statistics
# using CairoMakie

# Generating Land Cover Map
edge = 500

size_tile = (edge, edge)

# # spatial autocorrelation ##
Random.seed!(232323)
spatial_auto = 0.7
midpoint_sim = rand(MidpointDisplacement(spatial_auto), size_tile)
# heatmap(midpoint_sim)

## Number of classes ##
n_classes = 3
classes_dist = NeutralLandscapes.classify(midpoint_sim, ones(n_classes))



# heatmap(classes_dist)

## defining lst per class ##

class_1_lst = 20.0

class_2_lst = 22.0

class_3_lst = 23.8

alllst = (class_1_lst, class_2_lst, class_3_lst)

# delta variables

delta_1_org = abs(class_1_lst - class_2_lst)
delta_2_org = abs(class_1_lst - class_3_lst)
delta_3_org = abs(class_2_lst - class_3_lst)

delta_original = [delta_1_org, delta_2_org, delta_3_org]

#empty matrix
lst = fill(NaN, size_tile)

for i in eachindex(alllst)
    lst[findall(==(i), classes_dist)] .= alllst[i]
end

## Spatial resampling ##

# New spatial resolution in meters
size_new_pixel = 10 

size_new_edge = Int(edge/size_pixel_new)

# we will use Tiled View to create a view of each one of the new pixels
a = TiledView(classes_dist, (size_new_pixel, size_new_pixel),(0, 0); keep_center = false)

# we will create an array with the new resolution to save the results
new_res_array_classes = fill(0., (size_new_edge, size_new_edge, n_classes))

for i in 1:size_new_edge
    for j in 1:size_new_edge
        for c in 1:n_classes
            new_res_array_classes[i,j,c] = count(==(c), a[:,:,i,j]) / (size_new_pixel^2)
        end
    end
end
#=
for i in 1:n_classes
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Class " * string(i))
    temp = heatmap!(new_res_array_classes[:,:,i], colormap = Reverse(:bamako))
    Colorbar(fig[1, 2], temp, label = "occurrence")
    display(fig)
end
=#

a = TiledView(lst, (size_new_pixel, size_new_pixel),(0, 0); keep_center = false)

# we will create an array with the new resolution to save the results
new_res_array_lst = fill(NaN, (size_new_edge, size_new_edge))

for i in 1:size_new_edge
    for j in 1:size_new_edge
        new_res_array_lst[i,j] = mean(a[:,:,i,j])
    end
end

# heatmap(new_res_array_lst)

## Organizing data on YAXArray objects

# Land Cover classes
 
# axis
axlist = [
    RangeAxis("x", 1:size(new_res_array_classes)[1]),
    RangeAxis("y", 1:size(new_res_array_classes)[2]),
    CategoricalAxis("Classes", ["class" * string(i) for i in 1:n_classes])]
# YAXArray Cube
lcc_cube = YAXArray(axlist, new_res_array_classes)

# savecube(lcc_cube, "data/lcc_cube_test.zarr")

# testing Neutral landscape and tiled views

test_lcc = Cube(open_dataset("data/lcc_cube_test.zarr"))

@test lcc_cube.data == test_lcc.data

# LST cube

# axis
axlist = [
    RangeAxis("x", 1:size(new_res_array_classes)[1]),
    RangeAxis("y", 1:size(new_res_array_classes)[2])]

# YAXArray Cube
lst_cube = YAXArray(axlist, new_res_array_lst)

# savecube(lst_cube, "data/lst_cube_test.zarr")

test_lst = Cube(open_dataset("data/lst_cube_test.zarr"))

@test lst_cube.data == test_lst.data

## Running Space 4 Time

results_space4time = space4time_proc(lst_cube, lcc_cube; time_axis_name = nothing, lon_axis_name = "x", lat_axis_name = "y", classes_var_name = "Classes", winsize = 5, minpxl = 25, minDiffPxlspercentage = 0., classes_vec = ["class" * string(i) for i in 1:n_classes], max_value = 1, showprog = true, max_cache = 1e8)

metrics_transitions_cube = results_space4time.metrics_for_transitions

masking_without_delta = masking_proc(results_space4time.metrics_for_transitions; 
cube_rsquared = results_space4time.SummaryStats[summary_stat = "rsquared"], rsquared_thr = 0.2, 
cube_co_occurrence = results_space4time.metrics_for_transitions[Differences = "coocurence"], co_occurence_thr = 0.5, 
cube_delta = nothing, time_dim = "time", showprog = true)

transitions = getAxis("transitions", masking_without_delta)


results_delta_space4time = [mean(filter(!isnan, masking_without_delta.data[1,i,1,:,:])) for i in eachindex(transitions)]


# Testing space4time results

@test round.(results_delta_space4time; digits = 3) == round.(delta_original; digits = 3)

