
function empty_model(npred = 3, nsample =100)
    resp = GLM.LmResp(randn(nsample),similar(1:nsample,0))
    chol = GLM.cholpred([ones(nsample) rand(nsample, npred)])
    GLM.LinearModel(resp,chol)
end

function reset_model!(model;X=nothing,y=nothing)
    
    model.rr.mu .= 0.0
    
    if y !== nothing
        model.rr.y .= y
    end
    
    p = model.pp
    
    if X !== nothing
        p.X[:,2:end] .= X
        X2 = p.X
        newchol = cholesky!(Hermitian(float(X2'X2)))
        p.chol.factors .= newchol.factors
    end
    
    p.beta0 .= 0
    
    p.delbeta .= 0
    
end




function fit_with_data!(model, X, y)
    reset_model!(model,X=X, y=y)
    fit!(model)
end




function coocufun(out, q1, q2, p1, p2, out_pmindist, denom)
    #replace!(q1, NaN => missing)
    #replace!(q2, NaN => missing)
    # Where q1 and q2 are the ratio of presence of pft1 and pft2 in the moving window.
    #@show typeof(out) 
    #@show typeof(q1), typeof(q2), typeof(out)
    #println(q1, q2)
    # To write unit test ....
    
    #@show length(q1) == length(q2) && (!any(isnan.(q1))) && !any(isnan.(q2))
    
    #println("here coo")
    vecq1 = vec(q1)
    vecq2 = vec(q2)
    #p1 = range(0, 1, length = length(vecq1))
    #p2 = reverse(p1)
    if length(q1) == length(q2) && isfinite(sum(q1)) && isfinite(sum(q2))
        #@show typeof(vecq1)
        #@show typeof(vecq2)
        #@show size(vecq1) size(vecq2)
        #pmindist = minimum(([i - j for i in vecq1, j in p1]).^2 + ([i - j for i in vecq2, j in p2]).^2, dims = 1)
        minimum!(out_pmindist, ([i - j for i in vecq1, j in p1]).^2 + ([i - j for i in vecq2, j in p2]).^2)
        
        out[1] = 1 - (sum(sqrt.(out_pmindist)) / denom)
        
        #return 1 - (sum(sqrt.(pmindist)) / denom)
    else
        out[1] = 0. # 0 because assume not enough data on the moving window.
        #return 0
        
    end
end






"""
space4time(climate_cube, pfts_cube, pft_list::Vector{String}, winsize = 5, minpxl = 100, minDiffPxlspercentage = 40)

Compute the space for time analysis for a given climate variable.
    ...
# Arguments
    
    climate_cube: YAXARRAY cube with dimenssions: lon, lat, time.
    pfts_cube: YAXARRAY cube with dimenssions: pfts, lat,lon, time.
    ...
# Output
    Three output cubes are generated.
    out1: Summary statistics. YAXARRAY cube where summary_stat axis contains: 
    * "rsquare": XXXX 
    * "cumulative_variance": XXXX
    * "predicted": Prediction of Z for the central pixel with its real PFT combination.
    out2: 
    
#Examples
    
    
"""
function s4time(out_1, out_2, out_3, clim_var_cube_in, pfts_cube_in, loopvars; empty_models,
    pft_list::Vector{String}, 
   time_n::Int,
   max_value::Int, p1_static, p2_static, sigma1_glob, prederr_glob, predres_glob, minDiffPxls, tran_check, half, localcomp_fix_glob, pftsvarmat_f_glob, winsize = 5, transitions_n, pftstrans_comb_names, nc, out_pmindist_global, denom, minpxl)
    #println(size(clim_var_cube_in))
    #println(size(pfts_cube_in))
   #println(size(out_3))
   sigma1 = sigma1_glob[Threads.threadid()]
   prederr = prederr_glob[Threads.threadid()]
   predres = predres_glob[Threads.threadid()]
   localcomp_fix = localcomp_fix_glob[Threads.threadid()]
   pftsvarmat_f = pftsvarmat_f_glob[Threads.threadid()]
   my_empty_models = empty_models[Threads.threadid()]
   out_pmindist = out_pmindist_global[Threads.threadid()]
    #println(Threads.threadid())
    #println(loopvars)

   if max_value == 100
       pfts_cube_in = pfts_cube_in ./ 100
   end



   #@show typeof(pfts_cube_in), typeof(clim_var_cube_in)    
   #println(size(pfts_cube_in))
   #println(size(clim_var_cube_in))

   #@show typeof(out_1) typeof(out_2) typeof(out_3)
   #println(size(out_1), size(out_2), size(out_3))
   out_1 .= NaN
   out_2 .= NaN
   out_3 .= NaN

   #replace!(pfts_cube_in, missing => NaN)
   #replace!(clim_var_cube_in, missing => NaN)

   # be sure that values are really NaN
   #=
   There are two possible cases:
   1. When pfts_cube is time dimension less, or when
   time is present as an axis.!!! The first case is as implemented in the original code
   
   the second one here
   =#

   #clim_var_cube_2 = permutedims(clim_var_cube_in, (3,2,1))
   
   climvarmat = reshape(clim_var_cube_in, ((winsize^2), time_n))

   
   local_pft1 = [findall(pft_list .== pftstrans_comb_names[comp][1]) for comp in 1:transitions_n]
   local_pft1 = reduce(vcat, local_pft1)
   local_pft2 = [findall(pft_list .== pftstrans_comb_names[comp][2]) for comp in 1:transitions_n]
   local_pft2 = reduce(vcat, local_pft2)

   #println("before if")
   #println(out_3[1,1,1])

   #println(local_pft1)
   #println(local_pft2)
   pfts_cube_in_1 = replace!(pfts_cube_in, NaN => 0.)

   if time_n == 1

        pfts_cube_in_1 = permutedims(reshape(pfts_cube_in_1, (winsize, winsize, nc, 1)), (4,1,2,3))
        
   end

   
   for it in 1:time_n
       # co-ocurrence estimation
       #println(it)
       #println(transitions_n)
       for comp in eachindex(local_pft1)
        #println(comp)
        #println(size(pfts_cube_in))
        #println("I'm here")
        #println(pfts_cube_in[it,:,:,local_pft1[comp]])
        #println(pfts_cube_in[it,:,:,local_pft2[comp]])
           # define the pfts to be processed
           out_3[it,comp,3] = coocufun([0.], pfts_cube_in_1[it,:,:,local_pft1[comp]], pfts_cube_in_1[it,:,:,local_pft2[comp]], p1_static, p2_static,out_pmindist, denom)
           #println("all good")
       end       

       

       pfts_cube_in_2 = pfts_cube_in_1[it,:,:,:]
       #println(all(isnan, pfts_cube_in_2))

       #println(size(pfts_cube_in_2))
       
       pftsvarmat = reshape(pfts_cube_in_2, (winsize^2, nc))
       
       # for debug from R -------
       
       # pftsvarmat = Matrix(CSV.read("/Net/Groups/BGI/people/dpabon/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
       # pftsvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
       # climvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_temperature_example_from_R.csv", DataFrame))
       # check that there are not NaN values on pfts and at least one pft is present
       # @show sum(pftsvarmat), sum(pftsvarmat)
       #println(pftsvarmat)
       
       if isfinite(sum(pftsvarmat)) && sum(pftsvarmat) > 0.
           #println("test")
           #println(any(isnan.(pftsvarmat)))
           # check if pftsvarmat is 0 to 1 or 0 to 100
           #println(maximum(vec(pftsvarmat)))
           
           # make sure compositions are really precisely right.
           
           #localcomp_fix = mapslices(x->1-sum(x), pftsvarmat, dims = 2)
           localcomp_fix = map(x->1-sum(x), eachslice(pftsvarmat, dims = 1))
           #map!(x->1-sum(x), localcomp_fix, eachslice(pftsvarmat, dims = 1))
           #println(size(localcomp_fix))
           #println(localcomp_fix)
           pftsvarmat_f = [pftsvarmat localcomp_fix]
           
           map!((x) -> round(x, digits = 4), pftsvarmat_f, pftsvarmat_f)
           
           # some PFTs might not be present in the 5*5 window
           # these must be identified and removed, as they cannot be predicted
           
           #pftpres_check = vec(mapslices(sum, pftsvarmat, dims = 1) .> 0)
           pftpres_check = vec(sum(pftsvarmat_f, dims =1).> 0)

           pftpres_check[nc+1] = 0

           # println(pftpres_check)
           # @show typeof(pftpres_check)
           # pftpos = pft_list[pftpres_check[1:length(pft_list)]]
           
           # check that at least XX percent of the pixels is different
           
           #uniquepixels_char = mapslices(x->string(x), pftsvarmat, dims = 2)
           #uniquepixels_char = string(eachrow(pftsvarmat))
           uniquepixels_char = unique(eachslice(pftsvarmat_f, dims = 1))

           uniquepixels = length(uniquepixels_char)
           
           # Sometimes there is only 1 PFT in all 25 gridcells,
           # making the problem 0-dimensional
           
           if uniquepixels > minDiffPxls && sum(pftpres_check) > 1
               # println("test")
               # avoid divided by 0
               #lc1 = mapslices(x -p1_static, p2_static> x ./ (sum(x) + 0.000001), pftsvarmat, dims=2)
               lc1 = map(x -> x / (sum(x) + 0.000001), eachslice(pftsvarmat_f, dims = 1))
               lc1 = reduce(vcat, lc1')
               # centre the columns (to be in the centre wrt new space)
               
               #lc2 = mapslices(x -> x .- mean(x), lc1, dims = 1)
               lc2 = map(x -> x .- mean(x), eachslice(lc1, dims = 2))
               lc2 = reduce(hcat, lc2)
               # remember col means for the subsequent predictions
               
               #lcm = mapslices(mean, lc2, dims = 1)
               lCm = mean(lc2, dims = 1)

               #@show size(lc2)
               # decompose the resulting table
               # println("before svd")
               lcsvd = svd(lc2)
               # println("after svd")
               
               # related to "enough PFTs", is there enough variability between observations?
               # if all obs have exactly the same composition, the regression is not possible
               # so only do the regression if there is some variability...
               # println(sum(lcsvd.S))
               if sum(lcsvd.S) > 0
                   # n. of dimmensions that explain 100 % of the variance
                   
                   # when there are only two pfts are in the matrix 
                   # cumsum(lcsvd.S) / sum(lcsvd.S.^2) sometimes can be lower than 1 in that case
                   
                   #println(cumsum(lcsvd.S) / sum(lcsvd.S.^2))
                   
                   temp = round.(cumsum(lcsvd.S.^2) ./ sum(lcsvd.S.^2), digits = 8)
                   ndim = minimum(findall(temp.>= 1))
                   
                   # store results to output object
                   # cumulative variance
                   out_1[it,2] = sum(lcsvd.S)
                   
                   #println(out7_cumulated_variance)
                   
                   # dimmensions that explain 100 % of the variance
                   
                   lr = lc2 * lcsvd.V[:,1:ndim]
                   

                   # create bogus composition dataset
                   
                   #boguscomp = zeros(Float64, nc+1, nc+1)
                   #boguscomp[diagind(boguscomp)] .= 1

                   boguscomp = I(nc+1)
                   # println(boguscomp)
                   
                   # remove absent pfts from bogus predictor compositions and close compositions.
                   #
                   bogusc1 = mapslices(x -> x/sum(x), boguscomp[pftpres_check, pftpres_check], dims = 2)'
                   
                   #println(size(bogusc1))
                   #println(pftpres_check)
                   
                   # center the columns as the training data were centered
                   
                   bogusc2 = (I(sum(pftpres_check)) .- lCm[pftpres_check])
                   
                   bogusc2 = (bogusc1' .- lCm[pftpres_check])

                   bogusc3 = bogusc2 * lcsvd.V[pftpres_check,1:ndim]
                   #println(climvarmat[:,it])
                   
                   if sum(!isnan,climvarmat[:,it]) >= minpxl
                       #println("test")
                       # data = hcat(DataFrame(lt = convert(Vector{Float64}, climvarmat[:,it])), DataFrame(lr, :auto))
                       
                       # compreg = GLM.lm(Term(:lt) ~ sum(Term.(Symbol.(names(data[:, Not(:lt)])))), data)

                       #println("before fail")
                        try
                            fit_with_data!(my_empty_models[ndim], lr, climvarmat[:,it])
                             #my_empty_models[ndim] = GLM.lm(Float32.(lr), convert(Array{Float32}, climvarmat[:,it]))

                        catch e
                            println("error at", loopvars)
                            error()
                            
                        end
                       
                       
                       # continue only if there are no NA in the estimated coefficients
                       
                       coef_reg = GLM.coef(my_empty_models[ndim])
                       #println("original coef $coef_reg")
                       #println("second estimation coef $(lr\climvarmat[:,it])")
                       
                       if isfinite(sum(coef_reg))
                           # then do predictions for the log-normal approach
                           if isa(bogusc2, Vector)
                               
                               boguspred = predict(my_empty_models[ndim], [ones(length(bogusc3)) bogusc3])
                               # boguspred = GLM.predict(compreg, DataFrame( x1 = bogusc3))
                               
                           else
                               # boguspred = GLM.predict(compreg, DataFrame(bogusc3, :auto))
                               boguspred = predict(my_empty_models[ndim], [ones(size(bogusc3, 1)) bogusc3])
                               
                           end
                           
                           
                           x2pred = [ones(size(bogusc3,1), 1) bogusc3]
                           
                           vcv = GLM.vcov(my_empty_models[ndim])
                           # vcv = GLM.vcov(compreg)
                           
                           sigma = x2pred * vcv * x2pred'
                           
                           # now store the target variables
                           # but make sure appropiate temperatures go back to appropiate
                           # pfts (as absent pfts were removed)
                           
                           # value of climatevar for pure ptfs
                           
                           predres .= NaN
                           #println(boguspred)
                           predres[view(pftpres_check,1:nc)] = boguspred
                           #println(it)
                           
                           
                           out_2[it,:,1] = predres
                           
                           
                           prederr .= NaN
                           
                           prederr[view(pftpres_check,1:nc)] = sqrt.(diag(sigma))

                               
                           out_2[it,:,2] = prederr

                           
                           # prediction of varclim for the central pixel with its real pft combination
                           
                           out_1[it,3] = StatsModels.predict(my_empty_models[ndim])[half]
                           
                           # Rsquare of the regression
                           
                           out_1[it,1] = StatsModels.r2(my_empty_models[ndim])
                           # println(out_1)
                           # println(r2(compreg))
                           
                           # and now for the transitions
                           # only the PFTs identified in the pftlist are to be used
                           
                           sigma1 .= NaN
                           sigma1[view(pftpres_check,1:nc), view(pftpres_check,1:nc)] = sigma
                           
                           # calculate the difference on climatevar caused by going from one pft to another.
                           
                           diff_clim_pft_pred = round.((predres .- predres'), digits = 10)
                           
                           diff_clim_pft_pred = diff_clim_pft_pred[tran_check]
                           
                           # propagate the error (as variances) taking into account
                           # the covariance terms
                           # in the original implementation diffclim_pft_var = dZvar
                           
                           diff_clim_pft_pred_var = round.((diag(sigma1) .+ diag(sigma1)') .- 2 * sigma1, digits = 10)
                           
                           #print(diff_clim_pft_pred)
                           # flag out those with zero error (may occur with identical compositions for 2 pfts)
                           
                           diff_clim_pft_pred_var = diff_clim_pft_pred_var[tran_check]
                           
                           #println(diff_clim_pft_pred)
                           
                           # flag out those with zero error (may occur with identical compositions for 2 pfts)
                           
                           if any(round.(diff_clim_pft_pred, digits = 8) .== 0.)
                               diff_clim_pft_pred[round.(diff_clim_pft_pred, digits = 8) .== 0] .= NaN
                           end
                           
                           if any(round.(diff_clim_pft_pred_var, digits = 8) .== 0.)
                               diff_clim_pft_pred_var[round.(diff_clim_pft_pred_var; digits = 8) .== 0] .= NaN
                           end
                           
                           # mask out low co-ocurrence ask to greg!!! This can be performed  masking the pixels
                           # before all estimations
                           
                           #println("inside loop")
                           #println(diff_clim_pft_pred)

                           if length(diff_clim_pft_pred) == 1

                               out_3[it,1,1] = diff_clim_pft_pred[1]
                               #println(out1_delta)
                               if diff_clim_pft_pred_var[1] .< 0
                                   out_3[it,1,2] = NaN
                               else
                                   out_3[it,1,2] = sqrt.(diff_clim_pft_pred_var)[1]

                               end
                               

                           else

                               out_3[it,1:transitions_n,1] = diff_clim_pft_pred
                               #println(out1_delta)
                               #if diff_clim_pft_pred_var .< 0
                                   #out_3[it,1:transitions_n,2] .= NaN

                               #else
                               diff_clim_pft_pred_var[diff_clim_pft_pred_var .< 0] .=  NaN
                                out_3[it,1:transitions_n,2] = sqrt.(diff_clim_pft_pred_var)

                               #end
                           end
                       end    
                   end
               end 
           end
       end
       #println(it)
   end
end


 
 """
 # Space for time processor
 
 ## Arguments:
- ```cube_con``` : YAXARRAY with the continous variable to be analyized.
 
- ```cube_classes```: YAXARRAY with the discrete classes to be used in the space4time.
 
- ```time_axis_name``` : String or nothing. Name of the time axis on the input cubes. By default ```time_axis_name = "time"```. if ```time_axis_name = nothing```, not time dimension considered.
 
- ```lon_axis_name``` : String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lon"```
 
- ```lat_axis_name``` :  String. Name of the longitude axis on the input cubes. By default ```lon_axis_name = "lat"```
 
- ```classes_var_name``` : String. Name of the Variable containing the discrete classes. By default ```classes_var_name = "classes"```.
 
- ```winsize```: Edge size of the moving window on pixels. By default winsize = 5. E.g. ```winsize = 5``` will produce a moving window with 5^2 pixels.
 
- ```minpxl``` : Minimum number of pixels in the moving window. By default minpxl = 25. Change accordindly to your ```winsize``` parameter.
 
- ```minDiffPxlspercentage```: Percentage of minimum number pixels in the moving window that must have different compositions. Must be any value in the interval 30-100. By default minDiffPxlspercentage = 40
 
- ```classes_vec```: A string vector with the names of the classes on ```cube_classes``` to be used. e.g. from MPI-BGC internal structure ```classes_vec = ["Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests", "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests", "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands", "Permanent_Wetlands", "Croplands", "Urban_and_Built-up_Lands", "Cropland/Natural_Vegetation_Mosaics", "Permanent_Snow_and_Ice", "Barren", "Water_Bodies"]```

- ```max_value```: Indicates if the scale of the presence of the discrete classes if from 0 to 1 or 0 to 100 if ```max_value = 100``` then the data is re-scaled from 0 to 1. By default ```max_value = 1```

- ```showprog```: Show progress bar. By default ```showprog = true```

- ```max_cache```: Size of the cache to allocate temporarily sections of the cubes. By default ```max_cache = 1e8```

 ## Output:
 The ```space4time_proc``` produces a YAXARRAY.Dataset with three cubes:
 - SummaryStats cube has one axis ```summary_stat```, and three variables:
    - ```rsquared```:  
    - ```cumulative_variance```:
    - ```predicted```:
 - ```metrics_for_classes``` cube has one axis ```Values of Z for pure classes```, and two variables:
    - ```estimated```:
    - ```estimated_error```:
 - metrics_for_transitions has two axis ```transitions``` (all the transitions by pairs between the different classes), and ```Differences``` with three variables:
    - ```delta```:
    - ```delta_error```:
    - ```coocurence```:
 """ 
 function space4time_proc(cube_con, cube_classes; time_axis_name = "time", lon_axis_name = "lon", lat_axis_name = "lat", classes_var_name = "classes", winsize = 5, minDiffPxlspercentage = 40, classes_vec = NaN, max_value = NaN, minpxl = 25, showprog = true, max_cache=1e8)
    
    # Checking that winsize is odd
    
    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) -1
        
        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end
    #println(size(cube_con))
    #println(size(cube_classes))
    # assuming the first dimmension is time.
    if !isnothing(time_axis_name)
        time_n = length(getAxis(time_axis_name,cube_con).values)
    else
        time_n = 1
    end
    # assuming that pfts presence change in time. last dimmension refer to the number of pfts.
    # pfts_cube = rand(Float32, (5,5, length(classes_vec)))
    
    # number of classes
    # assuming the last dimmension is PFTs
    nc = length(classes_vec)
    
    sigma1_glob = [fill(NaN, (nc, nc)) for i in 1:Threads.nthreads()]
    
    prederr_glob = [fill(NaN, nc) for i in 1:Threads.nthreads()]
    
    predres_glob = [fill(NaN, nc) for i in 1:Threads.nthreads()]
    # lower triangular matrix index use further on
    
    ltriindex = NamedArray(LowerTriangular(fill(1, (nc, nc))))
    
    tril!(ltriindex, -1)
    
    setnames!(ltriindex, classes_vec, 1)
    setnames!(ltriindex, classes_vec, 2)
    
    tran_check = findall(>(0), ltriindex)
    
    # set names of transition combinations
    pftstrans_comb_names = collect(combinations(classes_vec, 2))
    
    # number of transitions
    
    transitions_n = length(pftstrans_comb_names)
    #println("transitions _ n ", transitions_n)
    
    # linear regression
    
    empty_models = [[empty_model(i, winsize^2) for i in 1:nc] for j in 1:Threads.nthreads()]
    
    p1_static = range(0, 1, length = winsize^2)
    
    p2_static = reverse(p1_static)

    out_pmindist_global = [zeros(1, winsize^2) for i in 1:Threads.nthreads()]

    denom = sum(sqrt.(sum.(eachrow([p1_static   p2_static].^2))))
    
    half = floor(Int, ceil((winsize^2)/2))
    
    minDiffPxls = (winsize^2 * minDiffPxlspercentage / 100)
    
    localcomp_fix_glob = [rand(winsize^2) for i in 1:Threads.nthreads()]
    
    pftsvarmat_f_glob = [rand(winsize^2, nc+1) for i in 1:Threads.nthreads()]
    
    
    # 
    if !isnothing(time_axis_name)

        
        indims = InDims(time_axis_name, MovingWindow(lon_axis_name, pre_step,after_step),
        MovingWindow(lat_axis_name, pre_step, after_step),  window_oob_value=NaN)
        
        indims_classes = InDims(time_axis_name, MovingWindow(lon_axis_name, pre_step,after_step),
        MovingWindow(lat_axis_name, pre_step, after_step), 
        classes_var_name, 
        window_oob_value=NaN)
        
        out_1_dims = OutDims(time_axis_name, CategoricalAxis("summary_stat", ["rsquared", "cumulative_variance", "predicted"]))
        
        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(time_axis_name,CategoricalAxis("classes", classes_vec), CategoricalAxis("Values_of_Z_for_pure_classes", ["estimated", "estimated_error"]))
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(time_axis_name, CategoricalAxis("transitions", [join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]), CategoricalAxis("Differences", ["delta", "delta_error", "coocurence"]))
        
    else
        indims = InDims(MovingWindow(lon_axis_name, pre_step,after_step),
        MovingWindow(lat_axis_name, pre_step, after_step),  window_oob_value=NaN)
        
        indims_classes = InDims(MovingWindow(lon_axis_name, pre_step,after_step),
        MovingWindow(lat_axis_name, pre_step, after_step), 
        classes_var_name, 
        window_oob_value=NaN)
        
        out_1_dims = OutDims(RangeAxis("time", [1]),CategoricalAxis("summary_stat", ["rsquared", "cumulative_variance", "predicted"]))
        
        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(RangeAxis("time", [1]), CategoricalAxis("classes", classes_vec), CategoricalAxis("Values_of_Z_for_pure_classes", ["estimated", "estimated_error"]))
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(RangeAxis("time", [1]), CategoricalAxis("transitions", [join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]), CategoricalAxis("Differences", ["delta", "delta_error", "coocurence"]))
        
    end
    
    
    outdims = (out_1_dims, out_2_dims, out_3_dims)
    #println(out_3_dims)
    
    out_1, out_2, out_3 = mapCube(s4time, (cube_con, cube_classes), 
    indims = (indims, indims_classes), 
    outdims = outdims, max_cache=max_cache,  showprog = showprog, include_loopvars=true; empty_models, pft_list = classes_vec, time_n = time_n, max_value = max_value, p1_static, p2_static, sigma1_glob, prederr_glob, predres_glob, minDiffPxls,tran_check, half, localcomp_fix_glob, pftsvarmat_f_glob, winsize = winsize, transitions_n = transitions_n, pftstrans_comb_names = pftstrans_comb_names, nc = nc, out_pmindist_global = out_pmindist_global, denom = denom, minpxl = minpxl)

    return Dataset(;SummaryStats=out_1, metrics_for_classes=out_2, metrics_for_transitions=out_3)
    
end