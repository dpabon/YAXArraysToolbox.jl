function coocufun(out, q1, q2, p1, p2, denom)
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
    if length(vecq1) == length(vecq2) && isfinite(sum(vecq1)) && isfinite(sum(vecq2)) && sum(vecq1)> 0. && sum(vecq2) > 0.
        #@show typeof(vecq1)
        #@show typeof(vecq2)
        #@show size(vecq1) size(vecq2)
        pmindist = minimum(([i - j for i in vecq1, j in p1]).^2 + ([i - j for i in vecq2, j in p2]).^2, dims = 1)
        #out_pmindist = minimum(([i - j for i in vecq1, j in p1]) .^ 2 + ([i - j for i in vecq2, j in p2]) .^ 2)

        out[1] = 1 - (sum(sqrt.(pmindist)) / denom)

        #return 1 - (sum(sqrt.(pmindist)) / denom)
    else
        out[1] = 0.0 # 0 because assume not enough data on the moving window.
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
    * "predicted": Mean prediction of Z for moving window with the real combination of values.
    out2: 
    
#Examples
    
    
"""
function s4time(
    out_1,
    out_2,
    out_3,
    clim_var_cube_in,
    pfts_cube_in,
    altitude_cube_in,
    loopvars;
    pft_list::Vector{String},
    time_n,
    max_value::Int,
    p1_static,
    p2_static,
    #sigma1_glob,
    #prederr_glob,
    #predres_glob,
    minDiffPxls,
    tran_check,
    half,
    #localcomp_fix_glob,
    #pftsvarmat_f_glob,
    winsize = 5,
    transitions_n,
    pftstrans_comb_names,
    nc,
    #out_pmindist_global,
    denom,
    minpxl,
    minDiffPxls_alt,
    )
    #println(size(clim_var_cube_in))
    #println(size(pfts_cube_in))
    #println(size(out_3))
    #igma1 = sigma1_glob[Threads.threadid()]
    sigma1 = fill(NaN, (nc, nc))
    #prederr = prederr_glob[Threads.threadid()]
    prederr = fill(NaN, nc) 
    #predres = predres_glob[Threads.threadid()]
    predres = fill(NaN, nc)
    #localcomp_fix = localcomp_fix_glob[Threads.threadid()]
    #pftsvarmat_f = pftsvarmat_f_glob[Threads.threadid()]
    #out_pmindist = out_pmindist_global[Threads.threadid()]
    #out_pmindist = zeros(1, winsize^2)
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
    
    if isnothing(time_n)
        climvarmat = reshape(clim_var_cube_in, (winsize^2))
        climvarmat = convert(Array{Float64}, climvarmat)
        
    else
        climvarmat = reshape(clim_var_cube_in, ((winsize^2), time_n))
        climvarmat = convert(Array{Float64}, climvarmat)
    end
    
    
    climvarmat = replace!(climvarmat, missing => NaN)
    
    
    local_pft1 =
    [findall(pft_list .== pftstrans_comb_names[comp][1]) for comp = 1:transitions_n]
    local_pft1 = reduce(vcat, local_pft1)
    local_pft2 =
    [findall(pft_list .== pftstrans_comb_names[comp][2]) for comp = 1:transitions_n]
    local_pft2 = reduce(vcat, local_pft2)
    
    #println("before if")
    #println(out_3[1,1,1])
    
    #println(local_pft1)
    #println(local_pft2)
    #pfts_cube_in_1 = replace!(pfts_cube_in, NaN => 0.0)
    pfts_cube_in_1 = pfts_cube_in
    #replace!(pfts_cube_in_1, missing => 0.0)
    #replace!(pfts_cube_in_1, NaN32 => 0.0)
    #replace!(pfts_cube_in_1, NaN16 => 0.0)
    pfts_cube_in_1 = convert(Array{Float64}, pfts_cube_in_1)
    
    
    if isnothing(time_n)
        
        #pfts_cube_in_1 =
        #permutedims(reshape(pfts_cube_in_1, (winsize, winsize, nc, 1)), (4, 1, 2, 3))
        
        for comp in eachindex(local_pft1)
            #println(comp)
            #println(size(pfts_cube_in))
            #println("I'm here")
            #println(pfts_cube_in[it,:,:,local_pft1[comp]])
            #println(pfts_cube_in[it,:,:,local_pft2[comp]])
            # define the pfts to be processed
            out_3[comp, 3] = coocufun(
            [0.0],
            pfts_cube_in_1[:, :, local_pft1[comp]],
            pfts_cube_in_1[:, :, local_pft2[comp]],
            p1_static,
            p2_static,
            denom,
            )
            #println("all good")
        end
        
        
        
        pfts_cube_in_2 = pfts_cube_in_1[:, :, :]
        #println(all(isnan, pfts_cube_in_2))
        
        #println(size(pfts_cube_in_2))
        
        pftsvarmat = reshape(pfts_cube_in_2, (winsize^2, nc))

        if count(!isnan, climvarmat[:]) >= minpxl
            

            # altitude processing
            altitude_center_1 = altitude_cube_in[Int(round(winsize/2)+1), Int(round(winsize/2)+1), 1]
            altitude_center_2 = altitude_cube_in[Int(round(winsize/2)+1), Int(round(winsize/2)+1), 2]

            altitude_1 = reshape(altitude_cube_in[:,:,1], winsize^2) .- altitude_center_1

            altitude_2 = reshape(altitude_cube_in[:,:,2], winsize^2) .- altitude_center_2

            altitude_1 = Float64.(altitude_1)
            altitude_2 = Float64.(altitude_2)

            
            if count(isnan, climvarmat[:]) != 0
                
                pftsvarmat = pftsvarmat[findall(!isnan, climvarmat[:]), :]
                
                altitude_1 = altitude_1[findall(!isnan, climvarmat[:])]

                altitude_2 = altitude_2[findall(!isnan, climvarmat[:])]

                climvarmat = filter(!isnan, climvarmat)
                

            end
            
            if isfinite(sum(pftsvarmat)) && sum(sum(pftsvarmat, dims = 1) .> 0.) > 1
                #println("test")
                #println(any(isnan.(pftsvarmat)))
                # check if pftsvarmat is 0 to 1 or 0 to 100
                #println(maximum(vec(pftsvarmat)))
                
                # make sure compositions are really precisely right.
                
                #localcomp_fix_glob = mapslices(x->1-sum(x), pftsvarmat, dims = 2)
                localcomp_fix = map(x -> 1 - sum(x), eachslice(pftsvarmat, dims = 1))
                #map!(x->1-sum(x), localcomp_fix_glob, eachslice(pftsvarmat, dims = 1))
                #println(size(localcomp_fix_glob))
                #println(localcomp_fix_glob)
                pftsvarmat_f = [pftsvarmat localcomp_fix]
                
                map!((x) -> round(x, digits = 4), pftsvarmat_f, pftsvarmat_f)
                
                # some PFTs might not be present in the 5*5 window
                # these must be identified and removed, as they cannot be predicted
                
                #pftpres_check = vec(mapslices(sum, pftsvarmat, dims = 1) .> 0)
                pftpres_check = vec(sum(pftsvarmat_f, dims = 1) .> 0)
                
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
                
                if uniquepixels >= minDiffPxls && sum(pftpres_check) > 1
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
                        
                        temp = cumsum(lcsvd.S .^ 2) ./ sum(lcsvd.S .^ 2)
                        try
                            minimum(findall(temp .>= 1))
                        catch e
                            return
                        end
                        
                        ndim = minimum(findall(temp .>= 1))
                        
                        # store results to output object
                        # cumulative variance
                        out_1[2] = sum(lcsvd.S)
                        
                        #println(out7_cumulated_variance)
                        
                        # dimmensions that explain 100 % of the variance
                        
                        lr = lc2 * lcsvd.V[:, 1:ndim]
                        
                        
                        # create bogus composition dataset
                        
                        #boguscomp = zeros(Float64, nc+1, nc+1)
                        #boguscomp[diagind(boguscomp)] .= 1
                        
                        boguscomp = I(nc + 1)
                        # println(boguscomp)
                        
                        # remove absent pfts from bogus predictor compositions and close compositions.
                        #
                        bogusc1 =
                        mapslices(
                        x -> x / sum(x),
                        boguscomp[pftpres_check, pftpres_check],
                        dims = 2,
                        )'
                        
                        #println(size(bogusc1))
                        #println(pftpres_check)
                        
                        # center the columns as the training data were centered
                        
                        bogusc2 = (I(sum(pftpres_check)) .- lCm[pftpres_check])
                        
                        bogusc2 = (bogusc1' .- lCm[pftpres_check])
                        
                        bogusc3 = bogusc2 * lcsvd.V[pftpres_check, 1:ndim]
                        #println(climvarmat[:,it])
                        #println("test")
                        # data = hcat(DataFrame(lt = convert(Vector{Float64}, climvarmat[:,it])), DataFrame(lr, :auto))
                        
                        # compreg = GLM.lm(Term(:lt) ~ sum(Term.(Symbol.(names(data[:, Not(:lt)])))), data)
                        
                        #println("before fail")

                        uniquepixels_char_altitude_1 = unique(round.(altitude_1, digits = 6))
                        uniquepixels_char_altitude_2 = unique(round.(altitude_2, digits = 6))
                
                        uniquepixels_altitude_1 = length(uniquepixels_char_altitude_1)
                        uniquepixels_altitude_2 = length(uniquepixels_char_altitude_2)

                        n_altitude = NaN

                        #println("unique altitude 1 = ")
                        #println(uniquepixels_char_altitude_1)

                        #println("unique altitude 2 = ")
                        #println(uniquepixels_char_altitude_2)

                        if uniquepixels_altitude_1 >= minDiffPxls_alt && uniquepixels_altitude_2 >= minDiffPxls_alt
                            ols = lm([ones(size(lr, 1)) lr altitude_1 altitude_2], identity.(climvarmat[:]); method=:qr, dropcollinear = false)
                            n_altitude = 4
                            out_1[4] = coeftable(ols).cols[4][end-1]
                            out_1[5] = coeftable(ols).cols[4][end]
                            out_1[8] = 4
                            
                        else
                            if uniquepixels_altitude_1 >= minDiffPxls_alt
                                ols = lm([ones(size(lr, 1)) lr altitude_1], identity.(climvarmat[:]); method=:qr, dropcollinear = false)
                                n_altitude = 2
                                out_1[4] = coeftable(ols).cols[4][end]
                                out_1[5] = NaN
                                out_1[8] = 2
                                
                            else
                                if uniquepixels_altitude_2 >= minDiffPxls_alt
                                    ols = lm([ones(size(lr, 1)) lr altitude_2], identity.(climvarmat[:]); method=:qr, dropcollinear = false)
                                    n_altitude = 3
                                    out_1[4] = NaN
                                    out_1[5] = coeftable(ols).cols[4][end]
                                    out_1[8] = 3
                                    
                                else 
                                    if uniquepixels_altitude_1 < minDiffPxls_alt && uniquepixels_altitude_2 < minDiffPxls_alt
                                        ols = lm([ones(size(lr, 1)) lr], identity.(climvarmat[:]); method=:qr, dropcollinear = false)
                                        n_altitude = 1
                                        out_1[4] = NaN
                                        out_1[5] = NaN
                                        out_1[8] = 1
                                    end
                                end
                            end
                        end
                        coef_reg = GLM.coef(ols)

                        #println("n_altitude =" *string(n_altitude))

                        #println(ols)
                        # continue only if there are no NA in the estimated coefficients
                        #println("climate_var = ")

                        #println(climvarmat[:])

                        #println("PFTs_var = ")
                        #println(lr)

                        #println("alt_1_var = ")
                        #println(altitude_1)

                        #println("alt_2_var = ")
                        #println(altitude_2)
                        
                        
                        #println("original coef $coef_reg")
                        #println("second estimation coef $(lr\climvarmat[:,it])")
                        
                        if isfinite(sum(coef_reg))
                            # then do predictions for the log-normal approach
                            if isa(bogusc2, Vector)
                                if n_altitude == 4
                                    boguspred = predict(
                                    ols,
                                    [ones(length(bogusc3)) bogusc3 zeros(length(bogusc3)) zeros(length(bogusc3))],
                                    )
                                else
                                    if n_altitude == 2 || n_altitude == 3
                                        boguspred = predict(
                                        ols,
                                        [ones(length(bogusc3)) bogusc3 zeros(length(bogusc3))],
                                        )
                                    else
                                        if n_altitude == 1
                                            boguspred = predict(
                                            ols,
                                            [ones(length(bogusc3)) bogusc3],
                                            )
                                        end
                                    end
                                end
                                
                                # boguspred = GLM.predict(compreg, DataFrame( x1 = bogusc3))
                                
                            else
                                # boguspred = GLM.predict(compreg, DataFrame(bogusc3, :auto))
                                if n_altitude == 4
                                    boguspred = predict(
                                    ols,
                                    [ones(size(bogusc3, 1)) bogusc3 zeros(size(bogusc3, 1)) zeros(size(bogusc3, 1))],
                                    )
                                else
                                    if n_altitude == 2 || n_altitude == 3
                                        boguspred = predict(
                                        ols,
                                        [ones(size(bogusc3, 1)) bogusc3 zeros(size(bogusc3, 1))],
                                        )
                                    else
                                        if n_altitude == 1
                                            boguspred = predict(
                                            ols,
                                            [ones(size(bogusc3, 1)) bogusc3],
                                            )
                                        end
                                    end
                                end
                            end
                            
                            if n_altitude == 4
                                x2pred = [ones(size(bogusc3, 1), 1) bogusc3 zeros(size(bogusc3, 1), 1) zeros(size(bogusc3, 1), 1)]

                            else
                                if n_altitude == 2 || n_altitude == 3
                                    x2pred = [ones(size(bogusc3, 1), 1) bogusc3 zeros(size(bogusc3, 1), 1)]
                                    # println("n_altitude = " * string(x2pred))

                                else
                                    if n_altitude  == 1
                                        x2pred = [ones(size(bogusc3, 1), 1) bogusc3]
                                    end
                                end
                            end
                            # println("n_altitude = " * string(n_altitude))
                            # x2pred = [ones(size(bogusc3, 1), 1) bogusc3 zeros(size(bogusc3, 1), 1) zeros(size(bogusc3, 1), 1)]
                            
                            vcv = GLM.vcov(ols)
                            # vcv = GLM.vcov(compreg)
                            
                            sigma = x2pred * vcv * x2pred'
                            
                            # now store the target variables
                            # but make sure appropiate temperatures go back to appropiate
                            # pfts (as absent pfts were removed)
                            
                            # value of climatevar for pure ptfs
                            
                            predres .= NaN
                            #println(boguspred)
                            predres[view(pftpres_check, 1:nc)] = boguspred
                            #println(it)
                            
                            
                            out_2[:, 1] = predres
                            
                            
                            prederr .= NaN
                            
                            prederr[view(pftpres_check, 1:nc)] = sqrt.(diag(sigma))
                            
                            
                            out_2[:, 2] = prederr
                            
                            
                            # prediction of varclim for the central pixel with its real pft combination
                            
                            out_1[3] = mean(StatsModels.predict(ols))
                            
                            # Rsquare of the regression
                            
                            out_1[1] = adjr2(ols)

                            
                            out_1[6] = aic(ols)
                            out_1[7] = aic(ols)

                            
                            # println(out_1)
                            # println(r2(compreg))
                            
                            # and now for the transitions
                            # only the PFTs identified in the pftlist are to be used
                            
                            sigma1 .= NaN
                            sigma1[view(pftpres_check, 1:nc), view(pftpres_check, 1:nc)] =
                            sigma
                            
                            # calculate the difference on climatevar caused by going from one pft to another.
                            
                            diff_clim_pft_pred = round.((predres .- predres'), digits = 10)
                            
                            diff_clim_pft_pred = diff_clim_pft_pred[tran_check]
                            
                            # propagate the error (as variances) taking into account
                            # the covariance terms
                            # in the original implementation diffclim_pft_var = dZvar
                            
                            diff_clim_pft_pred_var =
                            round.(
                            (diag(sigma1) .+ diag(sigma1)') .- 2 * sigma1,
                            digits = 10,
                            )
                            
                            #print(diff_clim_pft_pred)
                            # flag out those with zero error (may occur with identical compositions for 2 pfts)
                            
                            diff_clim_pft_pred_var = diff_clim_pft_pred_var[tran_check]
                            
                            #println(diff_clim_pft_pred)
                            
                            # flag out those with zero error (may occur with identical compositions for 2 pfts)
                            
                            if any(round.(diff_clim_pft_pred, digits = 8) .== 0.0)
                                diff_clim_pft_pred[round.(
                                diff_clim_pft_pred,
                                digits = 8,
                                ).==0] .= NaN
                            end
                            
                            if any(round.(diff_clim_pft_pred_var, digits = 8) .== 0.0)
                                diff_clim_pft_pred_var[round.(
                                diff_clim_pft_pred_var;
                                digits = 8,
                                ).==0] .= NaN
                            end
                            
                            # mask out low co-ocurrence ask to greg!!! This can be performed  masking the pixels
                            # before all estimations
                            
                            #println("inside loop")
                            #println(diff_clim_pft_pred)
                            
                            if length(diff_clim_pft_pred) == 1
                                
                                out_3[1, 1] = diff_clim_pft_pred[1]
                                #println(out1_delta)
                                if diff_clim_pft_pred_var[1] .< 0
                                    out_3[1, 2] = NaN
                                else
                                    out_3[1, 2] = sqrt.(diff_clim_pft_pred_var)[1]
                                    
                                end
                                
                                
                            else
                                
                                out_3[1:transitions_n, 1] = diff_clim_pft_pred
                                #println(out1_delta)
                                #if diff_clim_pft_pred_var .< 0
                                #out_3[it,1:transitions_n,2] .= NaN
                                
                                #else
                                diff_clim_pft_pred_var[diff_clim_pft_pred_var.<0] .= NaN
                                out_3[1:transitions_n, 2] =
                                sqrt.(diff_clim_pft_pred_var)
                                
                                #end
                            end
                        end
                    end
                end
            end
        end
        # for debug from R -------
        
        # pftsvarmat = Matrix(CSV.read("/Net/Groups/BGI/people/dpabon/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
        # pftsvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
        # climvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_temperature_example_from_R.csv", DataFrame))
        # check that there are not NaN values on pfts and at least one pft is present
        # @show sum(pftsvarmat), sum(pftsvarmat)
        #println(pftsvarmat)
    else
        
        if time_n == 1
            
            pfts_cube_in_1 =
            permutedims(reshape(pfts_cube_in_1, (winsize, winsize, nc, 1)), (4, 1, 2, 3))
            
        end
        
        for it = 1:time_n
            # co-ocurrence estimation
            
            #println(it)
            #println(local_pft1)
            #println(local_pft2)
            #println(transitions_n)
            for comp in eachindex(local_pft1)
                #println(comp)
                #println(size(pfts_cube_in))
                #println("I'm here")
                #println(pfts_cube_in[it,:,:,local_pft1[comp]])
                #println(pfts_cube_in[it,:,:,local_pft2[comp]])
                # define the pfts to be processed
                out_3[it, comp, 3] = coocufun(
                [0.0],
                pfts_cube_in_1[it, :, :, local_pft1[comp]],
                pfts_cube_in_1[it, :, :, local_pft2[comp]],
                p1_static,
                p2_static,
                denom,
                )
                #println("all good")
            end
            
            
            
            pfts_cube_in_2 = pfts_cube_in_1[it, :, :, :]
            #println(all(isnan, pfts_cube_in_2))
            
            #println(size(pfts_cube_in_2))
            
            pftsvarmat = reshape(pfts_cube_in_2, (winsize^2, nc))
            
            
            if count(!isnan, climvarmat[:, it]) >= minpxl

                # altitude processing
                altitude_center_1 = altitude_cube_in[round(winsize/2)+1, round(winsize/2)+1, 1]
                altitude_center_2 = altitude_cube_in[round(winsize/2)+1, round(winsize/2)+1, 2]

                altitude_1 = altitude_center_1 .- reshape(altitude_cube_in[:,:,1], winsize^2)
                altitude_2 = altitude_center_2 .- reshape(altitude_cube_in[:,:,2], winsize^2)

                
                if count(!isnan, climvarmat[:, it]) != 0
                    
                    pftsvarmat = pftsvarmat[findall(!isnan, climvarmat[:, it]), :]
                    
                    climvarmat_it = filter(!isnan, climvarmat[:,it])

                    altitude_1 = altitude_1[findall(!isnan, climvarmat[:, it])]

                    altitude_2 = altitude_2[findall(!isnan, climvarmat[:, it])]

                end

                # for debug from R -------
            
                # pftsvarmat = Matrix(CSV.read("/Net/Groups/BGI/people/dpabon/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
                # pftsvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
                # climvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_temperature_example_from_R.csv", DataFrame))
                # check that there are not NaN values on pfts and at least one pft is present
                # @show sum(pftsvarmat), sum(pftsvarmat)
                #println(pftsvarmat)
            
            
                #println(it)
                
                if isfinite(sum(pftsvarmat)) && sum(sum(pftsvarmat, dims = 1) .> 0.) > 1
                    #println("test")
                    #println(any(isnan.(pftsvarmat)))
                    # check if pftsvarmat is 0 to 1 or 0 to 100
                    #println(maximum(vec(pftsvarmat)))
                    
                    # make sure compositions are really precisely right.
                    
                    #localcomp_fix_glob = mapslices(x->1-sum(x), pftsvarmat, dims = 2)
                    localcomp_fix = map(x -> 1 - sum(x), eachslice(pftsvarmat, dims = 1))
                    #map!(x->1-sum(x), localcomp_fix_glob, eachslice(pftsvarmat, dims = 1))
                    #println(size(localcomp_fix_glob))
                    #println(localcomp_fix_glob)
                    pftsvarmat_f = [pftsvarmat localcomp_fix]
                    
                    map!((x) -> round(x, digits = 4), pftsvarmat_f, pftsvarmat_f)
                    
                    # some PFTs might not be present in the 5*5 window
                    # these must be identified and removed, as they cannot be predicted
                    
                    #pftpres_check = vec(mapslices(sum, pftsvarmat, dims = 1) .> 0)
                    pftpres_check = vec(sum(pftsvarmat_f, dims = 1) .> 0)
                    
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
                    
                    if uniquepixels >= minDiffPxls && sum(pftpres_check) > 1
                        #println("test")
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
                            
                            temp = round.(cumsum(lcsvd.S .^ 2) ./ sum(lcsvd.S .^ 2), digits = 8)
                            ndim = minimum(findall(temp .>= 1))
                            
                            # store results to output object
                            # cumulative variance
                            out_1[it, 2] = sum(lcsvd.S)
                            
                            #println(out7_cumulated_variance)
                            
                            # dimmensions that explain 100 % of the variance
                            
                            lr = lc2 * lcsvd.V[:, 1:ndim]
                            
                            
                            # create bogus composition dataset
                            
                            #boguscomp = zeros(Float64, nc+1, nc+1)
                            #boguscomp[diagind(boguscomp)] .= 1
                            
                            boguscomp = I(nc + 1)
                            # println(boguscomp)
                            
                            # remove absent pfts from bogus predictor compositions and close compositions.
                            #
                            bogusc1 =
                            mapslices(
                            x -> x / sum(x),
                            boguscomp[pftpres_check, pftpres_check],
                            dims = 2,
                            )'
                            
                            #println(size(bogusc1))
                            #println(pftpres_check)
                            
                            # center the columns as the training data were centered
                            
                            bogusc2 = (I(sum(pftpres_check)) .- lCm[pftpres_check])
                            
                            bogusc2 = (bogusc1' .- lCm[pftpres_check])
                            
                            bogusc3 = bogusc2 * lcsvd.V[pftpres_check, 1:ndim]
                            #println(climvarmat[:,it])
                            
                            #println("test")
                            # data = hcat(DataFrame(lt = convert(Vector{Float64}, climvarmat[:,it])), DataFrame(lr, :auto))
                            
                            # compreg = GLM.lm(Term(:lt) ~ sum(Term.(Symbol.(names(data[:, Not(:lt)])))), data)
                            
                            #println("before fail")

                            uniquepixels_char_altitude_1 = unique(round.(altitude_1, digits = 6))
                            uniquepixels_char_altitude_2 = unique(round.(altitude_2, digits = 6))
                
                            uniquepixels_altitude_1 = length(uniquepixels_char_altitude_1)
                            uniquepixels_altitude_2 = length(uniquepixels_char_altitude_2)

                            n_altitude = NaN

                            if uniquepixels_altitude_1 >= minDiffPxls_alt && uniquepixels_altitude_2 >= minDiffPxls_alt
                                ols = lm([ones(size(lr, 1)) lr altitude_1 altitude_2], identity.(climvarmat_it[:]); method=:qr, dropcollinear = false)
                                n_altitude = 4
                                out_1[4] = coeftable(ols).cols[4][end-1]
                                out_1[5] = coeftable(ols).cols[4][end]
                                out_1[8] = 4

                            else
                                if uniquepixels_altitude_1 >= minDiffPxls_alt
                                    ols = lm([ones(size(lr, 1)) lr altitude_1], identity.(climvarmat_it[:]); method=:qr, dropcollinear = false)
                                    n_altitude = 2
                                    out_1[4] = coeftable(ols).cols[4][end]
                                    out_1[5] = NaN
                                    out_1[8] = 2

                                else
                                    if uniquepixels_altitude_2 >= minDiffPxls_alt
                                        ols = lm([ones(size(lr, 1)) lr altitude_2], identity.(climvarmat_it[:]); method=:qr, dropcollinear = false)
                                        n_altitude = 3
                                        out_1[4] = NaN
                                        out_1[5] = coeftable(ols).cols[4][end]
                                        out_1[8] = 3
                                    else 
                                        if uniquepixels_altitude_1 < minDiffPxls_alt && uniquepixels_altitude_2 < minDiffPxls_alt
                                            ols = lm([ones(size(lr, 1)) lr], identity.(climvarmat_it[:]); method=:qr, dropcollinear = false)
                                            n_altitude = 1
                                            out_1[4] = NaN
                                            out_1[5] = NaN
                                            out_1[8] = 1
                                        end
                                    end
                                end
                            end                            
                            # continue only if there are no NA in the estimated coefficients
                            
                            coef_reg = GLM.coef(ols)
                            #println("original coef $coef_reg")
                            #println("second estimation coef $(lr\climvarmat[:,it])")
                            
                            if isfinite(sum(coef_reg))
                                # then do predictions for the log-normal approach
                                if isa(bogusc2, Vector)

                                    if n_altitude == 4
                                        boguspred = predict(
                                        ols,
                                        [ones(length(bogusc3)) bogusc3 zeros(length(bogusc3)) zeros(length(bogusc3))],
                                        )
                                    
                                    else
                                        if n_altitude == 2 || n_altitude == 3
                                            boguspred = predict(
                                            ols,
                                            [ones(length(bogusc3)) bogusc3 zeros(length(bogusc3))],
                                            )
                                        else
                                            if n_altitude  == 1
                                            boguspred = predict(
                                            ols,
                                            [ones(length(bogusc3)) bogusc3],
                                            )
                                            end
                                        end                            
                                    end
                                    # boguspred = GLM.predict(compreg, DataFrame( x1 = bogusc3))
                                    
                                else
                                    if n_altitude == 4
                                        boguspred = predict(
                                        ols,
                                        [ones(size(bogusc3, 1)) bogusc3 zeros(size(bogusc3, 1)) zeros(size(bogusc3, 1))],
                                        )

                                    else
                                        if n_altitude == 2 || n_altitude == 3
                                        boguspred = predict(
                                        ols,
                                        [ones(size(bogusc3, 1)) bogusc3 zeros(size(bogusc3, 1))],
                                        )

                                        else
                                            if n_altitude  == 1
                                                boguspred = predict(
                                                ols,
                                                [ones(size(bogusc3, 1)) bogusc3],
                                                )
                                                                    
                                            end
                                        end
                                    end
                                    # boguspred = GLM.predict(compreg, DataFrame(bogusc3, :auto))   
                                
                                end
                                
                                if n_altitude == 4
                                    x2pred = [ones(size(bogusc3, 1), 1) bogusc3 zeros(size(bogusc3, 1), 1) zeros(size(bogusc3, 1), 1)]

                                else
                                    if n_altitude == 2 || n_altitude == 3
                                        x2pred = [ones(size(bogusc3, 1), 1) bogusc3 zeros(size(bogusc3, 1), 1)]

                                    else
                                        if n_altitude  == 0
                                            x2pred = [ones(size(bogusc3, 1), 1) bogusc3]     
                                        end
                                    end
                                end
                                
                                vcv = GLM.vcov(ols)
                                # vcv = GLM.vcov(compreg)
                                
                                sigma = x2pred * vcv * x2pred'
                                
                                # now store the target variables
                                # but make sure appropiate temperatures go back to appropiate
                                # pfts (as absent pfts were removed)
                                
                                # value of climatevar for pure ptfs
                                
                                predres .= NaN
                                #println(boguspred)
                                predres[view(pftpres_check, 1:nc)] = boguspred
                                #println(it)
                                
                                out_2[it, :, 1] = predres
                                
                                prederr .= NaN
                                
                                prederr[view(pftpres_check, 1:nc)] = sqrt.(diag(sigma))
                                
                                out_2[it, :, 2] = prederr
                                
                                # prediction of varclim for the central pixel with its real pft combination
                                # It is possible that the center pixel correspond to an observation with NaN that was removed then
                                # I redefined as just mean of the prediction
                                
                                out_1[it, 3] = mean(StatsModels.predict(ols))
                                
                                # Rsquare of the regression
                                
                                out_1[it, 1] = adjr2(ols)
                                
                                out_1[it, 6] = aic(ols)
                                out_1[it, 7] = aicc(ols)
         
                                # println(out_1)
                                # println(r2(compreg))
                                
                                # and now for the transitions
                                # only the PFTs identified in the pftlist are to be used
                                
                                sigma1 .= NaN
                                sigma1[view(pftpres_check, 1:nc), view(pftpres_check, 1:nc)] =
                                sigma
                                
                                # calculate the difference on climatevar caused by going from one pft to another.
                                
                                diff_clim_pft_pred = round.((predres .- predres'), digits = 10)
                                
                                diff_clim_pft_pred = diff_clim_pft_pred[tran_check]
                                
                                # propagate the error (as variances) taking into account
                                # the covariance terms
                                # in the original implementation diffclim_pft_var = dZvar
                                
                                diff_clim_pft_pred_var =
                                round.(
                                (diag(sigma1) .+ diag(sigma1)') .- 2 * sigma1,
                                digits = 10,
                                )
                                
                                #print(diff_clim_pft_pred)
                                # flag out those with zero error (may occur with identical compositions for 2 pfts)
                                
                                diff_clim_pft_pred_var = diff_clim_pft_pred_var[tran_check]
                                
                                #println(diff_clim_pft_pred)
                                
                                # flag out those with zero error (may occur with identical compositions for 2 pfts)
                                
                                if any(round.(diff_clim_pft_pred, digits = 8) .== 0.0)
                                    diff_clim_pft_pred[round.(
                                    diff_clim_pft_pred,
                                    digits = 8,
                                    ).==0] .= NaN
                                end
                                
                                if any(round.(diff_clim_pft_pred_var, digits = 8) .== 0.0)
                                    diff_clim_pft_pred_var[round.(
                                    diff_clim_pft_pred_var;
                                    digits = 8,
                                    ).==0] .= NaN
                                end
                                
                                # mask out low co-ocurrence ask to greg!!! This can be performed  masking the pixels
                                # before all estimations
                                
                                #println("inside loop")
                                #println(diff_clim_pft_pred)
                                
                                if length(diff_clim_pft_pred) == 1
                                    
                                    out_3[it, 1, 1] = diff_clim_pft_pred[1]
                                    #println(out1_delta)
                                    if diff_clim_pft_pred_var[1] .< 0
                                        out_3[it, 1, 2] = NaN
                                    else
                                        out_3[it, 1, 2] = sqrt.(diff_clim_pft_pred_var)[1]
                                        
                                    end
                                    
                                else
                                    
                                    out_3[it, 1:transitions_n, 1] = diff_clim_pft_pred
                                    #println(out1_delta)
                                    #if diff_clim_pft_pred_var .< 0
                                    #out_3[it,1:transitions_n,2] .= NaN
                                    
                                    #else
                                    diff_clim_pft_pred_var[diff_clim_pft_pred_var.<0] .= NaN
                                    out_3[it, 1:transitions_n, 2] =
                                    sqrt.(diff_clim_pft_pred_var)
                                    
                                    #end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    GC.gc()
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
 
- ```minDiffPxls```: Minimum number pixels in the moving window that must have different compositions. Must be any value in the interval 1 to winsize^2. By default minDiffPxls = 15.
 
- ```classes_vec```: A string vector with the names of the classes on ```cube_classes``` to be used. e.g. from MPI-BGC internal structure ```classes_vec = ["Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests", "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests", "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands", "Permanent_Wetlands", "Croplands", "Urban_and_Built-up_Lands", "Cropland/Natural_Vegetation_Mosaics", "Permanent_Snow_and_Ice", "Barren", "Water_Bodies"]```

- ```max_value```: Indicates if the scale of the presence of the discrete classes if from 0 to 1 or 0 to 100 if ```max_value = 100``` then the data is re-scaled from 0 to 1. By default ```max_value = 1```

- ```showprog```: Show progress bar. By default ```showprog = true```

- ```max_cache```: Size of the cache to allocate temporarily sections of the cubes. By default ```max_cache = 1e8```

 ## Output:
 The ```space4time_proc``` produces a YAXARRAY.Dataset with three cubes:
 - summary_mov_window cube has one axis ```summary_stat```, and three variables:
    - ```rsquared```:  
    - ```cumulative_variance```:
    - ```predicted```: Mean prediction of Z for moving window with the real combination of values.

 - ```metrics_for_classes``` cube has one axis ```Values of Z for pure classes```, and two variables:
    - ```estimated```:
    - ```estimated_error```:
 - metrics_for_transitions has two axis ```transitions``` (all the transitions by pairs between the different classes), and ```Differences``` with three variables:
    - ```delta```: delta of the biophysical produced of going from one class the another.
    - ```delta_error```:
    - ```coocurence```:
 """
function space4time_proc(
    cube_con,
    cube_classes,
    cube_altitude;
    time_axis_name = :Ti,
    lon_axis_name = :lon,
    lat_axis_name = :lat,
    classes_var_name = :classes,
    altitude_var_name = :variable,
    winsize = 5,
    minDiffPxls = 15,
    classes_vec = NaN,
    altitude_vec = NaN,
    minDiffPxls_alt = 15,
    max_value = 1,
    minpxl = 25,
    showprog = true,
    max_cache = 1e8,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end
    #println(size(cube_con))
    #println(size(cube_classes))
    # assuming the first dimmension is time.
    if !isnothing(time_axis_name)
        time_n = try
             length(lookup(cube_con, Dim{time_axis_name}).data)
        catch e
            length(lookup(cube_con, time_axis_name).data)
        end
        
        time_seq = try
            lookup(cube_con, Dim{time_axis_name}).data 
        catch e
            lookup(cube_con, time_axis_name).data
        end
    else
        time_n = nothing
    end
    # assuming that pfts presence change in time. last dimmension refer to the number of pfts.
    # pfts_cube = rand(Float32, (5,5, length(classes_vec)))

    # number of classes
    # assuming the last dimmension is PFTs
    nc = length(classes_vec)

    #sigma1_glob = [fill(NaN, (nc, nc)) for i = 1:Threads.nthreads()]

    #prederr_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]

    #predres_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]
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


    p1_static = range(0, 1, length = winsize^2)

    p2_static = reverse(p1_static)

    #out_pmindist_global = [zeros(1, winsize^2) for i = 1:Threads.nthreads()]

    denom = sum(sqrt.(sum.(eachrow([p1_static p2_static] .^ 2))))

    half = floor(Int, ceil((winsize^2) / 2))

    #localcomp_fix_glob = [rand(winsize^2) for i = 1:Threads.nthreads()]

    #pftsvarmat_f_glob = [rand(winsize^2, nc + 1) for i = 1:Threads.nthreads()]


    # 
    if !isnothing(time_axis_name)


        indims = InDims(
            time_axis_name,
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            time_axis_name,
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        indims_altitude = InDims(
            time_axis_name,
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            altitude_var_name,
            window_oob_value = NaN,
        
        )
        

        out_1_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "p-val_alt_"*altitude_vec[1], "p-val_alt_"*altitude_vec[2], "aic", "aicc", "model_used"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    else
        indims = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        indims_altitude = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            altitude_var_name,
            window_oob_value = NaN,
        
        )

        out_1_dims = OutDims(
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "p-val_"*altitude_vec[1], "p-val_"*altitude_vec[2], "aic", "aicc","model_used"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    end


    outdims = (out_1_dims, out_2_dims, out_3_dims)
    #println(out_3_dims)

    out_1, out_2, out_3 = mapCube(
        s4time,
        (cube_con, cube_classes, cube_altitude),
        indims = (indims, indims_classes, indims_altitude),
        outdims = outdims,
        max_cache = max_cache,
        showprog = showprog,
        include_loopvars = true;
        pft_list = classes_vec,
        time_n = time_n,
        max_value = max_value,
        p1_static,
        p2_static,
        #sigma1_glob,
        #prederr_glob,
        #predres_glob,
        minDiffPxls,
        tran_check,
        half,
        #localcomp_fix_glob,
        #pftsvarmat_f_glob,
        winsize = winsize,
        transitions_n = transitions_n,
        pftstrans_comb_names = pftstrans_comb_names,
        nc = nc,
        #out_pmindist_global = out_pmindist_global,
        denom = denom,
        minpxl = minpxl,
        minDiffPxls_alt = minDiffPxls_alt,
    )

    return Dataset(;
        summary_mov_window = out_1,
        metrics_for_classes = out_2,
        metrics_for_transitions = out_3,
    )

end


"""
 # Space for time processor (space chunks)
 
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
function space4time_proc_space(
    cube_con,
    cube_classes;
    time_axis_name = :Ti,
    lon_axis_name = :lon,
    lat_axis_name = :lat,
    classes_var_name = :classes,
    winsize = 5,
    minDiffPxlspercentage = 40,
    classes_vec = NaN,
    max_value = 1,
    minpxl = 25,
    showprog = true,
    max_cache = 1e8,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end
    #println(size(cube_con))
    #println(size(cube_classes))
    # assuming the first dimmension is time.
    if !isnothing(time_axis_name)
        time_n = try
             length(lookup(cube_con, Dim{time_axis_name}).data)
        catch e
            length(lookup(cube_con, time_axis_name).data)
        end
        
        time_seq = try
            lookup(cube_con, Dim{time_axis_name}).data 
        catch e
            lookup(cube_con, time_axis_name).data
        end
    else
        time_n = nothing
    end
    # assuming that pfts presence change in time. last dimmension refer to the number of pfts.
    # pfts_cube = rand(Float32, (5,5, length(classes_vec)))

    # number of classes
    # assuming the last dimmension is PFTs
    nc = length(classes_vec)

    #sigma1_glob = [fill(NaN, (nc, nc)) for i = 1:Threads.nthreads()]

    #prederr_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]

    #predres_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]
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


    p1_static = range(0, 1, length = winsize^2)

    p2_static = reverse(p1_static)

    #out_pmindist_global = [zeros(1, winsize^2) for i = 1:Threads.nthreads()]

    denom = sum(sqrt.(sum.(eachrow([p1_static p2_static] .^ 2))))

    half = floor(Int, ceil((winsize^2) / 2))

    minDiffPxls = (winsize^2 * minDiffPxlspercentage / 100)

    #localcomp_fix_glob = [rand(winsize^2) for i = 1:Threads.nthreads()]

    #pftsvarmat_f_glob = [rand(winsize^2, nc + 1) for i = 1:Threads.nthreads()]


    # 
    if !isnothing(time_axis_name)


        indims = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        out_1_dims = OutDims(
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "aic", "aicc"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    else
        indims = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        out_1_dims = OutDims(
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "aic", "aicc"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    end


    outdims = (out_1_dims, out_2_dims, out_3_dims)
    #println(out_3_dims)

    out_1, out_2, out_3 = mapCube(
        s4time_space,
        (cube_con, cube_classes),
        indims = (indims, indims_classes),
        outdims = outdims,
        max_cache = max_cache,
        showprog = showprog,
        include_loopvars = true;
        pft_list = classes_vec,
        max_value = max_value,
        p1_static,
        p2_static,
        #sigma1_glob,
        #prederr_glob,
        #predres_glob,
        minDiffPxls,
        tran_check,
        half,
        #localcomp_fix_glob,
        #pftsvarmat_f_glob,
        winsize = winsize,
        transitions_n = transitions_n,
        pftstrans_comb_names = pftstrans_comb_names,
        nc = nc,
        #out_pmindist_global = out_pmindist_global,
        denom = denom,
        minpxl = minpxl,
    )

    return Dataset(;
        summary_stats = out_1,
        metrics_for_classes = out_2,
        metrics_for_transitions = out_3,
    )

end



### previous function without elevation





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
    * "predicted": Mean prediction of Z for moving window with the real combination of values.
    out2: 
    
#Examples
    
    
"""
function s4time_old(
    out_1,
    out_2,
    out_3,
    clim_var_cube_in,
    pfts_cube_in,
    loopvars;
    pft_list::Vector{String},
    time_n,
    max_value::Int,
    p1_static,
    p2_static,
    #sigma1_glob,
    #prederr_glob,
    #predres_glob,
    minDiffPxls,
    tran_check,
    half,
    #localcomp_fix_glob,
    #pftsvarmat_f_glob,
    winsize = 5,
    transitions_n,
    pftstrans_comb_names,
    nc,
    #out_pmindist_global,
    denom,
    minpxl,
    
    )
    #println(size(clim_var_cube_in))
    #println(size(pfts_cube_in))
    #println(size(out_3))
    #igma1 = sigma1_glob[Threads.threadid()]
    sigma1 = fill(NaN, (nc, nc))
    #prederr = prederr_glob[Threads.threadid()]
    prederr = fill(NaN, nc) 
    #predres = predres_glob[Threads.threadid()]
    predres = fill(NaN, nc)
    #localcomp_fix = localcomp_fix_glob[Threads.threadid()]
    #pftsvarmat_f = pftsvarmat_f_glob[Threads.threadid()]
    #out_pmindist = out_pmindist_global[Threads.threadid()]
    #out_pmindist = zeros(1, winsize^2)
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
    
    if isnothing(time_n)
        climvarmat = reshape(clim_var_cube_in, (winsize^2))
        climvarmat = convert(Array{Float64}, climvarmat)
        
    else
        climvarmat = reshape(clim_var_cube_in, ((winsize^2), time_n))
        climvarmat = convert(Array{Float64}, climvarmat)
    end
    
    
    climvarmat = replace!(climvarmat, missing => NaN)
    
    
    local_pft1 =
    [findall(pft_list .== pftstrans_comb_names[comp][1]) for comp = 1:transitions_n]
    local_pft1 = reduce(vcat, local_pft1)
    local_pft2 =
    [findall(pft_list .== pftstrans_comb_names[comp][2]) for comp = 1:transitions_n]
    local_pft2 = reduce(vcat, local_pft2)
    
    #println("before if")
    #println(out_3[1,1,1])
    
    #println(local_pft1)
    #println(local_pft2)
    #pfts_cube_in_1 = replace!(pfts_cube_in, NaN => 0.0)
    pfts_cube_in_1 = pfts_cube_in
    #replace!(pfts_cube_in_1, missing => 0.0)
    #replace!(pfts_cube_in_1, NaN32 => 0.0)
    #replace!(pfts_cube_in_1, NaN16 => 0.0)
    pfts_cube_in_1 = convert(Array{Float64}, pfts_cube_in_1)
    
    
    if isnothing(time_n)
        
        #pfts_cube_in_1 =
        #permutedims(reshape(pfts_cube_in_1, (winsize, winsize, nc, 1)), (4, 1, 2, 3))
        
        for comp in eachindex(local_pft1)
            #println(comp)
            #println(size(pfts_cube_in))
            #println("I'm here")
            #println(pfts_cube_in[it,:,:,local_pft1[comp]])
            #println(pfts_cube_in[it,:,:,local_pft2[comp]])
            # define the pfts to be processed
            out_3[comp, 3] = coocufun(
            [0.0],
            pfts_cube_in_1[:, :, local_pft1[comp]],
            pfts_cube_in_1[:, :, local_pft2[comp]],
            p1_static,
            p2_static,
            denom,
            )
            #println("all good")
        end
        
        
        
        pfts_cube_in_2 = pfts_cube_in_1[:, :, :]
        #println(all(isnan, pfts_cube_in_2))
        
        #println(size(pfts_cube_in_2))
        
        pftsvarmat = reshape(pfts_cube_in_2, (winsize^2, nc))
        
        if count(!isnan, climvarmat[:]) >= minpxl
            
            if count(!isnan, climvarmat[:]) != 0
                
                pftsvarmat = pftsvarmat[findall(!isnan, climvarmat[:]), :]
                
                climvarmat = filter(!isnan, climvarmat)
            end
            
            if isfinite(sum(pftsvarmat)) && sum(sum(pftsvarmat, dims = 1) .> 0.) > 1
                #println("test")
                #println(any(isnan.(pftsvarmat)))
                # check if pftsvarmat is 0 to 1 or 0 to 100
                #println(maximum(vec(pftsvarmat)))
                
                # make sure compositions are really precisely right.
                
                #localcomp_fix_glob = mapslices(x->1-sum(x), pftsvarmat, dims = 2)
                localcomp_fix = map(x -> 1 - sum(x), eachslice(pftsvarmat, dims = 1))
                #map!(x->1-sum(x), localcomp_fix_glob, eachslice(pftsvarmat, dims = 1))
                #println(size(localcomp_fix_glob))
                #println(localcomp_fix_glob)
                pftsvarmat_f = [pftsvarmat localcomp_fix]
                
                map!((x) -> round(x, digits = 4), pftsvarmat_f, pftsvarmat_f)
                
                # some PFTs might not be present in the 5*5 window
                # these must be identified and removed, as they cannot be predicted
                
                #pftpres_check = vec(mapslices(sum, pftsvarmat, dims = 1) .> 0)
                pftpres_check = vec(sum(pftsvarmat_f, dims = 1) .> 0)
                
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
                        
                        temp = cumsum(lcsvd.S .^ 2) ./ sum(lcsvd.S .^ 2)
                        try
                            minimum(findall(temp .>= 1))
                        catch e
                            return
                        end
                        
                        ndim = minimum(findall(temp .>= 1))
                        
                        # store results to output object
                        # cumulative variance
                        out_1[2] = sum(lcsvd.S)
                        
                        #println(out7_cumulated_variance)
                        
                        # dimmensions that explain 100 % of the variance
                        
                        lr = lc2 * lcsvd.V[:, 1:ndim]
                        
                        
                        # create bogus composition dataset
                        
                        #boguscomp = zeros(Float64, nc+1, nc+1)
                        #boguscomp[diagind(boguscomp)] .= 1
                        
                        boguscomp = I(nc + 1)
                        # println(boguscomp)
                        
                        # remove absent pfts from bogus predictor compositions and close compositions.
                        #
                        bogusc1 =
                        mapslices(
                        x -> x / sum(x),
                        boguscomp[pftpres_check, pftpres_check],
                        dims = 2,
                        )'
                        
                        #println(size(bogusc1))
                        #println(pftpres_check)
                        
                        # center the columns as the training data were centered
                        
                        bogusc2 = (I(sum(pftpres_check)) .- lCm[pftpres_check])
                        
                        bogusc2 = (bogusc1' .- lCm[pftpres_check])
                        
                        bogusc3 = bogusc2 * lcsvd.V[pftpres_check, 1:ndim]
                        #println(climvarmat[:,it])
                        #println("test")
                        # data = hcat(DataFrame(lt = convert(Vector{Float64}, climvarmat[:,it])), DataFrame(lr, :auto))
                        
                        # compreg = GLM.lm(Term(:lt) ~ sum(Term.(Symbol.(names(data[:, Not(:lt)])))), data)
                        
                        #println("before fail")
                        
                        ols = lm([ones(size(lr, 1)) lr], identity.(climvarmat[:]); method=:qr, dropcollinear = false)
                        
                        # continue only if there are no NA in the estimated coefficients
                        
                        coef_reg = GLM.coef(ols)
                        #println("original coef $coef_reg")
                        #println("second estimation coef $(lr\climvarmat[:,it])")
                        
                        if isfinite(sum(coef_reg))
                            # then do predictions for the log-normal approach
                            if isa(bogusc2, Vector)
                                
                                boguspred = predict(
                                ols,
                                [ones(length(bogusc3)) bogusc3],
                                )
                                # boguspred = GLM.predict(compreg, DataFrame( x1 = bogusc3))
                                
                            else
                                # boguspred = GLM.predict(compreg, DataFrame(bogusc3, :auto))
                                boguspred = predict(
                                ols,
                                [ones(size(bogusc3, 1)) bogusc3],
                                )
                                
                            end
                            
                            
                            x2pred = [ones(size(bogusc3, 1), 1) bogusc3]
                            
                            vcv = GLM.vcov(ols)
                            # vcv = GLM.vcov(compreg)
                            
                            sigma = x2pred * vcv * x2pred'
                            
                            # now store the target variables
                            # but make sure appropiate temperatures go back to appropiate
                            # pfts (as absent pfts were removed)
                            
                            # value of climatevar for pure ptfs
                            
                            predres .= NaN
                            #println(boguspred)
                            predres[view(pftpres_check, 1:nc)] = boguspred
                            #println(it)
                            
                            
                            out_2[:, 1] = predres
                            
                            
                            prederr .= NaN
                            
                            prederr[view(pftpres_check, 1:nc)] = sqrt.(diag(sigma))
                            
                            
                            out_2[:, 2] = prederr
                            
                            
                            # prediction of varclim for the central pixel with its real pft combination
                            
                            out_1[3] = mean(StatsModels.predict(ols))
                            
                            # Rsquare of the regression
                            
                            out_1[1] = adjr2(ols)
                            
                            out_1[4] = aic(ols)
                            out_1[5] = aicc(ols)
                            # println(out_1)
                            # println(r2(compreg))
                            
                            # and now for the transitions
                            # only the PFTs identified in the pftlist are to be used
                            
                            sigma1 .= NaN
                            sigma1[view(pftpres_check, 1:nc), view(pftpres_check, 1:nc)] =
                            sigma
                            
                            # calculate the difference on climatevar caused by going from one pft to another.
                            
                            diff_clim_pft_pred = round.((predres .- predres'), digits = 10)
                            
                            diff_clim_pft_pred = diff_clim_pft_pred[tran_check]
                            
                            # propagate the error (as variances) taking into account
                            # the covariance terms
                            # in the original implementation diffclim_pft_var = dZvar
                            
                            diff_clim_pft_pred_var =
                            round.(
                            (diag(sigma1) .+ diag(sigma1)') .- 2 * sigma1,
                            digits = 10,
                            )
                            
                            #print(diff_clim_pft_pred)
                            # flag out those with zero error (may occur with identical compositions for 2 pfts)
                            
                            diff_clim_pft_pred_var = diff_clim_pft_pred_var[tran_check]
                            
                            #println(diff_clim_pft_pred)
                            
                            # flag out those with zero error (may occur with identical compositions for 2 pfts)
                            
                            if any(round.(diff_clim_pft_pred, digits = 8) .== 0.0)
                                diff_clim_pft_pred[round.(
                                diff_clim_pft_pred,
                                digits = 8,
                                ).==0] .= NaN
                            end
                            
                            if any(round.(diff_clim_pft_pred_var, digits = 8) .== 0.0)
                                diff_clim_pft_pred_var[round.(
                                diff_clim_pft_pred_var;
                                digits = 8,
                                ).==0] .= NaN
                            end
                            
                            # mask out low co-ocurrence ask to greg!!! This can be performed  masking the pixels
                            # before all estimations
                            
                            #println("inside loop")
                            #println(diff_clim_pft_pred)
                            
                            if length(diff_clim_pft_pred) == 1
                                
                                out_3[1, 1] = diff_clim_pft_pred[1]
                                #println(out1_delta)
                                if diff_clim_pft_pred_var[1] .< 0
                                    out_3[1, 2] = NaN
                                else
                                    out_3[1, 2] = sqrt.(diff_clim_pft_pred_var)[1]
                                    
                                end
                                
                                
                            else
                                
                                out_3[1:transitions_n, 1] = diff_clim_pft_pred
                                #println(out1_delta)
                                #if diff_clim_pft_pred_var .< 0
                                #out_3[it,1:transitions_n,2] .= NaN
                                
                                #else
                                diff_clim_pft_pred_var[diff_clim_pft_pred_var.<0] .= NaN
                                out_3[1:transitions_n, 2] =
                                sqrt.(diff_clim_pft_pred_var)
                                
                                #end
                            end
                        end
                    end
                end
            end
        end
        # for debug from R -------
        
        # pftsvarmat = Matrix(CSV.read("/Net/Groups/BGI/people/dpabon/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
        # pftsvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
        # climvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_temperature_example_from_R.csv", DataFrame))
        # check that there are not NaN values on pfts and at least one pft is present
        # @show sum(pftsvarmat), sum(pftsvarmat)
        #println(pftsvarmat)
    else
        
        if time_n == 1
            
            pfts_cube_in_1 =
            permutedims(reshape(pfts_cube_in_1, (winsize, winsize, nc, 1)), (4, 1, 2, 3))
            
        end
        
        for it = 1:time_n
            # co-ocurrence estimation
            
            #println(it)
            #println(local_pft1)
            #println(local_pft2)
            #println(transitions_n)
            for comp in eachindex(local_pft1)
                #println(comp)
                #println(size(pfts_cube_in))
                #println("I'm here")
                #println(pfts_cube_in[it,:,:,local_pft1[comp]])
                #println(pfts_cube_in[it,:,:,local_pft2[comp]])
                # define the pfts to be processed
                out_3[it, comp, 3] = coocufun(
                [0.0],
                pfts_cube_in_1[it, :, :, local_pft1[comp]],
                pfts_cube_in_1[it, :, :, local_pft2[comp]],
                p1_static,
                p2_static,
                denom,
                )
                #println("all good")
            end
            
            
            
            pfts_cube_in_2 = pfts_cube_in_1[it, :, :, :]
            #println(all(isnan, pfts_cube_in_2))
            
            #println(size(pfts_cube_in_2))
            
            pftsvarmat = reshape(pfts_cube_in_2, (winsize^2, nc))
            
            
            if count(!isnan, climvarmat[:, it]) >= minpxl
                
                if count(!isnan, climvarmat[:, it]) != 0
                    
                    pftsvarmat = pftsvarmat[findall(!isnan, climvarmat[:, it]), :]
                    
                    climvarmat_it = filter(!isnan, climvarmat[:,it])
                end

                # for debug from R -------
            
                # pftsvarmat = Matrix(CSV.read("/Net/Groups/BGI/people/dpabon/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
                # pftsvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_composition_example_from_R.csv", DataFrame))
                # climvarmat = Matrix(CSV.read("/home/dpabon/Nextcloud/nfdi4earth_oemc/data/local_temperature_example_from_R.csv", DataFrame))
                # check that there are not NaN values on pfts and at least one pft is present
                # @show sum(pftsvarmat), sum(pftsvarmat)
                #println(pftsvarmat)
            
            
                #println(it)
                
                if isfinite(sum(pftsvarmat)) && sum(sum(pftsvarmat, dims = 1) .> 0.) > 1
                    #println("test")
                    #println(any(isnan.(pftsvarmat)))
                    # check if pftsvarmat is 0 to 1 or 0 to 100
                    #println(maximum(vec(pftsvarmat)))
                    
                    # make sure compositions are really precisely right.
                    
                    #localcomp_fix_glob = mapslices(x->1-sum(x), pftsvarmat, dims = 2)
                    localcomp_fix = map(x -> 1 - sum(x), eachslice(pftsvarmat, dims = 1))
                    #map!(x->1-sum(x), localcomp_fix_glob, eachslice(pftsvarmat, dims = 1))
                    #println(size(localcomp_fix_glob))
                    #println(localcomp_fix_glob)
                    pftsvarmat_f = [pftsvarmat localcomp_fix]
                    
                    map!((x) -> round(x, digits = 4), pftsvarmat_f, pftsvarmat_f)
                    
                    # some PFTs might not be present in the 5*5 window
                    # these must be identified and removed, as they cannot be predicted
                    
                    #pftpres_check = vec(mapslices(sum, pftsvarmat, dims = 1) .> 0)
                    pftpres_check = vec(sum(pftsvarmat_f, dims = 1) .> 0)
                    
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
                        #println("test")
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
                            
                            temp = round.(cumsum(lcsvd.S .^ 2) ./ sum(lcsvd.S .^ 2), digits = 8)
                            ndim = minimum(findall(temp .>= 1))
                            
                            # store results to output object
                            # cumulative variance
                            out_1[it, 2] = sum(lcsvd.S)
                            
                            #println(out7_cumulated_variance)
                            
                            # dimmensions that explain 100 % of the variance
                            
                            lr = lc2 * lcsvd.V[:, 1:ndim]
                            
                            
                            # create bogus composition dataset
                            
                            #boguscomp = zeros(Float64, nc+1, nc+1)
                            #boguscomp[diagind(boguscomp)] .= 1
                            
                            boguscomp = I(nc + 1)
                            # println(boguscomp)
                            
                            # remove absent pfts from bogus predictor compositions and close compositions.
                            #
                            bogusc1 =
                            mapslices(
                            x -> x / sum(x),
                            boguscomp[pftpres_check, pftpres_check],
                            dims = 2,
                            )'
                            
                            #println(size(bogusc1))
                            #println(pftpres_check)
                            
                            # center the columns as the training data were centered
                            
                            bogusc2 = (I(sum(pftpres_check)) .- lCm[pftpres_check])
                            
                            bogusc2 = (bogusc1' .- lCm[pftpres_check])
                            
                            bogusc3 = bogusc2 * lcsvd.V[pftpres_check, 1:ndim]
                            #println(climvarmat[:,it])
                            
                            #println("test")
                            # data = hcat(DataFrame(lt = convert(Vector{Float64}, climvarmat[:,it])), DataFrame(lr, :auto))
                            
                            # compreg = GLM.lm(Term(:lt) ~ sum(Term.(Symbol.(names(data[:, Not(:lt)])))), data)
                            
                            #println("before fail")
                            
                            ols = lm([ones(size(lr, 1)) lr], identity.(climvarmat_it[:]); method=:qr, dropcollinear = false)
                            
                            # continue only if there are no NA in the estimated coefficients
                            
                            coef_reg = GLM.coef(ols)
                            #println("original coef $coef_reg")
                            #println("second estimation coef $(lr\climvarmat[:,it])")
                            
                            if isfinite(sum(coef_reg))
                                # then do predictions for the log-normal approach
                                if isa(bogusc2, Vector)
                                    
                                    boguspred = predict(
                                    ols,
                                    [ones(length(bogusc3)) bogusc3],
                                    )
                                    # boguspred = GLM.predict(compreg, DataFrame( x1 = bogusc3))
                                    
                                else
                                    # boguspred = GLM.predict(compreg, DataFrame(bogusc3, :auto))
                                    boguspred = predict(
                                    ols,
                                    [ones(size(bogusc3, 1)) bogusc3],
                                    )
                                    
                                end
                                
                                
                                x2pred = [ones(size(bogusc3, 1), 1) bogusc3]
                                
                                vcv = GLM.vcov(ols)
                                # vcv = GLM.vcov(compreg)
                                
                                sigma = x2pred * vcv * x2pred'
                                
                                # now store the target variables
                                # but make sure appropiate temperatures go back to appropiate
                                # pfts (as absent pfts were removed)
                                
                                # value of climatevar for pure ptfs
                                
                                predres .= NaN
                                #println(boguspred)
                                predres[view(pftpres_check, 1:nc)] = boguspred
                                #println(it)
                                
                                out_2[it, :, 1] = predres
                                
                                
                                prederr .= NaN
                                
                                prederr[view(pftpres_check, 1:nc)] = sqrt.(diag(sigma))
                                
                                
                                out_2[it, :, 2] = prederr
                                
                                
                                # prediction of varclim for the central pixel with its real pft combination
                                # It is possible that the center pixel correspond to an observation with NaN that was removed then
                                # I redefined as just mean of the prediction
                                
                                out_1[it, 3] = mean(StatsModels.predict(ols))
                                
                                # Rsquare of the regression
                                
                                out_1[it, 1] = adjr2(ols)
                                out_1[it,4] = aic(ols)
                                out_1[it,5] = aicc(ols)

                                # println(out_1)
                                # println(r2(compreg))
                                
                                # and now for the transitions
                                # only the PFTs identified in the pftlist are to be used
                                
                                sigma1 .= NaN
                                sigma1[view(pftpres_check, 1:nc), view(pftpres_check, 1:nc)] =
                                sigma
                                
                                # calculate the difference on climatevar caused by going from one pft to another.
                                
                                diff_clim_pft_pred = round.((predres .- predres'), digits = 10)
                                
                                diff_clim_pft_pred = diff_clim_pft_pred[tran_check]
                                
                                # propagate the error (as variances) taking into account
                                # the covariance terms
                                # in the original implementation diffclim_pft_var = dZvar
                                
                                diff_clim_pft_pred_var =
                                round.(
                                (diag(sigma1) .+ diag(sigma1)') .- 2 * sigma1,
                                digits = 10,
                                )
                                
                                #print(diff_clim_pft_pred)
                                # flag out those with zero error (may occur with identical compositions for 2 pfts)
                                
                                diff_clim_pft_pred_var = diff_clim_pft_pred_var[tran_check]
                                
                                #println(diff_clim_pft_pred)
                                
                                # flag out those with zero error (may occur with identical compositions for 2 pfts)
                                
                                if any(round.(diff_clim_pft_pred, digits = 8) .== 0.0)
                                    diff_clim_pft_pred[round.(
                                    diff_clim_pft_pred,
                                    digits = 8,
                                    ).==0] .= NaN
                                end
                                
                                if any(round.(diff_clim_pft_pred_var, digits = 8) .== 0.0)
                                    diff_clim_pft_pred_var[round.(
                                    diff_clim_pft_pred_var;
                                    digits = 8,
                                    ).==0] .= NaN
                                end
                                
                                # mask out low co-ocurrence ask to greg!!! This can be performed  masking the pixels
                                # before all estimations
                                
                                #println("inside loop")
                                #println(diff_clim_pft_pred)
                                
                                if length(diff_clim_pft_pred) == 1
                                    
                                    out_3[it, 1, 1] = diff_clim_pft_pred[1]
                                    #println(out1_delta)
                                    if diff_clim_pft_pred_var[1] .< 0
                                        out_3[it, 1, 2] = NaN
                                    else
                                        out_3[it, 1, 2] = sqrt.(diff_clim_pft_pred_var)[1]
                                        
                                    end
                                    
                                    
                                else
                                    
                                    out_3[it, 1:transitions_n, 1] = diff_clim_pft_pred
                                    #println(out1_delta)
                                    #if diff_clim_pft_pred_var .< 0
                                    #out_3[it,1:transitions_n,2] .= NaN
                                    
                                    #else
                                    diff_clim_pft_pred_var[diff_clim_pft_pred_var.<0] .= NaN
                                    out_3[it, 1:transitions_n, 2] =
                                    sqrt.(diff_clim_pft_pred_var)
                                    
                                    #end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

"""
space4time_space(climate_cube, pfts_cube, pft_list::Vector{String}, winsize = 5, minpxl = 100, minDiffPxlspercentage = 40)

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
    * "predicted": Mean prediction of Z for moving window with the real combination of values.
    out2: 
    
#Examples
    
    
"""

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
    - ```predicted```: Mean prediction of Z for moving window with the real combination of values.

 - ```metrics_for_classes``` cube has one axis ```Values of Z for pure classes```, and two variables:
    - ```estimated```:
    - ```estimated_error```:
 - metrics_for_transitions has two axis ```transitions``` (all the transitions by pairs between the different classes), and ```Differences``` with three variables:
    - ```delta```: delta of the biophysical produced of going from one class the another.
    - ```delta_error```:
    - ```coocurence```:
 """
function space4time_proc_old(
    cube_con,
    cube_classes;
    time_axis_name = :Ti,
    lon_axis_name = :lon,
    lat_axis_name = :lat,
    classes_var_name = :classes,
    winsize = 5,
    minDiffPxlspercentage = 40,
    classes_vec = NaN,
    max_value = 1,
    minpxl = 25,
    showprog = true,
    max_cache = 1e8,
)

    # Checking that winsize is odd

    if isodd(winsize)
        pre_step = after_step = floor(winsize / 2)
    else
        pre_step = after_step = floor(winsize / 2) - 1

        @warn "Window size is not odd. Going on however... windowsize = $(winsize - 1)"
    end
    #println(size(cube_con))
    #println(size(cube_classes))
    # assuming the first dimmension is time.
    if !isnothing(time_axis_name)
        time_n = try
             length(lookup(cube_con, Dim{time_axis_name}).data)
        catch e
            length(lookup(cube_con, time_axis_name).data)
        end
        
        time_seq = try
            lookup(cube_con, Dim{time_axis_name}).data 
        catch e
            lookup(cube_con, time_axis_name).data
        end
    else
        time_n = nothing
    end
    # assuming that pfts presence change in time. last dimmension refer to the number of pfts.
    # pfts_cube = rand(Float32, (5,5, length(classes_vec)))

    # number of classes
    # assuming the last dimmension is PFTs
    nc = length(classes_vec)

    #sigma1_glob = [fill(NaN, (nc, nc)) for i = 1:Threads.nthreads()]

    #prederr_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]

    #predres_glob = [fill(NaN, nc) for i = 1:Threads.nthreads()]
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


    p1_static = range(0, 1, length = winsize^2)

    p2_static = reverse(p1_static)

    #out_pmindist_global = [zeros(1, winsize^2) for i = 1:Threads.nthreads()]

    denom = sum(sqrt.(sum.(eachrow([p1_static p2_static] .^ 2))))

    half = floor(Int, ceil((winsize^2) / 2))

    minDiffPxls = (winsize^2 * minDiffPxlspercentage / 100)

    #localcomp_fix_glob = [rand(winsize^2) for i = 1:Threads.nthreads()]

    #pftsvarmat_f_glob = [rand(winsize^2, nc + 1) for i = 1:Threads.nthreads()]


    # 
    if !isnothing(time_axis_name)


        indims = InDims(
            time_axis_name,
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            time_axis_name,
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        out_1_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "aic", "aicc"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{time_axis_name}(time_seq),
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    else
        indims = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            window_oob_value = NaN,
        )

        indims_classes = InDims(
            MovingWindow(lon_axis_name, pre_step, after_step),
            MovingWindow(lat_axis_name, pre_step, after_step),
            classes_var_name,
            window_oob_value = NaN,
        )

        out_1_dims = OutDims(
            Dim{:summary_stat}(["rsquared_adjusted", "cumulative_variance", "predicted", "aic", "aicc"]),
        )

        # Values of clim_var (z) for pure PFTs
        out_2_dims = OutDims(
            Dim{:classes}(classes_vec),
            Dim{:values_of_Z_for_pure_classes}(["estimated", "estimated_error"]),
        )
        #println([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)])
        # delta of clim_var produced by the transitions between PFTs
        out_3_dims = OutDims(
            Dim{:transitions}([join(pftstrans_comb_names[i], " to ") for i in eachindex(pftstrans_comb_names)]),
            Dim{:differences}(["delta", "delta_error", "coocurence"]),
        )

    end


    outdims = (out_1_dims, out_2_dims, out_3_dims)
    #println(out_3_dims)

    out_1, out_2, out_3 = mapCube(
        s4time_old,
        (cube_con, cube_classes),
        indims = (indims, indims_classes),
        outdims = outdims,
        max_cache = max_cache,
        showprog = showprog,
        include_loopvars = true;
        pft_list = classes_vec,
        time_n = time_n,
        max_value = max_value,
        p1_static,
        p2_static,
        #sigma1_glob,
        #prederr_glob,
        #predres_glob,
        minDiffPxls,
        tran_check,
        half,
        #localcomp_fix_glob,
        #pftsvarmat_f_glob,
        winsize = winsize,
        transitions_n = transitions_n,
        pftstrans_comb_names = pftstrans_comb_names,
        nc = nc,
        #out_pmindist_global = out_pmindist_global,
        denom = denom,
        minpxl = minpxl,
    )

    return Dataset(;
        summary_stats = out_1,
        metrics_for_classes = out_2,
        metrics_for_transitions = out_3,
    )

end

