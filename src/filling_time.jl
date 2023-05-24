function filling_time(xout, xin)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(xin)
            if !ismissing(xin[i])
                xout[:, i] .= xin[i]
            end
        end
    end
end
