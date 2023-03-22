using Statistics

function median_by_index(xout, xin; index_list = time_to_index)
    #@show size(xin)
    #@show typeof(xin)
    xout .= NaN
    if !all(isnan, xin)
        for i in eachindex(index_list)
            if !all(isnan, xin[index_list[i]])
                xout[i] = median(filter(!isnan, xin[index_list[i]]))
            end
        end
    end    
end 