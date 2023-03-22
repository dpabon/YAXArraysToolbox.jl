function dates_builder(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2]))
    end

    return out
end
