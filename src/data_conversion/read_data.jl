module ConvertData

export normalize_data, read_data_1D

"""
    normalize_data(values::AbstractVector{<:Real}) -> Vector{Float64}

Normalize `values` linearly to the range `[0, 1]`.
If all entries are equal, returns a vector of zeros.
"""
function normalize_data(values::AbstractVector{<:Real})
    minval = minimum(values)
    maxval = maximum(values)
    maxval == minval && return zeros(Float64, length(values))
    return (Float64.(values) .- minval) ./ (maxval - minval)
end

"""
    read_data_1D(filepath::AbstractString) -> x_values, y_values

Read 1D data from a `.txt` file where x-values come first and y-values start after
the delimiter line `# y values`. The first 3 lines are ignored.
"""
function read_data_1D(filepath::AbstractString)
    lines = readlines(filepath)[4:end]
    split_index = findfirst(==("# y values"), lines)
    split_index === nothing && error("Missing '# y values' delimiter in data.")

    x_values = parse.(Float64, lines[1:split_index-1])
    y_values = parse.(Float64, lines[split_index+1:end])
    return x_values, y_values
end

end # module ConvertData
