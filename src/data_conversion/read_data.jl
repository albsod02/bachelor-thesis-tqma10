module ConvertData

export normalize_data, read_data_1D

"""
    normalize_data(values::AbstractVector{<:Real}) -> Vector{Float64}

Normalize `values` linearly to the range `[0, 1]`.

If all entries are equal, returns a vector of zeros.

# Arguments
- `values`: A vector of real numbers.

# Returns
- A `Vector{Float64}` scaled to `[0, 1]`.
"""
function normalize_data(values::AbstractVector{<:Real})
    minval = minimum(values)
    maxval = maximum(values)
    maxval == minval && return zeros(Float64, length(values))
    return (Float64.(values) .- minval) ./ (maxval - minval)
end

"""
    read_data_1D(filepath::AbstractString) -> Tuple{Vector{Float64}, Vector{Float64}}

Read 1D data from a `.txt` file containing x-values and y-values separated by the
delimiter line `# y values`.

File format:
- The first 3 lines are ignored (e.g. headers).
- Lines up to (but not including) `# y values` are x-values (one per line).
- Lines after `# y values` are y-values (one per line).

# Arguments
- `filepath`: Path to the `.txt` file.

# Returns
- `(x_values, y_values)` as `Vector{Float64}`.

# Throws
- An error if the delimiter line `# y values` is missing.
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
