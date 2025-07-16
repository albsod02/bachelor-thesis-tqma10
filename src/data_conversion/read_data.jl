module ConvertData

export normalize_data, read_data_1D, read_data_2D

"""
    normalize_data(values::AbstractVector{<:Real}) -> Vector{Float64}

Normalize the input numeric vector to the range [0, 1].

# Arguments
- `values`: A vector of real numbers.

# Returns
- A vector of `Float64` values scaled linearly to the interval [0, 1].
"""
function normalize_data(values)
    min = minimum(values)
    max = maximum(values)
    return (values .- min) ./ (max - min)
end

"""
    read_data_1D(filepath::AbstractString) -> Tuple{Vector{Float64}, Vector{Float64}}

Read and normalize 1D data from a `.txt` file. The file should follow this format:
- First 3 lines are ignored (e.g., headers).
- Lines until delimeter `# y values` contain x-values.
- Lines after delimeter `# y values` contain y-values.

# Arguments
- `filepath`: Path to the `.txt` file containing the 1D data.

# Returns
- A tuple `(x_values, y_values)`, where both are vectors of normalized `Float64` values.
"""
function read_data_1D(filepath)
    lines = readlines(filepath)[4:end] # Skip first 3 lines
    split_index = findfirst(==("# y values"), lines)

    if split_index === nothing
        error("Missing '# y values' delimiter in data.")
    end

    x_values = parse.(Float64, lines[1:split_index-1])
    y_values = parse.(Float64, lines[split_index+1:end])
    
    return x_values, y_values
end

"""
    read_data_2D(filepath::AbstractString) -> Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}}

Placeholder function to read and return 2D data from a `.txt` file, where the z-values represent a function `b(x, y)`.

# Arguments
- `filepath`: Path to the `.txt` file containing the 2D data.

# Returns
- A tuple `(x_values, y_values, z_values)` where:
  - `x_values`: Vector of `Float64` representing x-axis points.
  - `y_values`: Vector of `Float64` representing y-axis points.
  - `z_values`: 2D matrix of `Float64` representing b(x, y) values.

*Note: Implementation needed.*
"""
function read_data_2D(filepath)
    # Implementation placeholder
end

end # module ConvertData
