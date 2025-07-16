module SettingsParser

export parse_settings

function parse_settings(path)
    settings = Dict{String, Any}()
    for line in eachline(path)
        line = strip(line)
        isempty(line) || startswith(line, "#") && continue
        parts = split(line, "=", limit=2)
        length(parts) < 2 && continue  # Skip lines that don’t contain "="
        key, val = parts
        key = strip(key)
        val = strip(val)
        
        # Try to convert value to number or bool
        if val == "true"
            val = true
        elseif val == "false"
            val = false
        elseif occursin(".", val)
            val = parse(Float64, val)
        elseif all(isdigit, val)
            val = parse(Int, val)
        end

        settings[key] = val
    end
    return settings
end


end # module