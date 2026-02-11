module SettingsParser

export parse_settings

"""
    parse_settings(path) -> Dict{String, Any}

Parse a simple `key = value` settings file into a dictionary.

- Empty lines and lines starting with `#` are ignored.
- Lines without `=` are skipped.
- Values are converted when possible:
  - `"true"` / `"false"` → `Bool`
  - strings containing `.` → `Float64`
  - digit-only strings → `Int`
  - otherwise kept as `String`
"""
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
