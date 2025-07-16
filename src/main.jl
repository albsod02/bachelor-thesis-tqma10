module Script

include("spline_interpolation/BSplineInterpolation.jl")
include("synthetic_data_generation/write_bathymetric_data.jl")
include("write_report/write_report.jl")
include("parse_settings.jl")


using TrixiBottomTopography
using CairoMakie

using .SplineInterpolation
using .SyntheticData
using .SettingsParser
using .WriteReport


function main()

    # Parameters for calculation
    settings = parse_settings(joinpath("../config", "settings" * ".txt"))

    filename = settings["filename"]
    path = joinpath("../data", "synthetic", "1D", filename * ".txt")

    # 1D smooth test function
    f(x) = (1 - 1 / (1 + exp(-30 * (x - 0.5)))) + 0.05 * (1 - 1 / (1 + exp(-40 * (x - 0.55)))) * sin(2 * pi * 15 * x)

    # Piecewise smooth 1D test function at x = 0.4
    #f(x) = x < 0.4 ? 0.6x + 0.05 * sin(10 * pi * x) : 1 - 0.4 * (1 - x) + 0.03 * cos(8 * pi * x)

    # Piecewise smooth 1D test function with disc. at x = 0.3 and x = 0.7
    #f(x) = x < 0.3 ? 0.6x + 0.05 * sin(10π * x) : x < 0.7 ? 1 - 0.4 * (x - 0.3) + 0.03 * cos(8π * x) : 0.5 + 0.2 * sin(12π * x)



    # Divide interval [0, end_point] in n sub-intervals
    end_point = settings["end_point"]
    n = settings["n"]

    # Add optional noise to the data defined by f
    add_noise = settings["add_noise"]

    alpha = settings["alpha"]
    variance = 1 / 10 ^ alpha

    # Get data defined from synthetic test function
    x, b = define_test_function(f; end_point, n, add_noise, variance)

    # Write bathymetric data to .txt file determined by variable path, if path already exist locally it will be overwritten
    write_bathymetric_data(path, x, b)  

    # Settings for interpolation
    end_condition = settings["end_condition"]
    smoothing_factor = settings["smoothing_factor"]
    add_knots = settings["add_knots"]

    interpolation_function = plot_interpolation_with_original_function(path, f, x, b, n, end_point, smoothing_factor, end_condition; add_noise, variance, add_knots)

    filename = "report"

    # Write report with error for the interpolation
    write_report(filename, f, interpolation_function, end_point, n)

    plot_error_graph(path, f, x, b, n, end_point, smoothing_factor, end_condition, alpha)

    println("Done")
end


end # module 

# Only run if the script is executed directly, not included
if abspath(PROGRAM_FILE) == @__FILE__
    Script.main()
end