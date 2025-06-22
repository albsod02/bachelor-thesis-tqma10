module Script

include("spline_interpolation/BSplineInterpolation.jl")
include("synthetic_data_generation/write_bathymetric_data.jl")
include("parse_settings.jl")

using TrixiBottomTopography
using CairoMakie

using .ConstructBSpline
using .SyntheticData
using .SettingsParser


function write_report(filename, f, interpolation_function, end_point, n)
    path = joinpath(filename * ".txt")

    L1_norm = interpolation_error("1", f, interpolation_function, end_point, n)
    L2_norm = interpolation_error("2", f, interpolation_function, end_point, n)
    L_inf_norm = interpolation_error("inf", f, interpolation_function, end_point, n)

    open(path, "w") do io

        println(io, "# Report\n")

        # Error with respect to different norms
        println(io, "L1-norm: ", L1_norm)
        println(io, "L2-norm: ", L2_norm)
        println(io, "L_inf-norm: ", L_inf_norm)
    end
end


function main()

    # Parameters for calculation
    settings = parse_settings(joinpath("../config", "settings" * ".txt"))

    filename = settings["filename"]
    path = joinpath("../data", "synthetic", "1D", filename * ".txt")

    # Smooth test function for generating synthetic bathymetric data in 1D
    #f(x) = 50 * sin(0.2 * x) + 10 * sin(1.5 * x) - 5 * exp(-((x - 15)^2) / 5) + 20
    #f(x) = 40 * cos(0.1 * x) - 15 * sin(0.8 * x) - 10 * exp(-((x - 25)^2) / 8) + 30
    #f(x) = sin(0.5 * x) + 10 * exp(-2 * x) + 0.1 * x + 0.1 * exp(0.1 * x) + 2
    #f(x) = 50 * sin(0.1 * x) - 40 * exp(-((x - 20)^2) / 0.5) + 10
    #f(x) = 30 * tanh(5 * (x - 10)) - 20 * exp(-0.2 * x) + 15
    #f(x) = 25 * sin(0.2 * x) + 100 * exp(-((x - 30)^2) / 0.1) - 5
    #f(x) = 10 * sin(0.2 * x) + 15 * sin(5 * x) - 10 * exp(-((x - 25)^2) / 4) + 5
    f(x) = 60 / (1 + exp(-(x - 12))) - 60 / (1 + exp(-(x - 18))) + 25

    # Divide interval [0, end_point] in n sub-intervals
    end_point = settings["end_point"]
    n = settings["n"]

    # Add optional noise to the data defined by f
    add_noise = settings["add_noise"]
    variance = settings["variance"]

    # Get data defined from synthetic test function
    x, b = define_test_function(f; end_point, n, add_noise, variance)

    # Write bathymetric data to .txt file determined by variable path, if path already exist locally it will be overwritten
    write_bathymetric_data(path, x, b)  

    # Settings for interpolation
    end_condition = settings["end_condition"]
    smoothing_factor = settings["smoothing_factor"]
    add_knots = settings["add_knots"]

    interpolation_function = plot_interpolation_with_original_function(path, f, x, b, n, end_point, smoothing_factor, end_condition, add_noise, variance, add_knots)

    filename = "report"

    # Write report with error for the interpolation
    write_report(filename, f, interpolation_function, end_point, n)
end


end # module 

# Only run if the script is executed directly, not included
if abspath(PROGRAM_FILE) == @__FILE__
    Script.main()
end