module WriteReport

export write_report 

using ..Errors: compute_relative_error, compute_relative_derivative_error

"""
    write_report(filename, b, reconstructed_function, end_point, n; n_dense::Int = 5000)

Write a text report `<filename>.txt` with relative errors for `b` and `b'`.

Errors are evaluated on a uniform grid of size `n_dense` over `[0, end_point]`.
The value `n` is the experiment grid size and is
recorded in the report.
"""
function write_report(filename, b, reconstructed_function, end_point, n; n_dense::Int = 5000)
    path = joinpath(filename * ".txt")

    # Errors for b
    L1_norm     = compute_relative_error("1",   b, reconstructed_function, end_point, n_dense)
    L2_norm     = compute_relative_error("2",   b, reconstructed_function, end_point, n_dense)
    L_inf_norm  = compute_relative_error("inf", b, reconstructed_function, end_point, n_dense)

    # Errors for b'
    L1_der_norm    = compute_relative_derivative_error("1",   b, reconstructed_function, end_point, n_dense)
    L2_der_norm    = compute_relative_derivative_error("2",   b, reconstructed_function, end_point, n_dense)
    L_inf_der_norm = compute_relative_derivative_error("inf", b, reconstructed_function, end_point, n_dense)

    open(path, "w") do io
        println(io, "# Report\n")
        println(io, "n (experiment grid): ", n)
        println(io, "n_dense (error eval grid): ", n_dense)
        println(io, "end_point: ", end_point)
        println(io, "\n===============\n")

        println(io, "L1-norm: ", L1_norm)
        println(io, "L2-norm: ", L2_norm)
        println(io, "L_inf-norm: ", L_inf_norm)

        println(io, "\n===============\n")

        println(io, "L1-der-norm: ", L1_der_norm)
        println(io, "L2-der-norm: ", L2_der_norm)
        println(io, "L_inf-der-norm: ", L_inf_der_norm)
    end

    return path
end

end # module
