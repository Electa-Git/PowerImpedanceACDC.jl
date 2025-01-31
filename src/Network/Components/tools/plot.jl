export bode, save_plot

"""
    bode(impedance::Array{Any}; omega_range=(-3, 5, 100), 
         titles::Array{String}=[""], omega=[], axis_type=:loglog, 
         save_data=false)

Plots impedance or frequency-dependent data as a Bode plot. The function supports various axis configurations and options for saving data.

# Arguments
- `impedance::Array{Any}`: A multidimensional array containing the impedance data to plot. Each element corresponds to a frequency-dependent impedance value.
- `omega_range::Tuple{Real, Real, Int}` (default `(-3, 5, 100)`): A tuple specifying the range of logarithmic frequency values:
  - `min_ω`: Minimum log-frequency (base 10).
  - `max_ω`: Maximum log-frequency (base 10).
  - `n_ω`: Number of frequency points.
- `titles::Array{String}` (default `[""]`): Optional titles or labels for the impedance data in LaTeX math format. Labels are matched to the data in `impedance`.
- `omega::Vector{Real}` (default `[]`): Optional precomputed frequency points. If provided, this overrides `omega_range`.
- `axis_type::Symbol` (default `:loglog`): Specifies the type of plot axes. Options:
  - `:loglog`: Logarithmic frequency and logarithmic impedance (in dB).
  - `:linlog`: Linear frequency and logarithmic impedance (in dB).
  - `:loglin`: Logarithmic frequency and linear impedance (magnitude).
  - `:linlin`: Linear frequency and linear impedance (magnitude).
  - `:logrealimag`: Logarithmic frequency with real and imaginary parts.
  - `:linrealimag`: Linear frequency with real and imaginary parts.
- `save_data::Bool` (default `false`): Whether to save the computed data (`frequency`, `magnitude`, `phase`) to a text file.

# Behavior
- Computes frequency points using `omega_range` if `omega` is not provided.
- Generates Bode plots with two subplots: one for the magnitude and one for the phase (or real/imaginary components, depending on `axis_type`).
- The `save_data` option writes frequency, magnitude, and phase data to a CSV file, with filenames timestamped.

# Plot Customization
- Titles/labels for the plots can be provided via `titles`.
- Handles various axis types to support different visual representations.

# Exceptions
- Throws `ArgumentError` if an invalid `axis_type` is provided.

# Example
```julia
impedance = [rand(ComplexF64, 2, 2) for _ in 1:100]
bode(impedance; omega_range=(-2, 4, 50), 
     titles=["Z₁₁", "Z₁₂", "Z₂₁", "Z₂₂"], 
     axis_type=:loglog, save_data=true)
"""
function bode(impedance :: Array{Any}; omega_range = (-3, 5, 100),
    titles :: Array{String} = [""], omega = [], axis_type = :loglog,
    save_data = false)
    (min_ω, max_ω, n_ω) = omega_range

    # set the type of plot
    if axis_type == :loglin
        p1 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$|H(j \\omega)|\$ ", grid = "major",})
        p2 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\arg(H(j\\omega)) \\, [^\\circ]\$", grid = "major",})
    elseif axis_type == :loglog
        p1 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$20 \\log_{10}|H(j \\omega)|\$ ", grid = "major",})
        p2 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\arg(H(j\\omega)) \\, [^\\circ]\$", grid = "major",})
    elseif axis_type == :linlog
        p1 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$|H(j \\omega)|\$ ", grid = "major",})
        p2 = @pgf Axis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\arg(H(j\\omega)) \\, [^\\circ]\$", grid = "major",})
    elseif axis_type == :linlin
        p1 = @pgf Axis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$|H(j \\omega)|\$ ", grid = "major",})
        p2 = @pgf Axis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\arg(H(j\\omega)) \\, [^\\circ]\$", grid = "major",})
    elseif axis_type == :logrealimag
        p1 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\Re \\{H(j \\omega)\\} \$", grid = "major",})
        p2 = @pgf SemiLogXAxis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\Im \\{H(j \\omega)\\} \$", grid = "major",})
    elseif axis_type == :linrealimag
        p1 = @pgf Axis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\Re \\{H(j \\omega)\\} \$", grid = "major",})
        p2 = @pgf Axis({xlabel = "\$\\omega \\, [\\mathrm{rad}/\\mathrm{s}]\$",
                ylabel = "\$\\Im \\{H(j \\omega)\\} \$", grid = "major",})
    else
        throw(ArgumentError("There is no axis_type $(axis_type)."))
    end

    if isempty(omega)
        n = (max_ω - min_ω) / n_ω
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = omega
    end

    for i in 1:size(impedance[1],1)
        for j in 1:size(impedance[1],2)
            (length(titles[1,1]) == 0) ? p_title = LegendEntry(string("\$H_{", i, "," , j, "}\$")) :
                                         p_title = LegendEntry(string("\$", titles[i, j], "\$"))

            # check every frequency
            mag = []
            phase = []

            if in(axis_type, [:loglog, :linlog])
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, 20log10(abs(imp[i,j]))))
                    push!(phase, angle(imp[i,j]))
                end
                phase = rad2deg.(phase)
            elseif in(axis_type, [:loglin, :linlin])
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, abs(imp[i,j])))
                    push!(phase, angle(imp[i,j]))
                end
                phase = rad2deg.(phase)
            else
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, real(imp[i,j])))
                    push!(phase, convert(Float64, imag(imp[i,j])))
                end
            end
            push!(p1, PlotInc(Table(omegas, mag)), p_title)
            push!(p2, PlotInc(Table(omegas, phase)), p_title)

            save_data && open(string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM"), "_z", i, j, ".txt"), "w") do io
                writedlm(io, [omegas mag phase], ',')
            end;
        end
    end

    p = @pgf GroupPlot(
        { group_style = { group_size="1 by 2"},
          no_markers,
          legend_pos="north west",
          height = "8cm",
          width = "12cm"
        },
        p1, p2)
end

function save_plot(figure, file_name = string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM")))
    pgfsave(string(file_name, ".pdf"), figure)
end
