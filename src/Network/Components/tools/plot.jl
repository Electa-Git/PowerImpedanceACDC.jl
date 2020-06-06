export bode, save_plot

"""
    function bode(impedance :: Array{Any}; omega_range = (-3, 5, 100),
        titles :: Array{String} = [""], omega = [], axis_type = :loglog,
        save_data = false)
Used for plotting determined impedance or a frequency dependent data as a bode
plot. Function takes frequency points in the form of `omega_range` or as mapped
values `omega`. For a nice diagrams, the labels can be given as `titles`.

It can be specified how to plot data:
- :loglog - logarithmic frequency scale and logarithmic impedance in dB
- :linlog - linear frequency scale and logarithmic impedance in dB
- :loglin - logarithmic frequency scale and linear impedance (magnitude)
- :logrealimag - logarithmic frequency scale and real/imaginary part
- :linrealimag - linear frequency scale and real/imaginary part

Data can be saved by setting `save_data = true`.
"""
function bode(impedance :: Array{Any}; omega_range = (-3, 5, 100),
    titles :: Array{String} = [""], omega = [], axis_type = :loglog,
    save_data = false)
    (min_ω, max_ω, n_ω) = omega_range
    gr()
    p1 = plot(reuse = true)
    p2 = plot(reuse = true)
    p = plot(reuse = true)
    ylabel1 = 0
    ylabel2 = 0
    xaxis = 0

    if isempty(omega)
        n = (max_ω - min_ω) / n_ω
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = omega
    end

    l = @layout [a;b]
    for i in 1:size(impedance[1],1)
        for j in 1:size(impedance[1],2)
            # check every frequency
            mag = []
            phase = []

            if in(axis_type, [:loglog, :linlog])
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, 20log10(abs(imp[i,j]))))
                    push!(phase, angle(imp[i,j]))
                end
                ylabel1 = L"20\log_{10}|H(j\omega)| "
                ylabel2 = L"\arg(H(j\omega)) \, [^\circ]"
                phase = rad2deg.(phase)
            elseif in(axis_type, [:loglin, :linlin])
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, abs(imp[i,j])))
                    push!(phase, angle(imp[i,j]))
                end
                ylabel1 = L"|H(j\omega)| "
                ylabel2 = L"\arg(H(j\omega)) \, [^\circ]"
                phase = rad2deg.(phase)
            else
                for imp in impedance
                    # add magnitude and phase
                    push!(mag, convert(Float64, real(imp[i,j])))
                    push!(phase, convert(Float64, imag(imp[i,j])))
                end
                ylabel1 = L"\Re\{H(j\omega)\} "
                ylabel2 = L"\Im\{H(j\omega)\}"
            end

            if in(axis_type, [:linlog, :linlin, :linrealimag])
                xaxis = :lin
            else
                xaxis = :log
            end

            save_data && open(string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM"), "_z", i, j, ".txt"), "w") do io
                writedlm(io, [omegas mag phase], ',')
            end;

            (length(titles[1,1]) == 0) ? p_title = string("Z_{", i, "," , j, "}") : p_title = titles[i, j]
            a = latexstring(string(p_title))
            plot!(p1, omegas, mag, xaxis=xaxis,
                xtickfont = font(6, "Courier"), ytickfont = font(6, "Courier"),
                xlabel=L"\omega \, [\mathrm{rad}/\mathrm{s}]", ylabel=ylabel1, leg = false,
                label = a)
            plot!(p2, omegas, phase, xaxis=xaxis,
                xtickfont = font(6, "Courier"), ytickfont = font(6, "Courier"),
                xlabel=L"\omega \, [\mathrm{rad}/\mathrm{s}]", ylabel=ylabel2,
                label=a)
        end
    end
    p = plot(p1, p2, layout = l, legend = :left)
    display(p)
end

function save_plot(file_name = string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM")))
    savefig(string(file_name, ".pdf"))
end
