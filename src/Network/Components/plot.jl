export bode, save_plot

"""
    function bode(impedance :: Array{Any}; omega_range = (-3, 5, 100),
        titles :: Array{String} = [""], omega = [])
Used for plotting determined impedance or a frequency dependent data as a bode
plot. Function takes frequency points in the form of `omega_range` or as mapped
values `omega`. For a nice diagrams, the labels can be given as `titles`.
"""
function bode(impedance :: Array{Any}; omega_range = (-3, 5, 100),
    titles :: Array{String} = [""], omega = [])
    (min_ω, max_ω, n_ω) = omega_range
    gr()
    p1 = plot(reuse = true)
    p2 = plot(reuse = true)
    p = plot(reuse = true)

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

            for imp in impedance
                # add magnitude and phase
                push!(mag, convert(Float64, 20log10(abs(imp[i,j]))))
                push!(phase, angle(imp[i,j]))
            end
            # unwrap!(convert.(Float64,phase), dims = ndims(phase))
            phase = rad2deg.(phase)

            open(string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM"), "_z", i, j, ".txt"), "w") do io
                writedlm(io, [omegas mag phase], ',')
            end;

            (length(titles[1,1]) == 0) ? p_title = string("Z_{", i, "," , j, "}") : p_title = titles[i, j]
            plot!(p1, omegas, mag, xaxis=:log, xtickfont = font(6, "Courier"), ytickfont = font(6, "Courier"),
                xlabel=L"\omega \, [\mbox{rad}/\mbox{s}]", ylabel=L"20\log_{10}|H(j\omega)| ", leg = false, #\, [\Omega]
                label = latexstring(p_title))
            plot!(p2, omegas, phase, xaxis=:log, xtickfont = font(6, "Courier"), ytickfont = font(6, "Courier"),
                xlabel=L"\omega \, [\mbox{rad}/\mbox{s}]", ylabel=L"\arg(H(j\omega)) \, [^\circ]", #leg = false,
                label = latexstring(p_title))
        end
    end
    p = plot(p1, p2, layout = l, legend = :left)
    display(p)
end

function save_plot(file_name = string("files/", Dates.Date(Dates.now()), "_", Dates.format(now(), "HH_MM")))
    savefig(string(file_name, ".pdf"))
end
