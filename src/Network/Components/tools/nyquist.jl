export nyquist

function nyquist(st; title :: String = "")
    p_title = string("\$ H(j \\omega) = ", title, "\$")
    p = @pgf Axis({xlabel = "\$\\Im \\{H(j \\omega)\\} \$",
            ylabel = "\$\\Re \\{H(j \\omega)\\} \$", grid = "major",
            title = p_title, no_markers,})

    tf = []
    for i in 1:size(st,1)
        if (real(st[i][3]) < 1.5) && (real(st[i][3]) > -1.5) &&
            (imag(st[i][3]) < 1.5) && (imag(st[i][3]) > -1.5)
            push!(tf, st[i][3])
        end
    end
    push!(p, PlotInc(Table(real(tf), imag(tf))))
    #push!(p, PlotInc(Table(real(tf), -imag(tf))))
    θ = 0:0.01:2*π
    push!(p, PlotInc(Table(cos.(θ), sin.(θ))))

end
