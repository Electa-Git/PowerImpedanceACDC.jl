export nyquist, nyquist_multiplot

function nyquist(st; title :: String = "")
    p_title = string("\$ H(j \\omega) = ", title, "\$")
    p = @pgf Axis({xlabel = "\$\\Re \\{H(j \\omega)\\} \$",
            ylabel = "\$\\Im \\{H(j \\omega)\\} \$", grid = "major",
            only_marks, })

    tf = []
    for i in 1:size(st,1)
        if (real(st[i][3]) < 3.5) && (real(st[i][3]) > -3.5) &&
            (imag(st[i][3]) < 3.5) && (imag(st[i][3]) > -3.5)
            push!(tf, st[i][3])
        end
    end
    push!(p, PlotInc(Table(real(tf), imag(tf))))
    θ = 0:0.01:2*π
    push!(p, PlotInc(Table(cos.(θ), sin.(θ))))

end

function nyquist_multiplot(func; titles = [""])
    p = @pgf Axis({xlabel = "\$\\Re \\{H(j \\omega)\\} \$",
            ylabel = "\$\\Im \\{H(j \\omega)\\} \$", grid = "major",
            title = p_title, only_marks, })

    tf = []
    for j in 1:size(func,1)
        st = func[j]
        title = titles[j]
        (length(titles[1]) == 0) ? p_title = LegendEntry(string("\$H_{", j, "}\$")) :
                                     p_title = LegendEntry(string("\$", titles[j], "\$"))

        for i in 1:size(st,1)
            if (real(st[i][3]) < 3.5) && (real(st[i][3]) > -3.5) &&
                (imag(st[i][3]) < 3.5) && (imag(st[i][3]) > -3.5)
                push!(tf, st[i][3])
            end
        end
        push!(p, PlotInc(Table(real(tf), imag(tf))), p_title)
    end

    θ = 0:0.01:2*π
    push!(p, PlotInc(Table(cos.(θ), sin.(θ))))
    p = @pgf GroupPlot(
        { group_style = { group_size="1 by 1"},
          no_markers,
          legend_pos="north west",
          height = "8cm",
          width = "12cm"
        },
        p)
end
