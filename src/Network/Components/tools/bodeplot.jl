export bodeplot

function bodeplot(L, omega; legend = [""])

    f = omega/(2*π)
    L_mag = abs.(L)
    L_mag_dB = 20*log10.(L_mag)
    L_ph = angle.(L).*(180/π)
    if minimum(L_mag) != maximum(L_mag) 
        ylims_mag = (minimum(L_mag_dB) - 0.1*abs(minimum(L_mag_dB)), maximum(L_mag_dB) + 0.1*abs(minimum(L_mag_dB)))
        ylims_ph = (minimum(L_ph) - 0.1*abs(minimum(L_ph)), maximum(L_ph) + 0.1*abs(minimum(L_ph)))
    else
        ylims_mag = (maximum(L_mag_dB)-100, maximum(L_mag_dB)+100)
        ylims_ph = (maximum(L_ph)-90, maximum(L_ph)+90)
    end
    
    if legend == [""]
        p1 = plot(f, L_mag_dB, legend = :none, linewidth = 3, linestyle = :auto)  # For different line styles: linestyle = :auto
    else
        p1 = plot(f, L_mag_dB, label = legend, linewidth = 3, linestyle = :auto)  # For different line styles: linestyle = :auto   
    end
    plot!(ylabel = "Magnitude [dB]", framestyle = :box, xaxis = :log10)
    plot!(xlims = (minimum(f), maximum(f)), ylims = ylims_mag)
    p2 = plot(f, L_ph, linewidth = 3, linestyle = :auto)  # For different line styles: linestyle = :auto
    plot!(xlabel = "Frequency [Hz]", ylabel = "Phase [deg]", framestyle = :box, legend = :none, xaxis = :log10)
    plot!(xlims = (minimum(f), maximum(f)), ylims = ylims_ph)
    plot!(yticks = -360:90:360)
    plot(p1, p2, layout = (2, 1))

end
