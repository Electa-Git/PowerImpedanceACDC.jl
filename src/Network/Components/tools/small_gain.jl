export small_gain

function small_gain(G1, G2, omega, title :: String = "Small gain: SVD")

    f = real(omega)./(2*pi)

    # Determine maximum singular values of the matrices and their product at each frequency
    S1 = maximum.(svdvals.(G1))
    S2 = maximum.(svdvals.(G2))
    S12 = maximum.(svdvals.(G1.*G2))

    #plotly() # To activate interactive plot
    p_abs = plot(f, S12, linewidth = 3, c = :black, label = "G1(s)*G2(s)")  
    plot!(f, S1, linewidth = 3, line = :dash, c = :blue, label = "G1(s)")  
    plot!(f, S2, linewidth = 3, line = :dash, c = :red, label = "G2(s)")  
    plot!(xlabel = "Frequency [Hz]", ylabel = "max(sigma( . ))", title = title, framestyle = :box, xaxis = :log10, yaxis = :log10) 
    plot!(xlims = (minimum(f),maximum(f)), legend = :topleft)  #xticks = [10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6]
    plot!(size=(1000,1000))
    display(p_abs)
    return passivity_index

end