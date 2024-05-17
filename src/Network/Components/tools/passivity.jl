export passivity

function passivity(G, omega, title :: String = "Passivity assessment")

    f = real(omega)./(2*pi)

    # Determine minimum real part of the eigenvalues of the matrix at each frequency
    passivity_index = minimum.(real.(eigvals.(G)))
    # plotly() # To activate interactive plot
    p_plot = plot(f, passivity_index, linewidth = 3, c = :blue, label = "G(s)")  
    plot!([f[1],f[end]], [0,0], linewidth = 1.5, c = :red, line = :dash, label = "Zero")  
    plot!(xlabel = "Frequency [Hz]", ylabel = "Re( lambda( . ) )", title = title, framestyle = :box, xaxis = :log10) 
    plot!(xlims = (minimum(f),maximum(f)), legend = :topleft)  #xticks = [10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6] ,ylims = (minimum(passivity_index),maximum(passivity_index))
    plot!(size=(1000,1000))
    display(p_plot)
    return passivity_index

end