export EVD

function EVD(Zcl_bus, omega, fmin, fmax)

    f = real(omega)./(2*pi)

    # Determine eigenvalues of the closed-loop bus impedance matrix at each frequency point 
    Λ = eigvals.(Zcl_bus)
    # Determine the related eigenvectors 
    Φ = eigvecs.(Zcl_bus)

    # Number of frequency points
    omegaₙ = size(Λ, 1)
    # Number of eigenvalues for each frequency point
    λₙ = size(Λ[1], 1)
    # Zcl_bus has dimensions λₙxλₙxomegaₙ 

    # Sorting algorithm of eigenvalues to create continuous lines in EVD plot 
    # Create empty matrix for sorted eigenvalues 
    Λₛ = zeros(Complex{Float64}, omegaₙ, λₙ)
    # First row of eigenvalues remains the same 
    Λₛ[1,:] = Λ[1]

    for i in 2:omegaₙ
        
        # Row of eigenvalues at previous frequency point 
        Λ₁ = Λ[i-1]
        # Row of eigenvalues at current frequency point
        Λ₂ = Λ[i]

        # Algorithm: 
        # see Munkres.jl https://juliapackages.com/p/munkres
        # see eigenshuffle.m https://nl.mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle
        # Sorting eigenvalues of row i to obtain minimum distance between row i and row i-1

        Λ₁_Re = real(Λ₁)
        Λ₂_Re = real(Λ₂)

        Λ₁_Re_Grid = [i for i in Λ₁_Re, j in 1:length(Λ₂_Re)]
        Λ₂_Re_Grid = [j for i in 1:length(Λ₁_Re), j in Λ₂_Re]

        Λ_Re = abs.(Λ₁_Re_Grid - Λ₂_Re_Grid)

        Λ₁_Im = imag(Λ₁)
        Λ₂_Im = imag(Λ₂)

        Λ₁_Im_Grid = [i for i in Λ₁_Im, j in 1:length(Λ₂_Im)]
        Λ₂_Im_Grid = [j for i in 1:length(Λ₁_Im), j in Λ₂_Im]

        Λ_Im = abs.(Λ₁_Im_Grid - Λ₂_Im_Grid)

        # Possible distances between eigenvalues
        D = sqrt.(Λ_Re.^2 + Λ_Im.^2)
        # Sequence of sorted eigenvalues
        S = munkres(D)

        Λ[i] = Λ[i][S]
        # Insert sorted eigenvalues in matrix Λₛ
        Λₛ[i, :] = Λ[i]
        # Eigenvectors are sorted based on same sequence 
        Φ[i] = Φ[i][:,S]

    end

    Λ = Λₛ

    # Normalizing eigenvectors
    for i in 1:omegaₙ

        for j in 1:λₙ
            N = 1/sqrt(transpose(Φ[i][:,j]) * Φ[i][:,j])
            Φ[i][:,j] = N * Φ[i][:,j]
        end

    end

    Ψ = inv.(Φ)

    # Oscillation modes are determined based on peak value of eigenvalue magnitude 
    λ_mag = abs.(Λ[:,1])
    λ_real = real.(Λ[:,1])
    λ_mag_dB = 20*log10.(λ_mag)

    # Determining x-axis limits 
    index_fmin = findmin(abs.(f .- fmin))[2]
    index_fmax = findmin(abs.(f .- fmax))[2]

    abs_lambda = abs.(Λ)

    # Determining y-axis limits
    max_lambda = zeros(λₙ, 1)
    min_lambda = zeros(λₙ, 1)

    for i in 1:λₙ
        max_lambda[i] = maximum(abs_lambda[index_fmin:index_fmax, i])
        min_lambda[i] = minimum(abs_lambda[index_fmin:index_fmax, i])
    end

    max_mode = maximum(max_lambda)
    y_max = maximum(max_lambda)
    y_min = minimum(min_lambda)

    # Find mode/eigenvalue with highest magnitude 
    mode = findmin(abs.(max_lambda .- max_mode))[2][1]
    index_mode = findmin(abs.(abs_lambda[:, mode] .- max_mode))[2]

    #plotly() # To activate interactive plot
    # Plotting of magnitude eigenvalues 
    p_abs = plot(f, λ_mag_dB, linewidth = 3, c = :blue, label = "Lambda " * string(1))  
    plot!(xlabel = "Frequency [Hz]", ylabel = "abs(lambda) [dB]", framestyle = :box, xaxis = :log10) 
    plot!(xlims = (fmin, fmax), ylims = (20*log10(y_min)-10, 20*log10(y_max)+10), legend = :topleft) #xticks = [10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6]

    color = palette(:default) #palette(:tab10)

    for i in 2:λₙ

        λ_mag = abs.(Λ[:,i])
        λ_mag_dB = 20*log10.(λ_mag)

        plot!(f, λ_mag_dB, linewidth = 3, c = color[i], label = "Lambda " * string(i))  

    end
    plot!(size=(1000,1000))
    display(p_abs)

    # Positive damping criterion 

    # Determining x-axis limits 
    index_fmin = findmin(abs.(f .- fmin))[2]
    index_fmax = findmin(abs.(f .- fmax))[2]
    
    real_lambda = real.(Λ)

    # Determining y-axis limits
    max_lambda = zeros(λₙ, 1)
    min_lambda = zeros(λₙ, 1)

    for i in 1:λₙ
        max_lambda[i] = maximum(real_lambda[index_fmin:index_fmax, i])
        min_lambda[i] = minimum(real_lambda[index_fmin:index_fmax, i])
    end

    min_mode_PDC = minimum(min_lambda)
    y_max = maximum(max_lambda)
    y_min = minimum(min_lambda)

    mode_PDC = findmin(abs.(min_lambda .- min_mode_PDC))[2][1]
    index_mode_PDC = findmin(abs.(real_lambda[:, mode_PDC] .- min_mode_PDC))[2]

    #plotly() # To activate interactive plot
    # Plotting of magnitude eigenvalues 
    p_re = plot(f, λ_real, linewidth = 3, c = :blue, label = "Lambda " * string(1))  
    plot!(xlabel = "Frequency [Hz]", ylabel = "Re(lambda)", framestyle = :box, xaxis = :log10) 
    plot!(xlims = (fmin, fmax), ylims = (1.1*y_min, 1.1*y_max), legend = :topleft) 

    color = palette(:default) #palette(:tab10)

    for i in 2:λₙ

        λ_real= real.(Λ[:,i])

        plot!(f, λ_real, linewidth = 3, c = color[i], label = "Lambda " * string(i))  

    end
    plot!(size=(1000,1000))
    display(p_re)

    p = plot(p_abs, p_re, layout = (2, 1))
    display(p)

    # Dominant oscillation mode is determined based on the frequency at which eigenvalue reaches highest magnitude 
    println("The oscillation frequency is ", round(f[index_mode]; digits = 1), " Hz")
    println("The critical eigenvalue is ", mode)


    # Dominant oscillation mode is determined based on the frequency at which eigenvalue reaches highest magnitude 
    println("The oscillation frequency based on PDC is ", round(f[index_mode_PDC]; digits = 1), " Hz")
    println("The critical eigenvalue is ", mode_PDC)

    # Observability, controllability and participation factors can be determined for dominant oscillation mode at the related frequency point 

    # Observability eigenvalue decomposition

    O = zeros(λₙ, λₙ);

    # k (Bus) = row; i (Mode) = column
    for i = 1:λₙ
        for k = 1:λₙ
            O[k,i] = abs(Φ[index_mode][k, i]);
        end

        Tot = sum(O[:,i])

        for k = 1:λₙ
            O[k, i] = O[k, i]/Tot;   
        end  
    end

    println("The observability of the critical eigenvalue at the buses is ", round.(O[:,mode]; digits = 2)) 

    # Controllability eigenvalue decomposition

    C = zeros(λₙ, λₙ);

    # k (Bus) = row; i (Mode) = column
    for i = 1:λₙ
        for k = 1:λₙ
            C[k,i] = abs(Ψ[index_mode][i, k]);
        end

        Tot = sum(C[:,i])

        for k = 1:λₙ
            C[k, i] = C[k, i]/Tot;   
        end  
    end

    println("The controllability of the critical eigenvalue at the buses is ", round.(C[:,mode]; digits = 2))

    # Participation factors eigenvalue decomposition

    P = zeros(λₙ, λₙ);

    # k (Bus) = row; i (Mode) = column
    for i = 1:λₙ
        for k = 1:λₙ
            P[k,i] = abs(Φ[index_mode][k, i] * Ψ[index_mode][i, k]);
        end

        Tot = sum(P[:,i])

        for k = 1:λₙ
            P[k, i] = P[k, i]/Tot;   
        end  

    end

    println("The participation factors of the buses in the critical eigenvalue are ", round.(P[:,mode]; digits = 2))

    # Participation factors eigenvalue decomposition

    P = zeros(λₙ, λₙ);

    # k (Bus) = row; i (Mode) = column
    for i = 1:λₙ
        for k = 1:λₙ
            P[k,i] = abs(Φ[index_mode_PDC][k, i] * Ψ[index_mode_PDC][i, k]);
        end

        Tot = sum(P[:,i])

        for k = 1:λₙ
            P[k, i] = P[k, i]/Tot;   
        end  

    end

    println("The participation factors of the buses in the critical eigenvalue based on PDC are ", round.(P[:,mode_PDC]; digits = 2))

end