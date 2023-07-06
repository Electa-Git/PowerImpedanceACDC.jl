export EVD

function EVD(Zcl_bus, omega, fmin, fmax)

    f = real(omega)./(2*pi)

    Λ = eigvals.(Zcl_bus)
    Φ = eigvecs.(Zcl_bus)
    Λd = Diagonal.(eigvals.(Zcl_bus))
    #Zcl_bus = Φ .* Λd .* inv.(Φ)
    #Zcl_bus[1] = Φ[1] * Λd[1] * inv(Φ[1])

    # Eigenvalue sorting algorithm 
    omegaₙ = size(Λ, 1)
    λₙ = size(Λ[1], 1)
    Sequence = zeros(Int, λₙ)

    for i in 1:λₙ
        Sequence[i] = i
    end

    for i in 2:omegaₙ

        R₁ = Λ[i-1]
        per = permutations(Λ[i])
        R₂ = []

        for item in per
            push!(R₂, item)
        end

        D = zeros(length(per), 1)

        for j in 1:length(R₂)
            Dk = 0
            for k in 1:λₙ
                Dk += abs(R₂[j][k] - R₁[k]) 
            end
            D[j] = Dk
        end

        index = findmin(D)[2]
        index = index[1]

        Λ[i] = R₂[index]

        per = permutations(Sequence)
        S = []

        for item in per
            push!(S, item)
        end

        Mi = zeros(Complex{Float64}, λₙ, λₙ)

        for j in 1:λₙ 
            Mi[:,j] = Φ[i][:,S[index][j]]
        end

        Φ[i] = Mi

    end

    # Normalizing eigenvectors
    for i in 1:omegaₙ

        for j in 1:λₙ
            N = 1/sqrt(transpose(Φ[i][:,j]) * Φ[i][:,j])
            Φ[i][:,j] = N * Φ[i][:,j]
        end

    end

    Ψ = inv.(Φ)
    Λd = Diagonal.(Λ)
    # Zcl_bus[1] = Φ[1] * Λd[1] * inv(Φ[1])

    # Conversion array of arrays to array

    Λₛ = zeros(Complex{Float64}, omegaₙ, λₙ)

    for i in 1:omegaₙ
        Λₛ[i, :] = Λ[i]
    end

    Λ = Λₛ

    λ_mag = abs.(Λ[:,1])
    λ_mag_dB = 20*log10.(λ_mag)

    index_fmin = findmin(abs.(f .- fmin))[2]
    index_fmax = findmin(abs.(f .- fmax))[2]

    abs_lambda = abs.(Λ)

    max_lambda = zeros(λₙ, 1)
    min_lambda = zeros(λₙ, 1)

    for i in 1:λₙ
        max_lambda[i] = maximum(abs_lambda[index_fmin:index_fmax, i])
        min_lambda[i] = minimum(abs_lambda[index_fmin:index_fmax, i])
    end

    max_mode = maximum(max_lambda)
    y_max = maximum(max_lambda)
    y_min = minimum(min_lambda)

    mode = findmin(abs.(max_lambda .- max_mode))[2][1]
    index_mode = findmin(abs.(abs_lambda[:, mode] .- max_mode))[2]

    #plotly() # interactive plot
    p = plot(f, λ_mag_dB, linewidth = 3, c = :blue, label = "Lambda " * string(1))  
    plot!(xlabel = "Frequency [Hz]", ylabel = "abs(lambda) [dB]", framestyle = :box, xaxis = :log10) 
    plot!(xlims = (fmin, fmax), ylims = (20*log10(y_min)-10, 20*log10(y_max)+10), legend = :topleft) #xticks = [10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6]

    color = [:blue :green :orange :purple :cyan :yellow :pink :brown :gray :honeydew]

    for i in 2:λₙ

        λ_mag = abs.(Λ[:,i])
        λ_mag_dB = 20*log10.(λ_mag)

        plot!(f, λ_mag_dB, linewidth = 3, c = color[i], label = "Lambda " * string(i))  

    end

    display(p)

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

    println("The oscillation frequency is ", round(f[index_mode]; digits = 1), " Hz")
    println("The critical eigenvalue is ", mode)
    #println("The critical right eigenvector is ", round.(Φ[index_mode][:, mode]; digits = 2))
    #println("The critical left eigenvector is ", round.(Ψ[index_mode][mode, :]; digits = 2))
    println("The participation factors of the buses in the critical eigenvalue are ", round.(P[:,mode]; digits =2))

end