#To do:
#Automatic answer if system is stable or not based on real axis crossing of eigenloci
#Adapt and integrate stabiliymargin.jl to determine stability margins of eigenloci

export nyquistplot

function nyquistplot(L; zoom :: String = "")

    Λ = eigvals.(L)

    # Sorting algorithm based on minimum sum of differences between eigenvalues at frequency i-1 and at frequency i

    omegaₙ = size(Λ, 1)
    λₙ = size(Λ[1], 1)

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

    end

    # Conversion array of arrays to array

    Λₛ = zeros(Complex{Float64}, omegaₙ, λₙ)

    for i in 1:omegaₙ
        Λₛ[i, :] = Λ[i]
    end

    Λ = Λₛ

    x = real(Λ)
    y = imag(Λ)
    
    # plotly() # interactive plot
    x_min = minimum(x)
    x_max = maximum(x)
    y_abs = abs.(y)
    y_max = maximum(y_abs)
    if zoom == "yes" 
        xlims = (-2.6, 0.6)
        ylims = (-1.6, 1.6)
    else
        xlims = (x_min - 0.1*abs(x_min), x_max + 0.1*abs(x_max))
        ylims = (-1.1*y_max, 1.1*y_max)
    end
    θ = 0:0.01:2*π
    plot(cos.(θ), sin.(θ), c = :red, linewidth = 3, line = :dash)
    plot!(xlabel = "Real axis", ylabel = "Imaginary axis", xlims = xlims, ylims = ylims, framestyle = :box, legend = :none)
    vline!([0], c = :black, line = :dash)
    hline!([0], c = :black, line = :dash)
    scatter!([-1], [0]; markercolor = :red, markershape = :xcross, markersize = 10)
    for i in 1:size(Λ,2) 
            λ = Λ[:,i]
            x = real(λ)
            y = imag(λ)
            y_abs = abs.(y)
            y_max = maximum(y_abs)
            color = [:blue :green :orange :purple :purple :cyan :yellow :pink :brown :gray] #Now colors up to 10 eigenloci
            plot!(x, y, linewidth = 3, c = color[i])
            plot!(x, -y, linewidth = 3, c = color[i], linestyle = :dash)
            index = findall(x -> x == y_max, y_abs)
            index = index[1]
                    if index != length(λ)
                            if x[index+1] < x[index]
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10)
                            else
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10)
                            end
                    end 
    end
    plot!(size=(1000,1000))

end

