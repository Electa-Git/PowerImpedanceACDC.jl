export nyquistplot

function nyquistplot(L, omega; zoom :: String = "", SM :: String = "", title :: String = "Nyquist plot")

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
    
    #plotly() # interactive plot
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
    np = plot(cos.(θ), sin.(θ), c = :red, linewidth = 3, line = :dash, label = :none)
    plot!(xlabel = "Real axis", ylabel = "Imaginary axis", title = title, xlims = xlims, ylims = ylims, framestyle = :box, legend=:topleft) #legend = :none
    vline!([0], c = :black, line = :dash, label = :none)
    hline!([0], c = :black, line = :dash, label = :none)
    scatter!([-1], [0]; markercolor = :red, markershape = :xcross, markersize = 10, label = :none)
    for i in 1:size(Λ,2) 
            λ = Λ[:,i]
            x = real(λ)
            y = imag(λ)
            y_abs = abs.(y)
            y_max = maximum(y_abs)
            color = [:blue :green :orange :purple :cyan :yellow :pink :brown :gray :honeydew] #Now colors up to 10 eigenloci
            plot!(x, y, linewidth = 3, c = color[i], label = "Lambda " * string(i))
            plot!(x, -y, linewidth = 3, c = color[i], linestyle = :dash, label = :none)
            index = findall(x -> x == y_max, y_abs)
            index = index[1]
                    if index != length(λ)
                            if x[index+1] < x[index]
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10, label = :none)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10, label = :none)
                            else
                                    scatter!([x[index]], [y[index]]; markercolor = color[i], markershape = :rtriangle, markersize = 10, label = :none)
                                    scatter!([x[index]], [-y[index]]; markercolor = color[i], markershape = :ltriangle, markersize = 10, label = :none)
                            end
                    end 
    end
    plot!(size=(1000,1000))
    display(np)

    x = real(Λ)
    y = imag(Λ)

    cw = []
    ccw = []

    for i in 1:λₙ
        cwi = 0      #Number of clockwise encirclements
        ccwi = 0     #Number of counterclockwise encirclements  
        for j in 2:omegaₙ
            if (y[j-1, i] < 0 && y[j, i] > 0 && x[j, i] < -1)
                cwi += 1
            elseif (y[j-1, i] > 0 && y[j, i] < 0  && x[j-1, i] < -1)
                ccwi += 1
            end
        end
        push!(cw, cwi)
        push!(ccw, ccwi)
    end

    # Nyquist stability criterion
    # N: Net number of clockwise encirclements
    # P: Number of RHP poles of loop gain (matrix)
    # Z: number of RHP poles of closed-loop system
    # N = Z - P 
    # If Z = N + P > 0, system is unstable 
    # P is equal to or greater than zero => System is unstable if N > 0

    N = sum(cw) - sum(ccw) 

    println("")
    if N > 0
        println("Result stability assessment: Unstable system \n")
    elseif N < 0
        println("Result stability assessment: Unstable subsystem \n")
    else
        println("Result stability assessment: Stable system if subsystems are stable \n")
    end

    ##### ----- STABILITY MARGINS ----- #####

    # Possibilities for SM are "PM", "GM" and "VM"
    for i in 1:λₙ
        println("Stability margins for λ", i, ":")
        stabilitymargin(Λ[:,i], omega, SM = SM)
        println("")
    end

end

