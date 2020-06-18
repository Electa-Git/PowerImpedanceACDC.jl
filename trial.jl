i = 1 #re-assign the value 1 to i
for key in keys(c.insulators) #As described at row 154 for the conductors-> depending on how many keys are present in the variable insulators -> it makes the same process for all the insulating layers present in the variable c
    (rᵢ, rₒ, μ, ϵ) = (c.insulators[key].rᵢ, c.insulators[key].rₒ, c.insulators[key].μᵣ*μ₀, c.insulators[key].ϵᵣ*ϵ₀)
    Zⁱ = s*μ/(2π) * log(rₒ/rᵢ) #insulator layer impedance -> eq 41 pag 24
    Pⁱ = log(rₒ/rᵢ) / (2π*ϵ) #p expression pag25 2nd row simulator_tutorial

    Z[i,i] += Zⁱ #for the impedance matrix diagonal values add also the insulator impedance-> showed in 1st row eq 44 pag 24 simulator_tutorial
    P[1:i,1:i] += ones(i,i) * Pⁱ
    i += 1
end
