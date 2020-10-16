
rᵢi₁ = 31.75e-3
rₒi₁ = 60.85e-3
a = 33.75e-3
b = 59.55e-3
ϵᵣi₁ = 2.26
x = log(rₒi₁ / rᵢi₁) / log( b / a) # x = log(rₒI₁/rᵢI₁)/log(bI₁/aI₁)
ϵᵣi₁ *=  x # ϵᵣi₁ = ϵᵣi₁* x


rᵢi₂ = 61.05e-3
rₒi₂ = 65.95e-3
rₒc₂ = 61.05e-3
ϵᵣi₂ = 2.26

#if in(:I2, keys(c.insulators)) #if I2 is present in the keywords of the variable insulators -> enter the if cycle
x = log(rₒi₂ / rₒc₂) / log(rₒi₂ / rᵢi₂) # x = log(rₒI₂/rₒC₂)/log(rₒI₂/rᵢI₂)
rᵢi₂ = rₒc₂ # RᵢI₂ = RₒC₂ -> The internal radius of insulator I2 is equal to the outer radius of conductor C2.
ϵᵣi₂ *=  x #Same as c.insulators[:I₂].ϵᵣ =  c.insulators[:I₂].ϵᵣ*x ->  ϵᵣI₂ = ϵᵣI₂ * x
#c.insulators[:I2].μᵣ /= x #Same as c.insulators[:I₂].μᵣ =  c.insulators[:I₂].μᵣ/x ->  μᵣI₂ = μᵣI₂ / x

#if in(:I1, keys(c.insulators)) && (c.insulators[:I1].a != 0) #If in the keys of the variable c is defined I₁ AND I₁ has a semiconductor section (outer radius, inner semiconductor different from 0) -> enter
x = log(rₒi₁ / rᵢi₁) / log( b / a) # x = log(rₒI₁/rᵢI₁)/log(bI₁/aI₁)
ϵᵣi₁ *=  x #Same as c.insulators[:I₁].ϵᵣ =  c.insulators[:I₁].ϵᵣ*x ->  ϵᵣI1 = ϵᵣI₁ * x
#N = 1.4
#c.insulators[:I1].μᵣ = c.insulators[:I1].μᵣ * (1 + 2π^2 * N^2 * (c.insulators[:I1].rₒ^2 - c.insulators[:I1].rᵢ^2) / log(c.insulators[:I1].rₒ / c.insulators[:I1].rᵢ)) # μᵣI₁ = μᵣI₁ * (1+2π²N² * (rₒI₁²-rᵢI₁²) / log(rₒI₁-rᵢI₁))
#end

#for key in keys(c.insulators) #As described at row 154 for the conductors-> depending on how many keys are present in the variable insulators -> it makes the same process for all the insulating layers present in the variable c
#(rᵢ, rₒ, μ, ϵ) = (c.insulators[key].rᵢ, c.insulators[key].rₒ, c.insulators[key].μᵣ*μ₀, c.insulators[key].ϵᵣ*ϵ₀)
#Zⁱ = s*μ/(2π) * log(rₒ/rᵢ) #insulator layer impedance -> eq 41 pag 24
#Pⁱ = log(rₒ/rᵢ) / (2π*ϵ) #p expression pag25 2nd row simulator_tutorial
