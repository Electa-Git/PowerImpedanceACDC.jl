# PowerImpedanceACDC.jl Documentation

```@meta
CurrentModule = PowerImpedanceACDC
```

## Overview

PowerImpedanceACDC.jl is a Julia package, which focuses on the impedance-based stability. Simulator provides support for the power system construction using detailed components. Components and their interconnections are modeled as multiport networks using ABCD parameters.

Implemented components are:
- Impedances
- AC and DC grid equivalents (sources)
- Transformers
- Overhead lines
- Cables
- MMCs (using various controllers)

## Installation
The latest stable release of PowerImpedancdACDC can be installed using the Julia package manager with
```julia
] add PowerImpedanceACDC
```

Test that the package works by running
```julia
] test PowerImpedanceACDC
```

## Example
The following will show how PowerImpedanceACDC can be utilized to assess the stability of a Single-Machine Infinite Bus system (SMIB). This should be replaced with working example!
 ```julia
 using DelimitedFiles
using SymEngine
using Plots
using LinearAlgebra
using PowerImpedanceACDC

s = symbols(:s)
ω = 2π*50

net = @network begin
    ac = ac_source(V = 230)
    x_sys = impedance(pins = 1, z = 0.04 + s*0.3/ω)
    

    ac[1.1] ⟷ x_sys[1.1] ⟷ bus1
    ac[2.1] ⟷ x_sys[2.1] ⟷ gnd
    
end

omega_Y_tlc1 = collect(range(2*pi*0.1, stop=2*pi*5000, step=1))



Y_dd = []


for i in 1:length(omega_Y_tlc1)
    Y1 = eval_abcd(net.elements[:x_sys].element_value, 1im*omega_Y_tlc1[i]) 
    push!(Y_dd, Y1[2,2])

end

bode_Ydd = bodeplot(Y_dd, omega_Y_tlc1, legend = "Y_dd")
```

## Citation

