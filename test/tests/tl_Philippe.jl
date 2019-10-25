net = @network begin
    vs = ac_source_3Φ(; amplitude = 0, frequency = 50, phase = 0, impedance = 0)

    tl = transmission_line(length = 200e3,  # [m], transmission line length
    conductors = Conductors(organization = :flat,
                            nᵇ      = 3,       # number of bundles (phases)
                            nˢᵇ     = 1,       # number of subconductors per bundle
                            Rᵈᶜ     = 0.063,   # [Ω/m], DC resistance for the entire conductor
                            gᶜ      = 1e-11,   # [S/m] shunt conductance
                            μᵣᶜ     = 1,       # [], relative conductor permeability
                            rᶜ      = 0.015,   # [m], conductor radius
                            yᵇᶜ     = 30,      # [m], height above the ground of the lowest bundle
                            Δyᵇᶜ    = 0,       # [m], vertical offset between the bundles
                            Δxᵇᶜ    = 10,      # [m], horizontal offset between the lowest bundles
                            Δ̃xᵇᶜ    = 0,       # [m], horizontal offset in group of bundles
                            dˢᵇ     = 0,       # [m], *sub*conductor spacing (symmetric)
                            dˢᵃᵍ    = 10),     # [m], sag offset
    groundwires = Groundwires(nᵍ    = 2,       # number of groundwires (typically 0 or 2)
                            Rᵍᵈᶜ    = 0.92,    # [Ω/m], groundwire DC resistance
                            μᵣᵍ     = 1,       # [], relative groundwire permeability
                            rᵍ      = 0.0062,  # [m], ground wire radius
                            Δxᵍ     = 6.5,     # [m], horizontal ofsset between groundwires
                            Δyᵍ     = 7.5,     # [m], vertical offset between the lowest conductor and groundwires
                            dᵍˢᵃᵍ   = 10),     # [m], sag offset
    earth_parameters = (1,1,100))              # ([], [], [Ωm]) = (μᵣ_earth, ϵᵣ_earth, ρ_earth)

    vs[2.1] ⟷ vs[2.2] ⟷ vs[2.3] ⟷ gnd
    vs[1.1] ⟷ tl[1.1] ⟷ Node11
    vs[1.2] ⟷ tl[1.2] ⟷ Node12
    vs[1.3] ⟷ tl[1.3] ⟷ Node13
    tl[2.1] ⟷ gnd1 # name MUST contain "gnd" letters, be can add some additional letters
    tl[2.2] ⟷ gnd2
    tl[2.3] ⟷ gnd3
end

imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node11, :Node12, :Node13],
                            output_pins = Any[:gnd1, :gnd2, :gnd3], omega_range = (0, 5, 1000))
bode(imp, omega_range = (0, 5, 1000))
