using Documenter, HVDCstability

makedocs(
    modules = [HVDCstability],
    format = Documenter.HTML(),
    sitename = "HVDCstability.jl",
    authors = "Aleksandra Lekic",
    pages = [
        "Home" => "index.md"
        "Manual" => [
            "Results" => "result-data.md",
        ]
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Components" => [
                "Source" => "source.md",
                "Impedance" => "impedance.md",
                "Transformer" => "transformer.md",
                "Shunt reactor" => "shunt.md",
                "Transmission line" => "transmission_line.md",
                "MMC" => "MMC.md"
                ]
        ]
    ]
)

deploydocs(
    repo = "github.com/Electa-Git/HVDCstability.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    forcepush = true
)
