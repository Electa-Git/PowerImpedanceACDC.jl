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
                "Impedance" => "impedance.md",
                "Transmission line" => "transmission_line.md",
                "Source" => "source.md",
                "Transformer" => "transformer.md",
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
