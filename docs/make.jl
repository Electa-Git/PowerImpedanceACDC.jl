using Documenter
using PowerImpedanceACDC

makedocs(
    sitename = "PowerImpedanceACDC",
    format = Documenter.HTML(inventory_version="0.0"),
    modules = [PowerImpedanceACDC],
    doctest = false,
    repo = "gitlab.kuleuven.be/electa/controlgroup/hvdcstability.jl.git",
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

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "gitlab.kuleuven.be/electa/controlgroup/hvdcstability.jl.git",
    devbranch="main",
    branch="gl-pages",
    deploy_config = Documenter.GitLab()
)
