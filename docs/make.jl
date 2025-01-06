using Documenter
using PowerImpedanceACDC

makedocs(
    sitename = "PowerImpedanceACDC",
    format = Documenter.HTML(),
    modules = [PowerImpedanceACDC],
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

#DOCUMENTER_KEY = ENV["DOCUMENTER_KEY"]

deploydocs(
    repo = "git@gitlab.kuleuven.be/electa/controlgroup/hvdcstability.jl.git", # u0167553:$DOCUMENTER_KEY@
    devbranch="docs",
    # branch="gl-pages",
    #deploy_config = Documenter.GitLab()
)
