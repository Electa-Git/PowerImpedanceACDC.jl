using Documenter
using PowerImpedanceACDC

makedocs(
    sitename = "PowerImpedanceACDC",
    format = Documenter.HTML(),
    modules = [PowerImpedanceACDC],
    repo = "gitlab.kuleuven.be/electa/controlgroup/hvdcstability.jl.git",
    pages = [
        "Getting started" => "index.md"
        "Introduction" =>  "introduction.md"
        "Manual" => [
            "Network" => "network.md",
            "Initialization" => "initialization.md",
            "Impedance & Stability" => "results.md",
            "Examples" => "example.md",
        ]
        
        "Components" => [
            "Source" => "source.md",
            "Impedance" => "impedance.md",
            "Transformer" => "transformer.md",
            "Shunt reactor" => "shunt.md",
            "Transmission line" => "transmission_line.md",
            "MMC" => "MMC.md",
            "TLC" => "TLC.md"
            ]
        
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

#DOCUMENTER_KEY = ENV["DOCUMENTER_KEY"]

deploydocs(
    repo = "gitlab.kuleuven.be/electa/controlgroup/hvdcstability.jl.git", # u0167553:$DOCUMENTER_KEY@
    devbranch="docs",
    # branch="gl-pages",
    #deploy_config = Documenter.GitLab()
)
