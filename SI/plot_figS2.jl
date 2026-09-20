"""
Reproduces `ProofOfMarkovianCorrection.pdf`.

Shows the effective rate constant k(t) computed from a HEOM simulation for a
few different numbers of auxiliary density operators (ADOs), demonstrating
convergence to the Markovian value as more bath modes are included.

Required data:
  - data/BathsData/enhancement_bath_sweep_precise/etac0.jld2
  - SI/MarkovianityHEOMApproximations.jld2
"""

include(joinpath(@__DIR__, "..", "common.jl"))
bathdat = load(require_file(joinpath(@__DIR__, "..", "data", "BathsData",
    "enhancement_bath_sweep_precise", "etac0.jld2")))
markovianity = 8 .* abs.(bathdat["d"]) ./ abs2.(real.(bathdat["z"]))

data_path = require_file(joinpath(@__DIR__, "MarkovianityHEOMApproximations.jld2"))
res = load(data_path, "res")

labels = round.(
    [0, sum(markovianity[1:end-2]), sum(markovianity[1:end-1])],
    digits=4
)

get_k(f, pr; chi=1) = f ./ (chi .- (1 + chi) .* pr)

plt = plot()
for (i, label) in zip([14, 2, 1], labels)
    plot!(plt, res[i].time, get_k(res[i].f, res[i].pr) * 1e6, label=label, lw=2, alpha=0.6)
end
plot!(plt,
    ylim=(3.5, 4),
    xscale=:log10,
    xlim=(1e2, 1e4),
    legendfontsize=6,
    legendtitlefontsize=7,
    legend_background_color=:transparent,
    fg_legend=:transparent,
    legendtitle=raw"$\sum_{k\in M}\frac{8|d_k|}{\mathsf{Re}(z_k)^2}$",
    size=(300, 200),
    dpi=300,
    xlabel=raw"$\mathsf{time}\quad[\mathsf{fs}]$",
    ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$"
)

savefig(plt, joinpath(@__DIR__, "figS2.pdf"))
