"""
Reproduces `ProofOfMarkovianCorrection.pdf`.

Shows the effective rate constant k(t) computed from a HEOM simulation for a
few different numbers of auxiliary density operators (ADOs), demonstrating
convergence to the Markovian value as more bath modes are included.

Required data (relative to FLATIRON_ROOT, see common.jl):
  - BathsData/Data/enhancement_bath_sweep_precise/etac0.jld2   (keys "z", "d")
  - Experiments/MarkovianityHEOMApproximations2/Data/*_renormedjump.jld2
    (each with key "obs", a DataFrame with columns :observables, :times)
"""

include(joinpath(@__DIR__, "..", "common.jl"))
# using DataFrames
# using Suppressor: @suppress_err

bathdat = load(require_file(joinpath(@__DIR__, "..", "data", "BathsData",
    "enhancement_bath_sweep_precise", "etac0.jld2")))
markovianity = 8 .* abs.(bathdat["d"]) ./ abs2.(real.(bathdat["z"]))

# data_path = joinpath(FLATIRON_ROOT, "Experiments", "MarkovianityHEOMApproximations2", "Data")
data_path = joinpath(@__DIR__, "MarkovianityHEOMApproximations.jld2")

res = load(data_path, "res")
# res = Dict()
# @suppress_err for file in readdir(data_path, join=true)
#     if !endswith(file, "_renormedjump.jld2")
#         continue
#     end
#     num_ados = parse(Int, match(r"include(.*)_renormedjump\.jld2", basename(file))[1])

#     obs = load(file, "obs")
#     pr = [o["PR"] for o in obs[!, :observables]]
#     f = [o["F"] for o in obs[!, :observables]]
#     time = obs[!, :times]

#     res[num_ados] = (; time, f, pr)
# end

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

# savefig(plt, joinpath(@__DIR__, "ProofOfMarkovianCorrection.pdf"))
savefig(plt, joinpath(@__DIR__, "fig11.pdf"))
