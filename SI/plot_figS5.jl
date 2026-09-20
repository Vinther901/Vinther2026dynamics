"""
Reproduces `EquilibriumPopulations_DoubleWells.pdf`.

Equilibrium (Gibbs + bath-induced correction) left/right well populations vs.
cavity frequency, for four asymmetric double-well systems (ratios 2, 4, 6, 8),
compared against the naive forward/backward-rate branching ratio
k_f/(k_f+k_b).

Required data (relative to FLATIRON_ROOT):
  - SystemsData/Data/degenerate_asymmetric_double_well_ratio{2.0,4.0,6.0,8.0}.jld2
    (keys "eigvecs", "eigvals", "x")
  - BathsData/Data/enhancement_bath_sweep_precise/*.jld2   (keys "z", "d")
  - main-script/ReactionRates.jld2   (key "data")
"""

include(joinpath(@__DIR__, "..", "common.jl"))

folder = joinpath(@__DIR__, "..", "data", "BathsData", "enhancement_bath_sweep_precise")
files = sort(filter(endswith(".jld2"), readdir(folder; join=true)))
bathdata = [
    let dct = load(f), name = splitext(basename(f))[1]
        omegac =
            startswith(name, "omegac") ? parse(Float64, name[7:end]) :
            startswith(name, "etac0") ? 0.0 :
            error("Could not parse omegac from $name")
        (; omegac, z=dct["z"], d=dct["d"])
    end
    for f in files
]
sort!(bathdata, by=x -> x.omegac)

println("Beware, this will take some time to run, since it computes the Gibbs-corrected populations for each bath file. Expect ~60min? (rough estimate)")
println("It's possible to reduce the number of bath files to speed this up, but then the plot will be less smooth. See the commented-out line below.")
# bathdata = bathdata = bathdata[[1,2,length(bathdata)÷3,2*length(bathdata)÷3,end]]  # <-- uncomment this line to reduce the number of bath files (and speed up the computation ~1min)
omegacs = [dat.omegac for dat in bathdata]

gPl_dct, gPr_dct, Pls_dct, Prs_dct = Dict(), Dict(), Dict(), Dict()

data = load(require_file(joinpath(@__DIR__, "..", "main-script", "ReactionRates.jld2")), "data")

for ratio in [2.0, 4.0, 6.0, 8.0]
    sysdat = load(require_file(joinpath(FLATIRON_ROOT, "SystemsData",
        "degenerate_asymmetric_double_well_ratio$(ratio).jld2")))

    N = 20
    evecs = sysdat["eigvecs"][:, 1:N]
    evals = sysdat["eigvals"][1:N]
    X = evecs' * diagm(sysdat["x"]) * sysdat["eigvecs"]

    T = convert_unit(300, :K, :au)
    beta = 1 / T
    weights = exp.(-beta .* evals)
    weights ./= sum(weights)
    gibbs = diagm(weights)
    dE = evals .- evals'

    Pleft = zeros(length(sysdat["x"]))
    Pleft[sysdat["x"].<=0] .= 1
    Pleft = diagm(Pleft)
    Pright = I - Pleft
    Pl = evecs' * Pleft * evecs
    Pr = evecs' * Pright * evecs

    taus = Dict(dat.omegac => get_gibbs_correction(dat.z, dat.d, X, dE, gibbs, beta; N) for dat in bathdata)

    gPl_dct[ratio] = tr(gibbs * Pl)
    gPr_dct[ratio] = tr(gibbs * Pr)
    Pls_dct[ratio] = [tr(taus[omegac] * Pl) for omegac in omegacs]
    Prs_dct[ratio] = [tr(taus[omegac] * Pr) for omegac in omegacs]
end


function get_k_tilde(dataf, datab)
    kf, kb = dataf.k, datab.k
    omegacsf, omegacsb = dataf.omegacs, datab.omegacs
    common = intersect(omegacsf, omegacsb)
    idxf = Dict(v => i for (i, v) in pairs(omegacsf))
    idxb = Dict(v => i for (i, v) in pairs(omegacsb))
    kf = [kf[idxf[w]] for w in common]
    kb = [kb[idxb[w]] for w in common]
    return kf ./ (kf .+ kb), kb ./ (kf .+ kb), common
end

function plot_pops(gPl, gPr, Pls, Prs, dataf, datab)
    plt = plot(size=(300, 200), dpi=300, xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$", legend=nothing)
    kftilde, kbtilde, common_omegacs = get_k_tilde(dataf, datab)
    common_omegacs[1] == 0 || error("expected common_omegacs to start at 0")
    plot!(plt, [kftilde[1]], color=:red, ls=:dash, seriestype=:hline, alpha=0.5, label=nothing)
    plot!(plt, [kbtilde[1]], color=:blue, ls=:dash, seriestype=:hline, alpha=0.5, label=nothing)
    plot!(plt, common_omegacs[2:end], kftilde[2:end], color=:red, ls=:dash, label=nothing)
    plot!(plt, common_omegacs[2:end], kbtilde[2:end], color=:blue, ls=:dash, label=nothing)
    plot!(plt, [gPl], seriestype=:hline, color=:red, ls=:dot, alpha=0.5, label=nothing)
    plot!(plt, [gPr], seriestype=:hline, color=:blue, ls=:dot, alpha=0.5, label=nothing)
    plot!(plt, [gPl + Pls[1]], seriestype=:hline, color=:red, alpha=0.5, label=nothing)
    plot!(plt, [gPr + Prs[1]], seriestype=:hline, color=:blue, alpha=0.5, label=nothing)
    plot!(plt, omegacs[2:end], gPl .+ Pls[2:end], color=:red, ls=:solid, lw=1.5, label=nothing)
    plot!(plt, omegacs[2:end], gPr .+ Prs[2:end], color=:blue, ls=:solid, lw=1.5, label=nothing)
    return plt
end

plts = []
for (i, (ratio, inv_ratio)) in enumerate(zip([2.0, 4.0, 6.0, 8.0], [0.5, 0.25, 0.167, 0.125]))
    push!(plts, plot!(
        plot_pops(gPl_dct[ratio], gPr_dct[ratio], Pls_dct[ratio], Prs_dct[ratio], data[inv_ratio], data[ratio]),
        title=raw"$\omega_r\approx " * "$([950, 700, 570, 500][i])" * raw"\mathsf{cm}^{-1}$"
    ))
end
curr_ylims = ylims(plts[1])

plot!(plts[1], [0], fill_between=[0], label=raw"$P^{eq}_\mathsf{left}$", color=:red, seriestype=:hline, lw=0)
plot!(plts[1], [0], fill_between=[0], label=raw"$P^{eq}_\mathsf{right}$", color=:blue, seriestype=:hline, lw=0)
plot!(plts[1], [1000, 1001], [0, 0], label=raw"$\tau$", color=:black, ls=:dot, lw=1.5)
plot!(plts[1], [1000, 1001], [0, 0], label=raw"$\tau+\tau_\mathsf{MF}^{(2)}$", color=:black, ls=:dash, lw=1.5)
plot!(plts[1], [1000, 1001], [0, 0], label=raw"$\frac{k_\leftarrow}{k_\leftarrow + k_\rightarrow}$", color=:red, ls=:dash)
plot!(plts[1], [1000, 1001], [0, 0], label=raw"$\frac{k_\rightarrow}{k_\leftarrow + k_\rightarrow}$", color=:blue, ls=:dash)
plot!(plts[1], legend=:inside, legend_background_color=:transparent, fg_legend=:transparent, legendfontsize=9)

for i in [2, 3, 4]
    plot!(plts[i], left_margin=-12Plots.mm, ytickfontcolor=:transparent, ylabel=nothing)
end

plt = plot(plts..., size=(800, 400), layout=(1, 4), link=:y, dpi=300, bottom_margin=5Plots.mm, ylim=curr_ylims)

savefig(plt, joinpath(@__DIR__, "figS5.pdf"))
