"""
Reproduces `EquilibriumPopulations_TripleWells.pdf`.

Equilibrium (Gibbs + bath-induced correction) populations of the left/middle/
right wells of the near-degenerate triple-well system, as a function of
cavity frequency, with an inset zoom near the left-well population.

Required data (relative to FLATIRON_ROOT):
  - SystemsData/Data/degenerate_asymmetric_triple_wells/overtonic_morefinetuned.jld2
    (keys "eigvecs", "eigvals", "x", "V", "params")
  - BathsData/Data/enhancement_bath_sweep_precise/*.jld2   (keys "z", "d";
    filenames "etac0.jld2" or "omegac<value>.jld2")
"""

include(joinpath(@__DIR__, "..", "common.jl"))

# sysdat = load(require_file(joinpath(@__DIR__, "..", "data", "SystemsData", "Data",
#     "degenerate_asymmetric_triple_wells", "overtonic_morefinetuned.jld2")))
sysdat = load(joinpath(@__DIR__, "..", "data", "SystemsData", "overtonic_morefinetuned.jld2"))

N = 20
evecs = sysdat["eigvecs"][:, 1:N]
evals = sysdat["eigvals"][1:N]
# X = evecs' * diagm(sysdat["x"] .* exp.(-abs.(sysdat["x"]) ./ 2) .* 2) * sysdat["eigvecs"]
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

println("Beware, this can take some time to run, since it computes the Gibbs-corrected populations for each bath file. Expect ~15min (rough estimate)")

folder = joinpath(FLATIRON_ROOT, "BathsData", "enhancement_bath_sweep_precise")
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
omegacs = [dat.omegac for dat in bathdata]

taus = Dict(dat.omegac => get_gibbs_correction(dat.z, dat.d, X, dE, gibbs, beta; N) for dat in bathdata)

function get_well_definitions(params=sysdat["params"])
    ω_l = get_dimful_param(params["ω_l"])
    ω_r = get_dimful_param(params["ω_r"])
    ω_0 = get_dimful_param(params["ω_0"])
    E_0 = get_dimful_param(params["E_0"])
    E_l = get_dimful_param(params["E_l"])
    E_r = get_dimful_param(params["E_r"])
    x_l0 = -sqrt((ω_0 + ω_l) * (E_0 + E_l) / ω_0 / ω_l)
    x_r0 = sqrt((ω_0 + ω_r) * (E_0 + E_r) / ω_0 / ω_r)
    return (; x_l0, x_r0)
end

function get_proj_reacts(sysdat; NHilbert)
    x = sysdat["x"]
    evecs = sysdat["eigvecs"][:, 1:NHilbert]
    well_definitions = get_well_definitions()
    proj_left = ones(eltype(x), size(x))
    proj_left[x.>well_definitions.x_l0] .= 0
    proj_0 = ones(eltype(x), size(x))
    proj_0[x.<well_definitions.x_l0] .= 0
    proj_0[x.>well_definitions.x_r0] .= 0
    proj_right = ones(eltype(x), size(x))
    proj_right[x.<well_definitions.x_r0] .= 0
    proj_left = evecs' * diagm(proj_left) * evecs
    proj_0 = evecs' * diagm(proj_0) * evecs
    proj_right = evecs' * diagm(proj_right) * evecs
    return (; proj_left, proj_0, proj_right)
end
projs = get_proj_reacts(sysdat; NHilbert=N)

gPl = tr(gibbs * projs.proj_left)
gP0 = tr(gibbs * projs.proj_0)
gPr = tr(gibbs * projs.proj_right)

Pls = [tr(taus[omegac] * projs.proj_left) for omegac in omegacs]
P0s = [tr(taus[omegac] * projs.proj_0) for omegac in omegacs]
Prs = [tr(taus[omegac] * projs.proj_right) for omegac in omegacs]

plt = plot(
    legend=(0.17, 0.475),
    legend_background_color=:transparent,
    fg_legend=:transparent,
    legendfontsize=9,
    size=(300, 250),
    dpi=300,
    xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
    title=raw"$\mathsf{Triple}$ $\mathsf{well}$ $\mathsf{P\!ES}$"
)
plot!([gPl], seriestype=:hline, color=:red, ls=:dot, alpha=0.5, label=nothing)
plot!([gP0], seriestype=:hline, color=:black, ls=:dot, alpha=0.5, label=nothing)
plot!([gPr], seriestype=:hline, color=:blue, ls=:dot, alpha=0.5, label=nothing)
plot!([gPl + Pls[1]], seriestype=:hline, color=:red, alpha=0.5, label=nothing)
plot!([gP0 + P0s[1]], seriestype=:hline, color=:black, alpha=0.5, label=nothing)
plot!([gPr + Prs[1]], seriestype=:hline, color=:blue, alpha=0.5, label=nothing)
plot!(omegacs[2:end], gPl .+ Pls[2:end], color=:red, ls=:solid, lw=1.5, label=nothing)
plot!(omegacs[2:end], gP0 .+ P0s[2:end], color=:black, ls=:solid, lw=1.5, label=nothing)
plot!(omegacs[2:end], gPr .+ Prs[2:end], color=:blue, ls=:solid, lw=1.5, label=nothing)
curr_ylims = ylims(plt)
plot!([NaN], [NaN], color=:blue, label=raw"$P^{eq}_\mathsf{right}$")
plot!([NaN], [NaN], color=:black, label=raw"$P^{eq}_\mathsf{middle}$")
plot!([NaN], [NaN], color=:red, label=raw"$P^{eq}_\mathsf{left}$", ylim=curr_ylims)

yzoom = (0.3128, 0.3132)
plot!(plt, [gPl + Pls[1]]; seriestype=:hline, subplot=2,
    inset=(1, bbox(0.07, 0.27, 0.38, 0.45, :bottom, :right)),
    color=:red, alpha=0.5, label=nothing, legend=false, ylim=yzoom,
    framestyle=:box, xlabel="", ylabel="", tickfontsize=5)
plot!(plt, omegacs[2:end], gPl .+ Pls[2:end]; subplot=2, color=:red, ls=:solid, lw=1.5, label=nothing)
plot!(plt, [0, 2050, 2220], [0.313, 0.3255, 0.313], color=:black, lw=0.2, alpha=0.5, label=nothing)

savefig(plt, joinpath(@__DIR__, "fig13_right.pdf"))
