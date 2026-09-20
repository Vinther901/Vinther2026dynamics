"""
Reproduces `PhononDressedEffectiveTripleWellPES.pdf`.

Shows the "phonon-dressed" effective potential energy surface of the
near-degenerate triple-well system: eigenstate probability densities plotted
at their Gibbs-corrected (bath-renormalized) energies.

Required data (relative to FLATIRON_ROOT):
  - SystemsData/Data/degenerate_asymmetric_triple_wells/overtonic_morefinetuned.jld2
    (keys "eigvecs", "eigvals", "x", "V")
  - BathsData/Data/enhancement_bath_sweep_precise/etac0.jld2   (keys "z", "d")
  - BathsData/Data/distinguished_enhancement_bath/omegac1570.jld2   (keys "z", "d")
"""

include(joinpath(@__DIR__, "..", "common.jl"))

sysdat = load(require_file(joinpath(@__DIR__, "..", "data", "SystemsData", "overtonic_morefinetuned.jld2")))

bath0dat = load(require_file(joinpath(@__DIR__, "..", "data", "BathsData",
    "enhancement_bath_sweep_precise", "etac0.jld2")))
z0, d0 = bath0dat["z"], bath0dat["d"]

bathdat = load(require_file(joinpath(@__DIR__, "..", "data", "BathsData",
    "distinguished_enhancement_bath", "omegac1570.jld2")))
z1, d1 = bathdat["z"], bathdat["d"]
z = vcat(z0, z1)
d = vcat(d0, d1)

N = 30
evecs = sysdat["eigvecs"][:, 1:N]
evals = sysdat["eigvals"][1:N]

X = evecs' * diagm(sysdat["x"]) * sysdat["eigvecs"]
T = convert_unit(300, :K, :au)
beta = 1 / T
weights = exp.(-beta .* evals)
weights ./= sum(weights)
gibbs = diagm(weights)

dE = evals .- evals'

tau = get_gibbs_correction(z, d, X, dE, gibbs, beta; N)

tmp = eigen(gibbs + tau)
rate = 1
wslims = (-850, 2250)
sys_xlims = (-3.5, 5.1)
num_states = 10
sys_ratio = 0.05
sys_fig_size = (300, 300)

cmap = vcat(palette(:tab10)..., palette(:tab20c)[end-3:end]...)
states, energies = evecs * tmp.vectors[:, end:-1:1], -1 / beta .* log.(abs.(tmp.values[end:-1:1]))

states, energies = get_states_and_energies([[8, 9]], states, energies, sysdat["x"])
E_displacement = evals[1]
energies = convert_unit.(energies .- energies[1], :au, :invcm)

function plot_sys(plt)
    x = sysdat["x"][1:rate:end]
    V = copy(sysdat["V"][1:rate:end])
    V .-= E_displacement
    V = convert_unit.(V, :au, :invcm)

    for i in 1:num_states
        prob_density = copy(states[1:rate:end, i])
        prob_density ./= maximum(abs.(prob_density))
        prob_density *= sys_ratio * (wslims[2] - wslims[1])
        prob_density *= sign(prob_density[argmax(abs.(prob_density))])

        energy_offset = energies[i]
        y_vals = prob_density .+ energy_offset
        color = cmap[i]
        plot!(plt, x, y_vals, fillrange=energy_offset, fillalpha=0.3, color=color, label=(i <= 5) ? i : nothing)
        hline!(plt, [energy_offset], color=color, alpha=0.3, label=false)
    end

    plot!(plt, x, V,
        ylim=wslims,
        legend=nothing,
        legend_background_color=:transparent,
        fg_legend=:transparent,
        xlim=sys_xlims,
        color=:black,
        label=nothing,
        xlabel=raw"$R\quad[\mathsf{a.\!u.\!\!}]$",
        ylabel=raw"$\mathsf{Frequency}\quad [\mathsf{cm}^{-1}]$",
        lw=2, alpha=0.5
    )
end
yticks = [-500, 0, 530, 800, 1000, 1250, 1570, 2000]
plt_sys = plot(size=sys_fig_size, dpi=300, yticks=yticks)
plot_sys(plt_sys)
plot!(plt_sys,
    title=raw"$\mathsf{P\!honon}$ $\mathsf{dressed}$ $\mathsf{effective}$ $\mathsf{P\!ES}$",
    titlefontsize=10,
    ylabel=raw"$\mathsf{Frequency}\quad [\mathsf{cm}^{-1}]$",
    bottom_margin=-2Plots.mm
)

savefig(plt_sys, joinpath(@__DIR__, "figS4_left.pdf"))
