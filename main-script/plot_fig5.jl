"""
Reproduces `asymmetric_transition_rates2.pdf`.

Forward/backward reaction rate constants vs. cavity frequency for four
asymmetric double-well systems (well-frequency ratios 2, 4, 6, 8), each shown
alongside an inset of the mixed eigenstates on the potential energy surface.

Required data (relative to FLATIRON_ROOT):
  - HarryPlotter/PlottingTools/ReactionRates.jld2  (key "data": Dict keyed by
    ratio/inverse-ratio, each with fields `.omegacs`, `.k`)
  - SystemsData/Data/degenerate_asymmetric_double_well_ratio{2.0,4.0,6.0,8.0}.jld2
    (keys "eigvecs", "eigvals", "x", "V")
"""

include(joinpath(@__DIR__, "..", "common.jl"))

# data = load(require_file(joinpath(FLATIRON_ROOT, "HarryPlotter", "PlottingTools", "ReactionRates.jld2")), "data")
data = load(require_file(joinpath(@__DIR__, "ReactionRates.jld2")), "data")

ratios = reverse([2.0, 4.0, 6.0, 8.0])
sysdata = Dict(
    ratio => load(require_file(joinpath(FLATIRON_ROOT, "SystemsData",
        "degenerate_asymmetric_double_well_ratio$(ratio).jld2")))
    for ratio in ratios
)

plts_sys = Dict()
rate = 1
for (key, sysdat) in sysdata
    i = Dict(2.0 => 1, 4.0 => 2, 6.0 => 3, 8.0 => 4)[key]
    wslims = (-850, 2150)
    sys_xlims = (-2, 3.5)
    num_states = 7
    sys_ratio = 0.05
    sys_fig_size = (250, 300)
    mix_states = [
        [[1, 2], [3, 4]],
        [[1, 2], [4, 5]],
        [[1, 2], [4, 5]],
        [[1, 2], [4, 5]],
    ][i]

    cmap = palette(:tab10)
    cmap = [cmap, cmap, cmap[[1, 2, 3, 5, 4, 6, 7]], cmap[[1, 2, 3, 5, 4, 6, 7]]][i]

    states, energies = get_states_and_energies(mix_states, sysdat["eigvecs"], sysdat["eigvals"], sysdat["x"])
    E_displacement = energies[1]
    energies = convert_unit.(energies .- E_displacement, :au, :invcm)

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
            lw=2,
            alpha=0.5
        )
    end
    yticks = [-500, 0, 500, 1000, 1200, 1500, 2000]
    plts_sys[key] = plot(size=sys_fig_size, dpi=300, yticks=yticks)
    plot_sys(plts_sys[key])
    plot!(plts_sys[key], ytickfontcolor=:transparent, ylabel=nothing)
    plot!(plts_sys[key], left_margin=-12Plots.mm)
end

plts = []
pairs = reverse([[0.125, 8.0], [0.167, 6.0], [0.25, 4.0], [0.5, 2.0]])

plot!(plts_sys[2.0],
    left_margin=0Plots.mm,
    ylabel=raw"$\mathsf{Frequency}\quad [\mathsf{cm}^{-1}]$",
    ytickfontcolor=:black,
    xtickfontcolor=:black,
)

for (i, pair) in enumerate(pairs)
    plt = plot(
        size=(550, 300),
        dpi=300,
        xlim=(350, 2000),
        ylim=(1, 4.5),
        legend=nothing,
        xmirror=true,
        xticks=[600, 1200, 1800],
        xminorticks=3,
        yminorticks=1,
        minorgrid=true,
        xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$"
    )

    ob, kb = copy(data[pair[1]].omegacs), copy(data[pair[1]].k * 1e6)
    of, kf = copy(data[pair[2]].omegacs), copy(data[pair[2]].k * 1e6)

    # Interpolation for the scarce backward rates near the peak, which is based on the forward rates rather than linear interpolation.
    # This is purely an aesthetic choice to make the plot look nicer than an uninformed linear interpolation. No physical meaning is implied by this interpolation.
    ind = argmin(abs.(of .- ob[argmax(kb)]))
    kb2 = kf[ind-3:ind+3] * maximum(kb) / kf[ind]
    push!(ob, of[ind-3:ind+3]...)
    push!(kb, kb2...)
    perm = sortperm(ob)
    kb = kb[perm]
    ob = ob[perm]

    plot!(plt, kb[[1]], seriestype=:hline, color=:red, lw=0.5, label=nothing)
    plot!(plt, kf[[1]], seriestype=:hline, color=:blue, lw=0.5, label=nothing)
    plot!(plt, ob[2:end], kb[2:end], color=:red, lw=1.5, label=raw"$k\!\!\!_\leftarrow$")
    plot!(plt, of[2:end], kf[2:end], color=:blue, lw=1.5, label=raw"$k\!\!\!_\rightarrow$")

    if i > 1
        plot!(plt, left_margin=-13Plots.mm, ytickfontcolor=:transparent)
    end
    if i == 2
        plot!(plt, legend=(0.75, 0.83), legend_background_color=:transparent,
            legend_foreground_color=:transparent, legendfontsize=14)
    end
    if i == 1
        plot!(plt, ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$")
    end

    push!(plts, plot(plt, plts_sys[pair[2]], layout=(2, 1)))
end

plt = plot(plts..., layout=(1, 4), size=(550, 400), dpi=300,
    top_margin=-1Plots.mm,
)

# savefig(plt, joinpath(@__DIR__, "asymmetric_transition_rates2.pdf"))
savefig(plt, joinpath(@__DIR__, "fig5.pdf"))
