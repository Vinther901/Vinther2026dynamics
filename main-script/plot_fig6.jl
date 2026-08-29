#!/usr/bin/env julia

# Recreates the two-panel probability-current figure plt_tmp from
# "ProbCurrent Figure Maker_latest.ipynb" (cells 17, 23-25).

include(joinpath(@__DIR__, "..", "common.jl"))

# -----------------------------------------------------------------------------
# User-configurable paths
# -----------------------------------------------------------------------------
const FLATIRON_ROOT = "/mnt/home/jvinther/ceph/FlatironStudy"
const OUTPUT_FILE = joinpath(@__DIR__, "prob_current_figure.pdf")

# -----------------------------------------------------------------------------
# Notebook helper: extract the rate contributions for transition i -> j.
# -----------------------------------------------------------------------------
function get_rates_ij(data, i, j; t0=6000)
    ind = argmin(abs.(data.times .- t0))
    r_coh = data.rate_coh[i, j, ind] * 1e6
    r_nu  = (data.rate_Mnu[i, j, ind] + data.rate_nu[i, j, ind]) * 1e6
    r_c   = (data.rate_Mc[i, j, ind] + data.rate_c[i, j, ind]) * 1e6
    r = r_coh + r_nu + r_c
    return (; r, r_coh, r_nu, r_c)
end

# -----------------------------------------------------------------------------
# Notebook helper: stacked positive/negative bars for the three contributions.
# -----------------------------------------------------------------------------
function stacked_fillbar(p, rates;
    order = [:r_coh, :r_nu, :r_c],
    colors = Dict(
        :r_coh => :lightsalmon,
        :r_nu  => :purple,
        :r_c   => :mediumseagreen,
    ),
    labels = Dict(
        :r_coh => nothing,
        :r_nu  => nothing,
        :r_c   => nothing,
    ),
    xcenter = 1.0,
    width = 0.8,
    linecolor = :black,
    linealpha = 0.0,
    legend = :topright,
    kwargs...
)
    vals = Dict(
        :r_coh => rates.r_coh,
        :r_nu  => rates.r_nu,
        :r_c   => rates.r_c,
    )

    valid_keys = Set(keys(vals))
    if Set(order) != valid_keys
        error("`order` must contain exactly [:r_coh, :r_nu, :r_c] in some order.")
    end

    left  = xcenter - width / 2
    right = xcenter + width / 2

    pos_base = 0.0
    neg_base = 0.0

    for key in order
        v = float(vals[key])

        if v >= 0
            y0 = pos_base
            y1 = pos_base + v
            pos_base = y1
        else
            y0 = neg_base
            y1 = neg_base + v
            neg_base = y1
        end

        plot!(
            p,
            [left, right],
            [y1, y1],
            fillrange = [y0, y0],
            fillcolor = colors[key],
            fillalpha = 1.0,
            linecolor = linecolor,
            linealpha = linealpha,
            linewidth = 1,
            label = labels[key],
        )
    end

    return p
end

# -----------------------------------------------------------------------------
# Transition list and labels used by the notebook for the TripleWell figure.
# This is the active definition in cell 13 at the point where cells 23-25 run.
# -----------------------------------------------------------------------------
# Pairs = [
#     [9, 7],
#     [9, 8],
#     [7, 3],
#     [1, 5],
#     [5, 9],
#     [9, 10],
#     [10, 8],
#     [8, 6],
#     [6, 4],
#     [4, 2],
# ]

# Pair_labels = [
#     raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\leftarrow}$",
#     raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\rightarrow}$",
#     raw"$\mathsf{\textbf{L}}\!_\textbf{0}$",
#     raw"$\mathsf{\textbf{M}}\!_\textbf{0}$",
#     raw"$\mathsf{\textbf{M}}\!_\textbf{1}$",
#     raw"$\mathsf{\textbf{R}}\!_\textbf{0}$",
#     raw"$\mathsf{\textbf{R}}\!_\textbf{1}$",
#     raw"$\mathsf{\textbf{R}}\!_\textbf{2}$",
#     raw"$\mathsf{\textbf{R}}\!_\textbf{3}$",
#     raw"$\mathsf{\textbf{R}}\!_\textbf{4}$",
# ]
Pairs = [[1, 2], [4, 3], [1, 4], [3, 2]]
Pair_labels = [raw"$\mathsf{\textbf{T}}\!_\textbf{0}$", raw"$\mathsf{\textbf{T}}\!_\textbf{1}$",
    raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\uparrow}$", raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\downarrow}$"]

# The notebook uses this size when constructing plt_rates and plt_rates2.
rate_fig_size = (650, 150)#(650, 300)

# -----------------------------------------------------------------------------
# Top panel: ratio 2.0 data -> plt_rates
# -----------------------------------------------------------------------------
# EXP_top = "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0"
# RATE_TAG = "markovianity0.01_ktol1e-8_MixStates12and34_distinguished"

# rate_file_top = joinpath(
#     FLATIRON_ROOT,
#     "Experiments",
#     EXP_top,
#     "rate_data_tag" * RATE_TAG * ".jld2",
# )
rate_data_top = load("EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2", "data")

omegacs_top = sort([keys(rate_data_top)...])

plt_rates = plot(
    xticks = (
        collect(1:(length(Pairs) + 1) * length(omegacs_top) - 1),
        repeat(vcat(Pair_labels, ""), length(omegacs_top))[1:end-1],
    ),
    grid = false,
    size = rate_fig_size,
)
hline!(plt_rates, [0.0], color=:black, linewidth=0.2, alpha=1.0, label="")

for (n_omegac, omegac) in enumerate(omegacs_top)
    for (n, (i, j)) in enumerate(Pairs)
        rates = get_rates_ij(rate_data_top[omegac], i, j)

        if n != length(Pairs) || n_omegac != length(omegacs_top)
            stacked_fillbar(
                plt_rates,
                rates;
                xcenter = n + (n_omegac - 1) * (length(Pairs) + 1),
            )
        else
            stacked_fillbar(
                plt_rates,
                rates;
                xcenter = n + (n_omegac - 1) * (length(Pairs) + 1),
                labels = Dict(
                    :r_coh => "Unitary",
                    :r_nu  => "Phonons",
                    :r_c   => "Cavity",
                ),
            )
        end
    end

    if n_omegac != length(omegacs_top)
        plot!(
            plt_rates,
            [length(Pairs) + 1 + (n_omegac - 1) * (length(Pairs) + 1)],
            seriestype = :vline,
            color = :black,
            label = "",
        )
    end

    if omegac == 0
        annotate!(
            plt_rates,
            (length(Pairs) + 1) / 2 + (n_omegac - 1) * (length(Pairs) + 1),
            3.2,
            (raw"$\mathsf{Outside}\;\,\mathsf{cavity}$", 9, :black),
        )
    else
        annotate!(
            plt_rates,
            (length(Pairs) + 1) / 2 + (n_omegac - 1) * (length(Pairs) + 1),
            3.2,
            (raw"$\omega_c=" * "$(omegac)" * raw"\mathsf{cm}^{-1}$", 9, :black),
        )
    end
end

plot!(
    plt_rates,
    ylim = (-0.2, 3.4),
    left_margin = 3Plots.mm,
    legend = nothing,
    legend_background_color = :transparent,
    legend_foreground_color = :transparent,
    xtickfontcolor = :transparent,
    bottom_margin = -3Plots.mm,
    grid = true,
)

# -----------------------------------------------------------------------------
# Bottom panel: ratio 0.5 data -> plt_rates2
# -----------------------------------------------------------------------------
# EXP_bottom = "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio0.5"

# rate_file_bottom = joinpath(
#     FLATIRON_ROOT,
#     "Experiments",
#     EXP_bottom,
#     "rate_data_tag" * RATE_TAG * ".jld2",
# )
rate_data_bottom = load("EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio0.5_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2", "data")

omegacs_bottom = sort([keys(rate_data_bottom)...])

plt_rates2 = plot(
    ylabel = raw"                        $\mathsf{prob. curr.}\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
    xticks = (
        collect(1:(length(Pairs) + 1) * length(omegacs_bottom) - 1),
        repeat(vcat(Pair_labels, ""), length(omegacs_bottom))[1:end-1],
    ),
    grid = false,
    size = rate_fig_size,
)
hline!(plt_rates2, [0.0], color=:black, linewidth=0.2, alpha=1.0, label="")

for (n_omegac, omegac) in enumerate(omegacs_bottom)
    for (n, (i, j)) in enumerate(Pairs)
        rates = get_rates_ij(rate_data_bottom[omegac], i, j)

        if n != length(Pairs) || n_omegac != length(omegacs_bottom)
            stacked_fillbar(
                plt_rates2,
                rates;
                xcenter = n + (n_omegac - 1) * (length(Pairs) + 1),
            )
        else
            stacked_fillbar(
                plt_rates2,
                rates;
                xcenter = n + (n_omegac - 1) * (length(Pairs) + 1),
                labels = Dict(
                    :r_coh => raw"$\mathsf{Tunneling}$",
                    :r_nu  => raw"$\mathsf{Solvent}$",
                    :r_c   => raw"$\mathsf{Cavity}$",
                ),
            )
        end
    end

    if n_omegac != length(omegacs_bottom)
        plot!(
            plt_rates2,
            [length(Pairs) + 1 + (n_omegac - 1) * (length(Pairs) + 1)],
            seriestype = :vline,
            color = :black,
            label = "",
        )
    end
end

plot!(
    plt_rates2,
    ylim = (-3.4, 0.2),
    left_margin = 3Plots.mm,
    legend = (0.1, 0.82),
    legend_background_color = :transparent,
    legend_foreground_color = :transparent,
    top_margin = -3Plots.mm,
    yflip = true,
    grid = true,
)

annotate!(
    plt_rates2,
    length(omegacs_bottom) * (length(Pairs) + 1) / 2 + 0.3,
    -2.8,
    (raw"$\mathsf{Reverse}\ \mathsf{Reaction}$", 14, :black),
)

# -----------------------------------------------------------------------------
# Final combined figure: this is the notebook's plt_tmp (cell 25).
# -----------------------------------------------------------------------------
plt_tmp = plot(
    plt_rates,
    plt_rates2,
    link = :xy,
    layout = grid(2, 1),
    size=(rate_fig_size[1], rate_fig_size[2]*2),
    # left_margin = -3Plots.mm,

)

# Save the combined figure.
savefig(plt_tmp, joinpath(@__DIR__, "fig6.pdf"))