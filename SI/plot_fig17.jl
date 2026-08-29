#!/usr/bin/env julia

# ============================================================================
# make_prob_current_figure_full.jl
#
# Reconstructs the combined probability-current figure from the notebook:
#
#     [ system / potential panel ] | [ rate panels: forward + reverse ]
#
# Main entry point:
#
#     plt_out = make_probability_current_figure(
#         sysdata, rate_data1, rate_data2, Pairs, Pair_labels;
#         output_file="prob_current_figure.pdf",
#     )
#
# IMPORTANT DATA CONTRACT FOR sysdata
# ------------------------------------
# The notebook's original `sysdat` only contains the system data such as
# `sysdat["x"]` and `sysdat["V"]`.  The system plot also needs `states` and
# the already shifted/converted `energies`, which the notebook stores as
# separate variables.
#
# Therefore this script expects `sysdata` to be a NamedTuple containing:
#
#     sysdata.sysdat       # original notebook sysdat dictionary/object
#     sysdata.states       # state wavefunctions / probability-density data
#     sysdata.energies     # energies in cm^-1 after the notebook conversion
#
# For example:
#
#     sysdata = (
#         sysdat          = sysdat,
#         states          = states,
#         energies        = energies,
#         E_displacement  = E_displacement,
#     )
#
# The two rate_data arguments should be the same nested structures used by
# the notebook, i.e. a mapping from cavity frequency (`omegac`) to a rate-data
# object containing fields such as `times`, `rate_coh`, `rate_Mnu`, `rate_nu`,
# `rate_Mc`, and `rate_c`.
# ============================================================================

include(joinpath(@__DIR__, "..", "common.jl"))
using LinearAlgebra
using Plots

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""
Extract the rate contributions for transition i -> j at the time point nearest
`t0`.  This follows the notebook's scaling to 10^-6 fs^-1.
"""
function get_rates_ij(data, i, j; t0=6000)
    ind = argmin(abs.(data.times .- t0))

    r_coh = data.rate_coh[i, j, ind] * 1e6
    r_nu  = (data.rate_Mnu[i, j, ind] + data.rate_nu[i, j, ind]) * 1e6
    r_c   = (data.rate_Mc[i, j, ind] + data.rate_c[i, j, ind]) * 1e6

    r = r_coh + r_nu + r_c
    return (; r, r_coh, r_nu, r_c)
end


"""Draw a stacked bar for the three rate contributions."""
function stacked_fillbar(
    p,
    rates;
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
)
    vals = Dict(
        :r_coh => rates.r_coh,
        :r_nu  => rates.r_nu,
        :r_c   => rates.r_c,
    )

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
# Manual arrow implementation copied from the notebook/pasted source
# -----------------------------------------------------------------------------

"""
Draw a manual arrow from `(x1,y1)` to `(x2,y2)` using a fixed visual-size
arrowhead.  This is the notebook's GR-safe implementation.
"""
function draw_arrow!(
    plt,
    x1,
    y1,
    x2,
    y2;
    head_len=0.035,
    head_width=0.025,
    use_plot_size=true,
    kwargs...
)
    xl = Plots.xlims(plt)
    yl = Plots.ylims(plt)

    xr = xl[2] - xl[1]
    yr = yl[2] - yl[1]

    xr == 0 && return plt
    yr == 0 && return plt

    W, H = if use_plot_size
        try
            sz = plt.attr[:size]
            float(sz[1]), float(sz[2])
        catch
            1.0, 1.0
        end
    else
        1.0, 1.0
    end

    to_scaled(x, y) = (
        W * (x - xl[1]) / xr,
        H * (y - yl[1]) / yr,
    )

    from_scaled(xs, ys) = (
        xl[1] + xr * xs / W,
        yl[1] + yr * ys / H,
    )

    x1s, y1s = to_scaled(x1, y1)
    x2s, y2s = to_scaled(x2, y2)

    dx = x2s - x1s
    dy = y2s - y1s
    L = hypot(dx, dy)

    L == 0 && return plt

    ux = dx / L
    uy = dy / L

    px = -uy
    py = ux

    ref_size = min(W, H)
    hl = min(head_len * ref_size, 0.85L)
    hw = min(head_width * ref_size, 0.85L)

    bxs = x2s - hl * ux
    bys = y2s - hl * uy

    xtri_s = [
        x2s,
        bxs + hw / 2 * px,
        bxs - hw / 2 * px,
        x2s,
    ]

    ytri_s = [
        y2s,
        bys + hw / 2 * py,
        bys - hw / 2 * py,
        y2s,
    ]

    tri_data = from_scaled.(xtri_s, ytri_s)
    xtri = first.(tri_data)
    ytri = last.(tri_data)

    bx_data, by_data = from_scaled(bxs, bys)

    kw = Dict{Symbol,Any}(kwargs)

    arrow_color = get(kw, :color, get(kw, :linecolor, :black))
    arrow_alpha = get(kw, :alpha, get(kw, :linealpha, 1.0))

    shaft_kw = copy(kw)
    shaft_kw[:label] = get(shaft_kw, :label, false)

    plot!(plt, [x1, bx_data], [y1, by_data]; shaft_kw...)

    head_kw = copy(kw)
    head_kw[:seriestype] = :shape
    head_kw[:label] = false
    head_kw[:fillcolor] = get(head_kw, :fillcolor, arrow_color)
    head_kw[:linecolor] = get(head_kw, :linecolor, arrow_color)
    head_kw[:fillalpha] = get(head_kw, :fillalpha, arrow_alpha)
    head_kw[:linealpha] = get(head_kw, :linealpha, arrow_alpha)

    if !haskey(head_kw, :lw) && !haskey(head_kw, :linewidth)
        head_kw[:lw] = 0
    end

    plot!(plt, xtri, ytri; head_kw...)
    return plt
end


"""Offset the endpoints of an arrow by a fixed visual radius."""
function offset_arrow_endpoints(
    plt,
    x1,
    y1,
    x2,
    y2;
    r=0.02,
    use_plot_size=true,
)
    xl = Plots.xlims(plt)
    yl = Plots.ylims(plt)

    xr = xl[2] - xl[1]
    yr = yl[2] - yl[1]

    xr == 0 && return x1, y1, x2, y2
    yr == 0 && return x1, y1, x2, y2

    W, H = if use_plot_size
        try
            sz = plt.attr[:size]
            float(sz[1]), float(sz[2])
        catch
            1.0, 1.0
        end
    else
        1.0, 1.0
    end

    to_scaled(x, y) = (
        W * (x - xl[1]) / xr,
        H * (y - yl[1]) / yr,
    )

    from_scaled(xs, ys) = (
        xl[1] + xr * xs / W,
        yl[1] + yr * ys / H,
    )

    x1s, y1s = to_scaled(x1, y1)
    x2s, y2s = to_scaled(x2, y2)

    dx = x2s - x1s
    dy = y2s - y1s
    L = hypot(dx, dy)

    L == 0 && return x1, y1, x2, y2

    ux = dx / L
    uy = dy / L

    r_scaled = r * min(W, H)

    x1s2 = x1s + r_scaled * ux
    y1s2 = y1s + r_scaled * uy
    x2s2 = x2s - r_scaled * ux
    y2s2 = y2s - r_scaled * uy

    x1p, y1p = from_scaled(x1s2, y1s2)
    x2p, y2p = from_scaled(x2s2, y2s2)

    return x1p, y1p, x2p, y2p
end


"""Draw labeled arrows between indexed eigenstates."""
function draw_state_arrows!(
    plt,
    pairs,
    labels;
    states,
    energies,
    x,
    arrow_color=:grey,
    arrow_lw=1.3,
    arrow_alpha=1.0,
    arrow_head_len=0.28,
    arrow_head_width=0.25,
    label_fontsize=10,
    label_color=:black,
    label_dy=70.0,
    label_dx=0.0,
    loc=0.5,
    kwargs...
)
    @assert length(pairs) == length(labels) "pairs and labels must have the same length."

    state_xcenter(i) = begin
        ψ = states[:, i]
        ρ = abs2.(ψ)
        sum(x .* ρ) / sum(ρ)
    end

    for (pair, lab) in zip(pairs, labels)
        i, j = pair

        x1 = state_xcenter(i)
        x2 = state_xcenter(j)

        y1 = energies[i] + 50
        y2 = energies[j] + 50

        x1o, y1o, x2o, y2o = offset_arrow_endpoints(
            plt,
            x1,
            y1,
            x2,
            y2;
            r=0.03,
        )

        draw_arrow!(
            plt,
            x1o,
            y1o,
            x2o,
            y2o;
            head_len=arrow_head_len,
            head_width=arrow_head_width,
            color=arrow_color,
            lw=arrow_lw,
            alpha=arrow_alpha,
            label=false,
            kwargs...
        )

        xm = loc * x1 + (1 - loc) * x2 + label_dx
        ym = loc * y1 + (1 - loc) * y2 + label_dy

        annotate!(
            plt,
            xm,
            ym,
            text(lab, label_fontsize, label_color, :center),
        )
    end

    return plt
end


"""Notebook wrapper used to draw the arrows three times with different styles."""
function wrapper_draw_arrows(
    plt_sys,
    Pairs,
    Pair_labels;
    states,
    energies,
    x,
    arrow_color=:darkgoldenrod,
    arrow_color2=:black,
    arrow_lw=2,
    arrow_alpha=0.1,
    arrow_head_len=0.103,
    arrow_head_width=0.07,
    label_dy=100.0,
    label_dx=0.0,
    fillalpha=0.0,
    label_fontsize=10,
    label_color=:darkgoldenrod,
    label_color2=:black,
    loc=0.5,
)
    draw_state_arrows!(
        plt_sys,
        Pairs,
        Pair_labels;
        states,
        energies,
        x,
        arrow_color,
        arrow_lw,
        arrow_alpha,
        arrow_head_len,
        arrow_head_width,
        label_dy,
        label_dx,
        fillalpha,
        label_fontsize=label_fontsize + 1,
        label_color,
        loc,
    )

    draw_state_arrows!(
        plt_sys,
        Pairs,
        Pair_labels;
        states,
        energies,
        x,
        arrow_color,
        arrow_lw,
        arrow_alpha,
        arrow_head_len,
        arrow_head_width,
        label_dy,
        label_dx,
        fillalpha,
        label_fontsize=label_fontsize - 1,
        label_color,
        loc,
    )

    draw_state_arrows!(
        plt_sys,
        Pairs,
        Pair_labels;
        states,
        energies,
        x,
        arrow_color=arrow_color2,
        arrow_lw=arrow_lw / 2,
        arrow_alpha=1.0,
        arrow_head_len,
        arrow_head_width,
        label_dy,
        label_dx,
        fillalpha,
        label_fontsize,
        label_color=label_color2,
        loc,
    )

    return plt_sys
end


# -----------------------------------------------------------------------------
# System panel: direct reconstruction of the notebook's `plot_sys` cell
# -----------------------------------------------------------------------------

function make_system_plot(
    sysdat,
    Pairs,
    Pair_labels;
    rate=25,
    ws=LinRange(-850, 2250, 5000),
    num_states=10,
    sys_ratio=0.05,
    sys_xlims=(-3.5, 5.1),
    sys_fig_size=(300, 300),
    yticks=[-500, 0, 530, 800, 1000, 1250, 1570, 2000],
    cmap=palette(:tab10),
    mix_states,
)
    states, energies = get_states_and_energies(mix_states, sysdat["eigvecs"], sysdat["eigvals"], sysdat["x"])
    E_displacement = energies[1]
    energies = convert_unit.(energies .- E_displacement, :au, :invcm)

    wslims = (ws[1], ws[end])

    x_full = sysdat["x"]
    V_full = sysdat["V"]

    x = x_full[1:rate:end]
    V = copy(V_full[1:rate:end])

    V .-= E_displacement
    V = convert_unit.(V, :au, :invcm)
    
    plt_sys = plot(
        size=sys_fig_size,
        dpi=300,
        yticks=yticks,
    )

    for i in 1:num_states
        prob_density = copy(states[1:rate:end, i])
        maxabs = maximum(abs.(prob_density))
        maxabs == 0 && continue

        prob_density ./= maxabs
        prob_density *= sys_ratio * (wslims[2] - wslims[1])
        prob_density *= sign(prob_density[argmax(abs.(prob_density))])

        energy_offset = energies[i]
        y_vals = prob_density .+ energy_offset
        color = cmap[i]

        plot!(
            plt_sys,
            x,
            y_vals,
            fillrange=energy_offset,
            fillalpha=0.3,
            color=color,
            label=(i <= 5) ? i : nothing,
        )

        hline!(
            plt_sys,
            [energy_offset],
            color=color,
            alpha=0.3,
            label=false,
        )
    end

    plot!(
        plt_sys,
        x,
        V,
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
        alpha=0.5,
    )

    # Exact active arrow settings from the notebook.
    arrow_alpha = 0.5
    arrow_head_len = 0.02
    arrow_head_width = 0.014
    label_dy = 30
    fillalpha = 1.0
    label_color = :white

    n = length(Pairs)

    if n >= 2
        wrapper_draw_arrows(
            plt_sys,
            Pairs[1:2],
            Pair_labels[1:2];
            states,
            energies,
            x=x_full,
            arrow_color=:white,
            arrow_lw=3,
            arrow_alpha,
            arrow_head_len,
            arrow_head_width,
            label_dy=-160.0,
            fillalpha,
            label_fontsize=12,
            label_color,
        )
    end

    if n >= 3
        wrapper_draw_arrows(
            plt_sys,
            Pairs[[3]],
            Pair_labels[[3]];
            states,
            energies,
            x=x_full,
            arrow_color=:white,
            arrow_lw=3,
            arrow_alpha,
            arrow_head_len,
            arrow_head_width,
            label_dy=0.0,
            fillalpha,
            label_fontsize=12,
            label_color,
        )
    end

    if n >= 4
        hi = min(5, n)
        wrapper_draw_arrows(
            plt_sys,
            Pairs[4:hi],
            Pair_labels[4:hi];
            states,
            energies,
            x=x_full,
            arrow_color=:white,
            arrow_lw=3,
            arrow_alpha,
            arrow_head_len,
            arrow_head_width,
            label_dy=0.0,
            label_dx=0.1,
            fillalpha,
            label_fontsize=12,
            label_color,
        )
    end

    if n >= 6
        wrapper_draw_arrows(
            plt_sys,
            Pairs[6:end],
            Pair_labels[6:end];
            states,
            energies,
            x=x_full,
            arrow_color=:white,
            arrow_lw=3,
            arrow_alpha,
            arrow_head_len,
            arrow_head_width,
            label_dy,
            label_dx=-0.3,
            fillalpha,
            label_fontsize=12,
            label_color,
        )
    end

    return plt_sys
end


# -----------------------------------------------------------------------------
# Rate panels
# -----------------------------------------------------------------------------

function make_rate_plot(
    rate_data,
    Pairs,
    Pair_labels;
    rate_fig_size=(650, 300),
    reverse=false,
    t0=6000,
)
    omegacs = sort(collect(keys(rate_data)))

    n_pairs = length(Pairs)
    n_omegac = length(omegacs)

    ticks = collect(1:(n_pairs + 1) * n_omegac - 1)
    ticklabels = repeat(vcat(Pair_labels, ""), n_omegac)[1:end-1]

    plt = plot(
        xticks=(ticks, ticklabels),
        grid=false,
        size=rate_fig_size,
    )

    if reverse
        plot!(
            plt,
            ylabel=raw"                        $\mathsf{prob. curr.}\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
        )
    end

    hline!(plt, [0.0], color=:black, linewidth=0.2, alpha=1.0, label="")

    for (n_ω, ωc) in enumerate(omegacs)
        for (n_pair, (i, j)) in enumerate(Pairs)
            rates = get_rates_ij(rate_data[ωc], i, j; t0)
            xcenter = n_pair + (n_ω - 1) * (n_pairs + 1)

            if n_pair == n_pairs && n_ω == n_omegac
                if reverse
                    labels = Dict(
                        :r_coh => raw"$\mathsf{Tunneling}$",
                        :r_nu  => raw"$\mathsf{Solvent}$",
                        :r_c   => raw"$\mathsf{Cavity}$",
                    )
                else
                    labels = Dict(
                        :r_coh => "Unitary",
                        :r_nu  => "Phonons",
                        :r_c   => "Cavity",
                    )
                end
            else
                labels = Dict(
                    :r_coh => nothing,
                    :r_nu  => nothing,
                    :r_c   => nothing,
                )
            end

            stacked_fillbar(
                plt,
                rates;
                xcenter,
                labels,
            )
        end

        if n_ω != n_omegac
            plot!(
                plt,
                [n_pairs + 1 + (n_ω - 1) * (n_pairs + 1)],
                seriestype=:vline,
                color=:black,
                label="",
            )
        end

        # Cavity labels only appear on the forward panel in the notebook.
        if !reverse
            xmid = (n_pairs + 1) / 2 + (n_ω - 1) * (n_pairs + 1)
            if ωc == 0
                annotate!(
                    plt,
                    xmid,
                    3.2,
                    (raw"$\mathsf{Outside}\;\,\mathsf{cavity}$", 9, :black),
                )
            else
                annotate!(
                    plt,
                    xmid,
                    3.2,
                    (raw"$\omega_c=" * "$(ωc)" * raw"\mathsf{cm}^{-1}$", 9, :black),
                )
            end
        end
    end

    if reverse
        plot!(
            plt,
            ylim=(-3.4, 0.2),
            left_margin=3Plots.mm,
            legend=(0.1, 0.82),
            legend_background_color=:transparent,
            legend_foreground_color=:transparent,
            top_margin=-7Plots.mm,
            yflip=true,
            grid=true,
        )

        annotate!(
            plt,
            n_omegac * (n_pairs + 1) / 2 + 0.3,
            -2.8,
            (raw"$\mathsf{Reverse}\ \mathsf{Reaction}$", 14, :black),
        )
    else
        plot!(
            plt,
            ylim=(-0.2, 3.4),
            left_margin=3Plots.mm,
            legend=nothing,
            legend_background_color=:transparent,
            legend_foreground_color=:transparent,
            xtickfontcolor=:transparent,
            bottom_margin=-7Plots.mm,
            grid=true,
        )
    end

    return plt
end


# -----------------------------------------------------------------------------
# Main figure builder
# -----------------------------------------------------------------------------

"""
Build the full figure from:

    sysdata
    rate_data1
    rate_data2
    Pairs
    Pair_labels

The result is returned as a Plots.jl object and optionally written to
`output_file`.
"""
function make_probability_current_figure(
    sysdata,
    rate_data1,
    rate_data2,
    Pairs,
    Pair_labels,
    mix_states;
    output_file="prob_current_figure.pdf",
    rate=25,
    ws=LinRange(-850, 2250, 5000),
    num_states=10,
    sys_ratio=0.05,
    sys_xlims=(-3.5, 5.1),
    sys_fig_size=(300, 300),
    rate_fig_size=(650, 300),
    yticks=[-500, 0, 530, 800, 1000, 1250, 1570, 2000],
    cmap=palette(:tab10),
    t0=6000,
)
    @assert length(Pairs) == length(Pair_labels) "Pairs and Pair_labels must have the same length."

    plt_sys = make_system_plot(
        sysdata,
        Pairs,
        Pair_labels;
        rate,
        ws,
        num_states,
        sys_ratio,
        sys_xlims,
        sys_fig_size,
        yticks,
        cmap,
        mix_states
    )

    plt_rates = make_rate_plot(
        rate_data1,
        Pairs,
        Pair_labels;
        rate_fig_size,
        reverse=false,
        t0,
    )

    plt_rates2 = make_rate_plot(
        rate_data2,
        Pairs,
        Pair_labels;
        rate_fig_size,
        reverse=true,
        t0,
    )

    # Exact final composition from the notebook.
    w, h = sys_fig_size[1] + rate_fig_size[1], sys_fig_size[2]

    plt_tmp = plot(
        plt_rates,
        plt_rates2,
        link=:xy,
        layout=grid(2, 1),
        left_margin=-3Plots.mm,
    )

    plt_out = plot(
        plot(plt_sys, right_margin=-3Plots.mm),
        plt_tmp,
        layout=grid(
            1,
            2,
            widths=[
                sys_fig_size[1] / w,
                (w - sys_fig_size[1]) / w,
            ],
        ),
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        size=(w, h),
        dpi=300,
    )

    if output_file !== nothing
        savefig(plt_out, output_file)
        println("Saved figure to: ", abspath(output_file))
    end

    return plt_out
end


# ============================================================================
# Example call
# ============================================================================
#
# Keep this commented if using the file via `include(...)`.
#
# sysdata = (
#     sysdat          = sysdat,
#     states          = states,
#     energies        = energies,
#     E_displacement  = E_displacement,
# )
#
# plt_out = make_probability_current_figure(
#     sysdata,
#     rate_data1,
#     rate_data2,
#     Pairs,
#     Pair_labels;
#     output_file="prob_current_figure.pdf",
# )
# ============================================================================

ratios = reverse([2.0, 4.0, 6.0, 8.0])
inv_ratios = reverse([0.5, 0.25, 0.167, 0.125])
# ratios = [2.0]
# inv_ratios = [0.5]
sysdata = Dict(
    ratio => load(require_file(joinpath(@__DIR__, "..", "data", "SystemsData",
        "degenerate_asymmetric_double_well_ratio$(ratio).jld2")))
    for ratio in ratios
)

rate_datas1 = Dict(
    ratio => load(require_file(joinpath(@__DIR__, "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio$(ratio)_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and45_distinguished.jld2")), "data")
    for ratio in ratios[1:end-1]
)

rate_datas1[2.0] = load(require_file(joinpath(@__DIR__, "..", "main-script",
    "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2")), "data")

rate_datas2 = Dict(
    ratio => load(require_file(joinpath(@__DIR__,
        "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio$(inv_ratio)_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and45_distinguished.jld2")), "data")
    for (ratio, inv_ratio) in zip(ratios[1:end-1], inv_ratios[1:end-1])
)

rate_datas2[2.0] = load(require_file(joinpath(@__DIR__, "..", "main-script",
    "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio0.5_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2")), "data")


# rate_datas1 = Dict(2.0 => load(require_file(joinpath(@__DIR__, "..", "main-script",
#     "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2")), "data")
# )

# rate_datas2 = Dict(2.0 => load(require_file(joinpath(@__DIR__, "..", "main-script",
#     "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio0.5_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2")), "data")
# )

for (key, sysdat) in sysdata
    i = Dict(2.0 => 1, 4.0 => 2, 6.0 => 3, 8.0 => 4)[key]
    # wslims = (-850, 2150)
    # sys_xlims = (-2, 3.5)
    # num_states = 7
    # sys_ratio = 0.05
    # sys_fig_size = (250, 300)
    mix_states = [
        [[1, 2], [3, 4]],
        [[1, 2], [4, 5]],
        [[1, 2], [4, 5]],
        [[1, 2], [4, 5]],
    ][i]
    output_file = [
        "fig17a.pdf",
        "fig17b.pdf",
        "fig17c.pdf",
        "fig17d.pdf",
    ][i]
    
    rates = [get_rates_ij(rate_datas1[key][0],i,j;t0=6000).r for i in 1:10 for j in (i+1):10]

    Pairs = [[i,j] for i in 1:10 for j in (i+1):10]
    for i in 1:length(Pairs)
        if sign(rates[i]) < 0
            reverse!(Pairs[i])
        end
    end
    Pairs = Pairs[abs.(rates) .> 5e-2]
    Alphabet = [raw"$A$", raw"$B$", raw"$C$", raw"$D$", raw"$E$", raw"$F$", raw"$G$", raw"$H$", raw"$I$", raw"$J$", raw"$K$", raw"$L$", raw"$M$", raw"$N$"]
    Pair_labels = [Alphabet[i] for i in 1:length(Pairs)]

    make_probability_current_figure(
        sysdat,
        rate_datas1[key],
        rate_datas2[key],
        Pairs,
        Pair_labels,
        mix_states;
        output_file=output_file,
        rate=25,
        ws=LinRange(-850, 2250, 5000),
        num_states=7,
        sys_ratio=0.05,
        sys_xlims=(-2, 3.5),
        sys_fig_size=(300, 300),
        rate_fig_size=(650, 300),
        yticks=[-500, 0, 500, 1000, 1200, 1500, 2000],
        cmap=palette(:tab10),
        t0=6000,
    )
end

