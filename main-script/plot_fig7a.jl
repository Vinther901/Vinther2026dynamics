"""
Reproduce Figure 7(a) (`fig7a.pdf`), the triple-well potential and localized
eigenstates. Requires `data/SystemsData/overtonic_morefinetuned.jld2`.
"""

include(joinpath(@__DIR__, "..", "common.jl"))
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

########################################################################################

sysdata = load(require_file(joinpath(@__DIR__, "..", "data", "SystemsData","overtonic_morefinetuned.jld2")))

Pairs = [
    [9,7],
    [9,8],
    [7,3],
    [1,5],
    [5,9],
    [9,10],[10,8],[8,6],[6,4],[4,2]]
Pair_labels = [
    raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\leftarrow}$", 
    raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\rightarrow}$", 
    raw"$\mathsf{\textbf{L}}\!_\textbf{0}$", 
    raw"$\mathsf{\textbf{M}}\!_\textbf{0}$", 
    raw"$\mathsf{\textbf{M}}\!_\textbf{1}$",
    raw"$\mathsf{\textbf{R}}\!_\textbf{0}$",
    raw"$\mathsf{\textbf{R}}\!_\textbf{1}$",
    raw"$\mathsf{\textbf{R}}\!_\textbf{2}$",
    raw"$\mathsf{\textbf{R}}\!_\textbf{3}$",
    raw"$\mathsf{\textbf{R}}\!_\textbf{4}$",
        ]

plt_sys = make_system_plot(
        sysdata,
        Pairs,
        Pair_labels;
        rate=25,
        ws=LinRange(-850, 2250, 5000),
        num_states=10,
        sys_ratio=0.05,
        sys_xlims=(-3.5, 5.1),
        sys_fig_size=(300, 300),
        yticks = [-500,0,530,800,1000,1250,1570,2000],
        cmap=palette(:tab10),
        mix_states=[[1,2,3],[7,8,9]]
    )

savefig(plt_sys, joinpath(@__DIR__, "fig7a.pdf"))


    
