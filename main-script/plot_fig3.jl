"""
Reproduces `SimpleNoiseSpectrums.pdf`.

Compares the noise spectrum S_beta(omega) of the phonon bath alone against
the "lossy cavity" contribution, evaluated at three representative
frequencies (the two well frequencies and the barrier frequency) of the
degenerate double-well/triple-well reaction coordinate.

Required data: `main-script/etac0.yaml`.
"""

include(joinpath(@__DIR__, "..", "common.jl"))

params = load_file(joinpath(@__DIR__, "etac0.yaml"))
params["η_c"] = 0.1

ws = LinRange(-500, 1600, 5000)
wslims = (ws[1], ws[end])
omegac_color = :mediumseagreen
S_lims = (-0.01, 6)

# cmap = vcat(palette(:tab10)..., palette(:tab20c)[end-3:end]...)
η_c = params["η_c"]

omega_l, omega_0, omega_r = 1200, 800, 400

function plot_S(plt)
    for (i, (omega, alpha)) in enumerate(zip([omega_r, omega_0, omega_l], [0.3, 0.5, 1.0]))
        params["ω_c"] = Dict("value" => omega, "unit" => "invcm")
        S = get_bath_func(params)

        params["η_c"] = 0.0
        S_nu = get_bath_func(params)
        params["η_c"] = η_c
        plot!(plt, [0], seriestype=:hline, color=:black, label=nothing)
        plot!(plt, ws, S.(convert_unit.(ws, :invcm, :au)) * 1e3,
            ylim=S_lims, xlim=wslims,
            label=i == 3 ? raw"$\mathsf{Lossy}\quad\!\!\!\mathsf{cavity}$" : nothing,
            legend=:topleft, legend_background_color=:transparent, fg_legend=:transparent,
            ls=i == 3 ? :solid : :dot, lw=i == 3 ? 3 : 1, alpha=alpha, color=omegac_color)
        plot!(plt, ws, S.(convert_unit.(ws, :invcm, :au)) * 1e3,
            fillrange=S_nu.(convert_unit.(ws, :invcm, :au)) * 1e3,
            label=nothing, lw=0, alpha=0.4 * alpha, color=omegac_color)
        plot!(plt, ws, S_nu.(convert_unit.(ws, :invcm, :au)) * 1e3,
            label=i == 3 ? raw"$\mathsf{P\!honon}$ $\mathsf{bath}$" : nothing,
            lw=3, alpha=alpha, color=:purple)
        plot!(plt, ws, S_nu.(convert_unit.(ws, :invcm, :au)) * 1e3,
            fillrange=0, label=nothing, alpha=0.4 * alpha, color=:purple)
        if i == 3
            plot!(plt, [omega], seriestype=:vline, color=omegac_color, ls=:dash, label=nothing)
        end
    end
end

plt = plot(size=(340, 180), dpi=300, legend=nothing,
    xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
    ylabel=raw"$S\!\!_\beta(\omega)\quad[\,\!\!\!\times 10^{-3} \mathsf{a.u.}]$",
    xminorticks=5)
plot_S(plt)
annotate!(plt, 1280, 5.57, (raw"$\omega_c$", 10, omegac_color))
plot!(plt, right_margin=2Plots.mm)

savefig(plt, joinpath(@__DIR__, "fig3.pdf"))
