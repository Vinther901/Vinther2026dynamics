"""
Reproduces `NoiseSpectrumsWithOC_filled.pdf`.

Grid of panels (rows = theta, columns = two representative cavity
frequencies) showing the cavity-induced change in the effective noise
spectrum `S_beta(omega)` seen by the reaction coordinate.

Required data (relative to FLATIRON_ROOT):
  - BathsData/YAML_scripts/enhancement_bath_sweep_precise/etac0.yaml
  - BathsData/YAML_scripts/isotropic_baths_precise/theta_<theta>/omegac220.yaml,
    for theta in [0.05, 0.25, 0.5, 0.75, 0.95]
"""

include(joinpath(@__DIR__, "..", "common.jl"))
using Printf

thetas = [0.05, 0.25, 0.5, 0.75, 0.95]

bath0_params = load_file(require_file(joinpath(@__DIR__, "etac0.yaml")))
S0 = get_bath_func(bath0_params)

omegas = 800:0.5:1500
n = length(thetas)

p = plot(layout=(n, 2),
    link=:both,
    size=(400, 80n),
    dpi=300,
    legend=nothing,
    xlim=(900, 1400),
    xticks=[900, 1000, 1100, 1200, 1300],
    margin=-2Plots.mm)
ylim = (0, 0.00325 * 1e3)
colors = palette(:rainbow1, length(thetas))

for (i, theta) in enumerate(thetas)
    bath_params = load_file(require_file(joinpath(@__DIR__, "theta$(theta).yaml")))

    bath_params["ω_c"]["value"] = 1000
    S = get_bath_func(bath_params)
    DeltaS(omega) = (S(omega) - S0(omega)) * 1e3

    is_last = i == n
    label_left = is_last ? theta : @sprintf("%.2f", theta)
    xlabel_left = is_last ? raw"$\omega_c\quad [\mathsf{cm}^{-1}]$" : ""

    plot!(p[i, 1], omegas, DeltaS.(convert_unit(omegas, :invcm, :au)),
        lw=3, fillbetween=-0.1, alpha=0.3, color=colors[i], label=nothing)
    plot!(p[i, 1], omegas, DeltaS.(convert_unit(omegas, :invcm, :au)),
        lw=3, color=colors[i], label=label_left,
        xformatter=is_last ? :auto : (_ -> ""),
        legend=:topright, legend_background_color=:transparent, fg_legend=:transparent,
        ylim=ylim, xlabel=xlabel_left)

    bath_params["ω_c"]["value"] = 1180
    S = get_bath_func(bath_params)

    plot!(p[i, 2], omegas, DeltaS.(convert_unit(omegas, :invcm, :au)),
        lw=3, fillbetween=-0.1, alpha=0.3, color=colors[i], label=nothing)
    plot!(p[i, 2], omegas, DeltaS.(convert_unit(omegas, :invcm, :au)),
        lw=3, ylabel="", yformatter=(_ -> ""), color=colors[i],
        xformatter=is_last ? :auto : (_ -> ""), ylim=ylim, label=nothing, xlabel=xlabel_left)

    plot!(p[i, 1], seriestype=:vline, [1000], color=:mediumseagreen, ls=:dash, alpha=0.7, lw=2, label=nothing)
    plot!(p[i, 1], seriestype=:vline, [1170], color=:grey, ls=:dot, alpha=0.5, lw=2, label=nothing)

    plot!(p[i, 2], seriestype=:vline, [1180], color=:mediumseagreen, ls=:dash, alpha=0.7, label=raw"$\omega_c$", lw=2)
    plot!(p[i, 2], seriestype=:vline, [1170], color=:grey, ls=:dot, alpha=0.5, label=raw"$\omega_Q$", lw=2)
end

plot!(p[1, 2], legend=:topleft, legendfontsize=12, legend_background_color=:transparent, fg_legend=:transparent)
plot!(p[3, 1], ylabel=raw"$\left\langle S_\beta^{(c)}(\omega)\right\rangle\!\!\!_\phi\quad[\,\!\!\!\times 10^{-3} \mathsf{a.u.}]$")
plot!(p[1, 1], legend=:topright, legendtitle=raw"$\theta$ [$\pi$]",
    legend_background_color=:transparent, fg_legend=:transparent)

bath_params = load_file(require_file(joinpath(@__DIR__, "theta$(theta).yaml")))
bath_params["η_c"] = 1e-6
S = get_bath_func(bath_params)
DeltaS(omega) = (S(omega) - S0(omega)) * 1e3

plot!(p[3, 1], omegas, DeltaS.(convert_unit(omegas, :invcm, :au)), color=:black, alpha=0.5, lw=1.2, label="")
annotate!(p[3, 1], 1300, 1, text(raw"Outside\nCavity", 8))
annotate!(p[3, 1], 1210, 1.4, text(raw"$\leftarrow$", 9))

savefig(p, joinpath(@__DIR__, "fig8.pdf"))
