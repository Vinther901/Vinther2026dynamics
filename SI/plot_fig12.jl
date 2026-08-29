"""
Reproduces `DistributionOfRates.pdf`.

Two-row heatmap of the distribution of environment-driven rate constants
k(omega_c, phi) over the cavity-polarization angle phi, for several angles
theta between the cavity axis and the molecular dipole, compared for the
"isotropic-average" (top row) vs. "solid-angle-average" (bottom row)
projections.

Required data (relative to FLATIRON_ROOT):
  - HarryPlotter/EnvironmentDrivenDistributionOfRates.jld2
    (keys "res", "res_iso": Dict{theta => (; K::Matrix)})
"""

include(joinpath(@__DIR__, "..", "common.jl"))
using StatsBase

c = 0.38
Nphis = 200
omegacsFGR = LinRange(150, 2250, 150)
thetas = [0.05, 0.25, 0.5, 0.75, 0.95]
colors = Dict(theta => palette(:rainbow1, length(thetas))[i] for (i, theta) in enumerate(thetas))

# _loaded = load(require_file(joinpath(FLATIRON_ROOT, "HarryPlotter", "EnvironmentDrivenDistributionOfRates.jld2")))
_loaded = load(require_file(joinpath(@__DIR__, "EnvironmentDrivenDistributionOfRates.jld2")))
res = _loaded["res"]
res_iso = _loaded["res_iso"]

function plot_K(theta; res, label=nothing, label2=true)
    K = (convert_unit(0.214^2 .* res[theta].K, :au, :invfs) .* c .+ 3.794426570656991e-6) * 1e6
    kmin, kmax = extrema(K)
    nbins_k = Nphis ÷ 2
    k_edges = range(kmin, kmax; length=nbins_k + 1)
    k_centers = (k_edges[1:end-1] .+ k_edges[2:end]) ./ 2
    D = zeros(length(k_centers), length(omegacsFGR))
    for i in eachindex(omegacsFGR)
        h = fit(Histogram, K[i, :], k_edges)
        D[:, i] .= h.weights ./ sum(h.weights)
    end
    plt = plot(size=(300, 200), ylim=(3.8, 10), xlim=(125, 2275))
    heatmap!(plt, omegacsFGR, k_centers, D;
        xlabel="ω", ylabel="k", framestyle=:box, colorbar=nothing,
        color=cgrad([:white, colors[theta], colors[theta] .* 0.5]))
    plot!(plt, omegacsFGR, (D' * (1 ./ k_centers)) .^ (-1), color=:black, alpha=1, label=nothing)
    plot!(plt, omegacsFGR, mean(K, dims=2)[:], color=colors[theta], label=label,
        legendbackgroundcolor=:transparent, fg_legend=:transparent, legendfontsize=12, lw=1.5)
    plot!(plt, omegacsFGR, maximum(K, dims=2)[:], color=colors[theta], label=nothing, ls=:dash, lw=0.8)
    plot!(plt, omegacsFGR, minimum(K, dims=2)[:], color=colors[theta], label=nothing, ls=:dash, lw=0.8)
    plot!(plt, [NaN], [NaN], color=:black, alpha=1,
        label=label2 ? raw"$\langle k^{-1}\rangle_{\phi}^{-1}$" : nothing)
    for x in [500, 1000, 1500, 2000]
        vline!(plt, [x]; color=:gray80, linewidth=0.5, alpha=0.6, label=false)
    end
    for y in 4:1:10
        hline!(plt, [y]; color=:gray80, linewidth=0.5, alpha=0.6, label=false)
    end
    return plt
end

plots_top = [plot_K(theta; res=res, label=raw"$\frac{1}{\pi}\int_0^\pi k(\phi)d\phi$") for theta in thetas]
for (i, plt) in enumerate(plots_top[2:end])
    plot!(plt, left_margin=-12Plots.mm, ylabel=nothing, ytickfontcolor=:transparent,
        xtickfontcolor=:transparent, title=raw"$\theta = " * "$(thetas[i+1])" * raw" \pi$",
        top_margin=3Plots.mm, legend=(0.65, 0.8))
end
plot!(plots_top[1], ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$", left_margin=8Plots.mm,
    xtickfontcolor=:transparent, title=raw"$\theta = " * "$(thetas[1])" * raw" \pi$",
    top_margin=3Plots.mm, legend=(0.65, 0.8))
ptop = plot(plots_top..., layout=(1, 5), size=(300 * 5, 200), grid=false, link=:xy)

plots_bottom = [plot_K(theta; res=res_iso, label=raw"$\frac{1}{2}\int_0^\pi k(\phi)\sin(\phi)d\phi$", label2=false) for theta in thetas]
for plt in plots_bottom[2:end]
    plot!(plt, left_margin=-12Plots.mm, ylabel=nothing, ytickfontcolor=:transparent,
        xlabel=raw"$\omega_c\quad [\mathsf{cm}^{-1}]$", bottom_margin=8Plots.mm, legend=(0.5, 0.78))
end
plot!(plots_bottom[1], ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$", left_margin=8Plots.mm,
    xlabel=raw"$\omega_c\quad [\mathsf{cm}^{-1}]$", bottom_margin=8Plots.mm, legend=(0.5, 0.78))
pbottom = plot(plots_bottom..., layout=(1, 5), size=(300 * 5, 200), grid=false, link=:xy)

pltout = plot(ptop, plot(pbottom, top_margin=-6Plots.mm), layout=(2, 1), size=(300 * 5, 200 * 2))

# savefig(pltout, joinpath(@__DIR__, "DistributionOfRates.pdf"))
savefig(pltout, joinpath(@__DIR__, "fig12.pdf"))
