"""
Reproduce Fig. S3 with the corrected 3D orientation distribution in the bottom row.

Run from any directory:
    julia /path/to/Vinther2026dynamics/SI/plot_figS3_new.jl

The top row uses the original planar data. The bottom row recomputes the FGR
rates at each (φ, ψ), with uniform midpoint sampling of μ = cos(φ) and ψ.
θ is the fixed angle between the reaction-coordinate and spectator dipoles.
The full distribution is retained before taking arithmetic/harmonic means.

Outputs SI/figS3_new.pdf. The original scripts and data are not overwritten.
Set FIGS3_NMU and FIGS3_NPSI to adjust angular resolution (default: 200 each).
Frequency integration uses rtol=1e-8 and the original 0–0.1 a.u. interval.
"""

include(joinpath(@__DIR__, "..", "common.jl"))
using StatsBase

c = 0.38
Nphis = 200
omegacsFGR = LinRange(150, 2250, 150)
thetas = [0.05, 0.25, 0.5, 0.75, 0.95]
colors = Dict(theta => palette(:rainbow1, length(thetas))[i] for (i, theta) in enumerate(thetas))

_loaded = load(require_file(joinpath(@__DIR__, "EnvironmentDrivenDistributionOfRates.jld2")))
res = _loaded["res"]
nmu = parse(Int, get(ENV, "FIGS3_NMU", "200"))
npsi = parse(Int, get(ENV, "FIGS3_NPSI", "200"))
result_3d = rate_distribution_3d(omegacsFGR, thetas .* π; nmu, npsi, progress=true)
res_3d = Dict(theta => result_3d.rates[theta * π] for theta in thetas)

function plot_K(theta; res, label=nothing, label2=true, solid_angle=false)
    K = (convert_unit(0.214^2 .* res[theta].K, :au, :invfs) .* c .+ 3.794426570656991e-6) * 1e6
    kmin, kmax = extrema(K)
    nbins_k = Nphis ÷ 2
    # Include the largest sample in the final left-closed histogram bin.
    k_edges = range(kmin, solid_angle ? nextfloat(kmax) : kmax; length=nbins_k + 1)
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
    # For the new distribution, compute the harmonic mean from the rates,
    # avoiding any dependence on the histogram bin width.
    harmonic = solid_angle ? vec(1 ./ mean(1 ./ K, dims=2)) : (D' * (1 ./ k_centers)) .^ (-1)
    plot!(plt, omegacsFGR, harmonic, color=:black, alpha=1, label=nothing)
    plot!(plt, omegacsFGR, mean(K, dims=2)[:], color=colors[theta], label=label,
        legendbackgroundcolor=:transparent, fg_legend=:transparent, legendfontsize=12, lw=1.5)
    plot!(plt, omegacsFGR, maximum(K, dims=2)[:], color=colors[theta], label=nothing, ls=:dash, lw=0.8)
    plot!(plt, omegacsFGR, minimum(K, dims=2)[:], color=colors[theta], label=nothing, ls=:dash, lw=0.8)
    plot!(plt, [NaN], [NaN], color=:black, alpha=1,
        label=solid_angle ? raw"$\langle k^{-1}\rangle_{\Omega}^{-1}$" :
            (label2 ? raw"$\langle k^{-1}\rangle_{\phi}^{-1}$" : nothing))
    for x in [500, 1000, 1500, 2000]
        vline!(plt, [x]; color=:gray80, linewidth=0.5, alpha=0.6, label=false)
    end
    for y in 4:1:10
        hline!(plt, [y]; color=:gray80, linewidth=0.5, alpha=0.6, label=false)
    end
    return plt
end

plots_top = [plot_K(theta; res=res, label=raw"$\frac{1}{2\pi}\int_0^{2\pi} k(\phi)d\phi$") for theta in thetas]
for (i, plt) in enumerate(plots_top[2:end])
    plot!(plt, left_margin=-12Plots.mm, ylabel=nothing, ytickfontcolor=:transparent,
        xtickfontcolor=:transparent, title=raw"$\theta = " * "$(thetas[i+1])" * raw" \pi$",
        top_margin=3Plots.mm, legend=(0.63, 0.8))
end
plot!(plots_top[1], ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$", left_margin=8Plots.mm,
    xtickfontcolor=:transparent, title=raw"$\theta = " * "$(thetas[1])" * raw" \pi$",
    top_margin=3Plots.mm, legend=(0.63, 0.8))
ptop = plot(plots_top..., layout=(1, 5), size=(300 * 5, 200), grid=false, link=:xy)

plots_bottom = [plot_K(theta; res=res_3d, label=raw"$\langle k\rangle_{\Omega}$", label2=false, solid_angle=true) for theta in thetas]
for plt in plots_bottom[2:end]
    plot!(plt, left_margin=-12Plots.mm, ylabel=nothing, ytickfontcolor=:transparent,
        xlabel=raw"$\omega_c\quad [\mathsf{cm}^{-1}]$", bottom_margin=8Plots.mm, legend=(0.63, 0.8))
end
plot!(plots_bottom[1], ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$", left_margin=8Plots.mm,
    xlabel=raw"$\omega_c\quad [\mathsf{cm}^{-1}]$", bottom_margin=8Plots.mm, legend=(0.63, 0.8))
pbottom = plot(plots_bottom..., layout=(1, 5), size=(300 * 5, 200), grid=false, link=:xy)

pltout = plot(ptop, plot(pbottom, top_margin=-6Plots.mm), layout=(2, 1), size=(300 * 5, 200 * 2))

savefig(pltout, joinpath(@__DIR__, "figS3.pdf"))

@info "Saved corrected figure" path=joinpath(@__DIR__, "figS3.pdf")
