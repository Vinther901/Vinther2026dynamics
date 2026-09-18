"""
Reproduce Figure 16 (`fig16.pdf`), the reaction rate as a function of phonon
coupling strength. The plotted data are embedded below, so no input file is
required.
"""

include(joinpath(@__DIR__, "..", "common.jl"))

etanu = [0.002, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
k0 = [6.28749960237726e-5, 7.314401399124137e-6, 4.882930636072713e-6, 3.7471262387744043e-6, 3.650830697366973e-6, 4.60380478014773e-6, 6.032602212477521e-6, 6.800492991192163e-6, 6.98745326165169e-6, 6.786675324118696e-6, 6.268119574597424e-6, 5.5566057564870156e-6, 4.657020247222081e-6]

plt = plot(
    etanu,
    k0 .* 1e6;
    color=:black,
    marker=:o,
    ylabel=raw"$k_0$ [ $\times 10^{-6}$ fs$^{-1}$]",
    xlabel=raw"$\eta_{\nu}$",
    ylim=(3, 10),
    linewidth=2,
    markersize=4,
    size=(300, 200),
    dpi=300,
    legend=nothing,
)
savefig(plt, joinpath(@__DIR__, "fig16.pdf"))
