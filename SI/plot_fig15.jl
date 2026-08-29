"""
Reproduces `ExampleOfKineticVsThermodynamicControl.pdf`.

Self-contained toy 3-state (left/middle/right) kinetics model illustrating
the difference between kinetically- and thermodynamically-controlled yields.
No external data files are required.
"""

include(joinpath(@__DIR__, "..", "common.jl"))
using OrdinaryDiffEqTsit5

w_l, w_r = 1570, 530
E_l, E_0, E_r = 1570, 800, 530

T = 208.5
beta = 1 / T

k_0l, k_0r = 0.1, 0.3

# pi_l = 0.31284
# pi_r = 0.35184
# pi_0 = 1 - pi_l - pi_r

# k_l0 = k_0l * pi_0 / pi_l
# k_r0 = k_0r * pi_0 / pi_r
k_l0 = 0.004
k_r0 = 0.006

K_l = k_0l + k_l0
K_r = k_0r + k_r0
Delta = sqrt((K_l + K_r)^2 - 4 * (K_l * K_r - k_0l * k_0r))
gamma_plus = 0.5 * (K_l + K_r + Delta)
gamma_minus = 0.5 * (K_l + K_r - Delta)

Y_l = k_0l * (gamma_plus - k_r0) / (gamma_plus * (gamma_plus - gamma_minus))
Y_r = k_0r * (gamma_plus - k_l0) / (gamma_plus * (gamma_plus - gamma_minus))

K = [
    -k_l0 k_0l 0.0
    k_l0 -(k_0l+k_0r) k_r0
    0.0 k_0r -k_r0
]

function master_equation!(dP, P, p, t)
    dP .= K * P
end

P0 = [0.0, 1.0, 0.0]
tspan = (0.0, 50.0)
problem = ODEProblem(master_equation!, P0, tspan)
solution = solve(problem, Tsit5())

function build_panel(; legend, solution)
    plt = plot(solution, label=[raw"$P_l$" raw"$P_0$" raw"$P_r$"],
        color=[:red :black :blue], linewidth=2, ylim=(0, 1))

    plot!(plt, [Y_l], seriestype=:hline, color=:red, ls=:dash,
        label=raw"$\frac{k_{0l}(\gamma_+-k_{r0})}{\gamma_+(\gamma_+-\gamma_-)}$")
    plot!(plt, [Y_r], seriestype=:hline, color=:blue, ls=:dash,
        label=raw"$\frac{k_{0r}(\gamma_+-k_{l0})}{\gamma_+(\gamma_+-\gamma_-)}$")
    plot!(plt, [3 / gamma_plus], seriestype=:vline, label=raw"$\mathcal{O}(\gamma_+^{-1})$", color=:darkgoldenrod)

    Gamma = k_0l * k_r0 + k_l0 * k_r0 + k_0r * k_l0
    plot!(plt, [k_0l * k_r0] / Gamma, seriestype=:hline, ls=:dot, color=:red, label=nothing)
    plot!(plt, [k_0r * k_l0] / Gamma, seriestype=:hline, ls=:dot, color=:blue, label=nothing)
    plot!(plt, [k_l0 * k_r0] / Gamma, seriestype=:hline, ls=:dot, color=:black, label=raw"$P_{(l,0,r)}^{eq}$")

    plot!(plt,
        legend=legend,
        size=(400, 300),
        dpi=300,
        xlabel=raw"$\mathsf{Time}$",
        ylabel=raw"$\mathsf{Population}$",
        xticks=[0, 25, 50],
        xminorticks=5,
        title=raw"$k_{0l}=" * "$(k_0l)" * raw"\quad " * raw"k_{0r}=" * "$(k_0r)" * raw"\quad " *
              raw"k_{l0}=" * "$(k_l0)" * raw"\quad " * raw"k_{r0}=" * "$(k_r0)" * raw"$",
        titlefontsize=11,
    )
    return plt
end

p2 = build_panel(legend=false, solution=solution)

k_0l *= 2
k_l0 *= 2

K_l = k_0l + k_l0
K_r = k_0r + k_r0
Delta = sqrt((K_l + K_r)^2 - 4 * (K_l * K_r - k_0l * k_0r))
gamma_plus = 0.5 * (K_l + K_r + Delta)
gamma_minus = 0.5 * (K_l + K_r - Delta)

Y_l = k_0l * (gamma_plus - k_r0) / (gamma_plus * (gamma_plus - gamma_minus))
Y_r = k_0r * (gamma_plus - k_l0) / (gamma_plus * (gamma_plus - gamma_minus))

K = [
    -k_l0 k_0l 0.0
    k_l0 -(k_0l+k_0r) k_r0
    0.0 k_0r -k_r0
]

function master_equation!(dP, P, p, t)
    dP .= K * P
end
problem = ODEProblem(master_equation!, P0, tspan)
solution = solve(problem, Tsit5())

p1 = build_panel(legend=:outerright, solution=solution)

plt = plot(p2, p1, layout=grid(1, 2, widths=[0.4, 0.6]), size=(1000, 300), dpi=300,
    bottom_margin=4Plots.mm,
    left_margin=5Plots.mm
)

# savefig(plt, joinpath(@__DIR__, "ExampleOfKineticVsThermodynamicControl.pdf"))
savefig(plt, joinpath(@__DIR__, "fig15.pdf"))
