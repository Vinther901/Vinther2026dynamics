"""
Reproduces `AdditionToCompetingPathways_main_plot.pdf`.

Branching ratio k_l/(k_l+k_r), k_r/(k_l+k_r) vs. cavity frequency, for the
"cold" reaction channel (k0l/k0r) of the near-degenerate triple-well system.
The rate values below are pre-computed HEOM results copied verbatim from the
source notebook (`ProbCurrent Figure Maker_latest.ipynb`); no external data
files are required to reproduce this particular figure.
"""

include(joinpath(@__DIR__, "..", "common.jl"))

omegacs = [0, 220, 270, 370, 420, 520, 530, 570, 670, 720, 800, 820, 870, 970, 1020, 1120, 1170, 1270, 1320, 1420, 1470, 1570, 1620, 1720, 1770, 1870, 1920, 2020, 2070, 2170]
k0l = [8.337941973170922e-8, 8.569656306258228e-8, 8.415014789177302e-8, 8.547726090310633e-8, 8.519771266192297e-8, 8.466982230562982e-8, 8.190345663060776e-8, 8.300492022708703e-8, 9.079699499592884e-8, 9.516369911252656e-8, 9.917942007974515e-8, 1.0026850360351777e-7, 9.34103790828767e-8, 8.869917015014817e-8, 8.55603627190414e-8, 8.953569995886163e-8, 8.833435383284136e-8, 9.388314577424422e-8, 9.69608026745865e-8, 1.2493335512244837e-7, 1.51889557485011e-7, 2.1469735154131193e-7, 2.0407346480792858e-7, 1.4086853003543274e-7, 1.2201710734844227e-7, 1.0370172207747732e-7, 9.805335964654704e-8, 9.350906151333254e-8, 8.763332062248001e-8, 8.452789936751738e-8]
k0r = [2.802691032457006e-7, 2.911285205661006e-7, 2.850739454009751e-7, 2.95358399583289e-7, 2.9832330009175595e-7, 3.142113903382259e-7, 3.0222203036005377e-7, 3.082492923492649e-7, 3.2175274730969954e-7, 3.370270511117063e-7, 3.4794870301648575e-7, 3.4905404397156677e-7, 3.2828004955130234e-7, 3.0301200527926653e-7, 2.8871715554175895e-7, 2.909166276256557e-7, 2.7979412315084224e-7, 2.7689949773594095e-7, 2.6911578142323453e-7, 2.6640149429180573e-7, 2.586337060624796e-7, 2.4226100865948255e-7, 2.460919363000388e-7, 2.6080199079041263e-7, 2.6567415106034415e-7, 2.7118256472409915e-7, 2.6997349134660386e-7, 2.733831593865171e-7, 2.640461897924406e-7, 2.61711377169858e-7]

plt = plot(
    legend=false,
    size=(140, 200),
    xticks=[500, 1500],
    xminorticks=5,
    yticks=5,
    yminorticks=5,
    ylim=(0, 1),
    xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
    ylabel=raw"$k_\leftrightharpoons/(k\!\!\!_\rightarrow + k\!\!\!_\leftarrow\!\!)$"
)
plot!(plt, omegacs[2:end], (k0l ./ (k0l .+ k0r))[2:end], color=:red, lw=2)
plot!(plt, omegacs[2:end], (k0r ./ (k0l .+ k0r))[2:end], color=:blue, lw=2)
annotate!(plt, 1600, 0.125, text(raw"$k\!\!\!_\leftarrow$", color=:red, pointsize=14))
annotate!(plt, 1600, 0.90, text(raw"$k\!\!\!_\rightarrow$", color=:blue, pointsize=14))

# savefig(plt, joinpath(@__DIR__, "AdditionToCompetingPathways_main_plot.pdf"))
savefig(plt, joinpath(@__DIR__, "fig7e.pdf"))
