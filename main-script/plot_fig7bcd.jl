"""
Creates the combined, publication-width Figure 7(b-e).

Stacked bar chart of the coherent/phonon/cavity contributions to the
transition-rate "probability current" between pairs of eigenstates of the
ratio-2.0 asymmetric double well, at three cavity frequencies, with the
forward/backward rate constants k_l/k_r and branching ratios shown above.

Required data (relative to FLATIRON_ROOT), for
tag = "markovianity0.01_ktol1e-8_MixStates12and34_distinguished":
  - Experiments/EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0/
    rate_data_tag<tag>.jld2   (key "data": Dict{omegac => rate_data} with
    fields .rate_coh, .rate_nu, .rate_Mnu, .rate_c, .rate_Mc, .times)
  - Experiments/EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0/
    YAML_scripts/*<tag>.yaml

The forward/backward rate arrays (k0l, k0r) below are pre-computed HEOM
results copied verbatim from the source notebook.
"""

include(joinpath(@__DIR__, "..", "common.jl"))

# const EXPERIMENTS_DIR = joinpath(FLATIRON_ROOT, "Experiments")
# const BATH_YAML_DIR = joinpath(FLATIRON_ROOT, "BathsData", "YAML_scripts")

# EXP = "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0"
# rate_tag = "markovianity0.01_ktol1e-8_MixStates12and34_distinguished"

# function get_bath_and_sys_params(EXP, tag)
#     exp_dir = joinpath(EXPERIMENTS_DIR, EXP)
#     yaml_files = filter(f -> endswith(f, tag * ".yaml"), readdir(joinpath(exp_dir, "YAML_scripts"); join=true))
#     isempty(yaml_files) && error("No YAML config matching tag '$tag' found in $(joinpath(exp_dir, "YAML_scripts"))")
#     exp_params = load_file(yaml_files[1])
#     params = load_file(require_file(joinpath(BATH_YAML_DIR, exp_params["bath"] * ".yaml")))
#     sysdat = load(require_file(joinpath(FLATIRON_ROOT, "SystemsData", "Data", exp_params["system"] * ".jld2")))
#     return params, sysdat
# end

# rate_data = load(require_file(joinpath(EXPERIMENTS_DIR, EXP * "_rate_data_tag" * rate_tag * ".jld2")), "data")
rate_data = load(joinpath(@__DIR__, "EnhancementSweep_etanu0.1_gammanu200_etac0.1_TripleWell_NearDegenerate_4_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates123and789_distinguished.jld2"), "data")
# _, sysdat = get_bath_and_sys_params(EXP, rate_tag)
sysdat = load(joinpath(@__DIR__, "..", "data", "SystemsData", "overtonic_morefinetuned.jld2"))

# ratio2.0 double well: the two lowest tunneling doublets, mixed pairwise
states, energies = get_states_and_energies([[1, 2], [3, 4]], sysdat["eigvecs"], sysdat["eigvals"], sysdat["x"])
E_displacement = energies[1]
rate = 25
sys_xlims = (-3.5, 5.1)

# Pairs = [[1, 2], [4, 3], [1, 4], [3, 2]]
# Pair_labels = [raw"$\mathsf{\textbf{T}}\!_\textbf{0}$", raw"$\mathsf{\textbf{T}}\!_\textbf{1}$",
#     raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\uparrow}$", raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\downarrow}$"]
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

# Final PDF dimensions. Plots.jl uses points for vector output:
# 510 pt = 7.08 in = 179.9 mm, a standard two-column publication width.
const FIGURE_SIZE = (650, 350)
const TOP_HEIGHT = 165
const BOTTOM_HEIGHT = FIGURE_SIZE[2] - TOP_HEIGHT

# omegacs = sort(collect(keys(rate_data)))
omegacs = [0, 530, 800, 1570]
plt_rates = plot(
    xticks=(collect(1:(length(Pairs)+1)*length(omegacs)-1), repeat(vcat(Pair_labels, ""), length(omegacs))[1:end-1]),
    grid=false,
    xlim=(0, (length(Pairs)+1)*length(omegacs)-0.5),
    size=(FIGURE_SIZE[1], BOTTOM_HEIGHT),
    tickfontsize=7,
    guidefontsize=10,
    legendfontsize=8
)
hline!(plt_rates, [0.0], color=:black, linewidth=0.2, alpha=1.0, label="")

for (n_omegac, omegac) in enumerate(omegacs)
    for (n, (i, j)) in enumerate(Pairs)
        rates = get_rates_ij(rate_data[omegac], i, j)
        if n != length(Pairs) || n_omegac != length(omegacs)
            stacked_fillbar!(plt_rates, rates; xcenter=n + (n_omegac - 1) * (length(Pairs) + 1))
        else
            stacked_fillbar!(plt_rates, rates; xcenter=n + (n_omegac - 1) * (length(Pairs) + 1),
                labels = Dict(
                            :r_coh => raw"$\mathsf{Tunneling}$",
                            :r_nu  => raw"$\mathsf{Solvent}$",
                            :r_c   => raw"$\mathsf{Cavity}$",
                        ))
        end
    end
    if n_omegac != length(omegacs)
        plot!(plt_rates, [length(Pairs)+1+(n_omegac-1)*(length(Pairs)+1)], seriestype=:vline, color=:black, label="")
    end
    label_text = omegac == 0 ? raw"$\mathsf{Outside}\;\,\mathsf{cavity}$" : raw"$\omega_c=" * "$(omegac)" * raw"\mathsf{cm}^{-1}$"
    annotate!(plt_rates, (length(Pairs) + 1) / 2 + (n_omegac - 1) * (length(Pairs) + 1), 0.72, (label_text, 9, :black))
end

plot!(plt_rates,
    ylim=(-0.05, 0.76),
    legend=(0.19, 0.70),
    left_margin=2Plots.mm,
    right_margin=1Plots.mm,
    bottom_margin=2Plots.mm,
    xlabel=raw"$\mathsf{Transitions}$",
    ylabel = raw"$\mathsf{prob. curr.}\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
    # legend=nothing,
    legend_background_color=:transparent,
    legend_foreground_color=:transparent,
    # xtickfontcolor=:transparent,
    # bottom_margin=-7Plots.mm,
    grid=true
)

# function plot_V!(plt; subplot, inset)
#     x = sysdat["x"][1:rate:end]
#     V = copy(sysdat["V"][1:rate:end])
#     V .-= E_displacement
#     V = convert_unit.(V, :au, :invcm)
#     plot!(plt, x, V;
#         subplot, inset,
#         legend=false, ticks=false, frame=:box,
#         ylim=(-900, 3000), xlim=(-3.2, 5),
#         color=:grey, lw=2, alpha=1.0, grid=false,
#         yticks=nothing, xticks=nothing, showaxis=false)
# end
# plot_V!(plt_rates; subplot=2, inset=(1, bbox(0.205, 0.12, 0.1, 0.2)))
# plot_V!(plt_rates; subplot=3, inset=(1, bbox(0.53, 0.12, 0.1, 0.2)))
# plot_V!(plt_rates; subplot=4, inset=(1, bbox(0.855, 0.12, 0.1, 0.2)))

# Forward/backward rate constants vs. cavity frequency (pre-computed HEOM results)
k_omegacs = [0, 220, 270, 370, 420, 520, 530, 570, 670, 720, 800, 820, 870, 970, 1020, 1120, 1170, 1270, 1320, 1420, 1470, 1570, 1620, 1720, 1770, 1870, 1920, 2020, 2070, 2170]
k0l = [8.337941973170922e-8, 8.569656306258228e-8, 8.415014789177302e-8, 8.547726090310633e-8, 8.519771266192297e-8, 8.466982230562982e-8, 8.190345663060776e-8, 8.300492022708703e-8, 9.079699499592884e-8, 9.516369911252656e-8, 9.917942007974515e-8, 1.0026850360351777e-7, 9.34103790828767e-8, 8.869917015014817e-8, 8.55603627190414e-8, 8.953569995886163e-8, 8.833435383284136e-8, 9.388314577424422e-8, 9.69608026745865e-8, 1.2493335512244837e-7, 1.51889557485011e-7, 2.1469735154131193e-7, 2.0407346480792858e-7, 1.4086853003543274e-7, 1.2201710734844227e-7, 1.0370172207747732e-7, 9.805335964654704e-8, 9.350906151333254e-8, 8.763332062248001e-8, 8.452789936751738e-8]
k0r = [2.802691032457006e-7, 2.911285205661006e-7, 2.850739454009751e-7, 2.95358399583289e-7, 2.9832330009175595e-7, 3.142113903382259e-7, 3.0222203036005377e-7, 3.082492923492649e-7, 3.2175274730969954e-7, 3.370270511117063e-7, 3.4794870301648575e-7, 3.4905404397156677e-7, 3.2828004955130234e-7, 3.0301200527926653e-7, 2.8871715554175895e-7, 2.909166276256557e-7, 2.7979412315084224e-7, 2.7689949773594095e-7, 2.6911578142323453e-7, 2.6640149429180573e-7, 2.586337060624796e-7, 2.4226100865948255e-7, 2.460919363000388e-7, 2.6080199079041263e-7, 2.6567415106034415e-7, 2.7118256472409915e-7, 2.6997349134660386e-7, 2.733831593865171e-7, 2.640461897924406e-7, 2.61711377169858e-7]

plt_k = plot()
plot!(plt_k, [k0l[1] * 1e6], seriestype=:hline, lw=1.2, color=:red, alpha=0.5, xlim=(250, 2000), ylim=(0.05, 0.4), label=nothing)
plot!(plt_k, [k0r[1] * 1e6], seriestype=:hline, lw=1.2, color=:blue, alpha=0.5, label=nothing)
plot!(plt_k, k_omegacs[2:end], k0r[2:end] * 1e6, color=:blue, lw=3, label=raw"$k\!\!\!_\rightarrow$")
plot!(plt_k, k_omegacs[2:end], k0l[2:end] * 1e6, color=:red, lw=3, label=raw"$k\!\!\!_\leftarrow$")
plot!(plt_k, [530, 800, 1570], seriestype=:vline, color=:grey, ls=:dash, label=nothing)
plot!(plt_k,
    size=(393, TOP_HEIGHT),
    legend=(0.50, 0.50),
    legendfontsize=9,
    tickfontsize=8,
    guidefontsize=10,
    legend_background_color=:transparent,
    fg_legend=:transparent,
    left_margin=1Plots.mm,
    right_margin=1Plots.mm,
    bottom_margin=1Plots.mm,
    xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
    ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$")

# Fig. 7e: branching ratios. The zero-cavity datum is omitted as in the
# original standalone script.
branch_left = k0l ./ (k0l .+ k0r)
branch_right = k0r ./ (k0l .+ k0r)
plt_branch = plot(
    k_omegacs[2:end], branch_left[2:end],
    color=:red, lw=2, label=nothing,
    xlim=(250, 2200), ylim=(0, 1),
    xticks=[500, 1500], xminorticks=5,
    yticks=0:0.25:1, yminorticks=5,
    grid=true,
    size=(117, TOP_HEIGHT),
    tickfontsize=8,
    guidefontsize=10,
    left_margin=1Plots.mm,
    right_margin=0Plots.mm,
    bottom_margin=1Plots.mm,
    xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
    ylabel=raw"$k_{\leftrightharpoons}/(k_{\rightarrow}+k_{\leftarrow})$")
plot!(plt_branch, k_omegacs[2:end], branch_right[2:end], color=:blue, lw=2, label=nothing)
annotate!(plt_branch, 1625, 0.14, text(raw"$k_{\leftarrow}$", color=:red, pointsize=9))
annotate!(plt_branch, 1625, 0.88, text(raw"$k_{\rightarrow}$", color=:blue, pointsize=9))

# Preserve Fig. 7e's narrow aspect while keeping the adjacent rate plot wide.
# The panel order places Fig. 7e on the right, after the existing panels.
top_row = plot(plt_k, plt_branch; layout=grid(1, 2, widths=[0.77, 0.23]))

plt_out = plot(top_row, plt_rates;
    layout=grid(2, 1, heights=[TOP_HEIGHT / FIGURE_SIZE[2], BOTTOM_HEIGHT / FIGURE_SIZE[2]]),
    size=FIGURE_SIZE,
    dpi=300
)

savefig(plt_out, joinpath(@__DIR__, "fig7bcd.pdf"))
