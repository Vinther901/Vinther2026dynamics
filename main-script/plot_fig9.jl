"""
Reproduces `EnvironmentDrivenRates.pdf`.

Compares HEOM-computed reaction rate constants (as a function of cavity
frequency, for several angles theta between the cavity axis and the
molecular dipole) against a Fermi's-Golden-Rule estimate built from the
noise-spectrum enhancement due to the cavity.

Required data (relative to FLATIRON_ROOT):
  - HarryPlotter/EnvironmentDrivenRates_data.jld2
    (only the raw HEOM simulation traces: pr_data_by_theta, baseline_test_pr/times,
    baseline_enh_pr/times — see the notebook cell below)
  - BathsData/YAML_scripts/enhancement_bath_sweep_precise/etac0.yaml
  - BathsData/YAML_scripts/isotropic_baths_precise/theta_<theta>/omegac220.yaml,
    for theta in [0.05, 0.25, 0.5, 0.75, 0.95]

The consolidated jld2 file is produced by running the notebook cell added at
the bottom of `FlatironStudy/HarryPlotter/PlottingTools/Figure Maker.ipynb`
(bundles the raw `Experiments/...` simulation outputs this figure depends on
into one portable file) — copy it to
`data/HarryPlotter/EnvironmentDrivenRates_data.jld2` here. The bath YAML
configs are ordinary source-repo files, not simulation outputs.
"""

include(joinpath(@__DIR__, "..", "common.jl"))
using QuadGK

# const BATH_YAML_DIR = joinpath(@__DIR__, "..", "data", "BathsData", "YAML_scripts")
# const DATA_FILE = require_file(joinpath(FLATIRON_ROOT, "HarryPlotter", "EnvironmentDrivenRates_data.jld2"))
_data = load(require_file(joinpath(@__DIR__, "EnvironmentDrivenRates_data.jld2")))

thetas = _data["thetas"]
pr_data_by_theta = _data["pr_data_by_theta"]
bath0_params = load_file(require_file(joinpath(@__DIR__, "etac0.yaml")))
bath_params_by_theta = Dict(
    theta => load_file(require_file(joinpath(@__DIR__, "theta$(theta).yaml")))
    for theta in thetas
)

# manually curated indices of "outlier" omegac points to drop per theta.
# One can check that the corresponding BCF wasn't decomposed properly, i.e. the fit was not accurate enough to warrant inclusion in the final plot.
rind_by_theta = Dict(
    0.05 => [1, 9],
    0.25 => [1, 24],
    0.5 => [1, 9, 12, 23],
    0.75 => [1, 14, 16, 26, 29],
    0.95 => [1, 13, 14, 18, 19, 26],
)

function get_kfs(pr_data)
    omegacs = sort(collect(keys(pr_data)))

    t0s = LinRange(1000, 7000, 7)
    kfs_t0sweep = zeros(length(omegacs), length(t0s))
    kbs_t0sweep = zeros(length(omegacs), length(t0s))

    for (i, omegac) in enumerate(omegacs)
        pr = pr_data[omegac].pr
        times = pr_data[omegac].times
        for (j, t0) in enumerate(t0s)
            kf, kb = solve_rate_constant(pr, times; t0)
            kfs_t0sweep[i, j] = kf
            kbs_t0sweep[i, j] = kb
        end
    end
    return kfs_t0sweep[:, end], omegacs
end

res = Dict()
for (theta, rind) in zip(thetas, [rind_by_theta[t] for t in thetas])
    kHEOMPR, omegacs = get_kfs(pr_data_by_theta[theta])
    inds = setdiff(eachindex(omegacs), rind)
    res[theta] = (; kHEOMPR=kHEOMPR[inds], omegacs=omegacs[inds])
end

# baseline "outside cavity" rate, k0HEOM
pr = _data["baseline_test_pr"]
times = _data["baseline_test_times"]
dprdt = (pr[2:end] - pr[1:end-1]) ./ (times[2:end] - times[1:end-1])
pr_mean = (pr[2:end] + pr[1:end-1]) ./ 2
k = -dprdt ./ (2 .* pr_mean .- 1)
k0HEOM = (minimum(k[end-4000:end]) + maximum(k[end-4000:end])) / 2

function get_koft(omegac, pr_data)
    pr = pr_data[omegac].pr
    times = pr_data[omegac].times
    dprdt = (pr[2:end] - pr[1:end-1]) ./ (times[2:end] - times[1:end-1])
    pr_mean = (pr[2:end] + pr[1:end-1]) ./ 2
    k = -dprdt ./ (2 .* pr_mean .- 1)
    return k, (times[2:end] + times[1:end-1]) ./ 2
end

for theta in thetas
    pr_data = pr_data_by_theta[theta]
    ks = [begin
        k, _ = get_koft(omegac, pr_data)
        (minimum(k[end-4000:end]) + maximum(k[end-4000:end])) / 2
    end for omegac in res[theta].omegacs]
    res[theta] = (; res[theta]..., kHEOM=ks)
end

# baseline "no environment enhancement" rate, kD
pr = _data["baseline_enh_pr"]
times = _data["baseline_enh_times"]
dprdt = (pr[2:end] - pr[1:end-1]) ./ (times[2:end] - times[1:end-1])
pr_mean = (pr[2:end] + pr[1:end-1]) ./ 2
k = -dprdt ./ (2 .* pr_mean .- 1)
kD = (minimum(k[end-4000:end]) + maximum(k[end-4000:end])) / 2

# FGR estimate for the cavity-induced rate enhancement
S0 = get_bath_func(bath0_params)

Gamma = convert_unit(30, :invcm, :au)
omega0 = convert_unit(1205, :invcm, :au)
c = 0.38

T = convert_unit(bath0_params["T"]["value"], :K, :au)
pdf(omega) = Gamma / ((omega - omega0)^2 + Gamma^2) / pi

omegacsFGR = 150:10:2250
for theta in thetas
    bath_params = bath_params_by_theta[theta]

    ks = map(omegacsFGR) do omegac
        bath_params["ω_c"]["value"] = omegac
        S = get_bath_func(bath_params)
        DeltaS(omega) = S(omega) - S0(omega)
        f(omega) = pdf(omega) * DeltaS(omega) * exp(-omega / T)
        k, _ = quadgk(f, 0, 0.1)
        k
    end
    res[theta] = (; res[theta]..., kFGR=ks)
end

plt = plot([k0HEOM * 1e6], seriestype=:hline,
    label=nothing,
    legendtitle=raw"$\theta$ [$\pi$]",
    legend_background_color=:transparent,
    fg_legend=:transparent,
    ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
    xlabel=raw"$\omega_c\quad [\mathsf{cm}^{-1}]$",
    color=:black,
    lw=1,
    xminorticks=5,
)
colors = palette(:rainbow1, length(thetas))
for (theta, color) in zip(thetas, colors)
    plot!(plt, omegacsFGR, (convert_unit(0.214^2 .* res[theta].kFGR, :au, :invfs) .* c .+ kD) * 1e6, lw=2; color, label=theta)
end
for (theta, color) in zip(thetas, colors)
    plot!(plt, res[theta].omegacs, res[theta].kHEOM * 1e6; color, lw=1, ls=:dash, marker=:o, label=nothing)
end
scatter!(plt, [0], [0], color=:grey, label="HEOM")
plot!(plt, [0], [0], color=:grey, label="FGR", lw=2, ylim=(3.9, 6.5), xlim=(150, 2250), size=(450, 350), dpi=300)

# savefig(plt, joinpath(@__DIR__, "EnvironmentDrivenRates.pdf"))
savefig(plt, joinpath(@__DIR__, "fig9.pdf"))
