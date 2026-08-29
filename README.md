# Figure-reproduction data for Vinther et al. (2026)

This repository contains the Julia scripts and simulation data needed to
reproduce the figures accompanying the manuscript **[paper title and citation
to be added when available]**.

Each plotting script is independent: it activates the repository's Julia
environment, reads only files included here, and writes PDF output next to the
script. The scripts do not rerun the underlying HEOM simulations.

## Requirements

- [Git](https://git-scm.com/)
- [Git LFS](https://git-lfs.com/)
- [Julia](https://julialang.org/downloads/) 1.12 (the version recorded in
  `Manifest.toml`; other recent Julia 1.x releases may also work)

The complete checkout occupies approximately **2.1 GB**: `main-script/` is
about **1.1 GB**, `SI/` about **1.0 GB**, and `data/` about **7.4 MB**. Git's
local LFS objects can require roughly another 2.1 GB. Sizes can differ slightly
between filesystems.

## Setup

Clone the repository and download the Git LFS objects before running a script:

```bash
git clone <repository-url> Vinther2026dynamics
cd Vinther2026dynamics
git lfs pull
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If the repository was downloaded as a source archive rather than cloned with
Git, the large `.jld2` files may be missing. A Git clone followed by
`git lfs pull` is the supported installation method. You can confirm that the
data are present with `git lfs ls-files`.

## Reproducing a figure

Run a script from any working directory. For example:

```bash
julia main-script/plot_fig5.jl
julia SI/plot_fig14.jl
```

The first run may take longer while Julia precompiles packages. All paths are resolved relative to the
script, so no source-code edits or environment variables are required.

## Repository layout

```text
Vinther2026dynamics/
├── Project.toml              Julia dependencies
├── Manifest.toml             pinned dependency versions
├── common.jl                 shared setup, physics routines, and plot helpers
├── main-script/              main-text scripts, local inputs, and PDF outputs
├── SI/                       supplementary scripts, local inputs, and outputs
└── data/
    ├── BathsData/            shared bath decompositions
    └── SystemsData/          shared potentials and eigenstates
```

The `.jld2` files are precomputed simulation results stored with Git LFS.
Small YAML files contain bath parameters. `common.jl` activates the project and
provides the unit conversions, bath models, equilibrium corrections, and
rate-plot helpers used across figures; it is not intended to be run directly.

## Figures and data dependencies

Paths in the table are relative to the repository root. Wildcards denote a
collection of cavity-frequency files.

| Manuscript figure | Script | Output | Input data |
|---|---|---|---|
| Figure 3 | `main-script/plot_fig3.jl` | `main-script/fig3.pdf` | `main-script/etac0.yaml` |
| Figure 4(b,c) | `main-script/plot_fig4bc.jl` | `main-script/fig4bc.pdf` | `main-script/EnhancementSweep_…ratio1.0…jld2`; total rates embedded in script |
| Figure 5 | `main-script/plot_fig5.jl` | `main-script/fig5.pdf` | `main-script/ReactionRates.jld2`; `data/SystemsData/degenerate_asymmetric_double_well_ratio{2,4,6,8}.0.jld2` |
| Figure 6 | `main-script/plot_fig6.jl` | `main-script/fig6.pdf` | `main-script/EnhancementSweep_…ratio{2.0,0.5}…jld2` |
| Figure 7(a) | `main-script/plot_fig7a.jl` | `main-script/fig7a.pdf` | `data/SystemsData/overtonic_morefinetuned.jld2` |
| Figure 7(b–d) | `main-script/plot_fig7bcd.jl` | `main-script/fig7bcd.pdf` | `main-script/EnhancementSweep_…TripleWell…jld2`; `data/SystemsData/overtonic_morefinetuned.jld2` |
| Figure 8 | `main-script/plot_fig8.jl` | `main-script/fig8.pdf` | `main-script/etac0.yaml`; `main-script/theta*.yaml` |
| Figure 9 | `main-script/plot_fig9.jl` | `main-script/fig9.pdf` | `main-script/EnvironmentDrivenRates_data.jld2`; `main-script/etac0.yaml`; `main-script/theta*.yaml` |
| Figure 11 | `SI/plot_fig11.jl` | `SI/fig11.pdf` | `SI/MarkovianityHEOMApproximations.jld2`; `data/BathsData/enhancement_bath_sweep_precise/etac0.jld2` |
| Figure 12 | `SI/plot_fig12.jl` | `SI/fig12.pdf` | `SI/EnvironmentDrivenDistributionOfRates.jld2` |
| Figure 13 (left) | `SI/plot_fig13_left.jl` | `SI/fig13_left.pdf` | `data/SystemsData/overtonic_morefinetuned.jld2`; bath files `etac0.jld2` and `distinguished_enhancement_bath/omegac1570.jld2` |
| Figure 13 (right) | `SI/plot_fig13_right.jl` | `SI/fig13_right.pdf` | `data/SystemsData/overtonic_morefinetuned.jld2`; `data/BathsData/enhancement_bath_sweep_precise/*.jld2` |
| Figure 14 | `SI/plot_fig14.jl` | `SI/fig14.pdf` | `main-script/ReactionRates.jld2`; double-well system files; `data/BathsData/enhancement_bath_sweep_precise/*.jld2` |
| Figure 15 | `SI/plot_fig15.jl` | `SI/fig15.pdf` | none (model parameters are embedded in the script) |
| Figure 16 | `SI/plot_fig16.jl` | `SI/fig16.pdf` | none (rate data are embedded in the script) |
| Figure 17(a–d) | `SI/plot_fig17.jl` | `SI/fig17a.pdf`–`SI/fig17d.pdf` | double-well system files; SI `EnhancementSweep_…ratio{0.125,0.167,0.25,4.0,6.0,8.0}…jld2`; main-text ratio `0.5` and `2.0` sweep files |

The long `EnhancementSweep_…` names retain the simulation parameters used to
produce each dataset.

## Notes on reproducibility

- The numerical simulation outputs are supplied; these scripts reproduce the
  analysis and figures, not the original high-cost simulations.
- Plot appearance can vary slightly with Julia, package, font, or PDF backend
  versions. The committed `Manifest.toml` is provided to minimize such drift.
- A missing input raises an error containing its expected path. In most cases,
  this means `git lfs pull` has not completed.

## License

See [LICENSE](LICENSE).
