# Vinther2026dynamics

Standalone Julia scripts that reproduce the figures used in the manuscript.
Each script is self-contained and can be run independently:

```bash
julia path/to/script.jl
```

Each script writes its figure(s) as a PDF into the same folder that the script is in.

## Layout

- `SI/` — scripts and data for the figures found in the SI.
- `data/` — data files that is shared by the `main-script/` and `SI/`
- `main-script/` — scripts and data for the figures found in the main manuscript.
- `common.jl` — shared imports, environment activation, and helper functions
  (Gibbs-state bath corrections, eigenstate mixing, rate-bar plotting, etc.)
  used by several of the figure scripts. Not meant to be run directly.

TODO:
- Include in readme: 
  - an instruction to instantiate the environment and to load the files pushed with git lfs
  - the size of the folders.
  - a complete overview of the files and folders and how the scripts and data-files go together.

<!-- - One script per figure (or closely related group of figures), listed below. -->

<!-- This folder is self-contained for data — only the analysis source code and the shared Julia environment are
pulled from sibling folders:

```
JuliaEnvs/MyJuliaEnv/         Julia project used by all scripts
FlatironStudy/src/            analysis code (ITensorHEOM module)
FlatironReproduceFigures/     this folder
  data/                       all data, self-contained
  common.jl
  *.jl                        one script per figure
``` -->

<!-- ## Scripts and required data

Paths below are relative to `FlatironReproduceFigures/data/`.

| Script | Figure(s) | Data required | Status |
|---|---|---|---|
| `kinetic_vs_thermodynamic_control.jl` | `ExampleOfKineticVsThermodynamicControl.pdf` | none (self-contained toy model) | ✅ ready |
| `addition_to_competing_pathways.jl` | `AdditionToCompetingPathways_main_plot.pdf` | none (rate values hard-coded from prior HEOM runs) | ✅ ready |
| `asymmetric_transition_rates.jl` | `asymmetric_transition_rates2.pdf` | `HarryPlotter/PlottingTools/ReactionRates.jld2`, `SystemsData/Data/degenerate_asymmetric_double_well_ratio{2,4,6,8}.0.jld2` | ✅ ready |
| `equilibrium_populations_double_wells.jl` | `EquilibriumPopulations_DoubleWells.pdf` | `SystemsData/Data/degenerate_asymmetric_double_well_ratio{2,4,6,8}.0.jld2`, `BathsData/Data/enhancement_bath_sweep_precise/*.jld2`, `HarryPlotter/PlottingTools/ReactionRates.jld2` | ✅ ready |
| `distribution_of_rates.jl` | `DistributionOfRates.pdf` | `HarryPlotter/EnvironmentDrivenDistributionOfRates.jld2` | ✅ ready |
| `phonon_dressed_triple_well_pes.jl` | `PhononDressedEffectiveTripleWellPES.pdf` | `SystemsData/Data/degenerate_asymmetric_triple_wells/overtonic_morefinetuned.jld2`, `BathsData/Data/enhancement_bath_sweep_precise/etac0.jld2`, `BathsData/Data/distinguished_enhancement_bath/omegac1570.jld2` | ✅ ready |
| `equilibrium_populations_triple_wells.jl` | `EquilibriumPopulations_TripleWells.pdf` | `SystemsData/Data/degenerate_asymmetric_triple_wells/overtonic_morefinetuned.jld2`, `BathsData/Data/enhancement_bath_sweep_precise/*.jld2` | ✅ ready |
| `noise_spectrums_with_oc.jl` | `NoiseSpectrumsWithOC_filled.pdf` | `BathsData/YAML_scripts/enhancement_bath_sweep_precise/etac0.yaml`, `BathsData/YAML_scripts/isotropic_baths_precise/theta_*/omegac220.yaml` | ✅ ready |
| `markovian_heom_approximation.jl` | `ProofOfMarkovianCorrection.pdf` | `BathsData/Data/enhancement_bath_sweep_precise/etac0.jld2` (present), `Experiments/MarkovianityHEOMApproximations2/Data/*_renormedjump.jld2` | ⏳ needs `Experiments/...` |
| `environment_driven_rates.jl` | `EnvironmentDrivenRates.pdf` | `HarryPlotter/EnvironmentDrivenRates_data.jld2` (consolidated raw simulation traces, see below); bath YAML configs already present | ⏳ needs consolidated file |
| `simple_noise_spectrums.jl` | `SimpleNoiseSpectrums.pdf` | `Experiments/EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio0.25` | ⏳ needs `Experiments/...` |
| `symmetric_transition_rates2.jl` | `symmetric_transition_rates2.pdf`, `CompetingPathways_transition_rates.pdf` | `Experiments/EnhancementSweep_etanu0.1_gammanu200_etac0.1_TripleWell_NearDegenerate_4` | ⏳ needs `Experiments/...` |
| `competing_pathways_rate_diff_omegac.jl` | `CompetingPathways_rate_diff_omegac.pdf` | `Experiments/EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio2.0` | ⏳ needs `Experiments/...` |

Each script has a docstring-style header comment listing its exact data
dependencies. If a required file is missing, the script raises a clear error
naming the missing path (see `require_file` in `common.jl`) rather than
failing with a generic I/O error.

Nine of thirteen scripts are fully self-contained and verified to run
end-to-end against the data in this folder. The remaining four need genuine
HEOM/rate-equation simulation outputs (`Experiments/...`) that only exist on
the cluster — copying bath/system config files alone cannot substitute for
those. `environment_driven_rates.jl`'s bath-config dependency has already
been resolved locally (via `BathsData/YAML_scripts/isotropic_baths_precise/`);
it only still needs the raw simulation traces.

### Consolidating cluster-only data into a single file

For `environment_driven_rates.jl`, run the notebook cell appended to the
bottom of `FlatironStudy/HarryPlotter/PlottingTools/Figure Maker.ipynb` (on
the machine where the `Experiments/...` data lives). It bundles the raw HEOM
simulation traces into one portable `EnvironmentDrivenRates_data.jld2` file —
copy just that one file to `data/HarryPlotter/EnvironmentDrivenRates_data.jld2`
here and the script will load it directly. -->

<!-- ## Notes on fidelity to the original notebooks

The original analysis notebooks (`FlatironStudy/HarryPlotter/**/*.ipynb`)
were exploratory and not written to be run top-to-bottom: `savefig` calls
were commented out, and in a few places later cells depended on global
variables set by cells for a *different* figure (rerun interactively with
different parameters). Each script here was reconstructed as a clean,
linear, top-to-bottom pipeline for exactly one figure; where the original
notebook's dependency chain was ambiguous, the choice made is noted in the
script's header comment.

Three additional exploratory variants from `ProbCurrent Figure Maker_latest.ipynb`
(`rate_plot3_ratio2.0_and_ratio0.5.pdf`, `symmetric_rate_diff_omegac2.pdf`) were
not curated here, since they depend on manual, non-linear reruns of the
notebook that could not be reliably reconstructed; `competing_pathways_rate_diff_omegac.jl`
and `symmetric_transition_rates2.jl` cover the same underlying data/analysis. -->
