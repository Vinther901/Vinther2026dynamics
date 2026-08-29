"""
common.jl — shared setup and helper functions for the figure-reproduction scripts
in this folder. Each `*.jl` script `include`s this file first.

The repository is self-contained: shared inputs live under `data/`, while
figure-specific inputs live next to their plotting scripts. See `README.md`
for the complete script-to-data map.
"""

const FLATIRON_ROOT = normpath(joinpath(@__DIR__, "data"))

using Pkg
Pkg.activate(@__DIR__)

using Plots
using JLD2: load
using YAML: load_file
using LinearAlgebra

"""
    require_file(path)

Throw a helpful error if a required data file is missing, instead of a raw
`SystemError` from `JLD2`/`YAML`.
"""
function require_file(path)
    if !isfile(path)
        error("Required data file not found: $path\n" *
              "Copy it into place under FLATIRON_ROOT ($FLATIRON_ROOT) before running this script.")
    end
    return path
end

"""
    get_states_and_energies(pairs, states, energies)

Diagonalize the dipole operator restricted to (near-)degenerate `pairs` of
eigenstate indices, to obtain "localized" states/energies. `pairs` is a
vector of index vectors, e.g. `[[1,2],[3,4]]`.
"""
function get_states_and_energies(pairs, states, energies, x)
    if !isempty(pairs)
        rot_mat = Matrix{Float64}(I, size(states, 2), size(states, 2))
        for pair in pairs
            xmat = states[:, pair]' * diagm(x) * states[:, pair]
            rot_mat[pair, pair] = eigen(xmat).vectors
        end
        energies = diag(rot_mat' * diagm(energies) * rot_mat)
        states = states * rot_mat
    end
    perm = sortperm(energies)
    return states[:, perm], energies[perm]
end

"""
    get_gibbs_correction(z, d, X, dE; N)

Second-order (in system-bath coupling) correction to the Gibbs state,
built from the bath decomposition `(z, d)` (poles/residues of the bath
correlation function) and the system dipole matrix `X` / energy gaps `dE`.
"""
function get_gibbs_correction(z, d, X, dE, gibbs, beta; N=size(X, 1))
    D = nu -> pi * nu * imag(sum(d ./ (z .* (nu .+ im .* z))))
    dD = nu -> pi * imag(sum(d ./ (z .* (nu .+ im .* z)))) -
               pi * nu * imag(sum(d ./ (z .* (nu .+ im .* z) .^ 2)))

    tau = zeros(N, N)
    for i in 1:N
        for j in 1:N
            Xm = zeros(N, N)
            Xm[i, j] = X[i, j]
            wm = dE[i, j]

            tau += beta * gibbs * (Xm * Xm' - I * tr(gibbs * Xm * Xm')) * D(wm)
            tau += (Xm' * gibbs * Xm - gibbs * Xm * Xm') * dD(wm)

            for ip in 1:N
                for jp in 1:N
                    if i == ip && j == jp
                        continue
                    elseif i == j && ip == jp
                        continue
                    end
                    Xn = zeros(N, N)
                    Xn[ip, jp] = X[ip, jp]
                    wn = dE[ip, jp]

                    tmp = Xn * Xm' * gibbs - Xm' * gibbs * Xn
                    tau += (tmp + tmp') * D(wm) / (wn - wm)
                end
            end
        end
    end
    return tau
end

"""
    stacked_fillbar!(p, rates; order, colors, labels, xcenter, width)

Draw one stacked bar (coherent/phonon/cavity rate contributions) at `xcenter`
on plot `p`. `rates` is a named tuple with fields `r_coh`, `r_nu`, `r_c`.
"""
function stacked_fillbar!(p, rates;
    order=[:r_coh, :r_nu, :r_c],
    colors=Dict(:r_coh => :lightsalmon, :r_nu => :purple, :r_c => :mediumseagreen),
    labels=Dict(:r_coh => nothing, :r_nu => nothing, :r_c => nothing),
    xcenter=1.0,
    width=0.8,
)
    vals = Dict(:r_coh => rates.r_coh, :r_nu => rates.r_nu, :r_c => rates.r_c)
    left, right = xcenter - width / 2, xcenter + width / 2
    pos_base, neg_base = 0.0, 0.0
    for key in order
        v = float(vals[key])
        y0, y1 = v >= 0 ? (pos_base, pos_base + v) : (neg_base, neg_base + v)
        v >= 0 ? (pos_base = y1) : (neg_base = y1)
        plot!(p, [left, right], [y1, y1];
            fillrange=[y0, y0], fillcolor=colors[key], fillalpha=1.0,
            linecolor=:black, linealpha=0.0, linewidth=1, label=labels[key])
    end
    return p
end

"""
    get_rates_ij(data, i, j; t0=6000)

Extract the coherent/phonon/cavity contributions (in units of ×1e-6 fs⁻¹)
to the transition rate between eigenstates `i` and `j` at time `t0`, from a
loaded `rate_data[omegac]` entry.
"""
function get_rates_ij(data, i, j; t0=6000)
    ind = argmin(abs.(data.times .- t0))
    r_coh = data.rate_coh[i, j, ind] * 1e6
    r_nu = (data.rate_Mnu[i, j, ind] + data.rate_nu[i, j, ind]) * 1e6
    r_c = (data.rate_Mc[i, j, ind] + data.rate_c[i, j, ind]) * 1e6
    r = r_coh + r_nu + r_c
    return (; r, r_coh, r_nu, r_c)
end

"""
    solve_rate_constant(pr, times; t0=0, t1=nothing)

Fit forward/backward rate constants `(k, k_b)` to a reaction-probability
time series `pr(times)` using the integrated master-equation method.
"""
function solve_rate_constant(pr, times; t0=0, t1=nothing)
    ind = argmin(abs.(times .- t0))
    if isnothing(t1)
        dtimes = times[ind:end] .- times[ind-1:end-1]
        pr, times = pr[ind:end], times[ind:end]
    else
        ind1 = argmin(abs.(times .- t1))
        dtimes = times[ind:ind1] .- times[ind-1:ind1-1]
        pr, times = pr[ind:ind1], times[ind:ind1]
    end
    PR = vcat(0, cumsum(pr[1:end-1])) .* dtimes
    k, k_b = hcat(-PR, -PR + times .- times[1]) \ (pr .- pr[1])
    return k, k_b
end


# =============================================================================
# Unit conversion, bath selection, and reaction-rate calculation
#
# Required bath constructors for get_bath_func:
#   PhononBathPlusMarkovianCavity(params)
#   PhononBathPlusSpectatorModesWithCavity(params)
#   S_isotropic(params)
# =============================================================================


# -----------------------------------------------------------------------------
# Unit-conversion constants
# -----------------------------------------------------------------------------

const FS_TO_AU = 41.341374575751
const AU_TO_FS = 1 / FS_TO_AU

const AU_TO_INVFS = 1 / AU_TO_FS
const INVFS_TO_AU = 1 / AU_TO_INVFS

const AU_TO_INVCM = 219474.6305
const INVCM_TO_AU = 1 / AU_TO_INVCM

const AU_TO_K = 315790.83525179856
const K_TO_AU = 1 / AU_TO_K


# -----------------------------------------------------------------------------
# Unit conversion
# -----------------------------------------------------------------------------

"""
    convert_unit(value, from_unit, to_unit)

Convert a scalar or array between supported units.

Supported units:

- `"au"` or `:au`: atomic units
- `"fs"` or `:fs`: femtoseconds
- `"invfs"` or `:invfs`: inverse femtoseconds
- `"invcm"` or `:invcm`: inverse centimetres
- `"K"` or `:K`: kelvin

Both strings and symbols are accepted for the unit arguments.

# Examples

    omega_au = convert_unit(1170.0, "invcm", "au")
    time_au = convert_unit(100.0, :fs, :au)
    times_au = convert_unit([0.0, 10.0, 20.0], "fs", "au")
"""
function convert_unit(value, from_unit, to_unit)
    from = Symbol(from_unit)
    to = Symbol(to_unit)

    supported_units = (:au, :fs, :invfs, :invcm, :K)

    from in supported_units ||
        throw(ArgumentError(
            "Unsupported input unit: $from_unit. " *
            "Supported units are: $(join(string.(supported_units), ", "))."
        ))

    to in supported_units ||
        throw(ArgumentError(
            "Unsupported output unit: $to_unit. " *
            "Supported units are: $(join(string.(supported_units), ", "))."
        ))

    from == to && return value

    # Convert the input to atomic units.
    value_au =
        from == :invcm ? value .* INVCM_TO_AU :
        from == :fs    ? value .* FS_TO_AU :
        from == :invfs ? value .* INVFS_TO_AU :
        from == :K     ? value .* K_TO_AU :
        value

    # Convert from atomic units to the requested output unit.
    return (
        to == :invcm ? value_au .* AU_TO_INVCM :
        to == :fs    ? value_au .* AU_TO_FS :
        to == :invfs ? value_au .* AU_TO_INVFS :
        to == :K     ? value_au .* AU_TO_K :
        value_au
    )
end


# -----------------------------------------------------------------------------
# Dimensionful parameter extraction
# -----------------------------------------------------------------------------

"""
    get_dimful_param(param)

Extract `"value"` and `"unit"` from a dictionary and convert the value to
atomic units.

Both string and symbol dictionary keys are supported.

# Examples

    omega = get_dimful_param(
        Dict("value" => 1170.0, "unit" => "invcm")
    )

    times = get_dimful_param(
        Dict(:value => [0.0, 10.0, 20.0], :unit => :fs)
    )
"""
function get_dimful_param(param::AbstractDict)
    value_key =
        haskey(param, "value") ? "value" :
        haskey(param, :value)  ? :value :
        throw(ArgumentError(
            "The parameter dictionary must contain a \"value\" or :value key."
        ))

    unit_key =
        haskey(param, "unit") ? "unit" :
        haskey(param, :unit)  ? :unit :
        throw(ArgumentError(
            "The parameter dictionary must contain a \"unit\" or :unit key."
        ))

    value = param[value_key]
    unit = param[unit_key]

    return convert_unit(value, unit, :au)
end


# -----------------------------------------------------------------------------
# Bath-function selection
# -----------------------------------------------------------------------------

"""
    get_bath_func(params)

Construct and return the spectral-noise function specified by the `"bath"`
entry in `params`.

Supported bath names:

- `"PhononBathPlusMarkovianCavity"`
- `"PhononBathPlusSpectatorModesWithCavity"`
- `"S_isotropic"`

The corresponding bath constructors must already be defined or imported into
the current scope. Each constructor receives `params` and returns a callable
spectral-noise function `S(omega)`.

# Example

    params = Dict(
        "bath" => "PhononBathPlusMarkovianCavity",
        # Remaining bath parameters...
    )

    S = get_bath_func(params)
    noise = S(omega)
"""
function get_bath_func(params::AbstractDict)
    bath_key =
        haskey(params, "bath") ? "bath" :
        haskey(params, :bath)  ? :bath :
        throw(ArgumentError(
            "The parameter dictionary must contain a \"bath\" or :bath key."
        ))

    bath_name = String(params[bath_key])

    if bath_name == "PhononBathPlusMarkovianCavity"
        return PhononBathPlusMarkovianCavity(params)

    elseif bath_name == "PhononBathPlusSpectatorModesWithCavity"
        return PhononBathPlusSpectatorModesWithCavity(params)

    elseif bath_name == "S_isotropic"
        return S_isotropic(params)

    else
        throw(ArgumentError(
            "Unknown bath \"$bath_name\". Supported baths are: " *
            "\"PhononBathPlusMarkovianCavity\", " *
            "\"PhononBathPlusSpectatorModesWithCavity\", and " *
            "\"S_isotropic\"."
        ))
    end
end


# -----------------------------------------------------------------------------
# Time-dependent reaction-rate calculation
# -----------------------------------------------------------------------------

"""
    get_k(
        pr;
        f=nothing,
        times=nothing,
        chi=1.0,
        finite_difference=false,
    )

Calculate the forward reaction rate from the rate equation

    dP_R/dt = k(t) * [chi - (1 + chi) * P_R(t)]

where `chi` is the ratio of the backward rate to the forward rate.

# Arguments

- `pr`: Right-well population or projection.
- `f`: Optional flux correlation function, interpreted as `d(pr)/dt`.
- `times`: Times corresponding to the entries in `pr`.
- `chi`: Ratio between backward and forward rates. Default: `1.0`.
- `finite_difference`: If `true`, calculate the flux using finite
  differences. If `false`, use the integrated rate equation.

# Returns

- With `f`: returns the rate vector `k`.
- With `times` and `finite_difference=true`: returns
  `(k, midpoint_times)`.
- With `times` and `finite_difference=false`: returns the rate vector `k`.

The first element returned by the integrated method is `NaN`, because both
the numerator and denominator are zero at the initial time.

# Examples

    k = get_k(pr; f=flux, chi=1.0)

    k, midpoint_times = get_k(
        pr;
        times=times,
        chi=1.0,
        finite_difference=true,
    )

    k = get_k(pr; times=times, chi=1.0)
"""
function get_k(
    pr::AbstractVector;
    f=nothing,
    times=nothing,
    chi::Real=1.0,
    finite_difference::Bool=false,
)
    isempty(pr) &&
        throw(ArgumentError("`pr` must contain at least one value."))

    isfinite(chi) ||
        throw(ArgumentError("`chi` must be finite."))

    chi >= 0 ||
        throw(ArgumentError("`chi` must be non-negative."))

    # -------------------------------------------------------------------------
    # Mode 1: The flux f = d(pr)/dt is supplied directly.
    # -------------------------------------------------------------------------

    if !isnothing(f)
        length(f) == length(pr) ||
            throw(DimensionMismatch(
                "`f` and `pr` must have the same length. " *
                "Received lengths $(length(f)) and $(length(pr))."
            ))

        denominator = chi .- (1 + chi) .* pr

        return f ./ denominator
    end

    # If no flux was supplied, times are required.
    isnothing(times) &&
        throw(ArgumentError(
            "Either `f` or `times` must be supplied."
        ))

    length(times) == length(pr) ||
        throw(DimensionMismatch(
            "`times` and `pr` must have the same length. " *
            "Received lengths $(length(times)) and $(length(pr))."
        ))

    length(times) >= 2 ||
        throw(ArgumentError(
            "At least two time points are required."
        ))

    delta_times = diff(times)

    all(delta_times .> 0) ||
        throw(ArgumentError(
            "`times` must be strictly increasing."
        ))

    # -------------------------------------------------------------------------
    # Mode 2: Estimate the instantaneous rate using finite differences.
    # -------------------------------------------------------------------------

    if finite_difference
        flux = diff(pr) ./ delta_times

        midpoint_times =
            (times[2:end] .+ times[1:end-1]) ./ 2

        midpoint_pr =
            (pr[2:end] .+ pr[1:end-1]) ./ 2

        denominator =
            chi .- (1 + chi) .* midpoint_pr

        rate = flux ./ denominator

        return rate, midpoint_times
    end

    # -------------------------------------------------------------------------
    # Mode 3: Estimate the rate using the integrated rate equation:
    #
    # pr(t) - pr(t0) =
    #     k * [
    #         chi * (t - t0)
    #         - (1 + chi) * integral(pr(s), s=t0..t)
    #     ]
    #
    # The population integral is evaluated with the trapezoidal rule.
    # -------------------------------------------------------------------------

    numeric_type = promote_type(
        eltype(pr),
        eltype(times),
        typeof(float(chi)),
    )

    integrated_pr = zeros(numeric_type, length(pr))

    @inbounds for i in 2:length(pr)
        interval_width = times[i] - times[i - 1]
        interval_mean_pr = (pr[i - 1] + pr[i]) / 2

        integrated_pr[i] =
            integrated_pr[i - 1] +
            interval_mean_pr * interval_width
    end

    elapsed_time = times .- first(times)

    numerator = pr .- first(pr)

    denominator =
        chi .* elapsed_time .-
        (1 + chi) .* integrated_pr

    return numerator ./ denominator
end


"""
    Parameters

A mutable struct to hold all system parameters.
This replaces the Python dictionary-based approach.
"""
mutable struct Parameters
    ω_c::Dict
    T::Dict
    τ_c::Dict
    η_c::Float64
    c_0::Dict
    α::Dict
    tilde_η_S::Float64
    λ_ν::Dict
    γ_ν::Dict
    N::Int
    tilde_η_Q::Float64
    λ_Q::Dict
    γ_Q::Dict
    c_Q::Dict
    ω_Q::Dict
    # ω_min::Dict
    # ω_max::Dict
end

"""
    initialize_parameters()::Parameters

Initialize parameters with default values matching the Python notebook.
"""
function initialize_parameters(ω_c, η_c)::Parameters
    params = Parameters(
        # Dict("value" => ω_c, "unit" => "invcm"),    # ω_c
        ω_c, # ω_c
        Dict("value" => 300.0, "unit" => "K"),         # T
        Dict("value" => 100.0, "unit" => "fs"),        # τ_c
        η_c,                                             # η_c <-- should actually be unitful, e.g. a.u.u
        Dict("value" => 0.0, "unit" => "au"),          # c_0 (computed)
        Dict("value" => 1/100.0, "unit" => "invfs"),          # α (computed)
        1.0,                                             # tilde_η_S
        Dict("value" => 83.6, "unit" => "invcm"),      # λ_ν
        Dict("value" => 200.0, "unit" => "invcm"),     # γ_ν
        1,                                               # N
        0.1,                                             # tilde_η_Q
        Dict("value" => 0.147, "unit" => "invcm"),     # λ_Q
        Dict("value" => 6000.0, "unit" => "invcm"),    # γ_Q
        Dict("value" => 4.9, "unit" => "invcm"),       # c_Q
        Dict("value" => 1170.0, "unit" => "invcm"),    # ω_Q
        # Dict("value" => 750.0, "unit" => "invcm"),     # ω_min
        # Dict("value" => 1750.0, "unit" => "invcm"),    # ω_max
    )
    
    # Compute derived parameters
    ω_c_au = get_dimful_param(params.ω_c)
    params.c_0["value"] = sqrt(2 * params.η_c^2 * ω_c_au^3)
    
    # T_au = get_dimful_param(params.T)
    # τ_c_au = get_dimful_param(params.τ_c)
    # params.α["value"] = (1 - exp(-ω_c_au / T_au)) / τ_c_au
    
    return params
end

# ============================================================================
# Bath Response Functions
# ============================================================================

"""
    J_μ(omega::Union{Float64, AbstractArray}, params::Parameters)::Union{Float64, Array}

Compute the imaginary part of the ν-bath self-energy.
"""
function J_μ(omega, params)
    λ_Q = get_dimful_param(params.λ_Q)
    γ_Q = get_dimful_param(params.γ_Q)
    return 2 .* λ_Q .* γ_Q .* omega ./ (omega.^2 .+ γ_Q^2)
end

"""
    R_tilde_μ(omega::Union{Float64, AbstractArray}, params::Parameters)::Union{Float64, Array}

Compute the real part of the ν-bath self-energy.
"""
function R_tilde_μ(omega, params)
    λ_Q = get_dimful_param(params.λ_Q)
    γ_Q = get_dimful_param(params.γ_Q)
    return 2 .* λ_Q .* omega.^2 ./ (omega.^2 .+ γ_Q^2)
end

"""
    J_c(omega::Union{Float64, AbstractArray}, params::Parameters)::Union{Float64, Array}

Compute the imaginary part of the cavity self-energy.
"""
function J_c(omega, params)
    α = get_dimful_param(params.α)
    return α .* omega
end

"""
    R_tilde_c(omega::Union{Float64, AbstractArray}, params::Parameters)::Union{Float64, Array}

Compute the real part of the cavity self-energy.
Returns zero (no contribution).
"""
function R_tilde_c(omega, params)
    return 0.0
end

# ============================================================================
# Main Spectral Density Function
# ============================================================================

"""
    J_isotropic(omega::Union{Float64, AbstractArray}, theta::Union{Float64, AbstractArray}, 
                 params::Parameters)::Union{Float64, Array}

Compute J(ω) using the analytical formula with proper angle averaging over α.

Based on the analytical equations:
    <1/(a-b*cos²(α))>_α = 1/√(a(a-b))
    <η_S²/(a-b*cos²(α))>_α = tilde_η_S² * [sin²(θ)/√(a(a-b)) + cos(2θ)/b * (√(a/(a-b)) - 1)]
    <η_S*η_Q/(a-b*cos²(α))>_α = tilde_η_S*tilde_η_Q * cos(θ)/b * (√(a/(a-b)) - 1)
    <η_S²*η_Q²/(a-b*cos²(α))>_α = tilde_η_S²*tilde_η_Q² * [sin²(θ)/b * (√(a/(a-b)) - 1) + 
                                     cos(2θ) * (a/b² * (√(a/(a-b)) - 1) - 1/(2b))]

where:
    a = A_c*A_μ
    b = tilde_η_Q²*(2*A_c*ω_c*η_c² + c_c²)

Arguments:
- omega: Frequency/frequencies (scalar or array)
- theta: Angle parameter in radians (scalar)
- params: Parameters structure

Returns:
- Imaginary part of the spectral density
"""
function J_isotropic(omega, theta,params)
    
    # Extract parameters
    ω_c = get_dimful_param(params.ω_c)
    η_c = params.η_c
    c_Q = get_dimful_param(params.c_Q)
    ω_Q = get_dimful_param(params.ω_Q)
    
    # Tilde parameters (angle-independent coupling constants)
    tilde_η_S = params.tilde_η_S
    tilde_η_Q = params.tilde_η_Q
    
    # Compute c_c (cavity coupling strength)
    c_c = sqrt(2 * ω_c^3 * η_c^2)
    
    # Compute A_μ and A_c
    A_μ = omega.^2 .- ω_Q^2 .- R_tilde_μ(omega, params) .- 1im .* J_μ(omega, params)
    A_c = omega.^2 .- ω_c^2 .- R_tilde_c(omega, params) .- 1im .* J_c(omega, params)
    
    # Define a and b for the bracket formulas
    a = A_c .* A_μ
    b = tilde_η_Q^2 .* (2 .* A_c .* ω_c .* η_c^2 .+ c_c^2)
    
    # Pre-compute common factors
    # inv_sqrt_ab = 1.0 / (a*sqrt(1 - b/a))
    inv_sqrt_ab = 1.0 ./ (a .* sqrt.(1 .- b ./ a))
    sqrt_ratio = a .* inv_sqrt_ab
    common_factor = sqrt_ratio .- 1
    
    # Term 1: A_c * c_Q² * <1/(a-b*cos²(α))>_α
    term1 = A_c .* c_Q^2 .* inv_sqrt_ab
    
    # Term 2: η_S² * A_μ * c_c² * <η_S²/(a-b*cos²(α))>_α
    term2 = tilde_η_S^2 .* A_μ .* c_c^2 .* 
            (sin.(theta).^2 .* inv_sqrt_ab .+ cos.(2 .* theta) ./ b .* common_factor)
    
    # Term 3: η_Q * η_S * (4*A_c*c_Q*ω_c*η_c² + 2*c_c²*c_Q) * <η_S*η_Q/(a-b*cos²(α))>_α
    term3_coeff = (4 .* A_c .* c_Q .* ω_c .* η_c^2 .+ 2 .* c_c^2 .* c_Q)
    term3 = tilde_η_Q .* tilde_η_S .* term3_coeff .* cos.(theta) ./ b .* common_factor
    
    # Term 4: η_Q² * η_S² * (4*A_c*ω_c²*η_c² + 2*c_c²*η_c²*ω_c) * <η_S²*η_Q²/(a-b*cos²(α))>_α
    term4_coeff = (4 .* A_c .* ω_c^2 .* η_c^4 .+ 2 .* c_c^2 .* η_c^2 .* ω_c)
    term4 = term4_coeff .* tilde_η_S^2 .* tilde_η_Q^2 .* 
            (sin.(theta).^2 ./ b .* common_factor .+ 
             cos.(2 .* theta) .* (a ./ b.^2 .* common_factor .- 1 ./ (2 .* b)))
    
    result = term1 .+ term2 .+ term3 .+ term4
    
    return imag.(result)
end


function PhononBathPlusMarkovianCavity(params::Dict)
    # Following e.g. https://pubs.acs.org/doi/10.1021/acs.jpclett.3c02985
    ω_c = get_dimful_param(params["ω_c"])
    T = get_dimful_param(params["T"])
    τ_c = get_dimful_param(params["τ_c"])
    η_c = params["η_c"]
    γ_ν = get_dimful_param(params["γ_ν"])
    λ_ν = get(params,"λ_ν",nothing)

    if isnothing(λ_ν)
        M = get_dimful_param(params["M"])
        ω_b = get_dimful_param(params["ω_b"])
        η_ν = params["η_ν"]
        λ_ν = 0.5 * η_ν * M * γ_ν * ω_b
    else
        λ_ν = get_dimful_param(λ_ν)
    end
        
    α = (1 - exp(-ω_c / T)) / τ_c

    # Phonon bath spectral density
    J_ν = ω -> begin
        numerator = 2 * λ_ν * γ_ν * ω
        denominator = ω^2 + γ_ν^2
        numerator / denominator
    end

    # Effective cavity bath spectral density
    J_eff = ω -> begin
        numerator = 2 * α * η_c^2 * ω_c^3 * ω
        denominator = (ω_c^2 - ω^2)^2 + (α * ω)^2
        numerator / denominator
    end

    # Spectral noise function
    S = ω -> begin
        if isapprox(ω, 0.0; atol=1e-9, rtol=0)
            return 2 * T * ( 2 * λ_ν / γ_ν + 2 * α * η_c^2 / ω_c )
        end
        2 / (1 - exp(-ω / T)) * (J_ν(ω) + J_eff(ω))
    end

    return S
end

function eqS62_S64_VegaEtAl_SD(params::Dict)
    # Based on the uniform and markovian limits of eq. S55 in the Supplemental Material of https://pubs.acs.org/doi/10.1021/jacs.5c03182#_i126
    # (i.e. Section B in the Supplemental Material)

    # # Spectator mode phonon bath
    λ_Q = get_dimful_param(params["λ_Q"]) #6000 #invcm
    γ_Q = get_dimful_param(params["γ_Q"]) #0.147 #invcm
    ω_Q = get_dimful_param(params["ω_Q"]) #1370 #invcm
    N = params["N"]  # Number of molecules
    C = get_dimful_param(params["C"])# 0.15 #invcm
    χ = params["χ"] 

    Lambda = N * C^2 / (2 * ω_Q^2)

    # # Cavity parameters
    ω_c = get_dimful_param(params["ω_c"]) #invcm
    η_c = params["η_c"]
    τ_c = get_dimful_param(params["τ_c"]) #500 fs
    γ_c = 1 / τ_c

    # # "parameterize in terms of rabi-splitting?"
    # Ω_R_tilde = 2 * sqrt(N) * η_c * ω_c   

    Gamma = ω -> begin
        2 * λ_Q / γ_Q + 2 * N * χ^2 * ω_c^3 * η_c^2 * γ_c / ( (ω_c^2 - ω^2)^2 + (ω * γ_c)^2 )
    # 2 * λ_Q / γ_Q + 0.5 * Ω_R_tilde^2 * ω_c / τ_c / ( (ω_c^2 - ω^2)^2 + (ω / τ_c)^2 )
    end
    
    R = ω -> begin
       2 * N * χ^2 * ω_c * η_c^2 * ω^2 * (ω^2 - ω_c^2 + γ_c^2) / ((ω_c^2 - ω^2)^2 + (ω * γ_c)^2)
    # 0.5 * Ω_R_tilde^2 * ω^2/ω_c * (ω^2 - ω_c^2 + 1/τ_c^2) / ( (ω_c^2 - ω^2)^2 + (ω / τ_c)^2 )
    end
    
    Jeff = ω -> begin
        Lambda * ω_Q^2 * ω * Gamma(ω) / ( (ω_Q^2 - ω^2 + R(ω))^2 + (ω * Gamma(ω))^2)
    # 0.5 * N * C^2 * Gamma(ω)*ω / ( ω^2*Gamma(ω) + (ω_Q^2 - ω^2 + R(ω))^2 )
    end

    Jeff0 = Lambda * ω_Q^2 * Gamma(0) / ω_Q^4
    # Jeff0 = 0.5 * N * C^2 * Gamma(0) / ω_Q^4

    return Jeff, Jeff0
end

function SystemPhononBath_SD(params::Dict)
    # Simple phonon bath spectral density
    γ_ν = get_dimful_param(params["γ_ν"])
    λ_ν = get(params,"λ_ν",nothing)

    if isnothing(λ_ν)
        M = get_dimful_param(params["M"])
        ω_b = get_dimful_param(params["ω_b"])
        η_ν = params["η_ν"]
        λ_ν = 0.5 * η_ν * M * γ_ν * ω_b
    else
        λ_ν = get_dimful_param(λ_ν)
    end

    # Phonon bath spectral density
    J_ν = ω -> begin
        numerator = 2 * λ_ν * γ_ν * ω
        denominator = ω^2 + γ_ν^2
        numerator / denominator
    end

    J_ν0 = 2 * λ_ν / γ_ν #J_ν/ω evaluated at ω = 0, which is needed for S(0)

    return J_ν, J_ν0
end

function PhononBathPlusSpectatorModesWithCavity(params::Dict)
    J_ν, J_ν0 = SystemPhononBath_SD(params)
    Jeff, Jeff0 = eqS62_S64_VegaEtAl_SD(params)

    T = get_dimful_param(params["T"])
    if T == 0
        S = ω -> begin
            if ω > 0
                return 2 * (J_ν(ω) + Jeff(ω))
            else
                return 0
            end
        end
    else
        S = ω -> begin
            
            if isapprox(ω, 0.0; atol=1e-9, rtol=0)
                return 2 * T * ( J_ν0 + Jeff0 )
            end
            2 / (1 - exp(-ω / T)) * (J_ν(ω) + Jeff(ω))
        end
    end
    return S
end

function S_isotropic(params::Dict)
    J_ν, J_ν0 = SystemPhononBath_SD(params)
    
    theta = params["theta"] * π
    Jparams = initialize_parameters(params["ω_c"], params["η_c"])
    J_iso = w -> J_isotropic(w, theta, Jparams)
    J_iso0 = J_iso(1e-8) / 1e-8
    
    T = get_dimful_param(params["T"])
    if T == 0
        S = ω -> begin
            if ω > 0
                return 2 * (J_ν(ω) + J_iso(ω))
            else
                return 0
            end
        end
    else
        S = ω -> begin
            
            if isapprox(ω, 0.0; atol=1e-9, rtol=0)
                return 2 * T * ( J_ν0 + J_iso0 )
            end
            2 / (1 - exp(-ω / T)) * (J_ν(ω) + J_iso(ω))
        end
    end
    return S
end
