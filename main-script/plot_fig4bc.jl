#!/usr/bin/env julia

"""
Reproduce Figure 4(b,c) (`fig4bc.pdf`), combining the total reaction-rate
sweep with its probability-current decomposition. The sweep totals are
embedded below; the decomposition is loaded from the ratio-1.0
`EnhancementSweep_*.jld2` file next to this script.
"""

include(joinpath(@__DIR__, "..", "common.jl"))
# =============================================================================
# Input data for the upper total-rate panel
# =============================================================================

omegacs = [0, 220, 270, 320, 370, 420, 470, 520, 570, 620, 670, 720, 770, 820, 870, 920, 970, 1020, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310, 1320, 1330, 1340, 1350, 1360, 1370, 1420, 1470, 1520, 1570, 1620, 1670, 1720, 1770, 1820, 1870, 1920, 1970, 2020, 2070, 2120, 2170]
ks = [3.7940183321822514e-6, 3.742202385979781e-6, 3.740256490644683e-6, 3.73690596612497e-6, 3.734510673926734e-6, 3.732302183300758e-6, 3.7336103155946184e-6, 3.7367887844471846e-6, 3.7420338020967603e-6, 3.747806643018377e-6, 3.7568999601780974e-6, 3.7692227356388876e-6, 3.7815443721519727e-6, 3.812534829694876e-6, 3.8473987623405845e-6, 3.910949937390333e-6, 4.01842950400876e-6, 4.220769728529371e-6, 4.6260784099294115e-6, 4.74973643352586e-6, 4.902431279083065e-6, 5.0890517830511974e-6, 5.301126962931856e-6, 5.559356505877585e-6, 5.864624669754931e-6, 6.223217524555175e-6, 6.617338146678098e-6, 7.043442725746948e-6, 7.46504194056426e-6, 7.813261989399179e-6, 8.033430709277648e-6, 8.071263878436278e-6, 7.917881541535931e-6, 7.605595790572847e-6, 7.213444811027539e-6, 6.788838983610306e-6, 6.381387307174468e-6, 6.013840554082595e-6, 5.6938762283943705e-6, 5.423852636015131e-6, 5.195853403822794e-6, 5.004882938207369e-6, 4.85028644758397e-6, 4.710720374616329e-6, 4.598367738771701e-6, 4.503341686759241e-6, 4.422345660593679e-6, 4.352825606914841e-6, 4.293800105523416e-6, 4.09524257686234e-6, 3.989016298810292e-6, 3.927421545291134e-6, 3.889147078010774e-6, 3.864351072063588e-6, 3.852656932283726e-6, 3.8360266984134786e-6, 3.828299588920249e-6, 3.8220883278779834e-6, 3.817198278381275e-6, 3.8057374648301607e-6, 3.8045830469571943e-6, 3.799556702864121e-6, 3.795127459409312e-6, 3.791953513600718e-6, 3.7841839574255077e-6]

@assert length(omegacs) == length(ks)

# Frequencies shown in the lower stacked-rate panel.
selected_omegacs = [500, 1200, 1270, 1500, 2000]

# =============================================================================
# Probability-current helpers
# =============================================================================

"""
Extract the probability-current contributions for transition i -> j at the
time point nearest t0.

The multiplication by 1e6 displays rates in units of 10^-6 fs^-1.
"""
function get_rates_ij(data, i, j; t0=6000)
    ind = argmin(abs.(data.times .- t0))

    r_coh = data.rate_coh[i, j, ind] * 1e6
    r_nu = (
        data.rate_Mnu[i, j, ind] +
        data.rate_nu[i, j, ind]
    ) * 1e6
    r_c = (
        data.rate_Mc[i, j, ind] +
        data.rate_c[i, j, ind]
    ) * 1e6

    return (; r=r_coh + r_nu + r_c, r_coh, r_nu, r_c)
end


"""
Draw one stacked probability-current bar.

Positive and negative contributions are stacked separately, so this also
works when one of the contributions is negative.
"""
function stacked_fillbar!(
    plt,
    rates;
    xcenter,
    width=0.80,
    labels=Dict(
        :r_coh => nothing,
        :r_nu => nothing,
        :r_c => nothing,
    ),
    colors=Dict(
        :r_coh => :lightsalmon,
        :r_nu => :purple,
        :r_c => :mediumseagreen,
    ),
)
    values = (
        r_coh=rates.r_coh,
        r_nu=rates.r_nu,
        r_c=rates.r_c,
    )

    left = xcenter - width / 2
    right = xcenter + width / 2

    positive_base = 0.0
    negative_base = 0.0

    for key in (:r_coh, :r_nu, :r_c)
        value = float(getproperty(values, key))

        if value >= 0
            y0 = positive_base
            y1 = positive_base + value
            positive_base = y1
        else
            y0 = negative_base
            y1 = negative_base + value
            negative_base = y1
        end

        plot!(
            plt,
            [left, right],
            [y1, y1];
            fillrange=[y0, y0],
            fillcolor=colors[key],
            fillalpha=1.0,
            linealpha=0.0,
            linewidth=0,
            label=labels[key],
        )
    end

    return plt
end


# =============================================================================
# Upper panel
# =============================================================================

function make_total_rate_plot(
    omegacs,
    ks;
    selected_omegacs=[1200, 1270, 1500],
    # xlims=(700, 1700),
    xlims=(200, 2200),
    ylims=(3, 9),
)
    rates_scaled = ks .* 1e6

    plt = plot(
        size=(650, 200),
        dpi=300,
        legend=nothing,
        grid=true,
        xlim=xlims,
        ylim=ylims,
        # xticks=[800, 1000, 1200, 1270, 1500],
        # xticks = 400:200:2000,
        xticks = 200:200:2200,
        xminorticks=2,
        xlabel=raw"$\omega_c\quad[\mathsf{cm}^{-1}]$",
        ylabel=raw"$k\quad[\,\!\!\!\times 10^{-6}\mathsf{fs}^{-1}]$",
    )

    # The omega_c = 0 result is treated as the outside-cavity reference.
    hline!(
        plt,
        [rates_scaled[1]];
        color=:blue,
        alpha=0.5,
        lw=1.2,
        label=false,
    )

    # Exclude omega_c = 0 from the frequency-dependent curve.
    plot!(
        plt,
        omegacs[2:end],
        rates_scaled[2:end];
        color=:blue,
        lw=3,
        label=false,
    )

    vline!(
        plt,
        selected_omegacs;
        color=:grey,
        linestyle=:dash,
        lw=1.2,
        alpha=0.9,
        label=false,
    )

    return plt
end


# =============================================================================
# Lower panel
# =============================================================================

function make_rate_plot(
    rate_data,
    Pairs,
    Pair_labels;
    rate_fig_size=(650, 300),
    reverse=false,
    t0=6000,
)
    omegacs = sort(collect(keys(rate_data)))

    n_pairs = length(Pairs)
    n_omegac = length(omegacs)

    ticks = collect(1:(n_pairs + 1) * n_omegac - 1)
    ticklabels = repeat(vcat(Pair_labels, ""), n_omegac)[1:end-1]

    plt = plot(
        xticks=(ticks, ticklabels),
        xlim=(ticks[1]- 1, ticks[end] + 0.5),
        grid=false,
        size=rate_fig_size,
        ylabel=raw"$\mathsf{prob. curr.}\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
    )

    if reverse
        plot!(
            plt,
            ylabel=raw"                        $\mathsf{prob. curr.}\quad[\,\!\!\!\times 10^{-6} \mathsf{fs}^{-1}]$",
        )
    end

    hline!(plt, [0.0], color=:black, linewidth=0.2, alpha=1.0, label="")

    for (n_ω, ωc) in enumerate(omegacs)
        for (n_pair, (i, j)) in enumerate(Pairs)
            rates = get_rates_ij(rate_data[ωc], i, j; t0)
            xcenter = n_pair + (n_ω - 1) * (n_pairs + 1)

            if n_pair == n_pairs && n_ω == n_omegac
                if reverse
                    labels = Dict(
                        :r_coh => raw"$\mathsf{Tunneling}$",
                        :r_nu  => raw"$\mathsf{Solvent}$",
                        :r_c   => raw"$\mathsf{Cavity}$",
                    )
                else
                    labels = Dict(
                        :r_coh => raw"$\mathsf{Tunneling}$",
                        :r_nu  => raw"$\mathsf{Solvent}$",
                        :r_c   => raw"$\mathsf{Cavity}$",
                    )
                end
            else
                labels = Dict(
                    :r_coh => nothing,
                    :r_nu  => nothing,
                    :r_c   => nothing,
                )
            end

            stacked_fillbar!(
                plt,
                rates;
                xcenter,
                labels,
            )
        end

        if n_ω != n_omegac
            plot!(
                plt,
                [n_pairs + 1 + (n_ω - 1) * (n_pairs + 1)],
                seriestype=:vline,
                color=:black,
                label="",
            )
        end

        # Cavity labels only appear on the forward panel in the notebook.
        if !reverse
            xmid = (n_pairs + 1) / 2 + (n_ω - 1) * (n_pairs + 1)
            if ωc == 0
                annotate!(
                    plt,
                    xmid,
                    5.18,
                    (raw"$\mathsf{Outside}\;\,\mathsf{cavity}$", 9, :black),
                )
            else
                annotate!(
                    plt,
                    xmid,
                    5.18,
                    (raw"$\omega_c=" * "$(ωc)" * raw"\mathsf{cm}^{-1}$", 9, :black),
                )
            end
        end
    end

    if reverse
        plot!(
            plt,
            ylim=(-3.4, 0.2),
            left_margin=3Plots.mm,
            legend=(0.1, 0.82),
            legend_background_color=:transparent,
            legend_foreground_color=:transparent,
            top_margin=-7Plots.mm,
            yflip=true,
            grid=true,
        )

        annotate!(
            plt,
            n_omegac * (n_pairs + 1) / 2 + 0.3,
            -2.8,
            (raw"$\mathsf{Reverse}\ \mathsf{Reaction}$", 14, :black),
        )
    else
        plot!(
            plt,
            ylim=(-0.2, 5.35),
            left_margin=3Plots.mm,
            legend=(0.1, 0.82),
            legend_background_color=:transparent,
            legend_foreground_color=:transparent,
            # xtickfontcolor=:transparent,
            bottom_margin=-7Plots.mm,
            grid=true,
        )
    end

    return plt
end


# =============================================================================
# Full figure
# =============================================================================

function make_rate_figure(
    rate_data,
    Pairs,
    Pair_labels;
    selected_omegacs=[1200, 1270, 1500],
    output_file="rate_figure.pdf",
    t0=6000,
)
    plt_k = make_total_rate_plot(
        omegacs,
        ks;
        selected_omegacs=selected_omegacs,
    )

    plt_rates = make_rate_plot(
        rate_data,
        Pairs,
        Pair_labels;
        t0=t0,
    )

    plt_out = plot(
        plt_k,
        plt_rates;
        layout=grid(
            2,
            1,
            heights=[0.35, 0.65],
        ),
        size=(650, 350),
        dpi=300,
        left_margin=5Plots.mm,
        right_margin=3Plots.mm,
        bottom_margin=4Plots.mm,
    )

    output_path = isabspath(output_file) ? output_file : joinpath(@__DIR__, output_file)
    savefig(plt_out, output_path)

    return plt_out
end


Pairs = [[2,1],[3,4],[2,3], [4,1], [5,3], [4,5], [4,2], [1,3]]
Pair_labels = [raw"$\mathsf{\textbf{T}}\!_\textbf{0}$", raw"$\mathsf{\textbf{T}}\!_\textbf{1}$", raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\uparrow}$", raw"$\mathsf{\textbf{T}}\!\!\!_\textbf{\downarrow}$", raw"$\mathsf{\textbf{a}}$", raw"$\mathsf{\textbf{b}}$", raw"$\mathsf{\textbf{c}}$", raw"$\mathsf{\textbf{d}}$"]
rate_data = load(
    require_file(
        joinpath(
            @__DIR__,
            "EnhancementSweep_etanu0.1_gammanu200_etac0.1_ratio1.0_rate_data_tagmarkovianity0.01_ktol1e-8_MixStates12and34_distinguished.jld2"
        )
    ),
    "data",
)
plt_out = make_rate_figure(
    rate_data,
    Pairs,
    Pair_labels;
    selected_omegacs=selected_omegacs,
    output_file="fig4bc.pdf",
    t0=6000,
)


