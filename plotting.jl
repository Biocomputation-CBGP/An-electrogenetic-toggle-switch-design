using Plots
using ModelingToolkit
using LaTeXStrings
const MT = ModelingToolkit

const PLOT_OPTS = Dict(
    :guidefontsize => 8,
    :tickfontsize => 8,
    :legendfontsize => 8,
    :grid => true,
    :size => (250, 200),
    :dpi => 900
)

function IVPlot(V, I; kwargs...)
    opts = merge(PLOT_OPTS, Dict(
        :xlabel => L"$V$ (Volts)",
        :xguidefontcolor => RGBA(1, 97/255, 0, 1),
        :ylabel => L"$I$ (Amperes)",
        :yguidefontcolor => RGBA(1, 176/255, 0, 1),
        :seriescolor => RGBA(1, 176/255, 0, 1),
        :markersize => 3,
        :markershape => :circle,
        :legend => false,
        :size => (250, 200),
        :guidefontsize => 8,
        :tickfontsize => 8,
    ), kwargs)
    scatter(V, I; opts...)
end

function QVPlot(V, Q; kwargs...)
    opts = merge(PLOT_OPTS, Dict(
        :xlabel => L"$V$ (Volts)",
        :xguidefontcolor => RGBA(1, 97/255, 0, 1),
        :ylabel => L"$Q$ $\left(\frac{mol}{m^3}\right)$",
        :yguidefontcolor => RGBA(1, 176/255, 0, 1),
        :seriescolor => RGBA(1, 176/255, 0, 1),
        :markersize => 3,
        :markershape => :circle,
        :legend => false,
        :size => (250, 200),
        :guidefontsize => 8,
        :tickfontsize => 8,
    ), kwargs)
    scatter(V, Q; opts...)
end

function SwitchPlot(results; kwargs...)
    opts = merge(PLOT_OPTS, Dict(
        :ylabel => L"$I$ (Amperes)",
        :xlabel => "Time (Hours)",
        :legend => false,
        :xformatter => x -> string((x / 3600) |> round |> Int),
        :rightmargin => Plots.mm * 15
    ), kwargs)

    plt = plot(
        results[:, 1],
        results[:, 3];
        seriescolor=RGBA(1, 176/255, 0, 1),
        linewidth=4,
        yguidefontcolor=RGBA(1, 176/255, 0, 1),
        opts...,
        grid=true
    )
    
    voltage = twinx(plt)
    plot!(
        voltage,
        results[:,1],
        results[:,2];
        opts...,
        ylabel=L"$V$ (Volts)",
        seriescolor=RGBA(1, 96/255, 0, 1),
        linestyle=:dash,
        yguidefontcolor=RGBA(1, 96/255, 0, 1),
        xformatter=x->"",
        xlabel="",
        grid=false
    )
    plt
end

function ChargeGradient(model, solutions; kwargs...)
    u = collect(MT.getvar(model, :u))
    vars = [u[i, 1] for i ∈ 1:size(u, 1)]
    idxs = [findfirst(isequal(x), states(model)) for x ∈ vars]
    opts = merge(PLOT_OPTS, kwargs, Dict(
        :xlabel => L"Position in biofilm $x$ ($\mu m$)",
        :ylabel => L"$\textrm{Charge}$ $\left(\frac{mol}{m^3}\right)$",
        :linewidth => 3
    ))

    data = reduce(hcat, sol[idxs] for sol ∈ solutions)
    plot(0:length(idxs)-1, data; opts...)
end

function SubstrateGradient(model, solutions; kwargs...)
    u = collect(MT.getvar(model, :u))
    vars = [u[i, 2] for i ∈ 1:size(u, 1)]
    idxs = [findfirst(isequal(x), states(model)) for x ∈ vars]
    opts = merge(PLOT_OPTS, kwargs, Dict(
        :xlabel => L"Position in biofilm $x$ ($\mu m$)",
        :ylabel => L"Substrate $\left(\frac{mol}{m^3}\right)$",
        :linewidth => 3
    ))

    data = reduce(hcat, sol[idxs] for sol ∈ solutions)
    plot(0:length(idxs)-1, data; opts...)
end

function RepressorGradient(model, solutions; kwargs...)
    u = collect(MT.getvar(model, :u))
    vars = [u[i, 4] for i ∈ 1:size(u, 1)]
    idxs = [findfirst(isequal(x), states(model)) for x ∈ vars]
    opts = merge(PLOT_OPTS, kwargs, Dict(
        :xlabel => L"Position in biofilm $x$ ($\mu m$)",
        :ylabel => L"Repressor $\left(\frac{mol}{m^3}\right)$",
        :linewidth => 3
    ))

    data = reduce(hcat, sol[idxs] for sol ∈ solutions)
    plot(0:length(idxs)-1, data; opts...)
end
