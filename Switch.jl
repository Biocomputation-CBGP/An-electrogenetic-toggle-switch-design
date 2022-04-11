module Switch

using Symbolics
using ModelingToolkit
using Catalyst
using Plots
using DifferentialEquations
using BifurcationKit
using Setfield
using LinearAlgebra

const MT = ModelingToolkit
const C = 18.3
const F = 96485.3

include("models.jl")
export WT, Model, MassAction, WTMassAction
include("analysis.jl")
export fixedpoint, fixedcurrent, current, montecarlo, highlowsweep, bifurcate, switching
export qmap, smap, amap, xmap, highmap, lowmap
include("plotting.jl")
export IVPlot, SwitchPlot, ChargeGradient, SubstrateGradient, RepressorGradient

inhibition(x, K; β=2)::Num = K^β / (K^β + x^β)
activation(x, K; β=2)::Num = x^β / (K^β + x^β)

function gradient_animation(model, solution, ts, fname)
    q_idxs = [findfirst(isequal(x), states(model)) for x ∈ @nonamespace model.u[:,1]]
    s_idxs = [findfirst(isequal(x), states(model)) for x ∈ @nonamespace model.u[:,2]]

    anim = @animate for t ∈ ts
        series = hcat(solution(t)[q_idxs], solution(t)[s_idxs])
        plot(series, label=["Charge" "Substrate"], legend=:bottomright, ylims=(0.0, 1.0))
        plot!(title="t = $(t) seconds")
    end

    gif(anim, fname, fps=16)
end

end
