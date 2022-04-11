using ModelingToolkit
const MT = ModelingToolkit
using LinearAlgebra
using DifferentialEquations
using Plots
using LaTeXStrings
using Symbolics

if "/home/lewis/sauce/julia/An-electrogenetic-toggle-switch-design" ∉ LOAD_PATH
    push!(LOAD_PATH, "/home/lewis/sauce/julia/An-electrogenetic-toggle-switch-design")
end
using Switch

@parameters (
    D[1:4,1:4]=zeros(4,4),
    S=1,
    Q=1,
    j₀=1e-6,
    K₁=Q/10,
    K₂=S/10,
    α₁,
    α₂,
    β₁,
    β₂,
    d₃,
    K₃,
    K₄,
    η
)

L = 20e-3
p = merge(
    Dict(D[i,j] => 0 for i ∈ 1:4 for j ∈ 1:4),
    Dict(D[1,1] => 1, D[2,2] => 1),
    Dict(S => 1, Q => 1, K₁ => Q / 10, K₂ => S / 10),
    Dict(α₁ => 1e-5, α₂ => 1, d₃ => 1e3, β₁ => 2, β₂ => 2),
    Dict(K₃ => α₂ / d₃ / 5, K₄ => Q / 20),
    Dict(j₀ => 1e-7, η => 0)
)

function RegulationInset(K=1/2, β=2, X₀=1)
    f = x -> convert(Float32, Switch.inhibition(x, K; β=β).val)
    xs = 0.0:X₀/100:X₀
    opts = merge(Switch.PLOT_OPTS, Dict(
        :ylabel => "Activity",
        :xlabel => "Effector",
        :legend => false
    ))
    plot(xs, f.(xs); opts...)
end


"""Figure 1 panel B shows the wildtype IV response"""
function Fig1B()
    sys = WTMassAction(1)
    u0 = merge(qmap(p, sys), smap(p, sys), xmap(p, sys), amap(p, sys))
    Vs = -0.2:0.01:0.4
    
    results = reduce(
        hcat,
        [v, fixedcurrent(merge(p, Dict(η => v)), u0, sys)] for v ∈ Vs
    )
    IVPlot(results[1,:], results[2,:])
end

"""Figure 1 panel C shows the toggle switch IV response"""
function Fig1C()
    sys = MassAction(1)
    u0 = merge(qmap(p, sys), smap(p, sys), xmap(p, sys), amap(p, sys))
    results = highlowsweep(p, sys, -0.2:0.01:0.4)
    IVPlot(results[:,1], results[:,2:3])
end

function dimensionless_model()
    @parameters P[1:4] β₁ β₂
    @variables τ A(τ) B(τ)
    P = collect(P)
    eqs = [
        Differential(τ)(A) ~ P[1] / (1 + B^β₁) - P[2] * A
        Differential(τ)(B) ~ (1 - P[3]*B) / (1 + A^β₂) - P[4] * B
    ]
    structural_simplify(ODESystem(eqs; name=:dimensionless_model))
end

function Fig2A(;p1=5, p2=1, p3=1/20, p4=1/5, b1=2, b2=2, b0=5)
    @parameters P[1:4] β₁ β₂
    @variables τ A(τ) B(τ)
    p = Dict(P[1]=>p1, P[2]=>p2, P[3]=>p3, P[4]=>p4, β₁=>b1, β₂=>b2)
    a = P[1] / P[2] / (1 + B^β₁)
    b = 1 / (P[4] + P[4] * A^β₂ + P[3])

    nullclineA = substitute(a, p)
    nullclineB = substitute(b, p)
    f(b) = convert(Float32, substitute(nullclineA, Dict(B => b)).val)
    g(a) = convert(Float32, substitute(nullclineB, Dict(A => a)).val)
    Bs = 0.001:0.001:b0
    
    opts = Dict(
        :xlabel => "B (Dimensionless)",
        :ylabel => "A (Dimensionless)",
        :linewidth => 3
    )
    plt = plot(Bs, f.(Bs); label=L"$\frac{dA}{dt} = 0$", merge(Switch.PLOT_OPTS, opts)...)
    plot!(plt, g.(f.(Bs)), f.(Bs); label=L"$\frac{dB}{dt} = 0$", merge(Switch.PLOT_OPTS, opts)...)
    plt
end

function Fig2B()
    sys = dimensionless_model()
    @parameters β₁ β₂ P[1:4]
    @variables τ A(τ) B(τ)
    p = Dict(β₁ => 2, β₂ => 2, P[1] => 5 , P[2] => 1, P[3] => 1/20, P[4] => 1/5)
    
    N = 128
    V = [v for v ∈ 0.1:0.005:0.3]
    results = []
    for i ∈ 1:length(V)
        p[P[4]] = V[i]
        points = []
        for n ∈ 1:N
            u0 = Dict(A => 5*rand(), B => 5*rand())
            push!(results, [V[i], fixedpoint(p, u0, sys)[2]])
        end
    end
    results = reduce(hcat, results)'

    opts = Dict(:xlabel => L"$P_4$ (dimensionless)", :ylabel => L"$B$ (dimensionless)",
                :yguidefontcolor => RGBA(1, 176/255, 0, 1),
                :legend => false,
                :seriescolor => RGBA(1, 176/255, 0, 1),
                :markersize => 3)
    scatter(results[:,1], results[:,2]; merge(Switch.PLOT_OPTS, opts)...)
end


function Fig2C()
    sys = MassAction(1)
    u0 = merge(qmap(p, sys), smap(p, sys), xmap(p, sys), amap(p, sys))
    results = highlowsweep(p, sys, 0.0:0.005:0.4)
    IVPlot(results[:,1], results[:,2:3])
end

function Fig2D()
    sys = MassAction(1)
    steps = [(120*60.0, 0.4), (720*60.0, 0.3), (2400*60.0, -0.2), (3000*60.0, 0.3), (4680*60.0, 0.3)]
    q = copy(p)
    q[η] = 0.3
    u0 = lowmap(q, sys)
    start = fixedpoint(q, u0, sys)
    u0 = Dict(states(sys)[i] => start[i] for i ∈ 1:length(start))
    tspan = (first(first(steps)), first(last(steps)))
    u, i = MT.getvar(sys, :u), MT.getvar(sys, :i)
    results = switching(sys, u0, 1.0, q, steps, [i, u[1,4], u[1,1]])
    opts = Dict(:ylabel => "I (Amperes)", :xlabel => "Time (Hours)",
                :legend => false,
                :linewidth => 2,
                :grid => false,
                :leftmargin => Plots.mm * 5,
                :bottommargin => Plots.mm * 5,
                :xformatter => x -> string((x / 3600) |> round |> Int),
                :size => (750, 200))

    return SwitchPlot(results; opts...)
end


function Fig3B(N)
    q = copy(p)
    q[η] = 0.3
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    αs = [1e-7, 1e-3]
    solutions = [fixedpoint(merge(q, Dict(α₁=>v / N)),
                            merge(qmap(merge(q, Dict(α₁=>v / N)), sys),
                                  smap(merge(q, Dict(α₁=>v / N)), sys),
                                  xmap(merge(q, Dict(α₁=>v / N)), sys),
                                  amap(merge(q, Dict(α₁=>v / N)), sys)),
                            sys)
                 for v ∈ αs]
    ChargeGradient(
        sys,
        solutions,
        labels=["Low" "High"],
        legend=:right
    )
end

function Fig3C(N)
    q = copy(p)
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    q[α₁] = p[α₁] / N
    sys = MassAction(N)
    u0 = merge(qmap(q, sys), smap(q, sys), xmap(q, sys), amap(q, sys))
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(results[:,1], results[:,2:3]; size=(250, 200))
end

function Fig3D(N)
    q = copy(p)
    q[η] = 0.3
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    αs = [1e-7, 1e-3]

    solutions = [fixedpoint(merge(q, Dict(α₁=>v / N)),
                            merge(qmap(merge(q, Dict(α₁=>v / N)), sys),
                                  smap(merge(q, Dict(α₁=>v / N)), sys),
                                  xmap(merge(q, Dict(α₁=>v / N)), sys),
                                  amap(merge(q, Dict(α₁=>v / N)), sys)),
                            sys)
                 for v ∈ αs]
    SubstrateGradient(
        sys,
        solutions,
        labels=["Low" "High"],
        legend=:left
    )
end

function Fig3E(N)
    q = copy(p)
    q[η] = 0.3
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    αs = [1e-7, 1e-3]

    solutions = [fixedpoint(merge(q, Dict(α₁=>v / N)),
                            merge(qmap(merge(q, Dict(α₁=>v / N)), sys),
                                  smap(merge(q, Dict(α₁=>v / N)), sys),
                                  xmap(merge(q, Dict(α₁=>v / N)), sys),
                                  amap(merge(q, Dict(α₁=>v / N)), sys)),
                            sys)
                 for v ∈ αs]
    RepressorGradient(
        sys,
        solutions,
        labels=["Low" "High"],
        legend=:right
    )
end

function Fig4A(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1000 / (L / N)^2
    sys = MassAction(N)
    u0 = merge(qmap(q, sys), smap(q, sys), xmap(q, sys), amap(q, sys))
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.05Q$",
        legend=:right
    )
end


function Fig4B(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.15
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1000 / (L / N)^2
    sys = MassAction(N)
    u0 = merge(qmap(q, sys), smap(q, sys), xmap(q, sys), amap(q, sys))
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.15Q$",
        legend=:right
    )
end

function Fig4C(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.3
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1000 / (L / N)^2
    sys = MassAction(N)
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.3Q$",
        legend=:right
    )
end

function Fig4D(N)
    sys = MassAction(N)
    steps = [(120*60.0, 1.0), (60*1024*60.0, 0.3)]
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.15
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1000 / (L / N)^2
    q[η] = 0.3
    u0 = lowmap(q, sys)
    start = fixedpoint(q, u0, sys)
    u0 = Dict(states(sys)[i] => start[i] for i ∈ 1:length(start))
    tspan = (first(first(steps)), first(last(steps)))
    u, i = MT.getvar(sys, :u), MT.getvar(sys, :i)
    results = switching(sys, u0, tspan[2] / 1024, q, steps, [i])
    opts = Dict(:ylabel => "I (Amperes)", :xlabel => "Time (Hours)",
                :legend => false,
                :linewidth => 2,
                :grid => false,
                :ylims => (0.0, 1.0),
                :leftmargin => Plots.mm * 5,
                :bottommargin => Plots.mm * 5,
                :xformatter => x -> string((x / 3600) |> round |> Int),
                :size => (750/2, 200))

    return SwitchPlot(results; opts...)
end

function Fig4E()
    qs = 0.001:0.001:1
    K  = 0.05
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig4F()
    qs = 0.001:0.001:1
    K  = 0.15
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig4G()
    qs = 0.001:0.001:1
    K  = 0.3
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig5D(N)
    sys = MassAction(N)
    steps = [(120*60.0, 0.4), (120*1024*60.0, 0.3)]
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.2
    q[β₁] = 4
    q[β₂] = 4
    q[D[1,1]] = 1e-9 / (L / N)^2
    q[D[2,2]] = 1e-9 / (L / N)^2
    q[η] = 0.3
    u0 = lowmap(q, sys)
    start = fixedpoint(q, u0, sys)
    u0 = Dict(states(sys)[i] => start[i] for i ∈ 1:length(start))
    tspan = (first(first(steps)), first(last(steps)))
    u, i = MT.getvar(sys, :u), MT.getvar(sys, :i)
    results = switching(sys, u0, tspan[2] / 1024, q, steps, [i])
    opts = Dict(:ylabel => "I (Amperes)", :xlabel => "Time (Hours)",
                :legend => false,
                :linewidth => 2,
                :grid => false,
                :ylims => (0.0, 1.0),
                :leftmargin => Plots.mm * 5,
                :bottommargin => Plots.mm * 5,
                :xformatter => x -> string((x / 3600) |> round |> Int),
                :size => (750/2, 200))

    return SwitchPlot(results; opts...)
end


function Fig5A(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    u0 = merge(qmap(q, sys), smap(q, sys), xmap(q, sys), amap(q, sys))
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.05Q$",
        legend=:right
    )
end


function Fig5B(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.1
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    u0 = merge(qmap(q, sys), smap(q, sys), xmap(q, sys), amap(q, sys))
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.1Q$",
        legend=:right
    )
end

function Fig5C(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = Q * 0.2
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$K_4 = 0.2Q$",
        legend=:right
    )
end

function Fig5E(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = 0.2
    q[β₁] = 2
    q[β₂] = 2
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$\beta_{1,2} = 2$",
        legend=:right
    )
end

function Fig5F(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = 0.2
    q[β₁] = 3
    q[β₂] = 3
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$\beta_{1,2} = 3$",
        legend=:right
    )
end

function Fig5G(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[K₄] = 0.2
    q[β₁] = 4
    q[β₂] = 4
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        label=L"$\beta_{1,2} = 4$",
        legend=:right
    )
end

function Fig4E()
    qs = 0.001:0.001:1
    K  = 0.05
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig4F()
    qs = 0.001:0.001:1
    K  = 0.15
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig5I()
    qs = 0.001:0.001:1
    K = 0.2
    f(x) = convert(Float32, Switch.inhibition(x, K, β=2).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig5J()
    qs = 0.001:0.001:1
    K = 0.2
    f(x) = convert(Float32, Switch.inhibition(x, K, β=3).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig5K()
    qs = 0.001:0.001:1
    K = 0.2
    f(x) = convert(Float32, Switch.inhibition(x, K, β=4).val)
    plot(qs, f.(qs);
         xlabel=L"Charge $\left(\frac{mol}{m^3}\right)$",
         ylabel="Normalised Expression",
         xguidefontcolor = RGBA(1, 176/255, 0, 1),
         legend=false,
         Switch.PLOT_OPTS...
             )  
end

function Fig5H(N)
    q = copy(p)
    q[α₁] = 1e-5/N
    q[D[1,1]] = 1e-10 / (L / N)^2
    q[D[2,2]] = 1e-10 / (L / N)^2
    sys = MassAction(N)

    plots = []
    factors = [2, 3, 4, 5]
    for f ∈ factors
        q[β₁] = f
        for g ∈ factors
            q[β₂] = g
            results = highlowsweep(q, sys, 0.0:0.05:0.6)
            push!(plots,
                  IVPlot(
                      repeat(results[:,1], 2),
                      vcat(results[:,2], results[:,3]),
                      legend=:false
                  ))
        end
    end
    return plots
    plot(plots..., layout=(length(factors), length(factors)))
end

function Fig6A(N)
    q = copy(p)
    q[K₄] = 0.2
    q[β₁] = 4
    q[β₂] = 4
    q[α₁] = 1e-5/N
    sys = MassAction(N)
    q[D[1,1]] = 1e-10 / (10e-3 / N)^2
    q[D[2,2]] = 1e-10 / (10e-3 / N)^2
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        legend=:right,
        label=L"Depth $10\mu m$"
    )
end

function Fig6B(N)
    q = copy(p)
    q[K₄] = 0.2
    q[β₁] = 4
    q[β₂] = 4
    q[α₁] = 1e-5/N
    sys = MassAction(N)
    q[D[1,1]] = 1e-10 / (20e-3 / N)^2
    q[D[2,2]] = 1e-10 / (20e-3 / N)^2
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        legend=:right,
        label=L"Depth $20\mu m$"
    )
end

function Fig6C(N)
    q = copy(p)
    q[K₄] = 0.2
    q[β₁] = 4
    q[β₂] = 4
    q[α₁] = 1e-5/N
    sys = MassAction(N)
    q[D[1,1]] = 1e-10 / (40e-3 / N)^2
    q[D[2,2]] = 1e-10 / (40e-3 / N)^2
    results = highlowsweep(q, sys, 0.0:0.01:0.6)
    IVPlot(
        repeat(results[:,1], 2),
        vcat(results[:,2], results[:,3]),
        legend=:right,
        label=L"Depth $40\mu m$"
    )
end

function dofigures()
    print("Remake 1B")
    panel = Fig1B()
    savefig(panel, "/home/lewis/images/toggle/figure-1-B.svg")
    savefig(panel, "/home/lewis/images/toggle/figure-1-B.ps")

    print(", 1C")
    panel = Fig1C()
    savefig(panel, "/home/lewis/images/toggle/figure-1-C.svg")
    savefig(panel, "/home/lewis/images/toggle/figure-1-C.ps")

    print(", 2A")
    panel = Fig2A()
    savefig(panel, "/home/lewis/images/toggle/figure-2-A.svg")
    savefig(panel, "/home/lewis/images/toggle/figure-2-A.ps")

    print(", 2B")
    panel = Fig2B()
    savefig(panel, "/home/lewis/images/toggle/figure-2-B.svg")
    savefig(panel, "/home/lewis/images/toggle/figure-2-B.ps")

    print(", 2C")
    panel = Fig2C()
    savefig(panel, "/home/lewis/images/toggle/figure-2-C.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-2-C.svg")

    print(", 2D")
    panel = Fig2D()
    savefig(panel, "/home/lewis/images/toggle/figure-2-D.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-2-D.svg")

    print(", 3B")
    panel = Fig3B(21)
    savefig(panel, "/home/lewis/images/toggle/figure-3-B.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-3-B.svg")

    print(", 3C")
    panel = Fig3C(21)
    savefig(panel, "/home/lewis/images/toggle/figure-3-C.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-3-C.svg")
    
    print(", 3D")
    panel = Fig3D(21)
    savefig(panel, "/home/lewis/images/toggle/figure-3-D.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-3-D.svg")

    print(", 4A")
    panel = Fig4A(21)
    savefig(panel, "/home/lewis/images/toggle/figure-4-A.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-A.svg")

    print(", 4B")
    panel = Fig4B(21)
    savefig(panel, "/home/lewis/images/toggle/figure-4-B.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-B.svg")

    print(", 4C")
    panel = Fig4C(21)
    savefig(panel, "/home/lewis/images/toggle/figure-4-C.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-C.svg")

    print(", 4D")
    panel = Fig4D(21)
    savefig(panel, "/home/lewis/images/toggle/figure-4-D.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-D.svg")

    print(", 4E")
    panel = Fig4E()
    savefig(panel, "/home/lewis/images/toggle/figure-4-E.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-E.svg")

    print(", 4F")
    panel = Fig4F()
    savefig(panel, "/home/lewis/images/toggle/figure-4-F.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-F.svg")

    print(", 4G")
    panel = Fig4G()
    savefig(panel, "/home/lewis/images/toggle/figure-4-G.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-4-G.svg")

    print(", 5A")
    panel = Fig5A(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-A.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-A.svg")

    print(", 5B")
    panel = Fig5B(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-B.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-B.svg")

    print(", 5C")
    panel = Fig5C(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-C.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-C.svg")

    print(", 5E")
    panel = Fig5E(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-E.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-E.svg")

    print(", 5F")
    panel = Fig5F(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-F.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-F.svg")

    print(", 5G")
    panel = Fig5G(21)
    savefig(panel, "/home/lewis/images/toggle/figure-5-G.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-G.svg")

    print(", 5I")
    panel = Fig5I()
    savefig(panel, "/home/lewis/images/toggle/figure-5-I.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-I.svg")

    print(", 5J")
    panel = Fig5J()
    savefig(panel, "/home/lewis/images/toggle/figure-5-J.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-J.svg")

    print(", 5K")
    panel = Fig5K()
    savefig(panel, "/home/lewis/images/toggle/figure-5-K.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-5-K.svg")

    print(", 6A")
    panel = Fig6A(21)
    savefig(panel, "/home/lewis/images/toggle/figure-6-A.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-6-A.svg")

    print(", 6B")
    panel = Fig6B(21)
    savefig(panel, "/home/lewis/images/toggle/figure-6-B.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-6-B.svg")

    print(", 6C")
    panel = Fig6C(21)
    savefig(panel, "/home/lewis/images/toggle/figure-6-C.ps")
    savefig(panel, "/home/lewis/images/toggle/figure-6-C.svg")

    # print(", 5D")
    # panel = Fig5D(21)
    # savefig(panel, "/home/lewis/images/toggle/figure-5-D.ps")
    # savefig(panel, "/home/lewis/images/toggle/figure-5-D.svg")
end
