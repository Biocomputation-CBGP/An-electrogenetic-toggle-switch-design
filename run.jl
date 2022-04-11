using ModelingToolkit
using LinearAlgebra
using DifferentialEquations
using Plots

if "/home/lewis/sauce/julia/An-electrogenetic-toggle-switch-design" ∉ LOAD_PATH
    push!(LOAD_PATH, "/home/lewis/sauce/julia/An-electrogenetic-toggle-switch-design")
end
using Switch

N = 100
M = 4

@variables t u[1:N,1:M](t) R[1:N,1:M](t)
@parameters D[1:M,1:M]=zeros(M,M) S Q η j₀ α₁ α₂ d₃ K₁ K₂

u = collect(u)
R = collect(R)
D = collect(D)
eqs = ReactionDiffusion(u, R, D, t)

electrochemical = [Electrode(u[1,:], R[1,:], Q, η, j₀); Empties(R[2:end,:])]
reactions = Wildtype(u, R, Q, α₁)
substrate = [Empties(R[1:end-1,:]); Supply(u[end,:], R[end,:], D, S)]


rxs = combi.(electrochemical, reactions, substrate)
@named system = ODESystem([eqs; rxs], t, [u; R][:], [D[:]; Q; η; j₀; S; α₁])
simplified = structural_simplify(system)

p = Dict(D[1,1] => 1, D[2,2] => 1)
p = merge(p, Dict(
    Q => 1, η => 0.1, j₀ => 3/100,
    α₁ => 1, α₂ => 1, d₃ => 1, K₁ => 1/2, K₂ => 1/2,
    S => 1
))
u0 = Dict(v => rand() for v ∈ u[:])
for x ∈ u[3,:]
    u0[x] => 1
end

problem = ODEProblem(simplified, u0, (0, 360.0), p)
solution = solve(problem)


