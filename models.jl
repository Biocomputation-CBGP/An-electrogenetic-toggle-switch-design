function add_rhs(eqs::Vararg{Equation})
    @assert all(b -> isequal(eqs[1].lhs)(b.lhs), eqs) "The equations do not have the same RHS"
    first(eqs).lhs ~ sum(getproperty.(eqs, :rhs))
end

function make_equations(func::Function, args...)
    func((x -> x isa Symbolics.Arr ? collect(x) : x).(args)...)
end

function ReactionDiffusion(u::T, R::T, D::T, t::Num) where {T<:Matrix{Num}}
    N, M = size(u)
    d = Differential(t)

    function dudx(i)
        if N > 1
            if 1 < i < N
                u[i-1,:] - 2u[i,:] + u[i+1,:]
            elseif i == N
                u[end-1,:] - u[end,:]
            elseif i == 1
                u[2,:] - u[1,:]
            end
        else
            zeros(M)
        end
    end
    reshape(reduce(vcat, [d.(u[i,:]) .~ D * (dudx(i)) + R[i,:] for i ∈ 1:N]), N, M)
end

Empty(R::Matrix{Num}) = collect(R .~ 0)
function Electrode(u::Matrix{Num}, R::Matrix{Num}, Q, η, j₀)
    eqs = Empty(R)
    ieq = -j₀ * (u[1,1] * exp(η*C) + u[1,1] / exp(η*C) - Q / exp(η*C))
    eqs[1, 1] = add_rhs(eqs[1, 1], R[1, 1] ~ ieq)
    eqs
end

function Supply(u::Matrix{Num}, R::Matrix{Num}, D::Matrix{Num}, S)
    eqs = Empty(R)
    eqs[end, 2] = add_rhs(eqs[end, 2], R[end,2] ~ D[2, 2] * (S - u[end, 2]))
    eqs
end

function Wildtype(u::Matrix{Num}, R::Matrix{Num}, Q, α₁, K₁, K₂)
    eqs = Empty(R)
    N, _ = size(R)
    for i ∈ 1:N
        flux = u[i, 3] * α₁ / (K₁ / (Q - u[i, 1]) + 1) / (K₂ / u[i, 2] + 1)
        eqs[i, 1] = add_rhs(eqs[i, 1], R[i,1] ~ flux)
        eqs[i, 2] = add_rhs(eqs[i, 2], R[i,2] ~ -flux)
    end
    eqs
end

function Wildtype(u::Matrix{Num}, R::Matrix{Num}, Q, α₁)
    eqs = Empty(R)
    N, _ = size(R)
    for i ∈ 1:N
        flux = u[i, 3] * α₁ * (Q - u[i, 1]) * u[i, 2]
        eqs[i, 1] = add_rhs(eqs[i, 1], R[i,1] ~ flux)
        eqs[i, 2] = add_rhs(eqs[i, 2], R[i,2] ~ -flux)
    end
    eqs
end

function Toggle(u::T, R::T, Q, α₁, α₂, β₁, β₂, d₃, K₃, K₄) where {T<:Matrix{Num}}
    eqs = Empty(R)
    N, _ = size(R)
    for i ∈ 1:N
        flux = u[i,3] * α₁ * (Q-u[i,1]) * u[i,2] * inhibition(u[i,4],K₃,β=β₁)
        eqs[i, 1] = add_rhs(eqs[i, 1], R[i, 1] ~ flux)
        eqs[i, 2] = add_rhs(eqs[i, 2], R[i, 2] ~ -flux)
        eqs[i, 4] = add_rhs(eqs[i, 4], R[i, 4] ~ α₂*inhibition(u[i,1],K₄,β=β₂)-d₃*u[i,4])
    end
    eqs
end

function Toggle(u::T, R::T, Q, α₁, α₂, β₁, β₂, d₃, K₁, K₂, K₃, K₄) where {T<:Matrix{Num}}
    eqs = Empty(R)
    N, _ = size(R)
    for i ∈ 1:N
        flux = u[i,3] * α₁ / (K₁/(Q-u[i,1])+1) / (K₂/u[i,2]+1) * inhibition(u[i,4],K₃,β=β₁)
        eqs[i, 1] = add_rhs(eqs[i, 1], R[i, 1] ~ flux)
        eqs[i, 2] = add_rhs(eqs[i, 2], R[i, 2] ~ -flux)
        eqs[i, 4] = add_rhs(eqs[i, 4], R[i, 4] ~ α₂*inhibition(u[i,1],K₄,β=β₂)-d₃*u[i,4])
    end
    eqs
end

j0 = 1e-7
S0 = 1
Q0 = 1
K1 = Q0 / 10
K2 = S0 / 10

function WT(N)
    @parameters (
        D[1:4,1:4] = zeros(4,4),
        S = S0,
        Q = Q0,
        j₀ = j0,
        K₁ = K1,
        K₂ = K2
    )  
    @parameters α₁ η
    @variables t u[1:N,1:4](t) R[1:N,1:4](t) i(t)
    u, R, D = collect(u), collect(R), collect(D)

    eqs = make_equations(ReactionDiffusion, u, R, D, t)
    electrochemical = make_equations(Electrode, u, R, Q, η, j₀)
    reactions = make_equations(Wildtype, u, R, Q, α₁, K₁, K₂)
    substrate = make_equations(Supply, u, R, D, S)

    rxs = add_rhs.(electrochemical, reactions, substrate)
    obs = [i ~ -electrochemical[1].rhs * Switch.F]
    structural_simplify(ODESystem(
        [
            eqs[:];
            rxs[:];
            obs
        ],
        t,
        [[u; R][:]; i],
        [D[:]; Q; j₀; S; α₁; K₁; K₂; η];
        name=:model
    ))
end

function WTMassAction(N)
    @parameters (
        D[1:4,1:4] = zeros(4,4),
        S = S0,
        Q = Q0,
        j₀ = j0,
        K₁ = K1,
        K₂ = K2
    )  
    @parameters α₁ η
    @variables t u[1:N,1:4](t) R[1:N,1:4](t) i(t)
    u, R, D = collect(u), collect(R), collect(D)

    eqs = make_equations(ReactionDiffusion, u, R, D, t)
    electrochemical = make_equations(Electrode, u, R, Q, η, j₀)
    reactions = make_equations(Wildtype, u, R, Q, α₁)
    substrate = make_equations(Supply, u, R, D, S)

    rxs = add_rhs.(electrochemical, reactions, substrate)
    obs = [i ~ -electrochemical[1].rhs * Switch.F]
    structural_simplify(ODESystem(
        [
            eqs[:];
            rxs[:];
            obs
        ],
        t,
        [[u; R][:]; i],
        [D[:]; Q; j₀; S; α₁; η];
        name=:model
    ))
end

function Model(N)
    @parameters (
        D[1:4,1:4] = zeros(4,4),
        S = S0,
        Q = Q0,
        j₀ = j0,
        K₁ = K1,
        K₂ = K2
    )  
    @parameters α₁ α₂ β₁ β₂ d₃ K₃ K₄ η
    @variables t u[1:N,1:4](t) R[1:N,1:4](t) i(t)
    u, R, D = collect(u), collect(R), collect(D)

    eqs = make_equations(ReactionDiffusion, u, R, D, t)
    electrochemical = make_equations(Electrode, u, R, Q, η, j₀)
    reactions = make_equations(Toggle, u, R, Q, α₁, α₂, β₁, β₂, d₃, K₁, K₂, K₃, K₄)
    substrate = make_equations(Supply, u, R, D, S)
    rxs = add_rhs.(electrochemical, reactions, substrate)

    obs = [i ~ -electrochemical[1].rhs * Switch.F]
    structural_simplify(ODESystem(
        [
            eqs[:];
            rxs[:];
            obs
        ],
        t,
        [[u; R][:]; i],
        [D[:]; Q; j₀; S; α₁; α₂; β₁; β₂; d₃; K₁; K₂; K₃; K₄; η];
        name=:model
    ))
end

function MassAction(N)
    @parameters (
        D[1:4,1:4] = zeros(4,4),
        S = S0,
        Q = Q0,
        j₀ = j0,
        K₁ = K1,
        K₂ = K2
    )
    @parameters α₁ α₂ β₁ β₂ d₃ K₃ K₄ η
    @variables t u[1:N,1:4](t) R[1:N,1:4](t) i(t)
    u, R, D = collect(u), collect(R), collect(D)

    eqs = make_equations(ReactionDiffusion, u, R, D, t)
    electrochemical = make_equations(Electrode, u, R, Q, η, j₀)
    reactions = make_equations(Toggle, u, R, Q, α₁, α₂, β₁, β₂, d₃, K₃, K₄)
    substrate = make_equations(Supply, u, R, D, S)
    rxs = add_rhs.(electrochemical, reactions, substrate)

    obs = [i ~ -electrochemical[1].rhs * Switch.F]
    structural_simplify(ODESystem(
        [
            eqs[:];
            rxs[:];
            obs
        ],
        t,
        [[u; R][:]; i],
        [D[:]; Q; j₀; S; α₁; α₂; β₁; β₂; d₃; K₃; K₄; η];
        name=:model
    ))
end


