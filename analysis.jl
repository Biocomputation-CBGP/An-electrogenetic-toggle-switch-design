

function switching(sys::ODESystem, x0, dt, p, steps, vars)
    @parameters η
    time = 0
    results = []
    q = copy(p)
    u0 = MT.varmap_to_vars(x0, states(sys); defaults=sys.defaults)
    for step ∈ steps
        problem = ODEProblem(sys, u0, (time, step[1]), q)
        solution = solve(problem, Rodas5(), abstol=1e-15, reltol=1e-15, saveat=dt)
        push!(results, hcat(
            solution.t,
            fill(convert(Float64, q[η].val), length(solution.t)),
            [solution[x] for x ∈ vars]...))
        time = step[1]
        q[η] = step[2]
        u0 = solution[end]
    end
    reduce(vcat, results)
end

function current(sys::ODESystem)
    i = MT.build_explicit_observed_function(sys, MT.getvar(sys, :i))
    (u, p) -> i(u, p, 0)
end

function fixedcurrent(p::Dict, u0::Dict, sys::ODESystem)
    current(sys)(
        fixedpoint(p, u0, sys),
        MT.varmap_to_vars(p, parameters(sys); defaults=sys.defaults)
    )
end

function fixedpoint(p::Dict, u0::Dict, sys::ODESystem)
    df = ODEFunction(sys)
    f(x) = x[1] => substitute(x[2], p)
    u0 = MT.varmap_to_vars(Dict(f.(collect(u0))), states(sys); defaults=sys.defaults)
    p = MT.varmap_to_vars(Dict(f.(collect(p))), parameters(sys); defaults=sys.defaults)
    u = zeros(length(states(sys)))
    opts = Dict(:abstol => 1e-15, :reltol => 1e-15, :cb => PositiveDomain(u))
    problem = SteadyStateProblem(df, u0, p; opts...)
    solve(problem, DynamicSS(Rodas5(); abstol=1e-15, reltol=1e-15)).u
end

function qmap(p::Dict, sys::ODESystem)
    q = rand()
    u, Q = MT.getvar(sys, :u), MT.getvar(sys, :Q)
    Dict(u[i,1] => substitute(q * Q, p) for i ∈ 1:size(u, 1))
end

function smap(p::Dict, sys::ODESystem)
    s = rand()
    u, S = MT.getvar(sys, :u), MT.getvar(sys, :S)
    Dict(u[i,2] => substitute(s * S, p) for i ∈ 1:size(u, 1))
end

function xmap(p::Dict, sys::ODESystem)
    u = MT.getvar(sys, :u)
    Dict(u[i, 3] => 1 for i ∈ 1:size(u, 1))
end

function amap(p::Dict, sys::ODESystem)
    u = MT.getvar(sys, :u)
    @parameters α₂ d₃
    if α₂ ∈ keys(p) && d₃ ∈ keys(p)
        a = rand()
        Dict(u[i,4] => substitute(a * α₂ / d₃, p) for i ∈ 1:size(u, 1))
    else
        Dict(u[i,4] => 0 for i ∈ 1:size(u, 1))
    end
end

function highmap(p::Dict, sys::ODESystem)
    u = MT.getvar(sys, :u)
    @parameters Q S α₂ d₃
    if α₂ ∈ keys(p) && d₃ ∈ keys(p)
        merge(
            Dict(u[i,4] => substitute(α₂ / d₃, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,1] => substitute(0.0001 * Q, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,2] => substitute(S, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,3] => 1 for i ∈ 1:size(u, 1))
        )
    else
        merge(
            Dict(u[i,4] => 0 for i ∈ 1:size(u, 1)),
            Dict(u[i,1] => substitute(0.0001 * Q, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,2] => substitute(S, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,3] => 1 for i ∈ 1:size(u, 1))
        )
    end
end


function lowmap(p::Dict, sys::ODESystem)
    u = MT.getvar(sys, :u)
    @parameters Q S α₂ d₃
    if α₂ ∈ keys(p) && d₃ ∈ keys(p)
        merge(
            Dict(u[i,4] => substitute(0.0001 * α₂ / d₃, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,1] => substitute(0.9999 * Q, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,2] => substitute(S, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,3] => 1 for i ∈ 1:size(u, 1))
        )
    else
        merge(
            Dict(u[i,4] => 0 for i ∈ 1:size(u, 1)),
            Dict(u[i,1] => substitute(0.9999 * Q, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,2] => substitute(S, p) for i ∈ 1:size(u, 1)),
            Dict(u[i,3] => 1 for i ∈ 1:size(u, 1))
        )
    end
end

function montecarlo(p::Dict, sys::ODESystem, N::Integer)
    @parameters η
    results = [p[η]; zeros(N)]
    for n ∈ 1:N
        u0 = merge(qmap(p, sys), smap(p, sys), xmap(p, sys), amap(p, sys))
        results[n + 1] = fixedcurrent(p, u0, sys)
    end
    results
end

function highlowsweep(p::Dict, sys::ODESystem, Vs)
    @parameters η
    results = zeros(length(Vs), 3)
    for (i, v) ∈ enumerate(Vs)
        results[i, 1] = v
        results[i, 2] = fixedcurrent(merge(p, Dict(η => v)), highmap(p, sys), sys)
        results[i, 3] = fixedcurrent(merge(p, Dict(η => v)), lowmap(p, sys), sys)
    end
    results
end

function bifurcate(p::Dict, sys::ODESystem, which, pspan)
    df = ODEFunction(sys, jac=true)
    F(u, p) = df(u, p, 0)
    J(u, p) = df.jac(u, p, 0)

    idx = findfirst(isequal(which), parameters(sys))
    p0 = MT.varmap_to_vars(p, parameters(sys); defaults=sys.defaults)

    opts = ContinuationPar(
        dsmax=5e-2, dsmin=1e-4, ds=1e-3,
        maxSteps=100000,
        pMin=first(pspan), pMax=last(pspan),
        detectBifurcation=3,
        newtonOptions=NewtonPar(tol=1e-9, verbose=false, maxIter=15)
    )
    j = current(sys)
    u0 = merge(qmap(p, sys), smap(p, sys), xmap(p, sys), amap(p, sys))
    f(x) = x[1] => substitute(x[2], p)
    u0 = MT.varmap_to_vars(Dict(f.(collect(u0))), states(sys); defaults=sys.defaults)
    operator = DeflationOperator(2, dot, 1.0, [zeros(length(states(sys))), u0])
    p0[idx] = first(pspan)
    continuation(
        F, J, p0, (@lens _[idx]), opts, operator;
        plot=false,
        perturbSolution = (x, p, id) -> (x .+ 0.8 .* rand(length(x)))
    )
end
