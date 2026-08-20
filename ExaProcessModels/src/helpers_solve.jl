###############################################################################
# Automatic segment and scenario ramp-up for the CCS design problem.
#   example: 
#   res = solve_with_increasing_complexity(cases = 1:6, target_segments = 30, ...)
#   res.model     # final solved JuMP model -> use it to seed other variations
#
# The number of segments and scenarios (N_segments, cases) is successively 
# increased. Each optimization problem is solved using the solution of the
# previous problem as a seed. If a probelm fails, a rescue ladder is applied:
#   (1) retry with a more robust barrier strategy
#   (2) polish the incumbent, retry from the polished point
#   (3) retry from older / coarser seeds
#   (4) insert an easier intermediate target (smaller step in N, or the new
#       scenario first at lower resolution)
#   (5) skip a non-final target (the next one then adds two scenarios at once)
###############################################################################

const OK_STATUS = ("Solve_Succeeded", "Solved_To_Acceptable_Level", "Feasible_Point_Found")

# design variables of the capture plant
const DESIGN_VARS = (:abs_H, :abs_D, :des_H, :des_D, :NC,
                     :pump_maxvflow_lean, :pump_maxvflow_rich,
                     :A_reboiler_max, :A_condenser_max,
                     :A_intercooler_max, :A_cooler_max, :comp_P_max)

const BASE_POINT   = (capture_rate = 0.95, sensitivity_flow = 1.0,
                      sensitivity_flh = 1.0, sensitivity_CCF = 1.0)

"""
    case_config(; base, multi_cases, nom_case, target_N = 30)

Bundles the case-study settings that are defined in `main.jl`:
`base` (kwargs forwarded to `full_process_model`), the part-load scenarios,
the nominal scenario, and the target column resolution.
"""
case_config(; base, multi_cases, nom_case, target_N::Int = 30) =
    (base = base, multi_cases = collect(Int, multi_cases),
     nom_case = collect(Int, nom_case), target_N = target_N)

# containers
struct Step
    model::JuMP.Model
    N::Int
    cases::Vector{Int}
    obj::Float64
    time::Float64
    label::String
end

struct ModelResult
    model::JuMP.Model
    N::Int
    cases::Vector{Int}
    objective::Float64
    converged::Bool
    trace::Vector{Step}
end

# polymorphic getter for convenience
model_of(x::JuMP.Model)          = x
model_of(x::Step)                = x.model
model_of(x::ModelResult)         = x.model
model_of(::Nothing)              = nothing

solved(m::JuMP.Model)  = try raw_status(m) in OK_STATUS catch; false end
obj_or_nan(m::JuMP.Model) = try objective_value(m) catch; NaN end

# model factory to create model variations
function make_factory(kw, attrs)
    """
        make_factory(kw, attrs) -> (N, cases) -> JuMP.Model

    `kw`    : keyword arguments forwarded to `full_process_model` (case_study, capture_rate,
            sensitivity_*, fix_design, design, ...)
    `attrs` : Ipopt attributes applied to every model that is built.
    """
    return function (N::Int, cases)
        model, _, _, _, _ = full_process_model(; N_abs = N, N_des = N,
                                                cases = collect(Int, cases), kw...)
        for (k, v) in attrs
            set_optimizer_attribute(model, k, v)
        end
        return model
    end
end

# one build-seed-solve
function attempt(factory, N, cases, seed;
                 label = "ramp", attrs = Dict{String,Any}(), hook = nothing, verbose = true)
    """
        attempt(factory, N, cases, seed; ...) -> Step or nothing

    Builds the model, transfers start values from `seed` (Model / Step / ModelResult /
    nothing), optionally applies `hook(model)` (used for multistart perturbation), solves.
    Returns `nothing` if the solve did not succeed.
    """
    t0    = time()
    model = factory(N, cases)
    for (k, v) in attrs
        set_optimizer_attribute(model, k, v)
    end

    src = model_of(seed)
    if src !== nothing && solved(src)
        try
            transfer_solution_same_variables(src, model; verbose = false)
        catch err
            @warn "start-value transfer failed" err
        end
    end
    hook === nothing || hook(model)

    ok = try
        optimize!(model)
        solved(model)
    catch err
        @warn "solver raised an exception" err
        false
    end

    dt = time() - t0
    verbose && @printf("  N=%-3d S=%-2d %-36s %-7s %7.1f s\n",
                       N, length(cases), label, ok ? "ok" : "FAILED", dt)
    return ok ? Step(model, N, collect(Int, cases), obj_or_nan(model), dt, label) : nothing
end

# return ramp complexity plan in consistent format
function default_plan(cases, target_N::Int; start_N::Int = 10, dN::Int = 4)
    """
        default_plan(cases, target_N; start_N, dN)

    Stage 1: refine the column discretisation on the first scenario only.
    Stage 2: add scenarios one by one at the target resolution.
    """
    cases = collect(Int, cases)
    Ns = start_N >= target_N ? [target_N] :
         unique(push!(collect(start_N:dN:target_N), target_N))
    plan = [(N, cases[1:1]) for N in Ns]
    append!(plan, [(target_N, cases[1:k]) for k in 2:length(cases)])
    return plan
end

# function to run the ramp plan, applies rescue methods to recover is optimization fails
function run_plan(factory, plan; seed = nothing, min_N::Int = 8,
                  rescue_budget::Int = 25, verbose = true)
    queue    = Tuple{Int,Vector{Int}}[(N, collect(Int, c)) for (N, c) in plan]
    trace    = Step[]
    inserted = Set{Tuple{Int,Vector{Int}}}()
    budget   = rescue_budget

    while !isempty(queue)
        N, cases = popfirst!(queue)
        current  = isempty(trace) ? seed : trace[end]
        step     = attempt(factory, N, cases, current; verbose = verbose)

        if step === nothing
            if budget <= 0
                @warn "rescue budget exhausted — stopping the ramp"
                break
            end
            budget -= 1
            step = rescue!(queue, factory, trace, (N, cases), seed, inserted, min_N;
                           verbose = verbose)
        end
        step === nothing || push!(trace, step)
    end
    return trace
end

function rescue!(queue, factory, trace, target, seed, inserted, min_N; verbose = true)
    """
        rescue!(queue, factory, trace, target, seed, inserted, min_N)

    Returns a solved `Step` for `target`, or `nothing` after making the ramp easier
    (inserting intermediate targets into `queue`, or skipping a non-final target).
    """
    N, cases = target
    verbose && println("  -> rescue for N=$N, S=$(length(cases))")
    current = isempty(trace) ? seed : trace[end]

    # (1) same model, different IPOPT settings (barrier update / more iterations)
    step = attempt(factory, N, cases, current; label = "retry: adaptive barrier",
                   attrs = Dict("mu_strategy" => "adaptive", "max_iter" => 900),
                   verbose = verbose)
    step === nothing || return step

    if !isempty(trace)
        inc = trace[end]
        # (2) polish the incumbent, then retry from the polished point
        pol = attempt(factory, inc.N, inc.cases, inc; label = "polish incumbent",
                      verbose = verbose)
        if pol !== nothing
            trace[end] = pol
            step = attempt(factory, N, cases, pol; label = "retry from polished",
                           verbose = verbose)
            step === nothing || return step
        end
        # (3) older, coarser seeds
        for older in Iterators.reverse(trace[1:end-1])
            step = attempt(factory, N, cases, older;
                           label = "retry from seed N=$(older.N),S=$(length(older.cases))",
                           verbose = verbose)
            step === nothing || return step
        end
    end

    # (4) take a smaller step
    incN = isempty(trace) ? min_N : trace[end].N
    incS = isempty(trace) ? 0     : length(trace[end].cases)
    mid = if length(cases) > incS && N > min_N
              (max(min_N, N ÷ 2), cases)            # new scenario first at lower resolution
          elseif N - incN > 1
              (incN + (N - incN) ÷ 2, cases)        # halve the resolution step
          else
              nothing
          end
    if mid !== nothing && !(mid in inserted)
        push!(inserted, mid)
        prepend!(queue, [mid, target])
        verbose && println("  -> inserting intermediate target N=$(mid[1]), S=$(length(mid[2]))")
        return nothing
    end

    # (5) skip a non-final target: the next one then adds two scenarios / a larger N step
    if !isempty(queue)
        verbose && println("  -> skipping N=$N, S=$(length(cases))")
        return nothing
    end
    verbose && println("  -> no rescue left for the final target")
    return nothing
end

# multistart clamp variables
function clamp_start(x, v; eps = 1e-9)
    has_lower_bound(v) && (x = max(x, lower_bound(v) + eps))
    has_upper_bound(v) && (x = min(x, upper_bound(v) - eps))
    return x
end

"Multiply the start values of the design variables by lognormal noise."
function perturb_design!(model, rng, spread)
    for sym in DESIGN_VARS
        v = get(model.obj_dict, sym, nothing)
        v isa VariableRef || continue
        is_fixed(v) && continue
        x0 = start_value(v)
        x0 === nothing && continue
        set_start_value(v, clamp_start(x0 * exp(spread * randn(rng)), v))
    end
end

design_of(m::JuMP.Model) = Dict{Symbol,Float64}(
    s => value(m[s]) for s in DESIGN_VARS if haskey(m.obj_dict, s))
design_of(r::ModelResult) = design_of(r.model)

"Multistart at full complexity: re-solve from perturbed design start values, keep the best.
The hook function modifies the model before solving: The seeds are updated for the multistart."
function multistart_polish(factory, best::Step; trials = 8, spread = 0.15, time_limit = 300.0,
                           rng = MersenneTwister(42), verbose = true)
    champion = best
    verbose && println("-- multistart: $trials perturbed design starts --")
    for k in 1:trials
        s = attempt(factory, best.N, best.cases, best;
                    label = "multistart $k/$trials",
                    attrs = Dict("max_cpu_time" => time_limit),
                    hook  = m -> perturb_design!(m, rng, spread), verbose = verbose)
        if s !== nothing && s.obj < champion.obj * (1 - 1e-8)
            verbose && @printf("     improved: %.6g -> %.6g\n", champion.obj, s.obj)
            champion = s
        end
    end
    return champion
end

# main function to solve model successively increasing the number of segments and scenarios
function solve_with_increasing_complexity(; cases,
        target_segments::Int = 30,
        start_segments::Int  = 10,
        segment_step::Int    = 4,
        plan                 = nothing,
        seed_model           = nothing,
        final_overrides      = NamedTuple(),
        polish::Bool         = true,
        multistart::Int      = 0,
        multistart_spread    = 0.20,
        multistart_time_limit= 300.0,
        rng                  = MersenneTwister(20250817),
        solver_attributes    = Dict{String,Any}("print_level" => 5),
        verbose::Bool        = true,
        model_kwargs...)
    """
        solve_with_increasing_complexity(; cases, target_segments = 30, ...) -> ModelResult

    Automatically ramps up model complexity until `(target_segments, cases)` is solved.

    Key options
    cases              : target scenario list, e.g. 1:6 (gas) or [7] (nominal/deterministic)
    target_segments    : final column resolution (default 30)
    start_segments/dN  : coarse start and step of the resolution ramp
    plan               : override the ramp completely, e.g. [(30, collect(1:6))] for one shot
    seed_model         : Model / Step / ModelResult from a related case study
    final_overrides    : kwargs switched on only at full complexity,
                        e.g. (use_average_capture_rate = true,)
    multistart         : number of perturbed design restarts in the final polishing step
    model_kwargs...    : everything else forwarded to `full_process_model`
    """

    cases   = collect(Int, cases)
    base_kw = values(model_kwargs)
    factory = make_factory(base_kw, solver_attributes)
    plan    = plan === nothing ?
              default_plan(cases, target_segments; start_N = start_segments, dN = segment_step) :
              plan

    verbose && println("== complexity ramp: $(length(plan)) planned steps ==")
    trace = run_plan(factory, plan; seed = seed_model, min_N = start_segments, verbose = verbose)
    isempty(trace) && error("continuation failed on the first model, check seed / model kwargs")

    # switch the model variant at full complexity (e.g. weighted-average capture rate)
    if !isempty(final_overrides)
        verbose && println("== final variant: $final_overrides ==")
        factory = make_factory(merge(base_kw, final_overrides), solver_attributes)
        append!(trace, run_plan(factory, [(trace[end].N, trace[end].cases)];
                                seed = trace[end], min_N = start_segments, verbose = verbose))
    end
    best = trace[end]

    if polish
        p = attempt(factory, best.N, best.cases, best; label = "final polish", verbose = verbose)
        if p !== nothing && (isnan(best.obj) || p.obj <= best.obj * (1 + 1e-6))
            push!(trace, p); best = p
        end
    end

    # perform multistart with the final model by varying the initial values for the design
    if multistart > 0
        champ = multistart_polish(factory, best; trials = multistart,
                                  spread = multistart_spread, time_limit = multistart_time_limit,
                                  rng = rng, verbose = verbose)
        champ === best || (push!(trace, champ); best = champ)
    end

    reached_target = best.N == target_segments && sort(best.cases) == sort(cases)   # was: converged
    reached_target || @warn "target complexity not reached" reached = (best.N, best.cases) requested = (target_segments, cases)

    res = ModelResult(best.model, best.N, best.cases, best.obj, reached_target, trace)
    verbose && print_trace(res)
    return res
end

function print_trace(res::ModelResult)
    println("\n== continuation trace ==")
    for (i, s) in enumerate(res.trace)
        @printf("%3d  N=%-3d S=%-2d %-34s TAC=%.6g  (%.1f s)\n",
                i, s.N, length(s.cases), s.label, s.obj, s.time)
    end
    @printf("final: N=%d  cases=%s  TAC=%.6g  converged=%s  status=%s\n",
            res.N, string(res.cases), res.objective, res.converged, raw_status(res.model))
end

###################################################################################
# Sensitivity analysis functions
# sweep definitions: values reported, base value, mapping to model kwargs
const SWEEPS = (
    capture_rate = (values = [70, 75, 80, 85, 90, 92, 95, 98, 99], base = 95,
                    kw = v -> (capture_rate = v / 100,)),
    mass_flow    = (values = [-20, -10, 0, 10, 20], base = 0,
                    kw = v -> (sensitivity_flow = 1 + v / 100,)),
    flh          = (values = [-20, -10, 0, 10, 20], base = 0,
                    kw = v -> (sensitivity_flh  = 1 + v / 100,)),
    CCF          = (values = [-40, -20, 0, 20, 40], base = 0,
                    kw = v -> (sensitivity_CCF  = 1 + v / 100,)),
)

# save the case solution
"""
A solved case variant. `eval` is the result that is reported (part-load
evaluation for the nominal variants); `design` is the free-design solve that
produced the fixed design, or `nothing` for the multi-load variants.
"""
struct CaseSolution
    eval::ModelResult
    design::Union{Nothing,ModelResult}
end

model_of(x::CaseSolution)   = model_of(x.eval)      # so it can be used as a seed
converged(x::CaseSolution)  = x.eval.converged
coa_of(x::CaseSolution)     = try value(x.eval.model[:coa]) catch; NaN end

plan_for(prev, cases, target_N) =           # ramp only if there is no seed
    prev === nothing ? nothing : [(target_N, collect(Int, cases))]

"Multi-load design: all scenarios, design free."
function solve_multi(cfg, p, prev; multistart = 0, multistart_spread = 0.2, multistart_timelimit = 300.0, average = false)
    r = solve_with_increasing_complexity(;
            cases           = cfg.multi_cases,
            target_segments = cfg.target_N,
            plan            = plan_for(prev, cfg.multi_cases, cfg.target_N),
            seed_model      = prev === nothing ? nothing : prev.eval,
            multistart      = multistart,
            multistart_spread = multistart_spread,
            multistart_time_limit= multistart_timelimit,
            final_overrides = average ? (use_average_capture_rate = true,) : NamedTuple(),
            cfg.base..., p...)
    return CaseSolution(r, nothing)
end

function solve_nominal(cfg, p, prev; multistart = 0, multistart_spread = 0.20,
        multistart_time_limit= 300.0, average = false, design = nothing)
    """
    Nominal design + part-load evaluation with that design fixed.

    `design`: reuse a previously computed design instead of re-solving it.
    """
    design_seed = prev === nothing ? nothing :
                  (prev.design === nothing ? prev.eval : prev.design)

    nom = design === nothing ?
        solve_with_increasing_complexity(;
            cases           = cfg.nom_case,
            target_segments = cfg.target_N,
            plan            = plan_for(design_seed, cfg.nom_case, cfg.target_N),
            seed_model      = design_seed,
            multistart      = multistart,
            multistart_spread = multistart_spread,
            multistart_time_limit= multistart_time_limit,
            cfg.base..., merge(p, (sensitivity_flh = 1.0,))...) :  # design ignores FLH scaling
        nothing

    d = nom === nothing ? design : design_of(nom)

    pl = solve_with_increasing_complexity(;
            cases           = cfg.multi_cases,
            target_segments = cfg.target_N,
            plan            = plan_for(prev, cfg.multi_cases, cfg.target_N),
            seed_model      = prev === nothing ? nothing : prev.eval,
            fix_design      = true, design = d,
            final_overrides = average ? (use_average_capture_rate = true,) : NamedTuple(),
            cfg.base..., p...)

    kept = nom === nothing ? (prev === nothing ? nothing : prev.design) : nom
    return CaseSolution(pl, kept)
end

const CASE_VARIANTS = (
    multi       = (cfg, p, prev; kw...) -> solve_multi(cfg, p, prev; average = false, kw...),
    multi_avg   = (cfg, p, prev; kw...) -> solve_multi(cfg, p, prev; average = true,  kw...),
    nominal     = (cfg, p, prev; kw...) -> solve_nominal(cfg, p, prev; average = false, kw...),
    nominal_avg = (cfg, p, prev; kw...) -> solve_nominal(cfg, p, prev; average = true,  kw...),
)

is_nominal(case::Symbol) = occursin("nominal", String(case))

function solve_base_cases(cfg; multistart::Int = 0, start_segments = 8)
    """
    Solve the four variants at the base point, in dependency order:
    multi (full ramp) -> multi_avg -> nominal -> nominal_avg.
    These anchor every sensitivity sweep.
    """
    println("\n############ base point: $BASE_POINT ############")

    println("\n--- multi-load ---")
    multi = CaseSolution(
        solve_with_increasing_complexity(; cases = cfg.multi_cases, target_segments = cfg.target_N,
                                           start_segments = start_segments, segment_step = 4,
                                           multistart = multistart, cfg.base..., BASE_POINT...),
        nothing)

    println("\n--- multi-load, average capture ---")
    multi_avg = solve_multi(cfg, BASE_POINT, multi, multistart = multistart; average = true)

    println("\n--- nominal ---")
    nominal = solve_nominal(cfg, BASE_POINT, multi, multistart = multistart)

    # with a single scenario the average and per-scenario targets coincide,
    # so the nominal design is reused instead of re-solved
    println("\n--- nominal, average capture ---")
    nominal_avg = solve_nominal(cfg, BASE_POINT, nominal; average = true,
                                design = design_of(nominal.design))

    base = Dict(:multi => multi, :multi_avg => multi_avg,
                :nominal => nominal, :nominal_avg => nominal_avg)

    return base
end

# sensitivity sweeps
"Values above and below `base`, each ordered by increasing distance from it."
function outward_branches(values, base)
    v = sort(collect(values))
    i = findfirst(==(base), v)
    i === nothing && error("base value $base is not contained in $values")
    return (v[i+1:end], reverse(v[1:i-1]))
end

function run_sweep(name::Symbol, cfg, base_solutions; label = "sweep")
    """
    Sweeps one parameter for all four case variants.

    We use each case to seed the next (multi -> multi_avg -> nominal -> nominal_avg).
    Across values, each case is seeded by the nearest converged solve in the
    same direction (outward from base).
    """
    spec = getfield(SWEEPS, name)

    # initialise output with base solutions
    out = Dict{Symbol, Dict{Any, CaseSolution}}(
              case => Dict{Any, CaseSolution}(spec.base => base_solutions[case])
              for case in keys(CASE_VARIANTS))

    for branch in outward_branches(spec.values, spec.base)
        # reset per-case seeds to base at the start of each direction
        prev_per_case = Dict{Symbol, CaseSolution}(pairs(base_solutions)...)

        for v in branch
            println("\n=== $name = $v ===")
            p              = merge(BASE_POINT, spec.kw(v))
            solved_at_v    = Dict{Symbol, CaseSolution}()
            last_converged = nothing   # best seed from cases already solved at this v

            for (case, solve_fn) in pairs(CASE_VARIANTS)
                println("\n=== $name = $v, case = $case ===")
                # cross-case seed at the same v (e.g. multi seeds multi_avg);
                # fall back to the same-case nearest converged v if not yet available
                seed = last_converged !== nothing ? last_converged : prev_per_case[case]

                # design reuse:
                #   FLH sweep: nominal design is independent of FLH -> reuse base design
                #   nominal_avg: reuse nominal[v]'s design if already solved at this v
                reuse = if name === :flh && is_nominal(case)
                    design_of(base_solutions[case].design)
                elseif case === :nominal_avg && haskey(solved_at_v, :nominal)
                    design_of(solved_at_v[:nominal].design)
                else
                    nothing
                end

                sol = try  # just in case solve fails:
                    reuse === nothing ? solve_fn(cfg, p, seed) :
                                        solve_fn(cfg, p, seed; design = reuse)
                catch err
                    @warn "solve failed for $case | $name = $v" err
                    nothing
                end

                if sol !== nothing
                    out[case][v] = sol
                    if converged(sol)
                        prev_per_case[case] = sol   # seed same case at next v
                        solved_at_v[case]   = sol   # seed next case at same v
                        last_converged      = sol
                    end
                end
            end
        end
    end
    return out
end

# summary
"Summarize in one DataFrame: cost of CO2 avoided vs. sweep value)."
function sweep_summary(sweeps::Dict{Symbol,<:Any})
    rows = NamedTuple[]
    for (name, sweep) in sweeps, (case, res) in sweep
        for (v, sol) in sort(collect(res); by = first)
            push!(rows, (sweep = name, case = case, value = v,
                         coa = coa_of(sol), tac = sol.eval.objective,
                         converged = converged(sol)))
        end
    end
    return DataFrame(rows)
end

function save_summary(df, label = "sensitivity_summary"; results_dir = "results")
    mkpath(results_dir)
    path = joinpath(results_dir, "$(label).xlsx")
    XLSX.writetable(path, df; overwrite = true, sheetname = "summary")
    println("saved: $path")
    return path
end
