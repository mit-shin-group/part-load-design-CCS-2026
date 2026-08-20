# Run the case study cases and sensitivity analysis.

using ExaProcessModels
using Dates
include("save_results.jl")

# coal case study
# case-study settings
cfg = case_config(
    base = (recycle = true, as_optimization = true, case_study = "coal",
            abs_P = 110.0, des_P = 200.0,
            intercooler_active = true, export_log = false,
            sensitivity_price_reduction = 0.0),
    multi_cases = 1:7,      # gas: 1:6, coal: 1:7, gas intermediate: 8:13
    nom_case    = [8],      # gas: [7], coal: [8], gas intermediate: [14]
    target_N    = 30,
)

label = "$(Dates.today())_$(cfg.base.case_study)"

# 1) base point for all four variants (successively ramps up model size for first problem)
base_solutions = solve_base_cases(cfg; multistart = 20, start_segments = 8)
for (case, sol) in base_solutions
    save_results(sol.eval, "$(label)_base_$(case)"; abs_P = cfg.base.abs_P, des_P = cfg.base.des_P,
                    sensitivity_CCF = BASE_POINT.sensitivity_CCF,
                    add_column_properties = true)
end

# 2) sensitivity studies
sweeps = Dict(name => run_sweep(name, cfg, base_solutions; label = label)
              for name in keys(SWEEPS))

for name in keys(sweeps) #capture_rate, mass_flow, flh, CCF
    specs = getfield(SWEEPS, name)
    for case in keys(sweeps[name])
        for v in keys(sweeps[name][case])
            v == specs.base && continue # if v = base value, skip as base already saved above
            sol = sweeps[name][case][v]
            converged(sol) || continue # skip failed solves
            p = merge(BASE_POINT, specs.kw(v)) # recover the parameters used for this point
            save_results(sol.eval, "$(label)_$(case)_$(name)_$(v)";
                         abs_P           = cfg.base.abs_P,
                         des_P           = cfg.base.des_P,
                         sensitivity_CCF = p.sensitivity_CCF)
        end
    end
end

# 3) collected summary
summary = sweep_summary(sweeps)
save_summary(summary, "$(label)_sensitivity_summary")

####################################################################
# gas case study
# case-study settings
cfg = case_config(
    base = (recycle = true, as_optimization = true, case_study = "gas",
            abs_P = 110.0, des_P = 200.0,
            intercooler_active = true, export_log = false,
            sensitivity_price_reduction = 0.0),
    multi_cases = 1:6,      # gas: 1:6, coal: 1:7, gas intermediate: 8:13
    nom_case    = [7],      # gas: [7], coal: [8], gas intermediate: [14]
    target_N    = 30,
)

label = "$(Dates.today())_$(cfg.base.case_study)"

# 1) base point for all four variants (successively ramps up model size for first problem)
base_solutions = solve_base_cases(cfg; multistart = 20, start_segments = 8)
for (case, sol) in base_solutions
    save_results(sol.eval, "$(label)_base_$(case)"; abs_P = cfg.base.abs_P, des_P = cfg.base.des_P,
                    sensitivity_CCF = BASE_POINT.sensitivity_CCF,
                    add_column_properties = true)
end

# 2) sensitivity studies
sweeps = Dict(name => run_sweep(name, cfg, base_solutions; label = label)
              for name in keys(SWEEPS))

for name in keys(sweeps) #capture_rate, mass_flow, flh, CCF
    specs = getfield(SWEEPS, name)
    for case in keys(sweeps[name])
        for v in keys(sweeps[name][case])
            v == specs.base && continue # if v = base value, skip as base already saved above
            sol = sweeps[name][case][v]
            converged(sol) || continue # skip failed solves
            p = merge(BASE_POINT, specs.kw(v)) # recover the parameters used for this point
            save_results(sol.eval, "$(label)_$(case)_$(name)_$(v)";
                         abs_P           = cfg.base.abs_P,
                         des_P           = cfg.base.des_P,
                         sensitivity_CCF = p.sensitivity_CCF)
        end
    end
end

# 3) collected summary
summary = sweep_summary(sweeps)
save_summary(summary, "$(label)_sensitivity_summary")
