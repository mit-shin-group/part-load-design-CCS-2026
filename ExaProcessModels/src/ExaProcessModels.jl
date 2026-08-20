module ExaProcessModels

# Bring in external dependencies used anywhere in the package
using JuMP, Ipopt, DataFrames, Random, Printf, MathOptInterface, HSL, HSL_jll, XLSX

import MathOptInterface as MOI

# Include files of the module
include("cost_models.jl")
include("helpers.jl")
include("helpers_casestudy.jl")
include("full_process_industry.jl")
include("helpers_solve.jl")

# Export the functions for users
export  full_process_model,
        transfer_solution_same_variables,
        export_variable_values_split,
        print_flash,
        compute_optimal_jacobian,
        build_scenario_permutations,
        group_ranges_by_scenario,
        scenario_block_nnz,
        count_dof,
        collect_results,
        collect_results_all_scenarios,
        solve_with_increasing_complexity,
        solve_base_cases,
        run_sweep,
        sweep_summary,
        save_summary,
        case_config, 
        design_of,
        ModelResult,
        species_index,
        CaseSolution,
        BASE_POINT, 
        SWEEPS,
        converged
end # module
