# Optimal design of carbon capture processes under part-load operating conditions
## Overview

This project models and optimizes the design of a post-combustion CO₂ capture process using monoethanolamine (MEA) as a solvent. The model is formulated as a nonlinear program in Julia/JuMP and solved with Ipopt.

The design optimization accounts for multiple flue gas conditions (part-load operation) by including multiple operating scenarios in a single stochastic optimization. Two case studies are provided: a natural gas combined cycle (NGCC) power plant and a bituminous coal power plant.

The model description and methodology can be found in:
>Shu DY, Huang B, Kim Y, Gandhi R, Field R, and Shin S. *Optimal design of carbon capture processes under part-load operating conditions.* 2026.

---

## Contents

- [Reproducing results via Makefile](#reproducing-results-makefile)
- [Quick start](#quick-start-manual)
- [Case studies](#case-studies)
- [Sensitivity analysis](#sensitivity-analysis)

---

## Reproducing results via Makefile

To easily reproduce the results from the paper, you can use the provided `Makefile`. This will automatically handle package instantiation and run the main case studies and sensitivity sweeps.

Open your terminal and run:

```bash
git clone git@github.com:mit-shin-group/2026_part_load_design_CCS.git
cd 2026_part_load_design_CCS
make
```

*Note: We recommend using the **HSL MA57** linear solver. Alternatively, Ipopt can use MUMPS, which may be slower.*

---

## Quick start

If you prefer to run the code interactively and explore the models, follow these steps:

**1. Clone the repository and navigate to the project directory:**
```bash
git clone git@github.com:mit-shin-group/2026_part_load_design_CCS.git
cd 2026_part_load_design_CCS
```

**2. Start Julia in the project environment and install dependencies:**
```bash
julia --project=.
```
Inside the Julia REPL, run:
```julia
julia> using Pkg
julia> Pkg.instantiate()
julia> Pkg.develop(path = "ExaProcessModels")
```

**3. Run the case study:**
```julia
julia> include("main.jl")
```

Results are written to a `results/`.

---

## Case studies

In `main.jl`, the case study can be adjusted using `cfg = case_config(...)` to switch between natural gas (`"gas"`) and coal (`"coal"`).

The script solves four design variants in dependency order:

| Variant | Design | Capture Target |
|---|---|---|
| `multi` | Free, multi-load | Per scenario |
| `multi_avg` | Free, multi-load | Weighted average |
| `nominal` | Fixed at nominal point, evaluated under part load | Per scenario |
| `nominal_avg` | Fixed at nominal point, evaluated under part load | Weighted average |


### Automated solution strategy

Solving the full 30-segment, multi-scenario problem directly is unreliable. Instead, we solve the model by successively ramping up the model complexity. This is handled automatically by `solve_with_increasing_complexity`:

1. Starts with a coarse segment resolution (as specified by the user), single-scenario model.
2. Successively increases the column resolution to 30 segments.
3. Adds part-load scenarios one by one, seeding each solve from the previous.
4. If a solve fails, a rescue ladder is applied automatically:
   - Retry with an adaptive Ipopt barrier strategy.
   - Resolve the last converged solution, and use the resolved model as seed.
   - Retry from earlier, coarser solutions.
   - Insert an easier intermediate step.
   - Skip a non-final step.

When seeding from a previously solved variant (e.g., using the multi-load solution to start the nominal solve), the full ramp can be skipped and a single seeded solve at full complexity is attempted instead to save some time.

### Multistart

An optional multistart polish can be enabled to improve the final design:

```julia
base_solutions = solve_base_cases(cfg; multistart = X)
```
Here, `X` is the number of multistarts performed.

Each trial re-solves the full model with the design variable start values perturbed by lognormal noise. The solution with the lowest overall objective function value (Total Annualized Cost) is kept.

---

## Sensitivity Analysis

`main.jl` runs four parameter sweeps after the base solutions are available:

| Sweep | Values |
|---|---|
| `capture_rate` | 70, 75, 80, 85, 90, 92, **95**, 98, 99 % |
| `mass_flow` | −20, −10, **0**, +10, +20 % |
| `full load hours` | −20, −10, **0**, +10, +20 % |
| `capital charge factor` | −40, −20, **0**, +20, +40 % |

*Bold values indicate the base point assumed in the case study.* 

Within each sweep, previous solutions are used as initial values to reduce solution time. 

*Note: For the full load hours sweep, the nominal plant design is independent of operating hours and is therefore solved once (at the base point) and reused across all sensitivity values.*

### Outputs
Results for each point are saved individually as:
- `date_case_sensitivities_summary.xlsx` (summarized results)
- `date_case_sensitivitity_variables.xlsx` (all optimal variable values)

A combined summary table tracking the cost of CO₂ avoided and total annualized cost vs. sensitivity values for all cases is saved as:
- `date_sensitivity_summary.xlsx`
