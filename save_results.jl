###############################################################################
# Result export: save_results saves a summary of the results
# export variables_to_excel saves all variables
#   save_results(res, label; abs_P, des_P, ...) 
#   export_variables_to_excel(model, filename)
###############################################################################
using XLSX, DataFrames, JuMP
using ExaProcessModels: ModelResult, collect_results_all_scenarios, species_index

# Excel-safe accessors
"Make a value writable by XLSX: `nothing`, `NaN` and `Inf` become `missing`."
excel_safe(x)                = x
excel_safe(::Nothing)        = missing
excel_safe(x::AbstractFloat) = isfinite(x) ? x : missing

"`try_get(f, v)` -> `f(v)`, or `missing` if it is unavailable (unsolved model, no bound, ...)."
try_get(f, v) = try excel_safe(f(v)) catch; missing end

# name classification, current component types
const VAR_CATEGORIES = ("abs", "HX", "des", "reb", "cond", "cool")

"abs_L[1,2] -> (abs_L, [1, 2])"
function parse_var_name(varname::AbstractString)
    i = findfirst(isequal('['), varname)
    i === nothing && return String(varname), Int[]

    base = String(varname[1:prevind(varname, i)])
    idx  = tryparse.(Int, strip.(split(strip(varname[i:end], ['[', ']']), ",")))
    return base, any(isnothing, idx) ? Int[] : Vector{Int}(idx)
end

"Unit a variable belongs to, taken from its name prefix (other if no prefix)."
function classify_variable(base_name::AbstractString)
    i = findfirst(p -> startswith(base_name, p), VAR_CATEGORIES)
    return i === nothing ? "other" : VAR_CATEGORIES[i]
end

# variable dump
idx_or_missing(idx, k) = k <= length(idx) ? idx[k] : missing

"One row per scalar variable: name, indices, value, start value, bounds, flags."
function variable_table(model::JuMP.Model)
    rows = NamedTuple[]
    for v in all_variables(model)
        nm        = String(name(v))
        base, idx = parse_var_name(nm)
        push!(rows, (
            variable_name = nm,
            base_name     = base,
            category      = classify_variable(base),
            ndims         = length(idx),
            idx1          = idx_or_missing(idx, 1),
            idx2          = idx_or_missing(idx, 2),
            idx3          = idx_or_missing(idx, 3),
            idx4          = idx_or_missing(idx, 4),
            indices       = join(idx, ","),
            value         = try_get(value, v),
            start_value   = try_get(start_value, v),
            fixed_value   = is_fixed(v) ? try_get(fix_value, v) : missing,
            lower_bound   = try_get(lower_bound, v),
            upper_bound   = try_get(upper_bound, v),
            is_fixed      = is_fixed(v),
            is_binary     = is_binary(v),
            is_integer    = is_integer(v),
        ))
    end
    df = DataFrame(rows)
    sort!(df, [:category, :base_name, :idx1, :idx2, :idx3, :idx4])
    return df
end

add_sheet!(xf, sheetname, df) =
    XLSX.writetable!(XLSX.addsheet!(xf, sheetname), collect(eachcol(df)), names(df))

function export_variables_to_excel(model::JuMP.Model, filename::AbstractString;
                                  separate_sheets::Bool = true,
                                  include_summary::Bool = true)
    """
        export_variables_to_excel(model, filename; separate_sheets = true, include_summary = true)

    Write all variables of `model` to `filename`: one sheet with everything, optionally
    one sheet per unit category, and optionally two small summary sheets.
    """
    df = variable_table(model)

    XLSX.openxlsx(filename, mode = "w") do xf
        sheet = xf[1]                                  # reuse the default sheet
        XLSX.rename!(sheet, "all_variables")
        XLSX.writetable!(sheet, collect(eachcol(df)), names(df))

        if separate_sheets
            for cat in unique(df.category)
                add_sheet!(xf, cat, df[df.category .== cat, :])
            end
        end

        if include_summary
            add_sheet!(xf, "summary", DataFrame(
                metric = ["num_variables", "num_with_value", "num_fixed"],
                value  = [nrow(df), count(!ismissing, df.value), count(df.is_fixed)]))
            add_sheet!(xf, "summary_by_category",
                       combine(groupby(df, :category), nrow => :n_variables))
        end
    end
    return df
end

# main export
"Fallback pressures for reporting; pressure drop is not modelled, so levels are fixed."
pressure_defaults(abs_P, des_P) = Dict{Symbol,Any}(
    :abs_P_in => abs_P * 1000, :abs_P => abs_P * 1000,   # Pa
    :des_P_in => des_P * 1000, :des_P => des_P * 1000,   # Pa
)

function save_results(res::ExaProcessModels.ModelResult, label::String;
                      abs_P::Float64,
                      des_P::Float64,
                      sensitivity_CCF::Float64    = 1.0,
                      add_column_properties::Bool = false,
                      results_dir::String         = "results")
    """
        save_results(res, label; abs_P, des_P, sensitivity_CCF, add_column_properties, results_dir)

    Write `<label>_summary.xlsx` (scenario results) and `<label>_variables.xlsx`
    (variable dump) for a solved `ModelResult`. Returns the summary DataFrame.
    """
    mkpath(results_dir)

    df = collect_results_all_scenarios(res.model;
        scenarios             = collect(1:length(res.cases)),
        N_abs                 = res.N,
        N_des                 = res.N,
        species_index         = species_index,
        defaults              = pressure_defaults(abs_P, des_P),
        add_column_properties = add_column_properties,
        sensitivity_CCF       = sensitivity_CCF,
    )

    stem = joinpath(results_dir, label)
    XLSX.writetable("$(stem)_summary.xlsx", df; overwrite = true, sheetname = "results")
    export_variables_to_excel(res.model, "$(stem)_variables.xlsx")

    println("saved: $(stem)_summary.xlsx\n       $(stem)_variables.xlsx")
    return df
end
