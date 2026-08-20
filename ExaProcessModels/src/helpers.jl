# This file contains constants and property calculation methods that are used in many units

#####################################################################################
# constants
# index species map
const species_index = Dict(:MEA => 1, :MEA2p => 2, :Hp => 3, :Carbamate => 4, :HCO3m => 5, :OHm => 6, :CO3_2m => 7, :CO2 => 8, :H2O => 9, :N2 => 10, :O2 => 11)

# list of species participating in reactions
const species_react = [:MEA, :MEA2p, :Hp, :Carbamate, :HCO3m, :OHm, :CO3_2m, :CO2, :H2O]

# molecular weights in g/mol
const MW = Dict(
    :MEA => 61.08,
    :MEA2p => 62.08,
    :Hp => 1.008,
    :Carbamate => 121.11,
    :HCO3m => 61.016,
    :OHm => 17.007,
    :CO3_2m => 60.01,
    :CO2 => 44.01,
    :H2O => 18.015,
    :N2 => 28.0134,
    :O2 => 31.9988
)

const density = 1000 # kg/m³, assume constant solvent density for simplicity
const R = 8.314 # J/(mol·K), ideal gas constant
const g = 9.81 # N/m2, standard acceleration gravity

#####################################################################################
# helpers for constraints declaration
function add_mole_balance!(model, x, L, y, V, x_Li, y_Vi, xi, S; species=species_react, N_seg=1)
    """
    Helper to add species mole balances for single- or multi-segment units.
    Arguments:
        model: JuMP model
        x: liquid mole fraction
        L: liquid flow rate
        y: vapor mole fraction
        V: vapor flow rate
        xi: reaction extent 
        S: scenario index set
        x_Li: inlet liquid molar flow rate
        y_Vi: inlet vapor molar flow rate
    Keyword args:
        species: vector of species symbols to include (default: all reactive species)
        N_seg: number of segments
    Returns: constraint container
    """
    n_reactions = size(stoich, 2)
    if N_seg == 1
        cons = @constraint(model, [s in S, c in species],
            x[s, 1, species_index[c]] * L[s, 1] + (y === nothing ? 0 : y[s, 1, species_index[c]] * V[s, 1])
            ==
            (y_Vi === nothing ? 0 : y_Vi[s, species_index[c]]) + (x_Li === nothing ? 0 : x_Li[s, species_index[c]]) +
            sum(stoich[species_index[c], r] * xi[s, 1, r] for r in 1:n_reactions))
        return cons
    else
        cons_1 = @constraint(model, [s in S, c in species],
            x[s, 1, species_index[c]] * L[s, 1] + (y === nothing ? 0 : y[s, 1, species_index[c]] * V[s, 1])
            ==
            (y_Vi === nothing ? 0 : y_Vi[s, species_index[c]]) + x[s, 2, species_index[c]] * L[s, 2] +
            sum(stoich[species_index[c], r] * xi[s, 1, r] for r in 1:n_reactions)) # bottom segment
        cons_N = @constraint(model, [s in S, c in species],
            x[s, N_seg, species_index[c]] * L[s, N_seg] + (y === nothing ? 0 : y[s, N_seg, species_index[c]] * V[s, N_seg])
            ==
            (y === nothing ? 0 : y[s, N_seg-1, species_index[c]] * V[s, N_seg-1]) + (x_Li === nothing ? 0 : x_Li[s, species_index[c]]) +
            sum(stoich[species_index[c], r] * xi[s, N_seg, r] for r in 1:n_reactions)) # top segment
        if N_seg > 2
            cons_n = @constraint(model, [s in S, n in 2:N_seg-1, c in species],
                x[s, n, species_index[c]] * L[s, n] + (y === nothing ? 0 : y[s, n, species_index[c]] * V[s, n])
                ==
                (y === nothing ? 0 : y[s, n-1, species_index[c]] * V[s, n-1]) + x[s, n+1, species_index[c]] * L[s, n+1] +
                sum(stoich[species_index[c], r] * xi[s, n, r] for r in 1:n_reactions))
            return (cons_1, cons_N, cons_n)
        end
        return (cons_1, cons_N)
    end
end

function add_energy_balance!(model, H_liq_in, H_liq_out, H_vap_in, H_vap_out, Q, S; N_seg=1)
    """
    Helper to add energy balances for single- or multi-segment units.
    Arguments:
        model: JuMP model
        H_liq_in: liquid enthalpy in
        H_liq_out: liquid enthalpy out
        H_vap_in: vapor enthalpy in
        H_vap_out: vapor enthalpy out
        Q: heat duty
        S: scenario index set
    Keyword args:
        N_seg: number of segments (default 1)
    Returns: constraint container
    """
    if N_seg == 1
        cons = @constraint(model, [s in S],
            (H_liq_in === nothing ? 0 : H_liq_in[s]) + (H_vap_in === nothing ? 0 : H_vap_in[s]) + (Q === nothing ? 0 : Q[s])
            ==
            H_liq_out[s, 1] + (H_vap_out === nothing ? 0 : H_vap_out[s, 1]))
        return cons
    else
        cons_1 = @constraint(model, [s in S],
            H_liq_out[s, 2] + (H_vap_in === nothing ? 0 : H_vap_in[s]) + (Q === nothing ? 0 : Q[s, 1])
            ==
            H_liq_out[s, 1] + (H_vap_out === nothing ? 0 : H_vap_out[s, 1])) # bottom segment
        cons_N = @constraint(model, [s in S],
            (H_liq_in === nothing ? 0 : H_liq_in[s]) + H_vap_out[s, N_seg-1] + (Q === nothing ? 0 : Q[s, N_seg])
            == H_liq_out[s, N_seg] + H_vap_out[s, N_seg]) # top segment
        if N_seg > 2
            cons_n = @constraint(model, [s in S, n in 2:N_seg-1],
                H_liq_out[s, n+1] + (H_vap_out === nothing ? 0 : H_vap_out[s, n-1]) + (Q === nothing ? 0 : Q[s, n]) 
                == H_liq_out[s, n] + H_vap_out[s, n])
            return (cons_1, cons_N, cons_n)
        end
        return (cons_1, cons_N)
    end
end

function add_VLE!(model, x, y, T, S, N_seg, density, P, y_in, calc_eff=false; 
                  beta=nothing,
                  mu=nothing,
                  sigma=nothing,
                  V_over_L=nothing,
                  T_in=nothing, 
                  alpha=nothing, 
                  x_MEA_composition_in=nothing, 
                  u_v=nothing, 
                  dz=nothing,
                  name=nothing)
    """
    Helper to add VLE constraints and expressions for single- or multi-segment units.
    Arguments:
        model: JuMP model
        x: liquid mole fraction, -
        y: vapor mole fraction, -
        T: temperature array, K
        S: scenario index set
        N_seg: number of segments
        density: liquid density, kg/m3
        P: pressure, atm
        y_in: inlet mole fractions in vapor phase, is set to 0 if nothing
        
        Inputs for approach to equilirbium model
        calc_eff: false for equilibrium assumption, true for approach to equilibrium model to approximate rate limitations
        beta: vector of regression coefficients for rate approximation
        mu: dict of mean values of each input variable for standardization in rate approximation
        sigma: dict of standard deviations of each input variable for standardization in rate approximation				   
        V_over_L: vapor to liquid mass flow ratio for each segment
        T_in: solvent inlet temperature to each segment in K
        alpha: loading for each segment
        x_MEA_composition_in: lean solvent MEA mass fraction (without CO2) kgMEA/(kgMEA+kgH2O)
        u_v: vapor superficial velocity in m/s
        Δz: segment height (required if calc_eff=true, can be set to nothing if not needed for rate approximation)
        
        Returns: tuple (cons_raoult_H2O, cons_Henry_efficiency_CO2_1, cons_Henry_efficiency_CO2_n) where the last one is only returned if N_seg > 1
    """
    # saturation pressure water
    P_H2O_sat_expr = get_Pvap.(T) / 101325 # atm
    # Raoult's law for H2O
    cons_raoult_H2O = @constraint(model, [s in S, n=1:N_seg], P * y[s, n, species_index[:H2O]] == x[s, n, species_index[:H2O]] * P_H2O_sat_expr[s, n])

    # Henry's law for CO2 with efficiency factor
    # apparent mole fractions
    MEA_app, H2O_app, CO2_app = compute_apparent_liquid(model, x, S, N_seg)
    # CO2 free solvent apparent molar fractions
    x_MEA_sf = @expression(model, [s in S, n in 1:N_seg], MEA_app[s,n] / (MEA_app[s,n] + H2O_app[s,n]))
    x_H2O_sf = @expression(model, [s in S, n in 1:N_seg], H2O_app[s,n] / (MEA_app[s,n] + H2O_app[s,n]))
    # Henry's coefficient
    He_CO2_expr = @expression(model, [s in S, n in 1:N_seg], get_H_hartono(x_MEA_sf[s,n], x_H2O_sf[s,n], CO2_app[s,n]/MEA_app[s,n], T[s,n]) / 101325) # atm·m^3/mol
    # average molecular weight
	MW_avg_expr = [sum(x[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_seg]  # kg/mol
    Py_CO2_equilibrium = @expression(model, [s in S, n=1:N_seg], He_CO2_expr[s, n] * (x[s, n, species_index[:CO2]] * density / MW_avg_expr[s, n]))
    
    # approximate rate limitations via approach to equilibrum model
    if calc_eff       
        #standardized
        ln_k = @expression(model, [s in S, n in 1:N_seg],
            beta[1] +
            beta[2] * (V_over_L[s,n] - mu[:V_over_L]) / sigma[:V_over_L] +
            beta[3] * (T_in[s,n] - 273.15-mu[:T_in]) / sigma[:T_in] +
            beta[4] * (P - mu[:P]) / sigma[:P] +
            beta[5] * (alpha[s,n] - mu[:alpha]) / sigma[:alpha] +
            beta[6] * (x_MEA_composition_in[s,n] - mu[:x_MEA_composition_in]) / sigma[:x_MEA_composition_in] +
            beta[7] * (u_v[s,n] - mu[:u_v]) / sigma[:u_v]
            +
            beta[8] * (V_over_L[s,n] * u_v[s,n]-mu[:V_over_L_u_v]) / sigma[:V_over_L_u_v] +
            beta[9] * ((T_in[s,n] - 273.15) * alpha[s,n] - mu[:T_in_alpha]) / sigma[:T_in_alpha] +
            beta[10] * (alpha[s,n] * x_MEA_composition_in[s,n] - mu[:alpha_x_MEA_composition_in]) / sigma[:alpha_x_MEA_composition_in]
        )
        k_CO2 = @expression(model, [s in S, n in 1:N_seg], exp(ln_k[s,n]))
        # segment efficiency η = 1 - exp(-k Δz)
        eff_CO2 = @expression(model, [s in S, n in 1:N_seg], 1 - exp(-k_CO2[s,n] * dz))
    else  # equilibrium assumption
        eff_CO2 = @expression(model, [s in S, n in 1:N_seg], 1.0)
    end

    # save expressions for results
    if name !== nothing && calc_eff
        key = Symbol(name * "_eff_CO2")
        model[key] = eff_CO2
        key = Symbol(name * "_V_over_L_inlet")
        model[key] = V_over_L
        key = Symbol(name * "_alpha_inlet")
        model[key] = alpha
        key = Symbol(name * "_x_MEA_composition_in")
        model[key] = x_MEA_composition_in
        key = Symbol(name * "_u_v_inlet")
        model[key] = u_v
        key = Symbol(name * "_T_inlet")
        model[key] = T_in
    end

    cons_Henry_efficiency_CO2_1 = @constraint(model, [s in S], y[s, 1, species_index[:CO2]] * P ==  (y_in === nothing ? 0 : y_in[s, species_index[:CO2]]) * P - ((y_in === nothing ? 0 : y_in[s, species_index[:CO2]]) * P - Py_CO2_equilibrium[s, 1]) * eff_CO2[s, 1])
    if N_seg > 1
        cons_Henry_efficiency_CO2_n = @constraint(model, [s in S, n=2:N_seg], y[s, n, species_index[:CO2]] * P ==  y[s, n-1, species_index[:CO2]] * P - (y[s, n-1, species_index[:CO2]] * P - Py_CO2_equilibrium[s, n]) * eff_CO2[s, n])
        return (cons_raoult_H2O, cons_Henry_efficiency_CO2_1, cons_Henry_efficiency_CO2_n)
    end
    return (cons_raoult_H2O, cons_Henry_efficiency_CO2_1)
end

function add_reaction_equilibria!(model, log_x, T_arr, S; N_seg=1)
    """
    Helper to add chemical reaction equilibrium constraints for single- or multi-segment units.
    Arguments:
        model: JuMP model
        log_x: log mole fraction of liquid species
        T_arr: temperature array 
        S: scenario index set
    Keyword args:
        N_seg: number of segments
    Returns: tuple of constraint containers (one for each reaction)
    """
    cons1 = @constraint(model, [s in S, n in 1:N_seg], get_lnK(T_arr[s, n])[1] + log_x[s, n, species_index[:MEA2p]] == log_x[s, n, species_index[:MEA]] + log_x[s, n, species_index[:Hp]])
    cons2 = @constraint(model, [s in S, n in 1:N_seg], get_lnK(T_arr[s, n])[2] + log_x[s, n, species_index[:HCO3m]] == log_x[s, n, species_index[:CO3_2m]] + log_x[s, n, species_index[:Hp]])
    cons3 = @constraint(model, [s in S, n in 1:N_seg], get_lnK(T_arr[s, n])[3] == log_x[s, n, species_index[:OHm]] + log_x[s, n, species_index[:Hp]])
    cons4 = @constraint(model, [s in S, n in 1:N_seg], get_lnK(T_arr[s, n])[4] + log_x[s, n, species_index[:Carbamate]] == log_x[s, n, species_index[:MEA]] + log_x[s, n, species_index[:HCO3m]])
    cons5 = @constraint(model, [s in S, n in 1:N_seg], get_lnK(T_arr[s, n])[5] + log_x[s, n, species_index[:CO2]] == log_x[s, n, species_index[:HCO3m]] + log_x[s, n, species_index[:Hp]])
    return (cons1, cons2, cons3, cons4, cons5)
end

function compute_apparent_liquid(model, x, S, N_seg=1)
    """
    Calculates apparent concentrations in liquid phase (overloaded)
    Arguments:
        model: JuMP model
        x: mole fraction or molar flow rate, expected to have two indices for scenario and species (s, i) or three indices for scenario, segment, and species (s, n, i)
        S: scenario index set
        N_seg: number of segments
    """
    mea_idx = [species_index[:MEA], species_index[:MEA2p], species_index[:Carbamate]]
    h2o_idx = [species_index[:H2O], species_index[:Hp], species_index[:OHm], species_index[:HCO3m]]
    co2_idx = [species_index[:CO2], species_index[:HCO3m], species_index[:CO3_2m], species_index[:Carbamate]]

    if ndims(x) == 2
        # x[s, i]
        mea = @expression(model, [s in S],
            sum(x[s, i] for i in mea_idx))
        h2o = @expression(model, [s in S],
            sum(x[s, i] for i in h2o_idx))
        co2 = @expression(model, [s in S],
            sum(x[s, i] for i in co2_idx))
    elseif ndims(x) == 3
        # x[s, n, i]
        mea = @expression(model, [s in S, n in 1:N_seg],
            sum(x[s, n, i] for i in mea_idx))
        h2o = @expression(model, [s in S, n in 1:N_seg],
            sum(x[s, n, i] for i in h2o_idx))
        co2 = @expression(model, [s in S, n in 1:N_seg],
            sum(x[s, n, i] for i in co2_idx))
    else
        throw(ArgumentError(
            "x must have two indices (s, i) or three indices (s, n, i)"
        ))
    end

    return (mea, h2o, co2)
end

function compute_liquid_enthalpy(model, app_Li, T_L_in, x_seg, T_seg, L_seg, S, N_seg)
    """
    Enthalpy calculation for liquid inlet, and the outlet of each segment
    Arguments:
        model: JuMP model
        app_Li: apparent liquid inlet molar flows in mol/s
        T_L_in: inlet temperature in K
        x_seg: mole fractions per segment
        T_seg: temperature in segment in K
        L_seg: molar flow in segment
        S: scenario index set
        N_seg: number of segments
    """
    # heat of absorption, assume constant
    h_absorption_CO2 = 84000 # J/mol
    # calculate apparent molar fractions
    mea_seg, h2o_seg, co2_seg = compute_apparent_liquid(model, x_seg, S, N_seg)
    # calculate liquid inlet enthalpy
    if app_Li !== nothing && T_L_in !== nothing
        H_liq_in = @expression(model, [s in S],
            (app_Li[1][s] * H_MEA_liq(T_L_in[s]) +
             app_Li[2][s] * H_H2O_liq(T_L_in[s]) +
             app_Li[3][s] * (H_CO2_vap(T_L_in[s]) - h_absorption_CO2)))
    else
        H_liq_in = nothing
    end

    H_liq_out = @expression(model, [s in S, n in 1:N_seg],
        (mea_seg[s, n] * L_seg[s, n] * H_MEA_liq(T_seg[s, n]) +
         h2o_seg[s, n] * L_seg[s, n] * H_H2O_liq(T_seg[s, n]) +
         co2_seg[s, n] * L_seg[s, n] * (H_CO2_vap(T_seg[s, n]) - h_absorption_CO2)))
    return (H_liq_in, H_liq_out)
end

# Enthalpy calculation for vapor inlet and segment enthalpy
function compute_vapor_enthalpy(model, y_Vi, T_V_in, y_seg, T_seg, V_seg, S, N_seg)
    """
    Enthalpy calculation for the vapor inlet, and the outlet of each segment
    Arguments:
        model: JuMP model
        y_Vi: vapor inlet molar flows in mol/s
        T_V_in: inlet temperature in K
        y_seg: mole fractions per segment
        T_seg: temperature in segment in K
        V_seg: molar flow in segment
        S: scenario index set
        N_seg: number of segments
    """
    # species indices for vapor phase
    h2o = species_index[:H2O]
    co2 = species_index[:CO2]
    n2 = species_index[:N2]
    o2 = species_index[:O2]

    if y_Vi !== nothing && T_V_in !== nothing
        H_vap_in = @expression(model, [s in S],
            y_Vi[s, h2o] * H_H2O_vap(T_V_in[s]) +
            y_Vi[s, co2] * H_CO2_vap(T_V_in[s]) +
            y_Vi[s, n2] * H_N2_vap(T_V_in[s]) +
            y_Vi[s, o2] * H_O2_vap(T_V_in[s]))
    else
        H_vap_in = nothing
    end

    H_vap_out = @expression(model, [s in S, n in 1:N_seg],
        y_seg[s, n, h2o] * V_seg[s, n] * H_H2O_vap(T_seg[s, n]) +
        y_seg[s, n, co2] * V_seg[s, n] * H_CO2_vap(T_seg[s, n]) +
        y_seg[s, n, n2] * V_seg[s, n] * H_N2_vap(T_seg[s, n]) +
        y_seg[s, n, o2] * V_seg[s, n] * H_O2_vap(T_seg[s, n]))

    return (H_vap_in, H_vap_out)
end

###############################################################################
# minimal column diameter calculation based on flooding
packing_factor = 12  # 1/ft, for Mellapak Plus 252Y, Kister HZ, Scherffius J, Afshar K, Abkar E. Chemical Engineering Progress (2007).
minimum_liquid_load = 0.2  # m3/m2h, minimum liquid load to ensure sufficient wetting, Sulzer Ltd. Structured Packings Enery-Efficient, Innovative & Profitable. (2024).
pressure_drop_flooding = 0.115 * packing_factor^0.7 # in. H2O/ft packing, pressure drop at flooding conditions, Kister HZ, Scherffius J, Afshar K, Abkar E. Chemical Engineering Progress (2007)

function calculate_CP_GPDC_flood(F_LV, pressure_drop_flooding=pressure_drop_flooding)
"""
Calculates the capacity factor based on the generalized pressure drop correlation (GPDC) for structured packing based on the flow parameter F_LV.
Tsai RE. Mass transfer area of structured packing. Dissertation. The University of Texas at Austin (2010).
"""
    C1 = 3.8617 
    C2 = 0.66090
    C3 = 6.3763
    C4 = 0.7206
    C5 = 0.2898
    C6 = -0.9093
    C7 = -0.6819
    return C1 * (pressure_drop_flooding)^C2 / (1+C3*pressure_drop_flooding^(C2/C4)*F_LV^C5)^C4 * (1-exp(C6*F_LV^C7))
end

function calculate_flow_parameter(L, G, rho_g, rho_l)
    """
    calculates the flow parameter

    L: liquid mass flow rate
    G: vapor mass flow rate
    rho_g: vapor density
    rho_l: liquid density
    units must be consistent for liquid and gas phases

    Kister HZ, Scherffius J, Afshar K, Abkar E. Chemical Engineering Progress (2007).
    """
    return L / G * sqrt(rho_g / rho_l)
end

function calculate_CP_superficial_velocity(U_s_ft_per_s, rho_g, rho_l, nu)
    """
    Calculates the capacity factor based on superficial gas velocity

    U_s_ft_per_s: superficial gas velocity in ft/s
    rho_g: vapor density
    rho_l: liquid density
    packing_factor: packing factor in 1/ft
    nu: kinematic viscosity of the liquid in cSt
    
    Kister HZ, Scherffius J, Afshar K, Abkar E. Chemical Engineering Progress (2007).
    """
    return U_s_ft_per_s * (rho_g / (rho_l-rho_g) * packing_factor)^0.5 * nu^0.05
end

function calculate_ideal_gas_density(P, MW, T)
    """
    Calculates ideal gas density in kg/m^3 based on ideal gas law
    P: pressure in atm
    MW: average molecular weight of the gas in kg/mol
    T: temperature in K
    """
    return P*101325 * MW / (R * T)  # density in kg/m3, pressure conversion with 101325 Pa/atm
end

###############################################################################
# helper for variable declarations
const INERTS = [:N2, :O2]
const NONVOLATILES = [:MEA, :MEA2p, :Hp, :Carbamate, :HCO3m, :OHm, :CO3_2m]

function fix_vapor(log_y, S, N_seg;species=NONVOLATILES)
    """ 
        Fixes non-volatile species vapor mole fraction to 0.
        Default non-volatile species are MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m
    """
    for s in S
        for n in 1:N_seg
            for spec in species
                fix(log_y[s, n, species_index[spec]], -50.0; force=true)
            end
        end
    end
end

function fix_liquid(log_x, S, N_seg; species = INERTS)
    """
    Fix inert species liquid mole fraction to a very small value
    Default inert species are N2, O2
    """
    for s in S
        for n in 1:N_seg
            for spec in species
                fix(log_x[s, n, species_index[spec]], -50.0; force=true)
            end
        end
    end
end

# helper for accessing initial values: check if a JuMP variable is fixed, if so return its fixed value, otherwise return its start value (from upstream)
check_val(v) = is_fixed(v) ? fix_value(v) : start_value(v)

function init_mole_frac(log_x, log_y, log_x0, log_y0, S, N_seg)
    """
    Helper to initialize mole fractions for single- or multi-segment units with the provided values.
    """
    # Initialize log mole fractions
    for s in S
        for n in 1:N_seg
            for c in 1:length(species_index)
                if log_x0 !== nothing
                    set_start_value(log_x[s, n, c], log_x0[s, n, c])
                end
                if log_y0 !== nothing
                    set_start_value(log_y[s, n, c], log_y0[s, n, c])
                end
            end
        end
    end
end

function init_LVT(L, V, T, L0, V0, T0, S, N_seg)
    """
    Helper to initialize L, V, T for single- or multi-segment units with the provided values.
    """
    # Initialize liquid flow rates, vapor flow rates, and temperatures
    for s in S
        for n in 1:N_seg
            if L0 !== nothing
                set_start_value(L[s, n], L0[s]) 
            end
            if V0 !== nothing
                set_start_value(V[s, n], V0[s]) 
            end
            if T0 !== nothing
                set_start_value(T[s, n], T0[s])  
            end
        end
    end
end

###############################################################################
# Stoichiometry matrix for reactions (9 liquid species x 5 reactions)
# Species order: MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m, CO2, H2O
# Reactions: 1: MEA protonation, 2: HCO3- dissociation, 3: H2O dissociation, 4: Carbamate formation, 5: CO2 hydration
stoich = [1  0  0  1  0;    # MEA
         -1  0  0  0  0;    # MEA2p (MEAH+)
          1  1  1  0  1;    # Hp
          0  0  0 -1  0;    # Carbamate (MEACOO-)
          0 -1  0  1  1;    # HCO3m
          0  0  1  0  0;    # OHm
          0  1  0  0  0;    # CO3_2m
          0  0  0  0 -1;    # CO2
          0  0 -1 -1 -1]    # H2O

function get_lnK(T)
    """
    Equilibrium constant calculation based on temperature.
    """
    A = [-3.038325, 216.049, 132.899, -0.52135, 231.465]
    B = [-7008.357, -12431.7, -13445.9, -2545.53, -12092.1]
    C = [0.0, -35.4819, -22.4773, 0.0, -36.7816]
    D = [-0.00313489, 0.0, 0.0, 0.0, 0.0]
    lnK = A .+ B ./ T .+ C .* log.(T) .+ D .* T
    return lnK
end

function get_Pvap(des_T)
    """
    Vapor pressure of water, in Pa
    from Perry's Engineer's Handbook, 8th edition
    """
    # Coefficients from correlation
    C1 = 73.649
    C2 = -7258.2
    C3 = -7.3037
    C4 = 4.1653e-6
    C5 = 2
    return exp(C1 + C2 / des_T + C3 * log(des_T) + C4 * des_T^C5)
end

function get_H_hartono(x_MEA, x_H2O, alpha, T_L)
    """
    Henry constant [Pa·m^3/mol] for physical CO2 solubility in loaded MEA, based on:
    Hartono A, Mba EO, Svendsen HF. Physical properties of partially CO2 loaded aqueous monoethanolamine (MEA). J Chem Eng Data 59(6):1808–1816 (2014).

    Arguments:
        x_MEA : solvent-basis mole fraction of MEA in the unloaded MEA/H2O solvent
        x_H2O : solvent-basis mole fraction of H2O in the unloaded MEA/H2O solvent
                (assumed x_MEA + x_H2O = 1)
        α     : CO2 loading, mol CO2 / mol MEA
        T_L   : liquid temperature [K]
    """
    M_MEA = MW[:MEA]/1000 # kg/mol
    M_H2O = MW[:H2O]/1000 # kg/mol

    # Pure-component Henry correlations
    H_N2O_MEA_kPa = 2.448e5 * exp(-1348.0 / T_L) / 1000 # kPa·m^3/mol
    H_CO2_H2O_kPa = 3.52e6 * exp(-2113.0 / T_L) / 1000 # kPa·m^3/mol
    H_N2O_H2O_kPa = 8.449e6 * exp(-2283.0 / T_L) / 1000 # kPa·m^3/mol

    # Pure-component density
    rho_MEA = 1000 # kg/m^3
    rho_H2O = 1000 # kg/m^3

    # Unloaded N2O model; volume fractions
    phi_MEA = (x_MEA * M_MEA / rho_MEA) / ((x_MEA * M_MEA / rho_MEA) + (x_H2O * M_H2O / rho_H2O))
    phi_H2O = 1.0 - phi_MEA

    # Excess molar volume V_E for unloaded MEA/H2O
    # temperatures in °C here
    t_C = T_L - 273.15
    V_E = (-1.9210 + 1.6792e-3 * t_C - 3.0951 * x_MEA + 3.4412 * x_MEA^2) * x_MEA * x_H2O * 1e-6   # m^3/mol

    # Eq. (11)
    lndH_unloaded = -1.0e6 * V_E * sqrt(1.0 - phi_MEA^2)
    # Eq. (10)
    lnH_N2O_unloaded = lndH_unloaded + phi_MEA * log(H_N2O_MEA_kPa) + phi_H2O * log(H_N2O_H2O_kPa)
    H_N2O_unloaded_kPa = exp(lnH_N2O_unloaded)

    # Loaded N2O model: Eqs. (12)-(13)
    b1 = 0.77
    b2 = 0.033

    lndH_loaded = b1 * x_MEA + b2 * x_MEA * alpha * T_L
    lnH_N2O_loaded = lndH_loaded + log(H_N2O_unloaded_kPa)
    H_N2O_loaded_kPa = exp(lnH_N2O_loaded)

    # N2O analogy for CO2
    H_CO2_loaded_kPa = H_N2O_loaded_kPa * (H_CO2_H2O_kPa / H_N2O_H2O_kPa)

    # Convert to Pa·m^3/mol
    He_CO2 = H_CO2_loaded_kPa * 1e3
    return He_CO2
end

###############################################################################
# Enthalpy calculation

function get_Hvap(T) 
    """
    Calculate heat of vaporization of water.
    Temperature in K
    """
    c = [5.66e4, 0.61204, -0.6257, 0.3988]
    Tc = 647.096
    Tr = T / Tc
    c[1] * (1 - Tr)^(c[2] + c[3] * Tr + c[4] * Tr^2)
end


# parameters for liquid and vapor enthalpy calculation
C_MEA = (2.6161, 3.706e-3, 3.787e-6, 0.0, 0.0)
C_H2O = (4.2107, -1.696e-3, 2.568e-5, -1.095e-7, 3.038e-10)
A_CO2 = (5.457, 1.045e-3, -1.157e5)
A_H2O = (3.47, 1.45e-3, 1.21e4)
A_N2 = (3.28, 5.93e-4, 4.0e3)
A_O2 = (3.639, 5.06e-4, -2.27e4)

function H_H2O_liq(T)
    """
    Calculates liquid enthalpy for water (J/mol)
    T: temperature in K
    """
    return MW[:H2O] * (
        C_H2O[1] * ((T - 273.15) - 25) +
        C_H2O[2] * ((T - 273.15)^2 - 25^2) / 2 +
        C_H2O[3] * ((T - 273.15)^3 - 25^3) / 3 +
        C_H2O[4] * ((T - 273.15)^4 - 25^4) / 4 +
        C_H2O[5] * ((T - 273.15)^5 - 25^5) / 5
    ) #25C is the reference temperature
end

function H_MEA_liq(T)
    """
    Calculates liquid enthalpy for MEA (J/mol)
    T: temperature in K
    """
    return MW[:MEA] * (
        C_MEA[1] * ((T - 273.15) - 25) +
        C_MEA[2] * ((T - 273.15)^2 - 25^2) / 2 +
        C_MEA[3] * ((T - 273.15)^3 - 25^3) / 3
    ) #25C is the reference temperature
end

H_H2O_vap(T) = R * (A_H2O[1] * (T - 298.15) + A_H2O[2] * ((T^2 - 298.15^2) / 2) - A_H2O[3] * (1 / T - 1 / 298.15)) + get_Hvap(298.15)
H_CO2_vap(T) = R * (A_CO2[1] * (T - 298.15) + A_CO2[2] * ((T^2 - 298.15^2) / 2) - A_CO2[3] * (1 / T - 1 / 298.15))
H_N2_vap(T) = R * (A_N2[1] * (T - 298.15) + A_N2[2] * ((T^2 - 298.15^2) / 2) - A_N2[3] * (1 / T - 1 / 298.15))
H_O2_vap(T) = R * (A_O2[1] * (T - 298.15) + A_O2[2] * ((T^2 - 298.15^2) / 2) - A_O2[3] * (1 / T - 1 / 298.15))

###############################################################################
# Helpers for heat exchanger property calculations
# Heat capacity of the liquid phase, 
# Hilliard MD. A predictive thermodynamic model for an aqueous blend of potassium carbonate, piperazine, and monoethanolamine for carbon dioxide capture from flue gas. Dissertation. University of Texas at Austin (2008).

C_arr_MEA = [2.6161, 3.706e-3, 3.787e-6, 0, 0]
C_arr_H2O = [4.2107, -1.696e-3, 2.568e-5, -1.095e-7, 3.038e-10]
function calc_Cp(mole_frac, t)
    """
    This function returns the specific heat capacity of liquid in each pass in J/(mol*K)
    mole_frac: [mole fraction MEA, mole fraction H2O] in the liquid phase
    t: temperature in K
    """
    Cp_MEA = MW[:MEA] * (C_arr_MEA[1] .+ C_arr_MEA[2] * (t .- 273.15) .+ C_arr_MEA[3] * (t .- 273.15) .^ 2 .+ C_arr_MEA[4] * (t .- 273.15) .^ 3 .+ C_arr_MEA[5] * (t .- 273.15) .^ 4)
    Cp_H2O = MW[:H2O] * (C_arr_H2O[1] .+ C_arr_H2O[2] * (t .- 273.15) .+ C_arr_H2O[3] * (t .- 273.15) .^ 2 .+ C_arr_H2O[4] * (t .- 273.15) .^ 3 .+ C_arr_H2O[5] * (t .- 273.15) .^ 4)

    return mole_frac[1] * Cp_MEA + mole_frac[2] * Cp_H2O
end

function conductivity_H2O(T)
    """
    This function estimates the thermal conductivity of pure water in W/(m*K) based on 
    Ramires, M. L. V. et al. Standard Reference Data for the Thermal Conductivity of Water. Journal of Physical and Chemical Reference Data 24, 1377–1381 (1995).
    T: temperature in K, 274K≤T≤370K
    """
    Tr = T / 293.15 # reduced temperature
    return 0.6065 * (-1.48445 + 4.12292 * Tr - 1.63866 * Tr^2)
end

function conductivity_MEA(T)
    """
    This function estimates the thermal conductivity of pure MEA in W/(m*K) based on 
    DiGuilio, R. M., McGregor, W. L. & Teja, A. S. Thermal conductivities of the ethanolamines. J. Chem. Eng. Data 37, 242–245 (1992).
    T: temperature in K, 298K≤T≤470K
    """
    return 0.27719 - 1.2509e-4 * T
end

function conductivity_Vredeveld(x, k)
    """
    The Vredeveld mixing rule estimates the thermal conductivity of a mixed liquid in W/(m*K) J/(s*m*K).
    x: vector of component fractions 
    k: vector of pure component thermal conductivities in W/(m*K)
    """
    return sum(x .* k .^ (-2))^(-1 / 2)
end

function calc_viscosity(t, mole_frac, alpha)
    """
    Dynamic viscosity, Pa s
    Akula P, Eslick J, Bhattacharyya D, Miller DC. Ind Eng Chem Res 60(14):5176–5193 (2021).
    """
    mu_param = [-0.0838, 2.8817, 33.651, 1817, 0.00847, 0.0103, -2.3890]
    r = MW[:MEA] * mole_frac[1] / (MW[:MEA] * mole_frac[1] + MW[:H2O] * mole_frac[2])
    omega = 100 * r

    mu_H2O = 1.002 * 10 .^ (1.3272 * (293.15 .- t .- 0.001053 .* (t .- 293.15) .^ 2) ./ (t .- 168.15))
    mu = mu_H2O ./ 1000 .* exp.(
        ((mu_param[1] * omega .+ mu_param[2]) .* t .+ mu_param[3] * omega .+ mu_param[4]) .*
        (alpha .* (mu_param[5] * omega .+ mu_param[6] .* t .+ mu_param[7]) .+ 1) .* omega ./ t .^ 2
    )
    return mu
end

###############################################################################
# region helpers for initialization based on previously solved models
function clamp_to_bounds(x::Float64, v::JuMP.VariableRef; eps=1e-12)
        hasLB = has_lower_bound(v)
        hasUB = has_upper_bound(v)

        lb = -Inf
        ub = Inf

        if hasLB
            lb = lower_bound(v)
        end
        if hasUB
            ub = upper_bound(v)
        end
    x1 = (lb != -Inf && x < lb) ? lb + eps : x
    x2 = (ub != Inf && x1 > ub) ? ub - eps : x1
    return x2
end

function _split_var_name(nm::AbstractString)
    """
    Splits the variable name from the indices. 
    E.g., x[1,2,3] --> "x", (1,2,3)
    """
    m = match(r"^([^\[]+)(?:\[(.*)\])?$", nm)
    m === nothing && error("Could not parse variable name: $nm")

    base = String(m.captures[1])
    idxs_str = m.captures[2]

    if idxs_str === nothing || isempty(idxs_str)
        return base, ()
    end

    idx = Tuple(parse(Int, strip(tok)) for tok in split(idxs_str, ","))
    return base, idx
end

function _ensure_length!(v::Vector{Int}, n::Int)
    while length(v) < n
        push!(v, 0)
    end
    return v
end

function _seed_value(v::JuMP.VariableRef)
    """
    gets the best seed available for a source variable.

    Preference order:
    1. fixed value            (if variable is fixed),
    2. solved value           (if model was solved),
    3. existing start value   (if no solved value is available).

    Returns:
    - Float64 seed if available
    - nothing if no usable value exists
    """
    # Prefer actual solution value if available; otherwise use fixed value or start value.
    if is_fixed(v)
        return Float64(fix_value(v))
    end

    val = try
        value(v)
    catch
        nothing
    end
    if val !== nothing && isfinite(val)
        return Float64(val)
    end

    sv = try
        start_value(v)
    catch
        nothing
    end
    if sv !== nothing && isfinite(sv)
        return Float64(sv)
    end

    return nothing
end

function _build_seed_dict(model::Model)
    """
    Collect initial values for all variables of the model.
    """
    exact = Dict{String,Float64}()
    by_base = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    extents = Dict{String,Vector{Int}}()

    for v in all_variables(model)
        seed = _seed_value(v)
        seed === nothing && continue

        nm = name(v)
        exact[nm] = seed

        base, idx = _split_var_name(nm)
        get!(by_base, base, Dict{Tuple{Vararg{Int}},Float64}())[idx] = seed

        ex = get!(extents, base, Int[])
        _ensure_length!(ex, length(idx))
        for d in 1:length(idx)
            ex[d] = max(ex[d], idx[d])
        end
    end

    return exact, by_base, extents
end

function _build_extents(model::Model)
    """
    builds only index extents for the target model.

    This tells us the target size of each variable family. For example,
    if target has abs_L[s=1:8, n=1:20], then:
    extents["abs_L"] = [8, 20]
    """
    extents = Dict{String,Vector{Int}}()

    for v in all_variables(model)
        base, idx = _split_var_name(name(v))
        ex = get!(extents, base, Int[])
        _ensure_length!(ex, length(idx))
        for d in 1:length(idx)
            ex[d] = max(ex[d], idx[d])
        end
    end

    return extents
end

function _scenario_source_index(
    s_target::Int,
    n_source::Int,
    scenario_map::AbstractDict{Int,Int},
    scenario_fallback::Symbol,
)
    """
        decides which source scenario should be used for a target
        scenario when the target has more scenarios than the source.

        Priority:
        1. explicit mapping from scenario_map
        2. If same number of scenarios, we assume they correspond directly and keep the same scenario index.
        3. fallback rule:
        - :first   -> always use scenario 1
        - :last    -> always use last source scenario
        - :nearest -> clamp target scenario to available source range
    """
    # 1) explicit mapping has highest priority
    if haskey(scenario_map, s_target)
        return clamp(scenario_map[s_target], 1, n_source)
    end

    # 2) if the target scenario already exists in the source, keep it
    if 1 <= s_target <= n_source
        return s_target
    end

    # 3) otherwise use fallback for genuinely new scenarios
    if scenario_fallback == :first
        return 1
    elseif scenario_fallback == :last
        return n_source
    elseif scenario_fallback == :nearest
        return clamp(s_target, 1, n_source)
    else
        error("Unsupported scenario_fallback = $scenario_fallback. Use :first, :last, or :nearest.")
    end
end


function _interp_bracket(i_target::Int, n_target::Int, n_source::Int)
    """
    interpolation between column segments, map a target segment index onto the
    source segment grid and return the two source brackets plus weights.

    Example:
    source has 10 segment, target has 20 segment.
    target segment 8 might map between source segment 4 and 5.

    Returns:
    i1, w1, i2, w2
    meaning:
    interpolated_value = w1 * value(i1) + w2 * value(i2)
    """
    # map target segment index to source segment coordinate on [1, n_source]
    if n_source <= 1
        return 1, 1.0, 1, 0.0
    elseif n_target <= 1
        return 1, 1.0, 1, 0.0
    end

    x = 1.0 + (i_target - 1) * (n_source - 1) / (n_target - 1)
    i1 = clamp(floor(Int, x), 1, n_source)
    i2 = clamp(ceil(Int, x), 1, n_source)

    w2 = x - i1
    w1 = 1.0 - w2
    return i1, w1, i2, w2
end

function _mapped_seed(
    base::String,
    idx_target::Tuple{Vararg{Int}},
    source_by_base::Dict{String,Dict{Tuple{Vararg{Int}},Float64}},
    source_extents::Dict{String,Vector{Int}},
    target_extents::Dict{String,Vector{Int}};
    interpolate_seg_dim::Bool = true,
    scenario_map::AbstractDict{Int,Int} = Dict{Int,Int}(),
    scenario_fallback::Symbol = :last,
)
    """
    obtains a seed for one target variable (with name base) by:
    - matching its base variable family,
    - mapping scenario index,
    - interpolating segment index if needed,
    - copying all other indices directly.

    Note: this assumes:
    - dimension 1 is scenario,
    - dimension 2 is segment/segment,
    - higher dimensions (species, reactions, etc.) are copied directly.
    """
    tbl = get(source_by_base, base, nothing)
    tbl === nothing && return nothing, :skipped

    src_ext = get(source_extents, base, Int[])
    tgt_ext = get(target_extents, base, Int[])

    # ranks should normally match for the same base variable name
    if length(idx_target) != length(src_ext)
        return nothing, :skipped
    end

    # scalars
    if isempty(idx_target)
        return get(tbl, (), nothing), :indexed
    end

    # interpolation indices (if no interpolation indices will be identical)
    idx_lo = Vector{Int}(undef, length(idx_target))
    idx_hi = Vector{Int}(undef, length(idx_target))

    #flags for reporting
    copied_scenario = false
    interpolated_seg = false
    w_hi = 0.0

    for d in 1:length(idx_target)  #loop over each dimension of the variable
        # check source and target extents for this dimension
        n_s = src_ext[d]
        n_t = d <= length(tgt_ext) ? tgt_ext[d] : n_s

        if d == 1  # number of scenarios
            # first dimension = scenario
            i_src = _scenario_source_index(idx_target[d], n_s, scenario_map, scenario_fallback)
            copied_scenario |= (i_src != idx_target[d])
            idx_lo[d] = i_src  # no interpolation needed, so both indices are the same
            idx_hi[d] = i_src

        elseif d == 2 && interpolate_seg_dim && n_s != n_t && n_s >= 1 && n_t >= 1  # number of column segments
            i1, _w1, i2, w2 = _interp_bracket(idx_target[d], n_t, n_s)  #get indices and weights for interpolation
            idx_lo[d] = i1
            idx_hi[d] = i2
            w_hi = w2
            interpolated_seg |= (i1 != i2 || n_s != n_t)

        else
            # all other dimensions copied directly (species, reactions, etc.)
            i_src = clamp(idx_target[d], 1, n_s)
            idx_lo[d] = i_src
            idx_hi[d] = i_src
        end
    end

    key_lo = Tuple(idx_lo)
    key_hi = Tuple(idx_hi)
    #e.g., key_lo = (5,6,8), key_hi = (5,7,8) --> copy from source scenario 5, interpolate between source segments 6 and 7, keep species index 8.

    val_lo = get(tbl, key_lo, nothing)  #lower source value
    val_lo === nothing && return nothing, :skipped

    seed =
        if interpolated_seg
            val_hi = get(tbl, key_hi, val_lo) # get higher source value
            (1.0 - w_hi) * val_lo + w_hi * val_hi # interpolate between lower and higher source values
        else
            val_lo # if no interpolation needed, just use the lower value (which is the only value in this case)
        end

    mode =
        if interpolated_seg && copied_scenario
            :copied_scenario_interpolated
        elseif interpolated_seg
            :interpolated
        elseif copied_scenario
            :copied_scenario
        else
            :indexed
        end

    return seed, mode
end

function transfer_solution_same_variables(
    source_model::Model,
    target_model::Model;
    interpolate_seg_dim::Bool = true,
    scenario_map::AbstractDict{Int,Int} = Dict{Int,Int}(),
    scenario_fallback::Symbol = :last,   # :first, :last, :nearest
    verbose::Bool = true,
    )
    """
    Transfers solution values or existing start values from `source_model` to
    `target_model`.

    Strategy:
    - Exact variable-name matches are transferred directly.
    - If an exact match does not exist, the function falls back to matching
    variables by base name and index structure.
    - For segment-like dimensions (assumed to be index 2), we interpolate
    from the coarse grid to the fine grid.
    - For additional scenarios (assumed to be index 1), we copy from
    an existing source scenario.

    Keyword arguments:
    - `interpolate_seg_dim=true`:
        interpolate along the second index (column segments) if source and target
        have different sizes in that dimension.
    - `scenario_map=Dict{Int,Int}()`:
        optional explicit mapping from target scenario => source scenario.
        Example: Dict(6 => 3, 7 => 3, 8 => 5)
    - `scenario_fallback=:last`:
        how to choose a source scenario if a target scenario is not in
        `scenario_map`. Allowed values: `:first`, `:last`, `:nearest`.
    - `verbose=true`:
        print summary output
    """
    source_exact, source_by_base, source_extents = _build_seed_dict(source_model)
    isempty(source_exact) && error("No usable values or start values found in source_model.")
    !(raw_status(source_model) in ["Solve_Succeeded", "Solved_To_Acceptable_Level", "Feasible_Point_Found"]) && error("Source model did not solve successfully. Status: $(raw_status(source_model))")
    
    target_extents = _build_extents(target_model)

    counts = Dict(
        :exact => 0,
        :indexed => 0,
        :copied_scenario => 0,
        :interpolated => 0,
        :copied_scenario_interpolated => 0,
        :skipped => 0,
    )

    for v_t in all_variables(target_model)
        # fixed vars do not need a start value
        is_fixed(v_t) && continue

        nm_t = name(v_t)

        seed = get(source_exact, nm_t, nothing)
        mode = :exact

        if seed === nothing
            base, idx_target = _split_var_name(nm_t)

            seed, mode = _mapped_seed(
                base,
                idx_target,
                source_by_base,
                source_extents,
                target_extents;
                interpolate_seg_dim = interpolate_seg_dim,
                scenario_map = scenario_map,
                scenario_fallback = scenario_fallback,
            )
        end

        if seed === nothing || !isfinite(seed)
            counts[:skipped] += 1
            continue
        end

        set_start_value(v_t, clamp_to_bounds(Float64(seed), v_t))
        counts[mode] += 1
    end

    if verbose
        @info(
            "Transferred start values",
            exact = counts[:exact],
            indexed = counts[:indexed],
            copied_scenario = counts[:copied_scenario],
            interpolated = counts[:interpolated],
            copied_scenario_interpolated = counts[:copied_scenario_interpolated],
            skipped = counts[:skipped],
        )
    end

    return nothing
end

###############################################################################
# helpers for result processing and output formatting
function export_variable_values_split(model)
    # initialize one dictionary per category
    abs = Dict{Symbol,Any}()
    HX = Dict{Symbol,Any}()
    des = Dict{Symbol,Any}()
    reb = Dict{Symbol,Any}()
    cond = Dict{Symbol,Any}()
    cool = Dict{Symbol,Any}()
    other = Dict{Symbol,Any}()

    # Temporary storage for indexed variables:
    # base_name => Dict{index_tuple => value}
    abs_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    HX_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    des_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    reb_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    cond_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    cool_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()
    other_temp = Dict{String,Dict{Tuple{Vararg{Int}},Float64}}()

    # choose which temp dict to use from name prefix
    choose_temp(prefix) =
        occursin("abs", prefix) ? abs_temp :
        occursin("HX", prefix) ? HX_temp :
        occursin("des", prefix) ? des_temp :
        occursin("reb", prefix) ? reb_temp :
        occursin("cond", prefix) ? cond_temp :
        occursin("cool", prefix) ? cool_temp :
        other_temp

    # parse "x[1,10]" -> (1,10)
    function parse_indices(nm_str::String)
        m = match(r"\[(.+)\]", nm_str)
        if m === nothing
            return ()
        end
        idx_str = m.captures[1]  # e.g. "1,10" or "2"
        index_parts = split(idx_str, ",")
        return Tuple(parse.(Int, index_parts))
    end

    # loop over all variables
    for v in all_variables(model)
        nm_str = name(v)
        val = try
            value(v)
        catch
            NaN
        end

        # get variable name without indices, then store value of variable under base_name and with index
        if occursin("[", nm_str)
            # Indexed variable
            base_name = split(nm_str, "[")[1]
            prefix = first(base_name, min(4, lastindex(base_name)))
            temp_dict = choose_temp(prefix)

            idx = parse_indices(nm_str)
            inner = get!(temp_dict, base_name, Dict{Tuple{Vararg{Int}},Float64}())
            inner[idx] = val
        else
            # Scalar variable (no indices), todo: probably useful for design variables, else delete
            nm_key = Symbol(nm_str)
            prefix = first(nm_str, min(4, lastindex(nm_str)))

            if occursin("abs", prefix)
                abs[nm_key] = val
            elseif occursin("HX", prefix)
                HX[nm_key] = val
            elseif occursin("des", prefix)
                des[nm_key] = val
            elseif occursin("reb", prefix)
                reb[nm_key] = val
            elseif occursin("cond", prefix)
                cond[nm_key] = val
            elseif occursin("cool", prefix)
                cool[nm_key] = val
            else
                other[nm_key] = val
            end
        end
    end

    function build_array(idxmap::Dict{Tuple{Vararg{Int}},Float64})
        """
        saves variable values in Dict{Tuple{Int,...},Float64} as an Array with one dim per index
        """
        if isempty(idxmap)
            return Float64[]
        end
        idxs = collect(keys(idxmap))  # e.g. [(1,1), (1,2), (2,1), (2,2)]
        nd = length(first(idxs))  # number of dimensions
        dims = ntuple(d -> maximum(t[d] for t in idxs), nd) # get maximum index in each dimension
        A = fill(NaN, dims...)
        for (I, val) in idxmap
            A[I...] = val
        end
        return A
    end

    # Convert indexed variables for each category
    for (base_name, idxmap) in abs_temp
        abs[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in HX_temp
        HX[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in des_temp
        des[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in reb_temp
        reb[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in cond_temp
        cond[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in cool_temp
        cool[Symbol(base_name)] = build_array(idxmap)
    end
    for (base_name, idxmap) in other_temp
        other[Symbol(base_name)] = build_array(idxmap)
    end

    return (abs, HX, des, reb, cond, cool, other)
end

# variables access helpers 
function base_name(var_container)
    """Given a JuMP variable, return the base name as a Symbol (without indices). E.g. for abs_L[1,1], return :abs_L."""
    # Take the first element's JuMP name and strip indices
    first_var = first(var_container)        # a VariableRef
    full = name(first_var)                 # "abs_L[1,1]"
    base = split(full, '[', limit=2)[1]    # "abs_L"
    return Symbol(base)
end

function slice_by_scenario(var, s)
    """From the DenseAxisArray or Array `var`, extract the slice corresponding to scenario index `s`.
    Assumes that the first dimension of `var` corresponds to scenarios."""
    if ndims(var) == 1
        return var[s]               # e.g. MU_H2O[s]
    elseif ndims(var) == 2
        return var[s, :]            # e.g. abs_L[s, :]
    else
        # More general slice for >2 dims: var[s, :, :, ...]
        return view(var, s, ntuple(_ -> Colon(), ndims(var) - 1)...)
    end
end

"""
    getcont(model, sym) -> container or nothing

Safe container lookup: returns `model[sym]` if present, else `nothing`.
"""
getcont(model::JuMP.Model, sym::Symbol) = (try model[sym] catch; nothing end)

function v(model::JuMP.Model, sym::Symbol, inds...;
           scenario::Int=1,
           defaults::Dict{Symbol,Any}=Dict{Symbol,Any}())
    """
    Generic scalar accessor.
    - If `defaults` has `sym`, return that.
    - Tries `cont[scenario, inds...]`, then `cont[inds...]`, then `cont` as scalar.
    - Returns `missing` if not available or not indexable that way.
    """
    cont = getcont(model, sym)
    cont === nothing && return get(defaults, sym, missing)
    cont isa Number && return cont

    # scenario-first
    try
        return value(cont[scenario, inds...])
    catch
    end

    # without scenario dimension
    try
        return value(cont[inds...])
    catch
    end

    # scalar container / single variable
    try
        return value(cont)
    catch
        return missing
    end
end

function vlog(model::JuMP.Model, logsym::Symbol, inds...;
              scenario::Int=1,
              defaults::Dict{Symbol,Any}=Dict{Symbol,Any}())
    """
    Like `v`, but exponentiates the retrieved scalar (for log-formulated variables).
    """
    x = v(model, logsym, inds...; scenario=scenario, defaults=defaults)
    ismissing(x) && return missing
    return exp(x)
end

function xlog_at(model::JuMP.Model, logsym::Symbol, segment::Int, sp::Symbol;
                 scenario::Int=1,
                 species_index, defaults)
    """
    Log-mole-fraction accessor for 3D containers `(scenario, segment, species)`.
    Name reflects that it is for log-formulated mole fractions.
    """
    k = get(species_index, sp, nothing)
    k === nothing && return missing
    return vlog(model, logsym, segment, k; scenario=scenario, defaults=defaults)
end

function loading_liq(model::JuMP.Model, logx_sym::Symbol, segment::Int, acid_species::AbstractVector{Symbol}, amine_species::AbstractVector{Symbol};
                     scenario::Int=1,
                     species_index, defaults::Dict{Symbol,Any}=Dict{Symbol,Any}())
    """
    Compute loading = x_acid / x_amine from log-liquid mole fractions.
    """
    x_acid  = sum(xlog_at(model, logx_sym, segment, sp;  scenario=scenario, species_index=species_index, defaults=defaults) for sp in acid_species)
    x_amine = sum(xlog_at(model, logx_sym, segment, sp;  scenario=scenario, species_index=species_index, defaults=defaults) for sp in amine_species)
    (ismissing(x_acid) || ismissing(x_amine) || x_amine == 0) && return missing
    return x_acid / x_amine
end

function apparent_frac(model::JuMP.Model,
                       sym::Symbol,
                       species_list::AbstractVector{Symbol};
                       scenario::Int=1,
                       segment::Union{Nothing,Int}=nothing,
                       is_log::Bool=false,
                       species_index)
    """
    Sum of mole fractions over `species_list` from container `sym`. Simplifies apparent species calculations.

    If `is_log=true`, applies exp() to each entry (log-formulated mole fractions).
    Returns `missing` if any species index or value is missing.
    """
    tot = 0.0
    for sp in species_list
        k = get(species_index, sp, nothing)
        k === nothing && return missing

        if segment === nothing
            is_log ? xi = vlog(model, sym, k; scenario=scenario) :
                     xi = v(model, sym, k; scenario=scenario)
        else
            is_log ? xi = vlog(model, sym, segment, k; scenario=scenario) :
                     xi = v(model, sym, segment, k; scenario=scenario)
        end

        ismissing(xi) && return missing
        tot += xi
    end
    return tot
end

function collect_results(model::JuMP.Model;
    scenario::Int=1,
    N_abs::Int,
    N_des::Int,
    species_index,
    defaults::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    add_column_properties::Bool=false,
    sensitivity_CCF::Float64
)
    """
    function to extract select model results in a dataframe
    """
    rows = NamedTuple[]
    add(stream, prop, val, units="") =
        push!(rows, (stream=stream, property=prop, value=val, units=units, scenario=scenario))

    # add results to be tabulated 
    # Sour gas (absorber vapor inlet)
    add("Sour gas", "Temperature", v(model, :abs_T_V_in; scenario=scenario, defaults=defaults), "K")
    add("Sour gas", "Pressure", v(model, :abs_P_in; scenario=scenario, defaults=defaults), "Pa")
    add("Sour gas", "Molar flow", v(model, :abs_V_in; scenario=scenario, defaults=defaults), "mol/s")
    add("Sour gas", "CO2 mole fraction", v(model, :abs_y_in, species_index[:CO2]; scenario=scenario, defaults=defaults), "-")

    # Lean amine to absorber
    add("Lean solvent, absorber feed", "Temperature", v(model, :abs_T_L_in; scenario=scenario, defaults=defaults), "K")
    add("Lean solvent, absorber feed", "Pressure", v(model, :abs_P_in; scenario=scenario, defaults=defaults), "Pa")
    add("Lean solvent, absorber feed", "Molar flow", v(model, :abs_L_in; scenario=scenario, defaults=defaults), "mol/s")

    MEA_app_lean_solvent_abs = apparent_frac(model, :abs_x_in, [:MEA, :MEA2p, :Carbamate]; scenario=scenario, species_index=species_index)
    H2O_app_lean_solvent_abs = apparent_frac(model, :abs_x_in, [:H2O, :OHm, :Hp, :HCO3m];scenario=scenario, species_index=species_index)
    CO2_app_lean_solvent_abs = apparent_frac(model, :abs_x_in, [:CO2, :HCO3m, :CO3_2m, :Carbamate]; scenario=scenario, species_index=species_index)
    add("Lean solvent, absorber feed", "MEA mole fraction", MEA_app_lean_solvent_abs, "-")
    add("Lean solvent, absorber feed", "H2O mole fraction", H2O_app_lean_solvent_abs, "-")
    add("Lean solvent, absorber feed", "CO2 mole fraction", CO2_app_lean_solvent_abs, "-")

    # Sweet gas (absorber vapor outlet)
    add("Sweet gas", "Temperature", v(model, :abs_T, N_abs; scenario=scenario, defaults=defaults), "K")
    add("Sweet gas", "Pressure", v(model, :abs_P, N_abs; scenario=scenario, defaults=defaults), "Pa")
    add("Sweet gas", "Molar flow", v(model, :abs_V, N_abs; scenario=scenario, defaults=defaults), "mol/s")
    add("Sweet gas", "CO2 mole fraction", xlog_at(model, :abs_log_y, N_abs, :CO2; scenario=scenario, species_index=species_index, defaults=defaults), "-")

    # Rich solvent (absorber liquid outlet)
    add("Rich solvent, absorber outlet", "Temperature", v(model, :abs_T, 1; scenario=scenario, defaults=defaults), "K")
    add("Rich solvent, absorber outlet", "Pressure", v(model, :abs_P, 1; scenario=scenario, defaults=defaults), "Pa")
    add("Rich solvent, absorber outlet", "Molar flow", v(model, :abs_L, 1; scenario=scenario, defaults=defaults), "mol/s")
    MEA_app_rich_solvent_abs =  apparent_frac(model, :abs_log_x, [:MEA, :MEA2p, :Carbamate]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    H2O_app_rich_solvent_abs =  apparent_frac(model, :abs_log_x, [:H2O, :OHm, :Hp, :HCO3m]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    CO2_app_rich_solvent_abs =  apparent_frac(model, :abs_log_x, [:CO2, :HCO3m, :CO3_2m, :Carbamate]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    add("Rich solvent, absorber outlet", "MEA mole fraction", MEA_app_rich_solvent_abs, "-")
    add("Rich solvent, absorber outlet", "H2O mole fraction", H2O_app_rich_solvent_abs, "-")
    add("Rich solvent, absorber outlet", "CO2 mole fraction", CO2_app_rich_solvent_abs, "-")

    # Stripper feed (may be missing)
    add("Rich solvent, stripper feed", "Temperature", v(model, :HX_T_C_passes_out, 1; scenario=scenario, defaults=defaults), "K")
    add("Rich solvent, stripper feed", "Pressure", v(model, :des_P_in; scenario=scenario, defaults=defaults), "Pa")
    add("Rich solvent, stripper feed", "Molar flow", v(model, :HX_C_L, 1; scenario=scenario, defaults=defaults), "mol/s")

    # Stripper bottoms (segment 1 liquid)
    add("Lean solvent, stripper outlet", "Temperature", v(model, :des_T, 1; scenario=scenario, defaults=defaults), "K")
    add("Lean solvent, stripper outlet", "Pressure", v(model, :des_P, 1; scenario=scenario, defaults=defaults), "Pa")
    add("Lean solvent, stripper outlet", "Molar flow", v(model, :des_L, 1; scenario=scenario, defaults=defaults), "mol/s")
    MEA_app_rich_solvent_des =  apparent_frac(model, :des_log_x, [:MEA, :MEA2p, :Carbamate]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    H2O_app_rich_solvent_des =  apparent_frac(model, :des_log_x, [:H2O, :OHm, :Hp, :HCO3m]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    CO2_app_rich_solvent_des =  apparent_frac(model, :des_log_x, [:CO2, :HCO3m, :CO3_2m, :Carbamate]; scenario=scenario, segment=1, is_log=true, species_index=species_index)
    add("Lean solvent, stripper outlet", "MEA mole fraction", MEA_app_rich_solvent_des, "-")
    add("Lean solvent, stripper outlet", "H2O mole fraction", H2O_app_rich_solvent_des, "-")
    add("Lean solvent, stripper outlet", "CO2 mole fraction", CO2_app_rich_solvent_des, "-")

    # CO2 production (stripper vapor outlet)
    add("Stripper vapor output", "Temperature", v(model, :des_T, N_des; scenario=scenario, defaults=defaults), "K")
    add("Stripper vapor output", "Pressure", v(model, :des_P, N_des; scenario=scenario, defaults=defaults), "Pa")
    add("Stripper vapor output", "Molar flow", v(model, :des_V, N_des; scenario=scenario, defaults=defaults), "mol/s")
    add("Stripper vapor output", "MEA mole fraction", xlog_at(model, :des_log_y, N_des, :MEA; scenario=scenario, species_index=species_index, defaults=defaults), "-")
    add("Stripper vapor output", "H2O mole fraction", xlog_at(model, :des_log_y, N_des, :H2O; scenario=scenario, species_index=species_index, defaults=defaults), "-")
    add("Stripper vapor output", "CO2 mole fraction", xlog_at(model, :des_log_y, N_des, :CO2; scenario=scenario, species_index=species_index, defaults=defaults), "-")

    # Duties
    add("Reboiler", "Heat flow", v(model, :Q_reb; scenario=scenario, defaults=defaults), "kW")
    add("Condenser", "Heat flow", v(model, :Q_cond; scenario=scenario, defaults=defaults), "kW")
	
    # Loadings
    add("Rich solvent, absorber outlet", "CO2 loading", loading_liq(model, :abs_log_x, 1, [:CO2, :HCO3m, :CO3_2m, :Carbamate], [:MEA, :MEA2p, :Carbamate]; scenario=scenario, species_index=species_index, defaults=defaults), "mol/mol")
    add("Lean solvent, desorber outlet", "CO2 loading", loading_liq(model, :des_log_x, 1, [:CO2, :HCO3m, :CO3_2m, :Carbamate], [:MEA, :MEA2p, :Carbamate]; scenario=scenario, species_index=species_index, defaults=defaults), "mol/mol")

    # capture rate
    CO2_FG_out = v(model, :abs_V, N_abs; scenario=scenario, defaults=defaults) * xlog_at(model, :abs_log_y, N_abs, :CO2; scenario=scenario, species_index=species_index, defaults=defaults) # absorber CO2 output rate in vapor in mol/s
    CO2_FG_in = v(model, :abs_V_in; scenario=scenario, defaults=defaults) * v(model, :abs_y_in, species_index[:CO2]; scenario=scenario, defaults=defaults) # absorber CO2 input rate in vapor in mol/s
    capture_rate = (CO2_FG_in - CO2_FG_out) / CO2_FG_in * 100
    # println("CO2 capture rate, scenario $scenario: ", capture_rate, "%")
    add("CO2 capture", "Capture rate", capture_rate, "%")
    add("CO2 capture", "Captured CO2 annual", v(model, :annual_CO2_captured; scenario=scenario, defaults=defaults), "t/year")					  
    # Condenser output
    add("Condenser output", "Temperature", v(model, :cond_T, 1; scenario=scenario, defaults=defaults), "K")
    add("Condenser output", "Molar flow", v(model, :cond_V, 1; scenario=scenario, defaults=defaults), "mol/s")
    add("Condenser output", "MEA mole fraction", xlog_at(model, :cond_log_y, 1, :MEA; scenario=scenario, species_index=species_index, defaults=defaults), "-")
    add("Condenser output", "H2O mole fraction", xlog_at(model, :cond_log_y, 1, :H2O; scenario=scenario, species_index=species_index, defaults=defaults), "-")
    add("Condenser output", "CO2 mole fraction", xlog_at(model, :cond_log_y, 1, :CO2; scenario=scenario, species_index=species_index, defaults=defaults), "-")

    # Absorber temperatures
    for segment in 1:N_abs
        add("Absorber", "Segment $segment temperature", v(model, :abs_T, segment; scenario=scenario, defaults=defaults), "K")
    end
    # Absorber efficiency (approach to equilibrium)
    for segment in 1:N_abs
        add("Absorber", "Segment $segment efficiency", v(model, :abs_eff_CO2, segment; scenario=scenario, defaults=defaults), "-")
    end
    # Absorber approach to flooding
    for segment in 1:N_abs
        add("Absorber", "Segment $segment approach to flooding", v(model, :abs_approach_to_flooding, segment; scenario=scenario, defaults=defaults), "-")
    end
    # Stripper temperatures
    for segment in 1:N_des
        add("Stripper", "Segment $segment temperature", v(model, :des_T, segment; scenario=scenario, defaults=defaults), "K")
    end
    # Stripper efficiencies (approach to equilibrium)
    for segment in 1:N_des
        add("Stripper", "Segment $segment efficiency", v(model, :des_eff_CO2, segment; scenario=scenario, defaults=defaults), "-")
    end
    # Stripper approach to flooding
    for segment in 1:N_des
        add("Stripper", "Segment $segment approach to flooding", v(model, :des_approach_to_flooding, segment; scenario=scenario, defaults=defaults), "-")
    end

    # Reboiler temperature
    add("Reboiler", "Temperature", v(model, :reb_T, 1; scenario=scenario, defaults=defaults), "K")

    # Makeup
    add("Makeup H2O", "Molar flow", v(model, :MU_H2O; scenario=scenario, defaults=defaults), "mol/s")
    add("Makeup MEA", "Molar flow", v(model, :MU_MEA; scenario=scenario, defaults=defaults), "mol/s")

    # reboiler and cooler duties
    add("Reboiler", "Specific reboiler duty", v(model, :Q_reb; scenario=scenario, defaults=defaults) / ((CO2_FG_in - CO2_FG_out) * ExaProcessModels.MW[:CO2]), "GJ/t")

    # Absorber properties per segment
    if add_column_properties
        for segment in 1:N_abs
            add("Absorber", "Segment $segment inlet V/L", v(model, :abs_V_over_L_inlet, segment; scenario=scenario, defaults=defaults), "kg/kg")
            add("Absorber", "Segment $segment inlet alpha", v(model, :abs_alpha_inlet, segment; scenario=scenario, defaults=defaults), "molCO2/molMEA")
            add("Absorber", "Segment $segment inlet superficial velocity", v(model, :abs_u_v_inlet, segment; scenario=scenario, defaults=defaults), "m/s")
            add("Absorber", "Segment $segment liquid loading", v(model, :abs_L_load_superficial, segment; scenario=scenario, defaults=defaults), "m3/m2h")
            add("Absorber", "Segment $segment inlet temperature", v(model, :abs_T_inlet, segment; scenario=scenario, defaults=defaults), "K")
        end
    end
    
    # Desorber properties per segment
    if add_column_properties
        for segment in 1:N_abs
            add("Desorber", "Segment $segment inlet V/L", v(model, :des_V_over_L_inlet, segment; scenario=scenario, defaults=defaults), "kg/kg")
            add("Desorber", "Segment $segment inlet alpha", v(model, :des_alpha_inlet, segment; scenario=scenario, defaults=defaults), "molCO2/molMEA")
            add("Desorber", "Segment $segment inlet superficial velocity", v(model, :des_u_v_inlet, segment; scenario=scenario, defaults=defaults), "m/s")
            add("Desorber", "Segment $segment liquid loading", v(model, :des_L_load_superficial, segment; scenario=scenario, defaults=defaults), "m3/m2h")
            add("Desorber", "Segment $segment inlet temperature", v(model, :des_T_inlet, segment; scenario=scenario, defaults=defaults), "K")
        end
    end
	
    # design variables
    add("Absorber design", "Packing height", v(model, :abs_H; scenario=scenario, defaults=defaults), "m")
    add("Absorber design", "Packing diameter", v(model, :abs_D; scenario=scenario, defaults=defaults), "m")
    add("Desorber design", "Packing height", v(model, :des_H; scenario=scenario, defaults=defaults), "m")
    add("Desorber design", "Packing diameter", v(model, :des_D; scenario=scenario, defaults=defaults), "m")
    add("Heat exchanger design", "number of channels per pass", v(model, :NC; scenario=scenario, defaults=defaults), "-")
    add("Heat exchanger design", "Area", v(model, :HX_length; scenario=scenario, defaults=defaults)*v(model, :HX_width; scenario=scenario, defaults=defaults)*1.2 * ft_per_m^2 * ((2 * v(model, :NC; scenario=scenario, defaults=defaults) * 4 - (1 + 3))), "ft2")
    add("Heat exchanger design, reboiler", "Area", v(model, :A_reboiler_max; scenario=scenario, defaults=defaults), "ft2")
    add("Heat exchanger design, condenser", "Area", v(model, :A_condenser_max; scenario=scenario, defaults=defaults), "ft2")
    add("Heat exchanger design, cooler", "Area", v(model, :A_cooler_max; scenario=scenario, defaults=defaults), "ft2")
    add("Heat exchanger design, intercooler", "Area", v(model, :A_intercooler_max; scenario=scenario, defaults=defaults), "ft2")
    add("Pump design", "Maximum volumetric flow, lean", v(model, :pump_maxvflow_lean; scenario=scenario, defaults=defaults), "gpm")
    add("Pump design", "Maximum volumetric flow, rich", v(model, :pump_maxvflow_rich; scenario=scenario, defaults=defaults), "gpm")
    add("Compressor design", "Maximum power", v(model, :comp_P_max; scenario=scenario, defaults=defaults), "kW")

    # objective
    add("Objective", "TPC absorber", v(model, :BEC_absorber; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC absorber packing", v(model, :BEC_absorber_packing; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC desorber", v(model, :BEC_desorber; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC desorber packing", v(model, :BEC_desorber_packing; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "2 x TPC pump lean", 2*v(model, :BEC_pump_lean ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "2 x TPC pump rich", 2*v(model, :BEC_pump_rich ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC HEX", v(model, :BEC_hex ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC HEX intercooler", v(model, :BEC_intercooler ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC HEX cooler", v(model, :BEC_cooler ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC HEX condenser", v(model, :BEC_condenser ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC HEX reboiler", v(model, :BEC_reboiler ; scenario=scenario, defaults=defaults) * (1+epc_factor+(contingency_process + contingency_project)) / retrofit_factor, "USD")
    add("Objective", "TPC compressor", v(model, :BEC_compressor ; scenario=scenario, defaults=defaults)*(1+tpc_factor_compressor) / retrofit_factor, "USD")
    add("Objective", "TPC", v(model, :tpc_CCS; scenario=scenario, defaults=defaults)+v(model, :tpc_compressor; scenario=scenario, defaults=defaults), "USD")
    add("Objective", "steam consumption", v(model, :steam_consumption ; scenario=scenario, defaults=defaults), "kWh/year")
    add("Objective", "electricity consumption", v(model, :electricity_consumption ; scenario=scenario, defaults=defaults), "MWh/year")
    add("Objective", "cooling  consumption", v(model, :coolingwater_consumption ; scenario=scenario, defaults=defaults), "m3/year")
    add("Objective", "TOC", v(model, :toc ; scenario=scenario, defaults=defaults), "USD")
    add("Objective", "annualized CAPEX", CCF*sensitivity_CCF * v(model, :toc ; scenario=scenario, defaults=defaults), "USD")
    add("Objective", "annual OPEX", v(model, :annual_opex ; scenario=scenario, defaults=defaults), "USD")
    add("Objective", "Cost of CO2 avoided", v(model, :coa ; scenario=scenario, defaults=defaults), "USD/t")

    # create dataframe
    df = DataFrame(rows)

    return df
end

function collect_results_all_scenarios(model::JuMP.Model;
    scenarios::AbstractVector{Int},
    N_abs::Int,
    N_des::Int,
    species_index,
    defaults::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    add_column_properties::Bool=false,
    sensitivity_CCF::Float64
)
    dfs = DataFrame[]
    for s in scenarios
        df_s = collect_results(model;
            scenario=s,
            N_abs=N_abs,
            N_des=N_des,
            species_index=species_index,
            defaults=defaults,
            add_column_properties=add_column_properties,
            sensitivity_CCF=sensitivity_CCF
        )
        push!(dfs, df_s)
    end

    df = vcat(dfs...; cols=:union)  # cols=:union if some scenarios yield missings/extra rows

    return df
end
