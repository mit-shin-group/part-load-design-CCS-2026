#############################################################################
# The variables are declared in the order of the material flow since the 
# downstream units depends on upstream variables.
# Units after the recycle are considered downstream.
# References
# [1] Akula P, Eslick J, Bhattacharyya D, Miller DC. Model development, validation, and optimization of an MEA-based post-combustion CO2 capture process under part-load and variable capture operations. Ind Eng Chem Res 60(14):5176–5193 (2021).
# [2] Bogaard RH. Thermal conductivity of selected stainless steels. In: Thermal Conductivity 18. Eds: Ashworth T, Smith DR pp. 175–185 Springer US (1985).
# [3] Morgan JC, Chinen AS, Anderson-Cook C, Tong C, Carroll J, Saha C, Omell B, Bhattacharyya D, Matuszewski M, Bhat KS, Miller DC. Development of a Framework for Sequential Bayesian Design of Experiments: Application to a Pilot-Scale Solvent-Based CO2 Capture Process. (2021).
# [4] Sulzer Ltd. Structured Packings Enery-Efficient, Innovative & Profitable. (2024).
# [5] Martorell, J.L., Rochelle, G.T., Baldea, M., Elliott, W., Bauer, C., 2023. Lessons Learned: Comparing Two Detailed Capital Cost Estimates for Carbon Capture by Amine Scrubbing. Ind. Eng. Chem. Res. 62, 4433–4443.
# [6]  Kim S, Léonard G. CO2 capture technologies and shortcut cost correlations for different inlet CO2 concentrations and flow rates. part 1: chemical absorption. International Journal of Greenhouse Gas Control 145:104391 (2025).
# [7] Fisher, K.S., Searcy, K., Rochelle, G.T., Ziaii, S., Schubert, C., 2007. Advanced Amine Solvent Formulations and Process Integration for Near-Term CO2 Capture Success (No. DOE/ER/84625-1 Final Report, 945367).
# [8]  Fisher KS, Searcy K, Rochelle GT, Ziaii S, Schubert C. Advanced Amine Solvent Formulations and Process Integration for Near-Term CO2 Capture Success. (2007).
#############################################################################

function full_process_model(;recycle::Bool=false, as_optimization::Bool=false, N_abs::Int=10, N_des::Int=10, 
                            abs_P::Float64=110.0, des_P::Float64=200.0, case_study::String="coal", 
                            capture_rate::Float64=0.95, use_average_capture_rate::Bool=false, 
                            cases::AbstractVector{<:Integer}=[1], intercooler_active::Bool=false, 
                            fix_design::Bool=false, design::Dict{Symbol,Float64}=Dict{Symbol,Float64}(), 
                            export_log::Bool=false, 
                            sensitivity_flow::Float64=1.0, sensitivity_flh::Float64=1.0, sensitivity_CCF::Float64=1.0, sensitivity_price_reduction::Float64=0.0)
    """
    recycle: set to true to close recycle loop
    as_optimization: set to true to solve design optimization problem, else simulation based on initial values
    N_abs, N_des: number of segments for the discretization of the absorber and desorber columns
    abs_P, des_P: pressure in columns in kPa
    case_study: select "coal" or "gas"
    capture_rate: set the capture rate target
    use_average_capture_rate: set to true to set weighted average capture rate target, else the target is enforced in each scenario
    cases: list of part-load operation cases selected (see definition in casestudy_helpers.jl)
    intercooler_active: set to true if an intercooler is to be included in the center of the absorber column
    fix_design: set to true to fix the design of the capture plant
    design: Dict of all components and their respective sizes, specify if the design is fixed
    export_log: set to true to save the solver log
    sensitivity_flow: vary the flue gas mass flow, 1.0 corresponds to the base case, 1.1 corresponds to a 10% increase
    sensitivity_flh: vary the full load hours
    sensitivity_CCF: vary the capital charge factor
    sensitivity_price_reduction: imposes a linear correlation between electricity price and part-load operation 
    """

    # initialize model and solver settings
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    set_optimizer_attribute(model, "max_iter", 900)
    if LIBHSL_isfunctional()
        set_optimizer_attribute(model, "linear_solver", "ma57")
    end
    if export_log
        set_optimizer_attribute(model, "output_file", "ipopt_out.txt")
        set_optimizer_attribute(model, "file_print_level", 6)
    end

    # initialize variable and constraint lists
    design_vars_list = []  #for design/optimization variables
    operation_vars_list = []  #for scenario dependent variables
    operation_cons_list = []  #for scenario dependent constraints
    design_cons_list = []  #for design-dependent constraints

    # case-study-specific parameters
    if case_study == "gas"
        println("Gas case study: using fitted parameters for low-CO2 flue gases.")
        beta_abs_selected = beta[:abs_gas] 
        mu_abs_selected = mu[:abs_gas]
        sigma_abs_selected = sigma[:abs_gas]
        n_trains = 4
    elseif case_study == "coal"
        println("Coal case study: using fitted parameters for high-CO2 flue gases.")
        beta_abs_selected = beta[:abs_coal]
        mu_abs_selected = mu[:abs_coal]
        sigma_abs_selected = sigma[:abs_coal]
        n_trains = 2
    else
        error("Invalid case study. Please choose either 'gas' or 'coal'.")
    end

    # load flue gas inlet data and some inlet and outlet data for initialization
    inout_data_s = make_inout_data(cases; case_study=case_study)
    S = 1:length(inout_data_s)  # scenario index


    # check if the variation of the full load hours is valid for sensitivity analysis
    if sensitivity_flh != 1.0
        total_weight = sum(inout_data_s[s].weight for s in S)
        if total_weight * sensitivity_flh > 1.0
            error("The scaled total weight of scenarios exceeds 1.0. Reduce the sensitivity_flh parameter.")
        end
    end

    #############################################################################
    # Absorber variables
    #############################################################################
    # Get absorber data
    yV = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           inout_data_s[s].vapor.y[:CO2], inout_data_s[s].vapor.y[:H2O], inout_data_s[s].vapor.y[:N2], inout_data_s[s].vapor.y[:O2]]
           for s in S] # molar composition of feed
    @variable(model, abs_y_in[s in S, i=1:length(yV[1])] == yV[s][i])
    @variable(model, abs_V_in[s in S] == inout_data_s[s].vapor.flow*sensitivity_flow) # Vapor feed flowrate, mol/s
    @variable(model, abs_T_V_in[s in S] == inout_data_s[s].vapor.T) # Vapor feed temperature, K
    append!(operation_vars_list, (abs_y_in, abs_V_in, abs_T_V_in))

    abs_P_atm = abs_P / 101.325 # atm, pressure of the absorber
    solvent_composition = 0.3 # 30wt% MEA solvent

    # define variables
    @variables(model, begin
        abs_L[s in S, n=1:N_abs] >= 0           # Liquid flowrate
        abs_V[s in S, n=1:N_abs] >= 0           # Vapor flowrate
        293.15 <= abs_T[s in S, 1:N_abs] <= 450 # Temperature, K, temperature limited to avoid thermal degradation (Davis, J., Rochelle, G., 2009. Thermal degradation of monoethanolamine at stripper conditions. Energy Procedia 1, 327–333.)
        abs_U_s_flood[s in S, n=1:N_abs] >= 0   # Superficial velocity for flooding constraint, ft/s
        abs_Q_intercool[s in S, n=1:N_abs] <= 0 # Intercooling duty for each segment, kW
        0 >= abs_log_x[s in S, n=1:N_abs, c=1:length(species_index)] >= -50 #log transform liquid mole fraction  
        0 >= abs_log_y[s in S, n=1:N_abs, c=1:length(species_index)] >= -50 #log transform vapor mole fractions 
        abs_xi[s in S, n=1:N_abs, r=1:5]        # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (abs_L, abs_V, abs_T, abs_U_s_flood, abs_Q_intercool, abs_log_x, abs_log_y, abs_xi))

    # fix locations of intercooler
    N_intercool = round(Int, N_abs/2)  # assume intercooler is placed in the middle of the absorber
    cons_intercooler_T = @constraint(model, cons_intercooler_T[s in S], abs_T[s, N_intercool] <= 273.15 + inout_data_s[s].T_intercool)
    push!(operation_cons_list, cons_intercooler_T)
    for s in S
        for n in 1:N_abs
            if n != N_intercool # 1 intercooler
                fix(abs_Q_intercool[s, n], 0.0, force=true) # no intercooling duty in other segments
            end
        end
    end

    # Fix non-volatile species vapor fraction and inert species liquid fraction to exp(-100)
    fix_vapor(abs_log_y, S, N_abs;species=NONVOLATILES)
    fix_liquid(abs_log_x, S, N_abs;species=INERTS)
    abs_x = exp.(abs_log_x)
    abs_y = exp.(abs_log_y)


    #############################################################################
    # Heat exchanger variables
    #############################################################################
    # Parameters
    # length and width of the plate and frame heat exchanger, [1]
    # introduced as expression for easier post processing and scenario generation
    @expression(model, HX_length, 1.657) # length of heat exchanger, m
    @expression(model, HX_width, 0.849) # Width of heat exchanger, m
    thickness = 0.0006 # Thickness of plate, m
    plate_gap = 0.0038 # Plate gap, m
    corrugation = 1.2 # Corrugation or enlargement factor
    d_e = 2 * plate_gap / corrugation # Equivalent channel diameter, m
    # conductivity assuming AISI 304 SS, [2]
    Lambda_P = 15.71 # Conductivity of plate in W/(m*K), value for 350 K and assumed constant for range of 300-400K.
    A_p = HX_length * HX_width * corrugation # heat transfer area of per channel, m^2
    N_passes = 4 # Number of passes, assumed fixed, [1]

    # Variables
    # number of channels per pass
    if !recycle || !as_optimization
        @variable(model, NC == 500) # initial guess for simulation
    else
        if fix_design
            @variable(model, NC == design[:NC])
        else
            @variable(model, 6000 >= NC >= 10)
            # note: NC is relaxed to a continuous variable; for large NC, the error is small
        end
    end
    push!(design_vars_list, NC)  #number of channels per pass

    # inlet and outlet temperatures of hot and cold side
    @variable(model, 450 >= HX_T_H_passes_in[s in S, 1:N_passes] >= 0)
    @variable(model, 450 >= HX_T_H_passes_out[s in S, 1:N_passes] >= 0)
    @variable(model, 450 >= HX_T_C_passes_in[s in S, 1:N_passes] >= 0)
    @variable(model, 450 >= HX_T_C_passes_out[s in S, 1:N_passes] >= 0)
    append!(operation_vars_list, (HX_T_H_passes_in, HX_T_H_passes_out, HX_T_C_passes_in, HX_T_C_passes_out))

    # HX cold stream liquid-only isothermal flash variables
    @variables(model, begin
        HX_C_L[s in S, n=1:1] >= 0
        0 >= HX_C_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        HX_C_xi[s in S, n = 1:1, 1:5] # reaction extents  for reactions 1-5
    end)
    append!(operation_vars_list, (HX_C_L, HX_C_log_x, HX_C_xi))
    HX_C_x = exp.(HX_C_log_x)

    # Fix inert species liquid fraction to 0
    fix_liquid(HX_C_log_x, S, 1;species=INERTS)

    #############################################################################
    # Desorber variables
    #############################################################################
    des_P_atm = des_P / 101.325 # atm , pressure of the desorber
    @variables(model, begin
        des_L[s in S, n=1:N_des] >= 0               # Liquid flowrate
        des_V[s in S, n=1:N_des] >= 0               # Vapor flowrate
        293.15 <= des_T[s in S, n=1:N_des] <= 450   # Temperature, K
        des_U_s_flood[s in S, n=1:N_des] >= 0       # Superficial velocity for flooding constraint, ft/s
        0 >= des_log_x[s in S, n=1:N_des, c=1:length(species_index)] >= -50 # log liquid mole fraction
        0 >= des_log_y[s in S, 1:N_des, 1:length(species_index)] >= -50     # log vapor mole fraction
        des_xi[s in S, n=1:N_des, 1:5]              # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (des_L, des_V, des_T, des_U_s_flood, des_log_x, des_log_y, des_xi))
    des_x = exp.(des_log_x)
    des_y = exp.(des_log_y)

    # Fix non-volatile and inert species vapor fraction and inert species liquid fraction to 0
    fix_vapor(des_log_y, S, N_des;species=[NONVOLATILES; INERTS])
    fix_liquid(des_log_x, S, N_des;species=INERTS)

    #############################################################################
    # Reboiler variables
    #############################################################################
    P_reb = des_P_atm   # Total pressure, atm
    @variables(model, begin
        reb_L[s in S, n = 1:1] >= 0             # Liquid flowrate
        reb_V[s in S, n = 1:1] >= 0             # Vapor flowrate
        450 >= reb_T[s in S, n = 1:1] >= 293.15 # Temperature, K
        0 >= reb_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        0 >= reb_log_y[s in S, n=1:1, c=1:length(species_index)] >= -50
        reb_xi[s in S,n = 1:1, r=1:5]           # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (reb_L, reb_V, reb_T, reb_log_x, reb_log_y, reb_xi))
    reb_x = exp.(reb_log_x)
    reb_y = exp.(reb_log_y)

    # Fix non-volatile and inert species vapor fraction and inert species liquid fraction to 0
    fix_vapor(reb_log_y, S, 1;species=[NONVOLATILES; INERTS])
    fix_liquid(reb_log_x, S, 1;species=INERTS)

    #############################################################################
    # Condenser variables
    #############################################################################
    P_cond = des_P_atm   # Total pressure, atm 
    @variables(model, begin
        cond_L[s in S, n = 1:1] >= 0                # Liquid flowrate
        cond_V[s in S, n = 1:1] >= 0                # Vapor flowrate
        293.15 <= cond_T[s in S, n = 1:1] <= 450    # Temperature, K
        0 >= cond_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        0 >= cond_log_y[s in S, n=1:1, c=1:length(species_index)] >= -50
        cond_xi[s in S, n = 1:1, r=1:5]             # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (cond_L, cond_V, cond_T, cond_log_x, cond_log_y, cond_xi))
    cond_x = exp.(cond_log_x)
    cond_y = exp.(cond_log_y)

    for s in S
        fix(cond_T[s,1], 40+273.15; force=true) # fix condenser temperature, [3]
    end

    # Fix non-volatile and inert species vapor fraction and inert species liquid fraction to 0
    fix_vapor(cond_log_y, S, 1;species=[NONVOLATILES; INERTS])
    fix_liquid(cond_log_x, S, 1;species=INERTS)

    #############################################################################
    # HX hot stream liquid-only isothermal flash variables
    #############################################################################
    @variables(model, begin
        HX_H_L[s in S, n = 1:1] >= 0
        0 >= HX_H_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        HX_H_xi[s in S, n = 1:1, r=1:5] # reaction extents
    end)
    append!(operation_vars_list, (HX_H_L, HX_H_log_x, HX_H_xi))
    HX_H_x = exp.(HX_H_log_x)

    # Fix inert species liquid fraction to 0
    fix_liquid(HX_H_log_x, S, 1;species=INERTS)

    #############################################################################
    # Cooler variables
    #############################################################################
    @variables(model, begin
        cool_L[s in S, n = 1:1] >= 0                # Liquid flowrate
        293.15 <= cool_T[s in S, n = 1:1] <= 450    # Temperature, K
        0 >= cool_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        cool_xi[s in S, n = 1:1, r=1:5]             # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (cool_L, cool_T, cool_log_x, cool_xi))
    cool_x = exp.(cool_log_x)

    # Fix inert species liquid fraction to 0
    fix_liquid(cool_log_x, S, 1;species=INERTS)

    #############################################################################
    # Mixer variables
    #############################################################################
    @variables(model, begin
        mix_L[s in S, n = 1:1] >= 0             # Liquid flowrate
        293.15 <= mix_T[s in S, n = 1:1] <= 450 # Temperature, K
        0 >= mix_log_x[s in S, n=1:1, c=1:length(species_index)] >= -50
        mix_xi[s in S, n = 1:1, r=1:5]          # Extent of reaction for reaction 1-5
    end)
    append!(operation_vars_list, (mix_L, mix_T, mix_log_x, mix_xi))
    mix_x = exp.(mix_log_x)

    # Fix inert species liquid fraction to 0
    fix_liquid(mix_log_x, S, 1;species=INERTS)

    #############################################################################
    # Compressor variables
    #############################################################################
    if fix_design
        @variable(model, comp_P_max == design[:comp_P_max])
    else
        @variable(model, comp_P_max >= 0) # maximum compressor power, kW)
    end
    push!(design_vars_list, comp_P_max)
    set_start_value(comp_P_max, 5100)


    #############################################################################
    # Stream connection, flow expression and variable definition
    #############################################################################
    #############################################################################
    # Liquid feed to absorber 
    #############################################################################
    L_in_s = [calc_inlet_species(inout_data_s[s])[1] for s in S]
    xL_s = [calc_inlet_species(inout_data_s[s])[2] for s in S]
    xL = [[x[:MEA], x[:MEA2p], x[:Hp], x[:Carbamate], x[:HCO3m], x[:OHm], x[:CO3_2m], x[:CO2], x[:H2O], 0.0, 0.0] for x in xL_s]  # last two entries for insoluble N2 and O2

    @variable(model, MU_H2O[s in S] <= inout_data_s[s].liquid_inlet.flow) # make-up water flowrate, mol/s
    # note: we allow negative makeup values in case of water buildup. However, water buildup does not seem to be an issue. 
    @variable(model, MU_MEA[s in S] <= inout_data_s[s].liquid_inlet.flow) # make-up MEA flowrate, mol/s
    append!(operation_vars_list, (MU_H2O, MU_MEA))

    for s in S
        set_start_value(MU_H2O[s], inout_data_s[s].initial_values.MU_H2O)
        set_start_value(MU_MEA[s], 0.0)
    end

    @variable(model, 1.5 * L_in_s[s] >= abs_L_in[s in S] >= 0.5 * L_in_s[s]) # Liquid feed flowrate, mol/s
    @variable(model, 1 >= abs_x_in[s in S, c=1:length(species_index)] >= 0) # Liquid feed mole fractions
    @variable(model, 50+273.15 >= abs_T_L_in[s in S] >= 20+273.15) # Liquid feed temperature, constrainted to 20-50C
    append!(operation_vars_list, (abs_L_in, abs_x_in, abs_T_L_in))

    # Make cooling duty, reboiler duty, and condenser duty decision variables and enforce liquid inlet input
    @variable(model, Q_cool[s in S] <= 0) # kW
    @variable(model, Q_reb[s in S] >= 0) # kW
    @variable(model, Q_cond[s in S] <= 0) # kW
    append!(operation_vars_list, (Q_cool, Q_reb, Q_cond))

    if recycle
        if as_optimization
            # fix solvent composition to 30 wt%
            cons_solvent_composition = @constraint(model, cons_solvent_composition[s in S], (abs_x_in[s, species_index[:MEA]] + abs_x_in[s, species_index[:MEA2p]] + abs_x_in[s, species_index[:Carbamate]]) * MW[:MEA] == solvent_composition * (MW[:H2O] * (abs_x_in[s, species_index[:H2O]] + abs_x_in[s, species_index[:Hp]] + abs_x_in[s, species_index[:OHm]] + abs_x_in[s, species_index[:HCO3m]]) + MW[:MEA] * (abs_x_in[s, species_index[:MEA]] + abs_x_in[s, species_index[:MEA2p]] + abs_x_in[s, species_index[:Carbamate]]))) # 30 wt% MEA
            push!(operation_cons_list, cons_solvent_composition)
        else # !as_optimization==true
            # fix lean solvent CO2 loading to initial value provided in scenario definition
            cons_solvent_composition_MEA = @constraint(model, cons_solvent_composition_MEA[s in S], (abs_x_in[s, species_index[:MEA]] + abs_x_in[s, species_index[:MEA2p]] + abs_x_in[s, species_index[:Carbamate]]) == inout_data_s[s].liquid_inlet.y[:MEA])
            cons_solvent_composition_CO2 = @constraint(model, cons_solvent_composition_CO2[s in S], (abs_x_in[s, species_index[:CO2]] + abs_x_in[s, species_index[:HCO3m]] + abs_x_in[s, species_index[:CO3_2m]] + abs_x_in[s, species_index[:Carbamate]]) == inout_data_s[s].liquid_inlet.y[:CO2])
            append!(operation_cons_list, (cons_solvent_composition_MEA, cons_solvent_composition_CO2))
            # fix inlet flowrate and temperature for simulation
            for s in S
                fix(abs_L_in[s], L_in_s[s]; force=true) # force the inlet flowrate
                fix(abs_T_L_in[s], inout_data_s[s].liquid_inlet.T; force=true)
            end
        end
    else # open-loop process system model
        for s in S
            # fix variables to initial values         
            fix(abs_L_in[s], L_in_s[s]; force=true)
            fix(abs_T_L_in[s], inout_data_s[s].liquid_inlet.T; force=true)
            for c in 1:length(species_index)
                fix(abs_x_in[s, c], xL[s][c]; force=true)
            end
            # fix cooler and reboiler duties, and make ups to initial values
            fix(Q_cool[s], inout_data_s[s].initial_values.Q_cool; force=true)
            fix(Q_reb[s], inout_data_s[s].Q_reboiler; force=true)
            fix(MU_H2O[s], inout_data_s[s].initial_values.MU_H2O; force=true)
            fix(MU_MEA[s], 0; force=true)
        end
    end

    #############################################################################
    # Absorber to HX_C
    # solvent exiting at absorber bottom is the cold fluid entering the HX
    #############################################################################
    F_C = [abs_L[s, 1] for s in S] #
    y_C = [[abs_x[s, 1, species_index[:MEA]], abs_x[s, 1, species_index[:MEA2p]], abs_x[s, 1, species_index[:Hp]], abs_x[s, 1, species_index[:Carbamate]],
        abs_x[s, 1, species_index[:HCO3m]], abs_x[s, 1, species_index[:OHm]], abs_x[s, 1, species_index[:CO3_2m]], abs_x[s, 1, species_index[:CO2]],
        abs_x[s, 1, species_index[:H2O]], abs_x[s, 1, species_index[:N2]], abs_x[s, 1, species_index[:O2]]] for s in S]
    @expression(model, T_C_IN[s in S], abs_T[s, 1]) # inlet temperature of cold side of heat exchanger

    #############################################################################
    # Desorber stream connections (heat exchanger, condenser, and reboiler)
    #############################################################################
    # Desorber liquid feed from HX_C
    des_L_in = HX_C_L
    des_T_L_in = [HX_T_C_passes_out[s, 1] for s in S]
    des_log_x_in = [[HX_C_log_x[s, 1, species_index[:MEA]], HX_C_log_x[s, 1, species_index[:MEA2p]], HX_C_log_x[s, 1, species_index[:Hp]], HX_C_log_x[s, 1, species_index[:Carbamate]], HX_C_log_x[s, 1, species_index[:HCO3m]], HX_C_log_x[s, 1, species_index[:OHm]], HX_C_log_x[s, 1, species_index[:CO3_2m]], HX_C_log_x[s, 1, species_index[:CO2]], HX_C_log_x[s, 1, species_index[:H2O]], HX_C_log_x[s, 1, species_index[:N2]], HX_C_log_x[s, 1, species_index[:O2]]] for s in S]
    des_x_in = [exp.(des_log_x_in[s]) for s in S]
    # Desorber liquid feed from condenser 
    x_cond = [[cond_x[s, 1, species_index[:MEA]], cond_x[s, 1, species_index[:MEA2p]], cond_x[s, 1, species_index[:Hp]], cond_x[s, 1, species_index[:Carbamate]], cond_x[s, 1, species_index[:HCO3m]], cond_x[s, 1, species_index[:OHm]], cond_x[s, 1, species_index[:CO3_2m]], cond_x[s, 1, species_index[:CO2]], cond_x[s, 1, species_index[:H2O]], cond_x[s, 1, species_index[:N2]], cond_x[s, 1, species_index[:O2]]] for s in S]
    # Desorber vapor feed from reboiler
    des_V_in = reb_V
    des_T_V_in = reb_T
    @expression(model, des_y_in[s in S, c=1:length(species_index)], reb_y[s, 1, c]) # Use the vapor output of the reboiler as the desorber vapor feed

    #############################################################################
    # Desorber to reboiler
    # liquid leaving the bottom segment of desorber is the liquid entering the reboiler
    #############################################################################
    reb_log_x_in = [[des_log_x[s, 1, species_index[:MEA]], des_log_x[s, 1, species_index[:MEA2p]], des_log_x[s, 1, species_index[:Hp]], des_log_x[s, 1, species_index[:Carbamate]], des_log_x[s, 1, species_index[:HCO3m]], des_log_x[s, 1, species_index[:OHm]], des_log_x[s, 1, species_index[:CO3_2m]], des_log_x[s, 1, species_index[:CO2]], des_log_x[s, 1, species_index[:H2O]], des_log_x[s, 1, species_index[:N2]], des_log_x[s, 1, species_index[:O2]]] for s in S]
    reb_x_in = [exp.(reb_log_x_in[s]) for s in S]
    reb_T_in = Dict()
    reb_L_in = Dict()
    for s in S
        reb_T_in[s] = des_T[s, 1] # Liquid feed temperature, K
        reb_L_in[s] = des_L[s, 1] # Liquid feed flowrate, mol/s
    end

    #############################################################################
    # Desorber vapor to condenser 
    # vapor leaving the top segment of desorber is the vapor entering the condenser
    #############################################################################
    cond_y_in = [[des_y[s, N_des, species_index[:MEA]], des_y[s, N_des, species_index[:MEA2p]], des_y[s, N_des, species_index[:Hp]], des_y[s, N_des, species_index[:Carbamate]],
        des_y[s, N_des, species_index[:HCO3m]], des_y[s, N_des, species_index[:OHm]], des_y[s, N_des, species_index[:CO3_2m]], des_y[s, N_des, species_index[:CO2]], des_y[s, N_des, species_index[:H2O]], des_y[s, N_des, species_index[:N2]], des_y[s, N_des, species_index[:O2]]] for s in S] # Use the vapor output of the desorber as the condenser feed
    cond_T_in = Dict()
    cond_V_in = Dict()
    for s in S
        cond_T_in[s] = des_T[s, N_des] # Vapor feed temperature, K
        cond_V_in[s] = des_V[s, N_des] # Vapor feed flowrate, mol/s
    end

    #############################################################################
    # Reboiler liquid to hot side heat exchanger
    #############################################################################
    F_H = [reb_L[s, 1] for s in S]
    y_H = [[reb_x[s, 1, species_index[:MEA]], reb_x[s, 1, species_index[:MEA2p]], reb_x[s, 1, species_index[:Hp]], reb_x[s, 1, species_index[:Carbamate]], reb_x[s, 1, species_index[:HCO3m]], reb_x[s, 1, species_index[:OHm]], reb_x[s, 1, species_index[:CO3_2m]], reb_x[s, 1, species_index[:CO2]], reb_x[s, 1, species_index[:H2O]], reb_x[s, 1, species_index[:N2]], reb_x[s, 1, species_index[:O2]]] for s in S]
    @expression(model, HX_T_H_IN[s in S], reb_T[s, 1])

    #############################################################################
    # Hot side heat exchanger to cooler
    #############################################################################
    cool_T_in = [HX_T_H_passes_out[s, N_passes] for s in S]
    cool_L_in = HX_H_L
    cool_x_in = [[HX_H_x[s, 1, species_index[:MEA]], HX_H_x[s, 1, species_index[:MEA2p]], HX_H_x[s, 1, species_index[:Hp]], HX_H_x[s, 1, species_index[:Carbamate]], HX_H_x[s, 1, species_index[:HCO3m]], HX_H_x[s, 1, species_index[:OHm]], HX_H_x[s, 1, species_index[:CO3_2m]], HX_H_x[s, 1, species_index[:CO2]], HX_H_x[s, 1, species_index[:H2O]], HX_H_x[s, 1, species_index[:N2]], HX_H_x[s, 1, species_index[:O2]]] for s in S]

    #############################################################################
    # Model constraints
    #############################################################################
    #############################################################################
    # Absorber
    #############################################################################
    # Species composition in solvent and vapor inlet
    @expression(model, abs_x_Li[s in S, c=1:length(species_index)], abs_L_in[s] * abs_x_in[s, c])
    @expression(model, abs_y_Vi[s in S, c=1:length(species_index)], abs_V_in[s] * abs_y_in[s, c])

    # Mole balance constraints
    if N_abs == 1 # one segment
        cons_abs_MB = add_mole_balance!(model, abs_x, abs_L, abs_y, abs_V, abs_x_Li, abs_y_Vi, abs_xi, S; N_seg=N_abs)
        push!(operation_cons_list, cons_abs_MB)
    elseif N_abs == 2 # two segments
        cons_abs_MB_1, cons_abs_MB_N = add_mole_balance!(model, abs_x, abs_L, abs_y, abs_V, abs_x_Li, abs_y_Vi, abs_xi, S; N_seg=N_abs)
        append!(operation_cons_list, [cons_abs_MB_1,cons_abs_MB_N])
    else # >2 segments
        cons_abs_MB_1, cons_abs_MB_N, cons_abs_MB_n = add_mole_balance!(model, abs_x, abs_L, abs_y, abs_V, abs_x_Li, abs_y_Vi, abs_xi, S; N_seg=N_abs)
        append!(operation_cons_list, [cons_abs_MB_1, cons_abs_MB_n, cons_abs_MB_N])
    end

    # N2 and O2 vapor mole balance constraints, assume no dissolution in solvent
    cons_abs_MB_N2 = @constraint(model, cons_abs_MB_N2[s in S, n=1:N_abs], abs_y[s, n, species_index[:N2]] * abs_V[s, n] == abs_y_Vi[s, species_index[:N2]])
    cons_abs_MB_O2 = @constraint(model, cons_abs_MB_O2[s in S, n=1:N_abs], abs_y[s, n, species_index[:O2]] * abs_V[s, n] == abs_y_Vi[s, species_index[:O2]])
    append!(operation_cons_list, [cons_abs_MB_N2, cons_abs_MB_O2])

    # Flooding and approach-to-equilibrium
    # design variables
    Height_init = 15 # initial value height of the column, m
    Diameter_init = 11 # initial value diameter of the column, m
    if !recycle && !as_optimization # fix height and diameter to initial values
        @variable(model, abs_H == Height_init)
        @variable(model, abs_D == Diameter_init)
    else
        if fix_design
            @variable(model, abs_H == design[:abs_H])
            @variable(model, abs_D == design[:abs_D])
        else
            @variable(model, 0.25*Height_init <= abs_H <= 3*Height_init)
            @variable(model, 0.25*Diameter_init <= abs_D <= 3*Diameter_init)
            set_start_value(abs_H, Height_init)
            set_start_value(abs_D, Diameter_init)
        end
    end
    append!(design_vars_list, [abs_H, abs_D])

    # Parameters/variables for flooding and approach-to-equilibrium model
    (abs_app_L_MEA, abs_app_L_H2O, abs_app_L_CO2) = compute_apparent_liquid(model, abs_x, S, N_abs) # apparent species amount in column solvent flow
    # CO2 loading in each segment, mol/mol
    @expression(model, abs_alpha[s in S, n in 1:N_abs], abs_app_L_CO2[s, n] / abs_app_L_MEA[s, n])
    # vapor superficial velocity in each segment, m/s
    @expression(model, abs_u_v[s in S, n in 1:N_abs], abs_V[s, n] * (R*abs_T[s,n]/(abs_P_atm*101325)) / (abs_D^2/4 * π))
    # Solvent and vapor average molar weights for every segment, kg/mol
    abs_x_MW_avg = [sum(abs_x[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_abs]
    abs_y_MW_avg = [sum(abs_y[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_abs]
    # vapor density in every segment, kg/m^3
    @expression(model, abs_rho_g[s in S, n in 1:N_abs], calculate_ideal_gas_density(abs_P_atm, abs_y_MW_avg[s,n], abs_T[s, n]))
    # liquid kinematic viscosity in every segment, m2/(s) / 10^6 = cSt
    @expression(model, abs_kinematic_viscosity[s in S, n in 1:N_abs], calc_viscosity(abs_T[s, n], [abs_app_L_MEA[s,n], abs_app_L_H2O[s,n]], abs_alpha[s,n])/density*1e6 )
    # Inlet vapor composition and flowrate  for each segment (from segment below or feed)
    @expression(model, abs_y_inlet[s in S, n in 1:N_abs, c in 1:length(species_index)], n == 1 ? abs_y_in[s, c] : abs_y[s, n-1, c])
    @expression(model, abs_V_inlet[s in S, n in 1:N_abs], n == 1 ? abs_V_in[s] : abs_V[s, n-1])
    # Inlet liquid composition and flowrate for each segment (from segment above or feed)
    @expression(model, abs_x_inlet[s in S, n in 1:N_abs, c in 1:length(species_index)], n == N_abs ? abs_x_in[s, c] : abs_x[s, n+1, c])
    @expression(model, abs_L_inlet[s in S, n in 1:N_abs], n == N_abs ? abs_L_in[s] : abs_L[s, n+1])
    # Inlet solvent composition without CO2 for each segment (from segment above or feed)
    @expression(model, abs_x_MEA_composition_in[s in S, n in 1:N_abs], n == N_abs ? solvent_composition : abs_app_L_MEA[s, n]*MW[:MEA] / (abs_app_L_MEA[s, n]*MW[:MEA] + abs_app_L_H2O[s, n]*MW[:H2O]))
    # Inlet apparent species mole fraction
    @expression(model, abs_app_L_MEA_inlet[s in S, n in 1:N_abs], abs_x_inlet[s, n, species_index[:MEA]] + abs_x_inlet[s, n, species_index[:MEA2p]] + abs_x_inlet[s, n, species_index[:Carbamate]])
    @expression(model, abs_app_L_CO2_inlet[s in S, n in 1:N_abs], abs_x_inlet[s, n, species_index[:CO2]] + abs_x_inlet[s, n, species_index[:HCO3m]] + abs_x_inlet[s, n, species_index[:CO3_2m]] + abs_x_inlet[s, n, species_index[:Carbamate]])
    # Inlet solvent CO2 loading, mol/mol
    @expression(model, abs_alpha_inlet[s in S, n in 1:N_abs], abs_app_L_CO2_inlet[s, n] / abs_app_L_MEA_inlet[s, n])
    # Inlet solvent and vapor average molar weights, kg/mol
    abs_x_MW_avg_inlet = [sum(abs_x_inlet[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_abs]
    abs_y_MW_avg_inlet = [sum(abs_y_inlet[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_abs]
    # Inlet vaport-to-liquid ratio, kg/kg
    @expression(model, abs_V_over_L_inlet[s in S, n in 1:N_abs], (abs_V_inlet[s, n] * abs_y_MW_avg_inlet[s, n]) / (abs_L_inlet[s, n] * abs_x_MW_avg_inlet[s, n]))
    # Inlet superficial velocity, m/s
    @expression(model, abs_u_v_inlet[s in S, n in 1:N_abs], abs_V_inlet[s, n] * (R*abs_T[s,n]/(abs_P_atm*101325)) / (abs_D^2/4 * π))
    # Inlet temperature in each segment, K
    @expression(model, abs_T_inlet[s in S, n in 1:N_abs], n == N_abs ? abs_T_L_in[s] : abs_T[s, n+1])

    # Diameter constraints to avoid flooding
    @expression(model, abs_flow_parameter[s in S, n in 1:N_abs], calculate_flow_parameter(abs_L[s,n]*abs_x_MW_avg[s,n], abs_V[s,n]*abs_y_MW_avg[s,n], abs_rho_g[s,n], density)) # vapor flow parameter
    @expression(model, abs_approach_to_flooding[s in S, n in 1:N_abs], abs_u_v[s, n]*ft_per_m/abs_U_s_flood[s, n]) # dimensionless, ratio of superficial velocity to flooding velocity
    # constraints coupling GPDC flooding correlation with superficial velocity flooding correlation to solve for diameter
    cons_abs_CP_at_flooding = @constraint(model, cons_abs_CP_at_flooding[s in S, n=1:N_abs], calculate_CP_GPDC_flood(abs_flow_parameter[s, n]) == calculate_CP_superficial_velocity(abs_U_s_flood[s, n], abs_rho_g[s, n], density, abs_kinematic_viscosity[s, n]))
    cons_abs_superficial_velocity_limit = @constraint(model, cons_abs_superficial_velocity_limit[s in S, n=1:N_abs], abs_u_v[s, n]*3.28084 <= 0.8 * abs_U_s_flood[s, n]) # superficial velocity is maximum 80% of the flooding velocity, note that abs_U_s_flood is in ft/s while u_v is in m/s, 3.28084 ft/m
    append!(operation_cons_list, [cons_abs_CP_at_flooding, cons_abs_superficial_velocity_limit])

    # minimum liquid load, [4]
    @expression(model, abs_L_load_superficial[s in S, n in 1:N_abs], abs_L[s,n] * abs_x_MW_avg[s,n] / density / (abs_D^2/4 * π) * 3600) # m3/m2h
    cons_abs_min_liquid_load = @constraint(model, cons_abs_min_liquid_load[s in S, n=1:N_abs], abs_L_load_superficial[s,n] >= minimum_liquid_load)
    push!(operation_cons_list, cons_abs_min_liquid_load)

    # VLE constraints
    if N_abs > 1
        (cons_abs_raoult_H2O, cons_abs_Henry_CO2_1, cons_abs_Henry_CO2_n) = add_VLE!(
            model, abs_x, abs_y, abs_T, S, N_abs, density, abs_P_atm, abs_y_in, true; 
            beta=beta_abs_selected, 
            mu=mu_abs_selected,
            sigma=sigma_abs_selected,
            V_over_L=abs_V_over_L_inlet, 
            T_in=abs_T_inlet,
            alpha=abs_alpha_inlet, 
            x_MEA_composition_in=abs_x_MEA_composition_in, 
            u_v=abs_u_v_inlet, 
            dz=abs_H/N_abs,
            name="abs")
        append!(operation_cons_list, [cons_abs_raoult_H2O, cons_abs_Henry_CO2_1, cons_abs_Henry_CO2_n])
    else
        (cons_abs_raoult_H2O, cons_abs_Henry_CO2_1) = add_VLE!(
            model, abs_x, abs_y, abs_T, S, N_abs, density, abs_P_atm, abs_y_in, false; 
            beta=nothing, 
            mu=nothing,
            sigma=nothing,
            V_over_L=nothing,
            T_in = nothing,
            alpha=nothing, 
            x_MEA_composition_in=nothing, 
            u_v=nothing, 
            dz=nothing,
            name="abs")
        append!(operation_cons_list, [cons_abs_raoult_H2O, cons_abs_Henry_CO2_1])
    end

    # Reaction equilibrium constraints
    (cons_abs_rxn1, cons_abs_rxn2, cons_abs_rxn3, cons_abs_rxn4, cons_abs_rxn5) = add_reaction_equilibria!(model, abs_log_x, abs_T, S; N_seg=N_abs)
    append!(operation_cons_list, [cons_abs_rxn1, cons_abs_rxn2, cons_abs_rxn3, cons_abs_rxn4, cons_abs_rxn5])
    # Summation constraint (mole fraction)
    cons_abs_sum_liq = @constraint(model, cons_abs_sum_liq[s in S, n=1:N_abs], sum(abs_x[s, n, c] for c in 1:length(species_index)) == 1)
    cons_abs_sum_vap = @constraint(model, cons_abs_sum_vap[s in S, n=1:N_abs], sum(abs_y[s, n, c] for c in 1:length(species_index)) == 1)
    append!(operation_cons_list, [cons_abs_sum_liq, cons_abs_sum_vap])
    # Liquid and vapor enthalpy
    abs_app_Li = compute_apparent_liquid(model, abs_x_Li, S) # apparent species amount in liquid feed   
    (abs_H_liq_in, abs_H_liq_out) = compute_liquid_enthalpy(model, abs_app_Li, abs_T_L_in, abs_x, abs_T, abs_L, S, N_abs)
    (abs_H_vap_in, abs_H_vap_out) = compute_vapor_enthalpy(model, abs_y_Vi, abs_T_V_in, abs_y, abs_T, abs_V, S, N_abs)
    # Energy balance constraints
    if N_abs == 1
        cons_abs_EB = add_energy_balance!(model, abs_H_liq_in, abs_H_liq_out, abs_H_vap_in, abs_H_vap_out, nothing, S; N_seg=N_abs)
        push!(operation_cons_list, cons_abs_EB)
    elseif N_abs == 2
        cons_abs_EB_1, cons_abs_EB_N = add_energy_balance!(model, abs_H_liq_in, abs_H_liq_out, abs_H_vap_in, abs_H_vap_out, nothing , S; N_seg=N_abs)
        append!(operation_cons_list, [cons_abs_EB_1, cons_abs_EB_N])
    else
        cons_abs_EB_1, cons_abs_EB_N, cons_abs_EB_n = add_energy_balance!(model, abs_H_liq_in, abs_H_liq_out, abs_H_vap_in, abs_H_vap_out, abs_Q_intercool*1000, S; N_seg=N_abs)
        append!(operation_cons_list, [cons_abs_EB_1, cons_abs_EB_n, cons_abs_EB_N])
    end

    # Initialization
    # Use fixed values or initial values of entering streams as initial values
    abs_L0 = [L_in_s[s] for s in S] # liquid flow to absorber
    abs_V0 = [check_val(abs_V_in[s]) for s in S] # vapor flow to absorber
    abs_TV0 = [check_val(abs_T_V_in[s]) for s in S] # vapor T to absorber
    init_LVT(abs_L, abs_V, abs_T, abs_L0, abs_V0, abs_TV0, S, N_abs)

    abs_x0 = xL # liquid mole fraction to absorber
    abs_y0 = [check_val.(abs_y_in[s, :]) for s in S] # vapor mole fraction to absorber
    abs_log_x0 = [log(abs_x0[s][c]) for s in 1:length(S), n in 1:N_abs, c in 1:length(species_index)]
    abs_log_y0 = [log(abs_y0[s][c]) for s in 1:length(S), n in 1:N_abs, c in 1:length(species_index)]
    init_mole_frac(abs_log_x, abs_log_y, abs_log_x0, abs_log_y0, S, N_abs)

    for s in S
        for n in 1:N_abs
            if n == N_intercool
                set_start_value(abs_Q_intercool[s, n], inout_data_s[s].initial_values.abs_Q_intercool) # initialize intercooling
            end
                set_start_value(abs_U_s_flood[s,n], 10.0)
        end
    end

    # rough initial guesses for abs_xi
    abs_xi_0 = Dict(
        1 => -15,
        2 => 0.03,
        3 => -0.001,
        4 => -10,
        5 => 15
    )
    for s in S
        for n in 1:N_abs
            for r in 1:5  #5 reactions
                set_start_value(abs_xi[s, n, r], abs_xi_0[r])
            end
        end
    end

    #############################################################################
    # Heat exchanger constraints
    #############################################################################
    # CO2 loading in cold and hot streams
    @expression(model, alpha_H[s in S], (reb_x[s, 1, species_index[:CO2]] + reb_x[s, 1, species_index[:CO3_2m]] + reb_x[s, 1, species_index[:HCO3m]] + reb_x[s, 1, species_index[:Carbamate]]) / (reb_x[s, 1, species_index[:MEA]] + reb_x[s, 1, species_index[:MEA2p]] + reb_x[s, 1, species_index[:Carbamate]]))
    @expression(model, alpha_C[s in S], (abs_x[s, 1, species_index[:CO2]] + abs_x[s, 1, species_index[:CO3_2m]] + abs_x[s, 1, species_index[:HCO3m]] + abs_x[s, 1, species_index[:Carbamate]]) / (abs_x[s, 1, species_index[:MEA]] + abs_x[s, 1, species_index[:MEA2p]] + abs_x[s, 1, species_index[:Carbamate]]))
    # Apparent species mole fractions in liquid feeds
    (mole_frac_H_MEA, mole_frac_H_H2O, mole_frac_H_CO2) = compute_apparent_liquid(model, reb_x, S, 1)
    (mole_frac_C_MEA, mole_frac_C_H2O, mole_frac_C_CO2) = compute_apparent_liquid(model, abs_x, S, 1)
    @expression(model, Tavg[s in S], 0.5 * (HX_T_H_IN[s] + T_C_IN[s]))  # K
    # Average molecular weight of hot and cold fluids
    @expression(model, MW_H[s in S], MW[:MEA] * mole_frac_H_MEA[s,1] + MW[:H2O] * mole_frac_H_H2O[s,1] + MW[:CO2] * mole_frac_H_CO2[s,1])
    @expression(model, MW_C[s in S], MW[:MEA] * mole_frac_C_MEA[s,1] + MW[:H2O] * mole_frac_C_H2O[s,1] + MW[:CO2] * mole_frac_C_CO2[s,1])
    # Heat exchanger correlations
    # Heat capacity for each pass (J/mol·K)
    @expression(model, Cp_H[s in S, p in 1:N_passes], calc_Cp([mole_frac_H_MEA[s,1], mole_frac_H_H2O[s,1]], HX_T_H_passes_in[s, p]))
    @expression(model, Cp_C[s in S, p in 1:N_passes], calc_Cp([mole_frac_C_MEA[s,1], mole_frac_C_H2O[s,1]], HX_T_C_passes_in[s, p]))
    # Capacitance rate, J/(s*K)
    @expression(model, C_rate_H[s in S, p in 1:N_passes], F_H[s,1] * Cp_H[s, p] / NC)
    @expression(model, C_rate_C[s in S, p in 1:N_passes], F_C[s,1] * Cp_C[s, p] / NC)
    # minimum and maximum capacitance rates
    delta1 = delta2 = 1e-9
    @expression(model, C_rate_min[s in S, p in 1:N_passes], 0.5 * (C_rate_C[s, p] + C_rate_H[s, p] - sqrt((C_rate_C[s, p] - C_rate_H[s, p])^2 + delta1))) # J/(s*K) 
    @expression(model, C_rate_max[s in S, p in 1:N_passes], 0.5 * (C_rate_C[s, p] + C_rate_H[s, p] + sqrt((C_rate_C[s, p] - C_rate_H[s, p])^2 + delta2)))
    @expression(model, C_rate_ratio[s in S, p in 1:N_passes], C_rate_min[s, p] / C_rate_max[s, p]) # capacitance ratio 
    # Dynamic viscosity of hot and cold fluids,  Pa·s
    @expression(model, mu_H[s in S, p in 1:N_passes], calc_viscosity(HX_T_H_passes_in[s, p], [mole_frac_H_MEA[s,1], mole_frac_H_H2O[s,1]], alpha_H[s]))
    @expression(model, mu_C[s in S, p in 1:N_passes], calc_viscosity(HX_T_C_passes_in[s, p], [mole_frac_C_MEA[s,1], mole_frac_C_H2O[s,1]], alpha_C[s]))
    # Apparent velocity of hot and cold fluids, m/s
    cs_Area = plate_gap * HX_width # Cross-sectional area, m^2
    @expression(model, v_H[s in S], (F_H[s,1] / NC) * MW_H[s] / (density * 1000 * cs_Area))
    @expression(model, v_C[s in S], (F_C[s,1] / NC) * MW_C[s] / (density * 1000 * cs_Area))
    # Reynolds numbers of hot and cold fluids
    @expression(model, Re_H[s in S, p in 1:N_passes], density * v_H[s] * d_e / mu_H[s, p])
    @expression(model, Re_C[s in S, p in 1:N_passes], density * v_C[s] * d_e / mu_C[s, p])
    # Thermal conductivity
    @expression(model, Lambda_H[s in S, p in 1:N_passes], conductivity_Vredeveld([mole_frac_H_MEA[s,1], mole_frac_H_H2O[s,1]], [conductivity_MEA((HX_T_H_passes_in[s, p] + HX_T_C_passes_in[s, p]) / 2), conductivity_H2O((HX_T_H_passes_in[s, p] + HX_T_C_passes_in[s, p]) / 2)]))
    @expression(model, Lambda_C[s in S, p in 1:N_passes], conductivity_Vredeveld([mole_frac_C_MEA[s,1], mole_frac_C_H2O[s,1]], [conductivity_MEA((HX_T_H_passes_in[s, p] + HX_T_C_passes_in[s, p]) / 2), conductivity_H2O((HX_T_H_passes_in[s, p] + HX_T_C_passes_in[s, p]) / 2)]))
    # Prandtl Number of hot and cold fluid
    @expression(model, Pr_H[s in S, p in 1:N_passes], mu_H[s, p] * Cp_H[s, p] / (MW_H[s] / 1000) / Lambda_H[s, p])
    @expression(model, Pr_C[s in S, p in 1:N_passes], mu_C[s, p] * Cp_C[s, p] / (MW_C[s] / 1000) / Lambda_C[s, p])
    # Heat transfer coefficient
    HT_param = [0.4, 0.5746, 0.333, 1.441, 0.206]
    @expression(model, h_H[s in S, p in 1:N_passes], HT_param[1] * Lambda_H[s, p] / d_e * (Re_H[s, p])^HT_param[2] * (Pr_H[s, p])^HT_param[3])
    @expression(model, h_C[s in S, p in 1:N_passes], HT_param[1] * Lambda_C[s, p] / d_e * (Re_C[s, p])^HT_param[2] * (Pr_C[s, p])^HT_param[3]) # W/(m^2*K)
    # Overall heat transfer coefficient
    @expression(model, U[s in S, p in 1:N_passes], 1 / (1 / h_H[s, p] + 1 / h_C[s, p] + thickness / Lambda_P)) # W/(m^2*K)
    # Number of transfer units (NTU) for each pass
    @expression(model, NTU[s in S, p in 1:N_passes], U[s, p] * A_p / C_rate_min[s, p])
    # Effectiveness for each pass
    @expression(model, eff[s in S, p in 1:N_passes], (1 - exp(-NTU[s, p] * (1 - C_rate_ratio[s, p]))) / (1 - C_rate_ratio[s, p] * exp(-NTU[s, p] * (1 - C_rate_ratio[s, p]))))

    # Constraints
    cons_HX_T_C_in_fix = @constraint(model, cons_HX_T_C_in_fix[s in S], HX_T_C_passes_in[s, N_passes] == T_C_IN[s]) # fix inlet temperature of cold fluid
    cons_HX_T_H_in_boundary = @constraint(model, cons_HX_T_H_in_boundary[s in S, i=1:N_passes-1], HX_T_C_passes_in[s, i] == HX_T_C_passes_out[s, i+1]) # pass boundary condition
    const_HX_T_C_out = @constraint(model, const_HX_T_C_out[s in S, i=1:N_passes], HX_T_C_passes_out[s, i] == HX_T_C_passes_in[s, i] + eff[s, i] * C_rate_min[s, i] * (HX_T_H_passes_in[s, i] - HX_T_C_passes_in[s, i]) / C_rate_C[s, i]) # temperature of cold fluid at each pass
    append!(operation_cons_list, [cons_HX_T_C_in_fix, cons_HX_T_H_in_boundary, const_HX_T_C_out])

    cons_HX_T_H_in_fix = @constraint(model, cons_HX_T_H_in_fix[s in S], HX_T_H_passes_in[s, 1] == HX_T_H_IN[s]) # fix inlet temperature of hot fluid
    cons_HX_T_H_out_boundary = @constraint(model, cons_HX_T_H_out_boundary[s in S, i=1:N_passes-1], HX_T_H_passes_out[s, i] == HX_T_H_passes_in[s, i+1]) # pass boundary condition
    cons_HX_T_H_out = @constraint(model, cons_HX_T_H_out[s in S, i=1:N_passes], HX_T_H_passes_out[s, i] == HX_T_H_passes_in[s, i] - eff[s, i] * C_rate_min[s, i] * (HX_T_H_passes_in[s, i] - HX_T_C_passes_in[s, i]) / C_rate_H[s, i])
    append!(operation_cons_list, [cons_HX_T_H_in_fix, cons_HX_T_H_out_boundary, cons_HX_T_H_out])

    # initialization
    T_H0 = 395 # initial guess for the stream temperature exiting the reboiler
    T_C0 = check_val.(T_C_IN)
    for s in S
        # Hot-side inlet temperature, initial guesses
        set_start_value(HX_T_H_passes_in[s, 1], T_H0)
        for i in 2:N_passes
            set_start_value(HX_T_H_passes_in[s, i], T_H0 - ((T_H0 - T_C0[s]) / (N_passes - 1)) * (i - 1))
        end
        for i in 1:N_passes
            set_start_value(HX_T_H_passes_out[s, i], start_value(HX_T_H_passes_in[s, i]))
            set_start_value(HX_T_C_passes_in[s, i], start_value(HX_T_H_passes_in[s, i]))
            set_start_value(HX_T_C_passes_out[s, i], start_value(HX_T_H_passes_in[s, i]) - 10)
        end
    end

    #############################################################################
    # HX cold stream liquid-only isothermal flash constraints
    #############################################################################
    @expression(model, HX_C_x_Li[s in S, c=1:length(species_index)], F_C[s] * y_C[s][c]) # cold stream inlet molar flowrate of each species
    # Mole balance constraints
    cons_HX_C_MB = add_mole_balance!(model, HX_C_x, HX_C_L, nothing, nothing, HX_C_x_Li, nothing, HX_C_xi, S; N_seg=1)
    # Reaction equilibrium constraints (use first pass out temperature)
    (cons_HX_C_rxn1, cons_HX_C_rxn2, cons_HX_C_rxn3, cons_HX_C_rxn4, cons_HX_C_rxn5) = add_reaction_equilibria!(model, HX_C_log_x, reshape(HX_T_C_passes_out[:, 1], length(S), 1), S; N_seg=1)
    # Summation constraint (mole fraction)
    cons_HX_C_sum_liq = @constraint(model, cons_HX_C_sum_liq[s in S], sum(HX_C_x[s, 1, c] for c in 1:length(species_index)) == 1)
    append!(operation_cons_list, (cons_HX_C_MB, cons_HX_C_rxn1, cons_HX_C_rxn2, cons_HX_C_rxn3, cons_HX_C_rxn4, cons_HX_C_rxn5, cons_HX_C_sum_liq))

    # Initialization
    HX_C_L0 = [check_val(F_C[s]) for s in S]
    init_LVT(HX_C_L, nothing, nothing, HX_C_L0, nothing, nothing, S, 1)
    HX_C_log_x0 = reshape([start_value(abs_log_x[s, 1, c]) for s in 1:length(S), c in 1:length(species_index)], (length(S), 1, length(species_index)))
    init_mole_frac(HX_C_log_x, nothing, HX_C_log_x0, nothing, S, 1)

    HX_C_xi_0 = Dict(
        1 => 20,
        2 => -1,
        3 => -0.001,
        4 => 80,
        5 => -20
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(HX_C_xi[s, 1, r], HX_C_xi_0[r])
        end
    end 

    #############################################################################
    # Desorber constraint
    #############################################################################
    # Initial amount of the 11 species, moles
    @expression(model, des_Lf[s in S, c=1:length(species_index)], des_L_in[s,1] * des_x_in[s][c]) # liquid feed from HX
    @expression(model, des_Lc[s in S, c=1:length(species_index)], cond_L[s,1] * x_cond[s][c]) # condenser liquid feed
    @expression(model, des_x_Li[s in S, c=1:length(species_index)], des_Lf[s, c] + des_Lc[s, c]) # total liquid feed to desorber
    @expression(model, des_y_Vi[s in S, c=1:length(species_index)], des_V_in[s,1] * des_y_in[s,c]) # vapor feed from reboiler

    # Apparent species amounts in liquids from HX and from condenser
    des_app_Lf = compute_apparent_liquid(model, des_Lf, S) # apparent species in liquid feed from HX
    des_app_Lc = compute_apparent_liquid(model, des_Lc, S) # apparent species in liquid feed from condenser

    # Mole balance constraints
    if N_des == 1
        cons_des_MB = add_mole_balance!(model, des_x, des_L, des_y, des_V, des_x_Li, des_y_Vi, des_xi, S; N_seg=N_des)
        push!(operation_cons_list, cons_des_MB)
    elseif N_des == 2
        cons_des_MB_1, cons_des_MB_N = add_mole_balance!(model, des_x, des_L, des_y, des_V, des_x_Li, des_y_Vi, des_xi, S; N_seg=N_des)
        append!(operation_cons_list, [cons_des_MB_1, cons_des_MB_N])
    else
        cons_des_MB_1, cons_des_MB_N, cons_des_MB_n = add_mole_balance!(model, des_x, des_L, des_y, des_V, des_x_Li, des_y_Vi, des_xi, S; N_seg=N_des)
        append!(operation_cons_list, [cons_des_MB_1, cons_des_MB_n, cons_des_MB_N])
    end

    # Flooding and efficiency calculations desorber
    # design variables
    Height_desorber_init = 7 # initial value height of the column, m
    Diameter_desorber_init = 5 # initial value diameter of the column, m
    if !recycle || !as_optimization
        @variable(model, des_H == Height_desorber_init)
        @variable(model, des_D == Diameter_desorber_init)
    else
        if fix_design
            @variable(model, des_H == design[:des_H])
            @variable(model, des_D == design[:des_D])
        else
            @variable(model, 0.25*Height_desorber_init <= des_H <= 2*Height_desorber_init)
            @variable(model, 0.25*Diameter_desorber_init <= des_D <= 2*Diameter_desorber_init)
            set_start_value(des_D, Diameter_desorber_init)
            set_start_value(des_H, Height_desorber_init)
        end
    end
    append!(design_vars_list, [des_H, des_D])

    # Parameters/variables for flooding and approach-to-equilibrium model
    (des_app_L_MEA, des_app_L_H2O, des_app_L_CO2) = compute_apparent_liquid(model, des_x, S, N_des) 
    # CO2 loading in each segment, mol/mol
    @expression(model, des_alpha[s in S, n in 1:N_des], des_app_L_CO2[s, n] / des_app_L_MEA[s, n])
    # vapor superficial velocity in each segment, m/s
    @expression(model, des_u_v[s in S, n in 1:N_des], des_V[s, n] * (R*des_T[s,n]/(des_P_atm*101325)) / (des_D^2/4 * π))
    # Solvent and vapor average molar weights for every segment, kg/mol
    des_x_MW_avg = [sum(des_x[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_des]
    des_y_MW_avg = [sum(des_y[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_des]
    # vapor density in every segment, kg/m^3
    @expression(model, des_rho_g[s in S, n in 1:N_des], calculate_ideal_gas_density(des_P_atm, des_y_MW_avg[s,n], des_T[s, n]))
    # liquid kinematic viscosity in every segment, m2/(s) / 10^6 = cSt
    @expression(model, des_kinematic_viscosity[s in S, n in 1:N_des], calc_viscosity(des_T[s, n], [des_app_L_MEA[s,n], des_app_L_H2O[s,n]], des_alpha[s,n])/density*1e6 )
    # Inlet vapor composition and flowrate  for each segment (from segment below or feed)
    @expression(model, des_y_inlet[s in S, n in 1:N_des, c in 1:length(species_index)], n == 1 ? des_y_in[s, c] : des_y[s, n-1, c])
    @expression(model, des_V_inlet[s in S, n in 1:N_des], n == 1 ? des_V_in[s, 1] : des_V[s, n-1])
    # Inlet liquid composition and flowrate for each segment (from segment above or feed)
    @expression(model, des_x_inlet[s in S, n in 1:N_des, c in 1:length(species_index)], n == N_des ? des_x_in[s][c] : des_x[s, n+1, c])
    @expression(model, des_L_inlet[s in S, n in 1:N_des], n == N_des ? des_L_in[s, 1] : des_L[s, n+1])
    # Inlet apparent species mole fraction
    @expression(model, des_app_L_MEA_inlet[s in S, n in 1:N_des], des_x_inlet[s, n, species_index[:MEA]] + des_x_inlet[s, n, species_index[:MEA2p]] + des_x_inlet[s, n, species_index[:Carbamate]])
    @expression(model, des_app_L_H2O_inlet[s in S, n in 1:N_des], des_x_inlet[s, n, species_index[:H2O]] + des_x_inlet[s, n, species_index[:OHm]] + des_x_inlet[s, n, species_index[:Hp]] + des_x_inlet[s, n, species_index[:HCO3m]])
    @expression(model, des_app_L_CO2_inlet[s in S, n in 1:N_des], des_x_inlet[s, n, species_index[:CO2]] + des_x_inlet[s, n, species_index[:HCO3m]] + des_x_inlet[s, n, species_index[:CO3_2m]] + des_x_inlet[s, n, species_index[:Carbamate]])
    # Inlet solvent composition without CO2 for each segment (from segment above or feed)
    @expression(model, des_x_MEA_composition_in[s in S, n in 1:N_des], des_app_L_MEA_inlet[s,n] * MW[:MEA] / (des_app_L_MEA_inlet[s,n] * MW[:MEA] + des_app_L_H2O_inlet[s,n] * MW[:H2O]))
    # Inlet solvent CO2 loading, mol/mol
    @expression(model, des_alpha_inlet[s in S, n in 1:N_des], des_app_L_CO2_inlet[s, n] / des_app_L_MEA_inlet[s, n])
    # Inlet solvent and vapor average molar weights, kg/mol
    des_x_MW_avg_inlet = [sum(des_x_inlet[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_des]
    des_y_MW_avg_inlet = [sum(des_y_inlet[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_des]
    # Inlet vaport-to-liquid ratio, kg/kg
    @expression(model, des_V_over_L_inlet[s in S, n in 1:N_des], (des_V_inlet[s, n] * des_y_MW_avg_inlet[s, n]) / (des_L_inlet[s, n] * des_x_MW_avg_inlet[s, n]))
    # Inlet superficial velocity, m/s
    @expression(model, des_u_v_inlet[s in S, n in 1:N_des], des_V_inlet[s, n] * (R*des_T[s,n]/(des_P_atm*101325)) / (des_D^2/4 * π))
    # Inlet temperature in each segment, K
    @expression(model, des_T_inlet[s in S, n in 1:N_des], n == 1 ? des_T_L_in[s] : des_T[s, n-1])

    # diameter constraints to avoid flooding
    @expression(model, des_flow_parameter[s in S, n in 1:N_des], calculate_flow_parameter(des_L[s,n]*des_x_MW_avg[s,n], des_V[s,n]*des_y_MW_avg[s,n], des_rho_g[s,n], density)) # vapor flow parameter
    #constraints coupling GPDC flooding correlation with superficial velocity flooding correlation to solve for diameter
    cons_des_CP_at_flooding = @constraint(model, cons_des_CP_at_flooding[s in S, n=1:N_des], calculate_CP_GPDC_flood(des_flow_parameter[s, n]) == calculate_CP_superficial_velocity(des_U_s_flood[s, n], des_rho_g[s, n], density, des_kinematic_viscosity[s, n]))
    cons_des_superficial_velocity_limit = @constraint(model, cons_des_superficial_velocity_limit[s in S, n=1:N_des], des_u_v[s, n]*3.28084 <= 0.8 * des_U_s_flood[s, n]) # superficial velocity is maximum 80% of the flooding velocity, note that des_U_s_flood is in ft/s while u_v is in m/s, 3.28084 ft/m
    append!(operation_cons_list, [cons_des_CP_at_flooding, cons_des_superficial_velocity_limit])
    @expression(model, des_approach_to_flooding[s in S, n in 1:N_des], des_u_v[s, n]*3.28084/des_U_s_flood[s, n]) # dimensionless, ratio of superficial velocity to flooding velocity, 3.28084 ft/m

    # minimum liquid load, [4]
    @expression(model, des_L_load_superficial[s in S, n in 1:N_des], des_L[s,n] * des_x_MW_avg[s,n] / density / (des_D^2/4 * π) * 3600) # m3/m2h
    cons_des_min_liquid_load = @constraint(model, cons_des_min_liquid_load[s in S, n=1:N_des], des_L_load_superficial[s,n] >= minimum_liquid_load)
    push!(operation_cons_list, cons_des_min_liquid_load)

    # VLE constraints
    if N_des > 1
        (cons_des_raoult_H2O, cons_des_Henry_CO2_1, cons_des_Henry_CO2_n) = add_VLE!(model, des_x, des_y, des_T, S, N_des, density, des_P_atm, des_y_in, true; 
                    beta=beta[:des],
                    mu=mu[:des],
                    sigma=sigma[:des], 
                    V_over_L=des_V_over_L_inlet, 
                    T_in=des_T_inlet,
                    alpha=des_alpha_inlet, 
                    x_MEA_composition_in=des_x_MEA_composition_in,
                    u_v=des_u_v_inlet,
                    dz=des_H/N_des,
                    name="des")
        append!(operation_cons_list, [cons_des_raoult_H2O,cons_des_Henry_CO2_1, cons_des_Henry_CO2_n])
    else
        (cons_des_raoult_H2O, cons_des_Henry_CO2_1) = add_VLE!(model, des_x, des_y, des_T, S, N_des, density, des_P_atm, des_y_in, false; 
                    beta=beta[:des],
                    mu=mu[:des],
                    sigma=sigma[:des],
                    V_over_L=nothing, 
                    T_in=nothing,
                    alpha=nothing, 
                    x_MEA_composition_in=nothing, 
                    u_v=nothing, 
                    dz=nothing,
                    name="des")
        append!(operation_cons_list, [cons_des_raoult_H2O, cons_des_Henry_CO2_1])
    end

    # Reaction equilibrium constraints
    (cons_des_rxn1, cons_des_rxn2, cons_des_rxn3, cons_des_rxn4, cons_des_rxn5) = add_reaction_equilibria!(model, des_log_x, des_T, S; N_seg=N_des)
    append!(operation_cons_list, [cons_des_rxn1, cons_des_rxn2, cons_des_rxn3, cons_des_rxn4, cons_des_rxn5])
    # Summation constraint (mole fraction)
    cons_des_sum_liq = @constraint(model, cons_des_sum_liq[s in S, n=1:N_des], sum(des_x[s, n, i] for i in 1:length(species_index)) == 1)
    cons_des_sum_vap = @constraint(model, cons_des_sum_vap[s in S, n=1:N_des], des_y[s, n, species_index[:CO2]] + des_y[s, n, species_index[:H2O]] == 1 - des_y[s, n, species_index[:O2]] - des_y[s, n, species_index[:N2]])
    append!(operation_cons_list, [cons_des_sum_liq, cons_des_sum_vap])
    # Liquid and vapor enthalpy
    (des_H_liq_feed, des_H_liq_out) = compute_liquid_enthalpy(model, des_app_Lf, des_T_L_in, des_x, des_T, des_L, S, N_des)
    (des_H_liq_cond, _) = compute_liquid_enthalpy(model, des_app_Lc, cond_T[:,1], des_x, des_T, des_L, S, N_des)
    (des_H_vap_in, des_H_vap_out) = compute_vapor_enthalpy(model, des_y_Vi, des_T_V_in[:,1], des_y, des_T, des_V, S, N_des)
    # Energy balance constraint
    if N_des == 1
        cons_des_EB = add_energy_balance!(model, des_H_liq_feed + des_H_liq_cond, des_H_liq_out, des_H_vap_in, des_H_vap_out, nothing, S; N_seg=N_des)
        push!(operation_cons_list, cons_des_EB)
    elseif N_des == 2
        cons_des_EB_1, cons_des_EB_N = add_energy_balance!(model, des_H_liq_feed + des_H_liq_cond, des_H_liq_out, des_H_vap_in, des_H_vap_out, nothing, S; N_seg=N_des)
        append!(operation_cons_list, [cons_des_EB_1, cons_des_EB_N])
    else
        cons_des_EB_1, cons_des_EB_N, cons_des_EB_n = add_energy_balance!(model, des_H_liq_feed + des_H_liq_cond, des_H_liq_out, des_H_vap_in, des_H_vap_out, nothing, S; N_seg=N_des)
        append!(operation_cons_list, [cons_des_EB_1, cons_des_EB_n, cons_des_EB_N])
    end

    # Initialization
    des_L0 = [isnothing(check_val(des_L_in[s,1])) ? inout_data_s[s].liquid_inlet.flow : check_val(des_L_in[s,1]) for s in S]
    des_TL0 = [isnothing(check_val(des_T_L_in[s])) ? inout_data_s[s].liquid_inlet.T : check_val(des_T_L_in[s]) for s in S]
    init_LVT(des_L, des_V, des_T, des_L0, 0.2 .* des_L0, des_TL0, S, N_des)

    des_log_x0 = [something(check_val(des_log_x_in[s][c]), -5) for s in 1:length(S), n in 1:N_des, c in 1:length(species_index)]
    des_log_y0 = fill(-50.0, (length(S), N_des, length(species_index)))
    des_log_y0[:, :, species_index[:CO2]] .= log(0.2)
    des_log_y0[:, :, species_index[:H2O]] .= log(0.8)
    init_mole_frac(des_log_x, des_log_y, des_log_x0, des_log_y0, S, N_des)
    
    for s in S
        for n in 1:N_des
            set_start_value(des_U_s_flood[s,n], 5.0)
        end
    end

    # initial guess of des_xi
    des_xi_0 = Dict(
        1 => -3,
        2 => -0.005,
        3 => -0.001,
        4 => -2,
        5 => 3
    )
    for s in S
        for n in 1:N_des
            for r in 1:5  #5 reactions
                set_start_value(des_xi[s, n, r], des_xi_0[r])
            end
        end
    end

    #############################################################################
    # Reboiler constraints
    #############################################################################
    # reboiler solvent inlet molar flows per species
    @expression(model, reb_x_Li[s in S, c=1:length(species_index)], reb_L_in[s] * reb_x_in[s][c])
    # Apparent species amounts in the liquid feed
    reb_app_Li = compute_apparent_liquid(model, reb_x_Li, S)
    # Mole balance constraint
    cons_reb_MB = add_mole_balance!(model, reb_x, reb_L, reb_y, reb_V, reb_x_Li, nothing, reb_xi, S; N_seg=1)
    push!(operation_cons_list, cons_reb_MB)
    # VLE constraints
    (cons_reb_raoult_H2O, cons_reb_Henry_CO2) = add_VLE!(model, reb_x, reb_y, reb_T, S, 1, density, P_reb, nothing, false)
    append!(operation_cons_list, [cons_reb_raoult_H2O,cons_reb_Henry_CO2])
    # Reaction equilibrium constraints
    (cons_reb_rxn1, cons_reb_rxn2, cons_reb_rxn3, cons_reb_rxn4, cons_reb_rxn5) = add_reaction_equilibria!(model, reb_log_x, reshape(reb_T, length(S), 1), S; N_seg=1)
    append!(operation_cons_list, [cons_reb_rxn1, cons_reb_rxn2, cons_reb_rxn3, cons_reb_rxn4, cons_reb_rxn5])
    # Summation constraint (mole fraction)
    cons_reb_sum_liq = @constraint(model, cons_reb_sum_liq[s in S], sum(reb_x[s, 1, c] for c in 1:length(species_index)) == 1)
    cons_reb_sum_vap = @constraint(model, cons_reb_sum_vap[s in S], sum(reb_y[s, 1, c] for c in 1:length(species_index)) == 1)
    append!(operation_cons_list, [cons_reb_sum_liq, cons_reb_sum_vap])
    # Liquid and vapor enthalpy
    (reb_H_liq_in, reb_H_liq_out) = compute_liquid_enthalpy(model, reb_app_Li, reb_T_in, reb_x, reb_T, reb_L, S, 1)
    (_, reb_H_vap_out) = compute_vapor_enthalpy(model, nothing, nothing, reb_y, reb_T, reb_V, S, 1)
    # Energy balance constraint
    cons_reb_EB = add_energy_balance!(model, reb_H_liq_in, reb_H_liq_out, nothing, reb_H_vap_out, 1000*Q_reb, S; N_seg=1)
    push!(operation_cons_list, cons_reb_EB)

    # Initialization
    reb_L0 = [check_val(reb_L_in[s]) for s in S]
    reb_T0 = [393 for s in S]
    init_LVT(reb_L, reb_V, reb_T, 0.8 .* reb_L0, 0.2 .* reb_L0, reb_T0, S, 1) #Initialize  reboiler segment L, V, T

    reb_log_x0 = reshape([something(check_val(des_log_x[s, 1, c]), -5) for s in 1:length(S), c in 1:length(species_index)], (length(S), 1, length(species_index)))
    reb_log_y0 = fill(-50.0, (length(S), 1, length(species_index)))
    reb_log_y0[:, 1, species_index[:CO2]] .= log(0.02)
    reb_log_y0[:, 1, species_index[:H2O]] .= log(0.98)
    init_mole_frac(reb_log_x, reb_log_y, reb_log_x0, reb_log_y0, S, 1)

    reb_xi_0 = Dict(
        1 => 60,
        2 => -0.007,
        3 => -0.001,
        4 => 50,
        5 => -60
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(reb_xi[s, 1, r], reb_xi_0[r])
        end
    end

    #############################################################################
    # Condenser constraints
    #############################################################################
    # condenser vapor inlet flows per species
    @expression(model, cond_y_Vi[s in S, c=1:length(species_index)], cond_V_in[s] * cond_y_in[s][c])
    # Mole balance constraint
    cons_cond_MB = add_mole_balance!(model, cond_x, cond_L, cond_y, cond_V, nothing, cond_y_Vi, cond_xi, S; N_seg=1)
    push!(operation_cons_list, cons_cond_MB)
    # VLE constraints
    (cons_cond_raoult_H2O, cons_cond_Henry_CO2) = add_VLE!(model, cond_x, cond_y, cond_T, S, 1, density, P_cond, nothing, false)
    append!(operation_cons_list, [cons_cond_raoult_H2O,cons_cond_Henry_CO2])
    # Reaction equilibrium constraints
    (cons_cond_rxn1, cons_cond_rxn2, cons_cond_rxn3, cons_cond_rxn4, cons_cond_rxn5) = add_reaction_equilibria!(model, cond_log_x, reshape(cond_T, length(S), 1), S; N_seg=1)
    append!(operation_cons_list, [cons_cond_rxn1, cons_cond_rxn2, cons_cond_rxn3, cons_cond_rxn4, cons_cond_rxn5])
    # Summation constraint (mole fraction)
    cons_cond_sum_liq = @constraint(model, cons_cond_sum_liq[s in S], sum(cond_x[s, 1, c] for c in 1:length(species_index)) == 1)
    cons_cond_sum_vap = @constraint(model, cons_cond_sum_vap[s in S], sum(cond_y[s, 1, c] for c in 1:length(species_index)) == 1) # cond_y_CO2 + cond_y_H2O = 1 (No cond_O2 and cond_N2 in the condenser section)
    append!(operation_cons_list, [cons_cond_sum_liq, cons_cond_sum_vap])
    # Liquid and vapor enthalpy
    (_, cond_H_liq_out) = compute_liquid_enthalpy(model, nothing, nothing, cond_x, cond_T, cond_L, S, 1)
    (cond_H_vap_in, cond_H_vap_out) = compute_vapor_enthalpy(model, cond_y_Vi, cond_T_in, cond_y, cond_T, cond_V, S, 1)
    # Energy balance constraint
    cons_cond_EB = add_energy_balance!(model, nothing, cond_H_liq_out, cond_H_vap_in, cond_H_vap_out, 1000*Q_cond, S; N_seg=1)
    push!(operation_cons_list, cons_cond_EB)

    # Initialization
    cond_L0 = [check_val(cond_V_in[s]) for s in S]
    cond_T0 = [check_val(cond_T[s,1]) for s in S]
    init_LVT(cond_L, cond_V, cond_T, cond_L0*0.9, 0.1 .* cond_L0, cond_T0, S, 1)

    cond_log_x0 = fill(-50.0, (length(S), 1, length(species_index)))
    cond_log_x0[:, 1, species_index[:CO2]] .= log(0.02)
    cond_log_x0[:, 1, species_index[:H2O]] .= log(0.98)
    
    cond_log_y0 = fill(-50.0, (length(S), 1, length(species_index)))
    cond_log_y0[:, 1, species_index[:CO2]] .= log(0.95)
    cond_log_y0[:, 1, species_index[:H2O]] .= log(0.05)
    init_mole_frac(cond_log_x, cond_log_y, cond_log_x0, cond_log_y0, S, 1)

    cond_xi_0 = Dict( # initial guess
        1 => -6,
        2 => 0.0,
        3 => 0.0,
        4 => 0.0,
        5 => 0.004
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(cond_xi[s, 1, r], cond_xi_0[r])
        end
        set_start_value(Q_cond[s], inout_data_s[s].initial_values.Q_cond)
    end 

    #############################################################################
    # HX hot stream liquid-only isothermal flash constraints
    #############################################################################
    # molar flows per species in inlet of heat exchanger
    @expression(model, HX_H_x_Li[s in S, c=1:length(species_index)], F_H[s] * y_H[s][c]) 
    # Mole balance constraints
    cons_HX_H_MB = add_mole_balance!(model, HX_H_x, HX_H_L, nothing, nothing, HX_H_x_Li, nothing, HX_H_xi, S; N_seg=1)
    push!(operation_cons_list, cons_HX_H_MB)
    # Reaction equilibrium constraints (use last pass out temperature)
    (cons_HX_H_rxn1, cons_HX_H_rxn2, cons_HX_H_rxn3, cons_HX_H_rxn4, cons_HX_H_rxn5) = add_reaction_equilibria!(model, HX_H_log_x, reshape(HX_T_H_passes_out[:, N_passes], length(S), 1), S; N_seg=1)
    append!(operation_cons_list, [cons_HX_H_rxn1, cons_HX_H_rxn2, cons_HX_H_rxn3, cons_HX_H_rxn4, cons_HX_H_rxn5])
    # Summation constraint (mole fraction)
    cons_HX_H_sum_liq = @constraint(model, cons_HX_H_sum_liq[s in S], sum(HX_H_x[s, 1, c] for c in 1:length(species_index)) == 1)
    push!(operation_cons_list, cons_HX_H_sum_liq)

    # Initialization
    HX_H_L0 = [check_val(F_H[s,1]) for s in S]
    init_LVT(HX_H_L, nothing, nothing, HX_H_L0, nothing, nothing, S, 1)

    HX_H_log_x0 = check_val.(reb_log_x)
    init_mole_frac(HX_H_log_x, nothing, HX_H_log_x0, nothing, S, 1)

    HX_H_xi_0 = Dict(
        1 => -1,
        2 => 0.5,
        3 => 0.05,
        4 => -10,
        5 => 1
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(HX_H_xi[s, 1, r], HX_H_xi_0[r])
        end
    end 
    #endregion

    #############################################################################
    # Cooler constraints
    #############################################################################
    # inlet mole flows per species
    @expression(model, cool_x_Li[s in S, c=1:length(species_index)], cool_L_in[s,1] * cool_x_in[s][c])
    # Apparent species amounts in the liquid feed
    cool_app_Li = compute_apparent_liquid(model, cool_x_Li, S)
    # Mole balance constraints
    cons_cool_MB = add_mole_balance!(model, cool_x, cool_L, nothing, nothing, cool_x_Li, nothing, cool_xi, S; N_seg=1)
    push!(operation_cons_list, cons_cool_MB)
    # Reaction equilibrium constraints
    (cons_cool_rxn1, cons_cool_rxn2, cons_cool_rxn3, cons_cool_rxn4, cons_cool_rxn5) = add_reaction_equilibria!(model, cool_log_x, reshape(cool_T, length(S), 1), S; N_seg=1)
    append!(operation_cons_list, [cons_cool_rxn1, cons_cool_rxn2, cons_cool_rxn3, cons_cool_rxn4, cons_cool_rxn5])
    # Summation constraint (mole fraction)
    cons_cool_sum_liq = @constraint(model, cons_cool_sum_liq[s in S], sum(cool_x[s, 1, c] for c in 1:length(species_index)) == 1)
    push!(operation_cons_list, cons_cool_sum_liq)
    # Liquid and vapor enthalpy
    (cool_H_liq_in, cool_H_liq_out) = compute_liquid_enthalpy(model, cool_app_Li, cool_T_in, cool_x, cool_T, cool_L, S, 1)
    # Energy balance constraint
    cons_cool_EB = add_energy_balance!(model, cool_H_liq_in, cool_H_liq_out, nothing, nothing, 1000*Q_cool, S; N_seg=1)
    push!(operation_cons_list, cons_cool_EB)

    # Initialization
    cool_L0 = [check_val(cool_L_in[s,1]) for s in S]
    cool_T0 = [check_val(cool_T_in[s,1]) for s in S]
    init_LVT(cool_L, nothing, cool_T, cool_L0, nothing, cool_T0, S, 1)
    cool_log_x0 = check_val.(HX_H_log_x)
    init_mole_frac(cool_log_x, nothing, cool_log_x0, nothing, S, 1)
    cool_xi_0 = Dict(
        1 => 0.0,
        2 => 0.0,
        3 => 0.0,
        4 => 0.0,
        5 => 0.0
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(cool_xi[s, 1, r], cool_xi_0[r])
        end
    end 

    #############################################################################
    # Mixer constraints 
    #############################################################################
    # equilibrium mixer, isothermal assuming make-ups have the same temperature as the liquid output of cooler
    @expression(model, mix_x_Li[s in S, c=1:length(species_index)], cool_L[s,1] * cool_x[s, 1, c])
    for s in S
        mix_x_Li[s, species_index[:H2O]] += MU_H2O[s]
        mix_x_Li[s, species_index[:MEA]] += MU_MEA[s]
    end
    cons_mix_MB = add_mole_balance!(model, mix_x, mix_L, nothing, nothing, mix_x_Li, nothing, mix_xi, S; N_seg=1)
    push!(operation_cons_list, cons_mix_MB)
    (cons_mix_rxn1, cons_mix_rxn2, cons_mix_rxn3, cons_mix_rxn4, cons_mix_rxn5) = add_reaction_equilibria!(model, mix_log_x, reshape(mix_T, length(S), 1), S; N_seg=1)
    append!(operation_cons_list, [cons_mix_rxn1, cons_mix_rxn2, cons_mix_rxn3, cons_mix_rxn4, cons_mix_rxn5])
    cons_mix_sum_liq = @constraint(model, cons_mix_sum_liq[s in S], sum(mix_x[s, 1, c] for c in 1:length(species_index)) == 1)
    push!(operation_cons_list, cons_mix_sum_liq)
    cons_mix_T_coupling = @constraint(model, cons_mix_T_coupling[s in S], mix_T[s,1] == cool_T[s,1])
    push!(operation_cons_list, cons_mix_T_coupling)
    if recycle # close the loop
        cons_abs_mix_L = @constraint(model, cons_abs_mix_L[s in S], abs_L_in[s] == mix_L[s,1])
        push!(operation_cons_list, cons_abs_mix_L)
        cons_abs_mix_x = @constraint(model, cons_abs_mix_x[s in S, c=1:length(species_index)], abs_x_in[s, c] == mix_x[s,1,c])
        push!(operation_cons_list, cons_abs_mix_x)
        cons_abs_mix_T = @constraint(model, cons_abs_mix_T[s in S], abs_T_L_in[s] == mix_T[s,1])
        push!(operation_cons_list, cons_abs_mix_T)
    end

    # initialization
    mix_L0 = [check_val(cool_L[s,1]) + check_val(MU_H2O[s]) + check_val(MU_MEA[s]) for s in S]
    mix_T0 = [check_val(cool_T[s,1]) for s in S]
    init_LVT(mix_L, nothing, mix_T, mix_L0, nothing, mix_T0, S, 1)
    mix_log_x0 = cool_log_x0  # use same initial values as for the cooler
    for s in S
        cool_L_start = cool_L0[s]
        amount = [cool_L_start * exp(cool_log_x0[s,1,c]) for c in 1:length(species_index)]
        amount[species_index[:H2O]] += something(check_val(MU_H2O[s]), 0.0)
        amount[species_index[:MEA]] += something(check_val(MU_MEA[s]), 0.0)
        total_amt = sum(amount)
        for c in 1:length(species_index)
            mix_log_x0[s,1,c] = log(max(amount[c] / total_amt, 1e-50))
        end
    end
    init_mole_frac(mix_log_x, nothing, mix_log_x0, nothing, S, 1)

    mix_xi_0 = Dict(
        1 => -0.00006,
        2 => 0.00005,
        3 => 0.000009,
        4 => 0.000007,
        5 => -0.0
    )
    for s in S
        for r in 1:5  #5 reactions
            set_start_value(mix_xi[s, 1, r], mix_xi_0[r])
        end
    end

    #############################################################################
    # Objective function
    #############################################################################
    if !as_optimization
        # The following variables are only needed if we perform a design optimization.
        # We define them here so recycle and non-recycle cases have the same variables.
        # This simplifies initialization of the recycle case with the non-recycle case.

        # initial value for maximum pump volume flow
        # average molar weight of solvent inlet to absorber
        @expression(model, MW_inlet_avg[s in S], sum(xL[s][species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3) # kg/mol
        max_flow = maximum(L_in_s[s] * MW_inlet_avg[s] for s in S) / density * gal_per_m3 * 60 # in gpm
        @variable(model, pump_maxvflow_rich == max_flow) # in gpm
        @variable(model, pump_maxvflow_lean == max_flow) # in gpm
        @variable(model, A_reboiler_max == 17000)
        @variable(model, A_condenser_max == 21000)
        @variable(model, T_coolingwater_out_intercooler[s in S] == (105-32)*5/9 + 273.15) 
        @variable(model, A_intercooler_max == 47000)
        append!(design_vars_list, [pump_maxvflow_rich, pump_maxvflow_lean, A_reboiler_max, A_condenser_max, A_intercooler_max]) #,A_cooler_max
        append!(operation_vars_list, [T_coolingwater_out_intercooler]) #T_coolingwater_out_cooler
        fix(comp_P_max, 3000, force=true)
    else
        #set capture target
        if use_average_capture_rate == true
            # set scenario weighted average CO2 capture target
            cons_capture_target = @constraint(model, cons_capture_target, sum(abs_V_in[s] * abs_y_in[s, species_index[:CO2]]* inout_data_s[s].weight*sensitivity_flh*8760 - abs_y[s, N_abs, species_index[:CO2]] * abs_V[s, N_abs]* inout_data_s[s].weight*sensitivity_flh*8760 for s in S) >= capture_rate * sum(abs_V_in[s] * abs_y_in[s, species_index[:CO2]]* inout_data_s[s].weight*sensitivity_flh*8760 for s in S))
            push!(design_cons_list, cons_capture_target)
        else
            # set capture target for every scenario
            cons_capture_target = @constraint(model, cons_capture_target[s in S], abs_V_in[s] * abs_y_in[s, species_index[:CO2]] - abs_y[s, N_abs, species_index[:CO2]] * abs_V[s, N_abs] >= capture_rate * abs_V_in[s] * abs_y_in[s, species_index[:CO2]]) # e.g., set 95% CO2 capture target in every scenario
            push!(operation_cons_list, cons_capture_target)
        end

        # Bare erected costs
        # columns
        # applies additional height factors, factors are based on [5]
        abs_add_h_factor = 1.269979508
        des_add_h_factor = 0.846567797
        @expression(model, BEC_absorber, bec_column(abs_H,abs_D, abs_add_h_factor*abs_H, maximum_allowable_stress=20000))
        @expression(model, BEC_absorber_packing, bec_packing(abs_H, abs_D))
        @expression(model, BEC_desorber, bec_column(des_H,des_D, des_add_h_factor*des_H, maximum_allowable_stress=15000))
        @expression(model, BEC_desorber_packing, bec_packing(des_H, des_D))
   
        # pump, rich
        abs_outlet_MW_avg = [sum(abs_x[s, 1, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S]
        @expression(model, volume_flow_pump_rich[s in S], abs_L[s,1] * abs_outlet_MW_avg[s] / density * gal_per_m3 * 60) # in gpm
        # maximum volume flow among scenarios
        if fix_design
            @variable(model, pump_maxvflow_rich == design[:pump_maxvflow_rich]) # in gpm
        else
            @variable(model, pump_maxvflow_rich >= 250)
            set_start_value(pump_maxvflow_rich, 5500)
        end
        push!(design_vars_list, pump_maxvflow_rich)
        cons_pump_maxvflow_rich = @constraint(model, cons_pump_maxvflow_rich[s in S], pump_maxvflow_rich >= volume_flow_pump_rich[s])
        push!(operation_cons_list, cons_pump_maxvflow_rich)
        pump_head_rich = (des_P_atm-abs_P_atm) * 101325 / g / density + (des_H+des_add_h_factor*des_H)  # in m, pressure change + hydrostatic pressure
        @expression(model, head_pump_rich, pump_head_rich * ft_per_m) # in ft
        # note: the maximum pump size is 5000 gpm. The solvent flow for large-scale CCS plants is much higher.
        # Therefore, we assume a pump network with parallel pumps. We assume pumps of 5000 gpm each and assume a continuous number of pumps to meet the required flow (integer relaxation)
        v_max_pump = 5000 # gpm
        bec_per_pump_rich = bec_pump(v_max_pump,head_pump_rich,density* lb_per_kg / gal_per_m3)
        @expression(model, BEC_pump_rich, bec_per_pump_rich * pump_maxvflow_rich / v_max_pump)
        
        # pump, lean
        reb_outlet_MW_avg = [sum(reb_x[s, 1, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S]
        @expression(model, volume_flow_pump_lean[s in S], reb_L[s, 1] * reb_outlet_MW_avg[s] / density * gal_per_m3 * 60) # in gpm
        # find the maximum volume flow among scenarios
        if fix_design
            @variable(model, pump_maxvflow_lean == design[:pump_maxvflow_lean]) # in gpm
        else
            @variable(model, pump_maxvflow_lean >= 250)
        end
        push!(design_vars_list, pump_maxvflow_lean)
        set_start_value(pump_maxvflow_lean, 5500)
        cons_pump_maxvflow_lean = @constraint(model, cons_pump_maxvflow_lean[s in S], pump_maxvflow_lean >= volume_flow_pump_lean[s])
        push!(operation_cons_list, cons_pump_maxvflow_lean)
        pump_head_lean = abs_H+abs_add_h_factor*abs_H  # in m
        @expression(model, head_pump_lean, pump_head_lean * ft_per_m) # in ft, pressure change + hydrostatic pressure
        # note: the maximum pump size is 5000 gpm. The solvent flow for large-scale CCS plants is much higher.
        # Therefore, we assume a pump network with parallel pumps. We assume pumps of 5000 gpm each and assume a continuous number of pumps to meet the required flow (integer relaxation)
        bec_per_pump_lean = bec_pump(v_max_pump,head_pump_lean,density* lb_per_kg / gal_per_m3)
        @expression(model, BEC_pump_lean, bec_per_pump_lean * pump_maxvflow_lean / v_max_pump)

        # Cross heat exchanger
        # 3 dividers, based on [1]
        @expression(model, BEC_hex, bec_heat_exchanger(A_p * ft_per_m^2 * (2 * NC * N_passes - (1 + 3)))) 

        # compressor
        # as an order of magnitude estimate, we assume constant electricity demand per t of CO2 compressed
        electricity_factor_compression = 87*3600 # kJ/t, [6]
        intercooling_compression = 0.6*1000000 # kJ/t, [6]
        cond_y_MW_avg = [sum(cond_y[s, 1, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S] # kg/mol
        @expression(model, comp_P[s in S], electricity_factor_compression * cond_y_MW_avg[s]*cond_V[s,1]/1000) # kW
        cons_comp_P_max = @constraint(model, cons_comp_P_max[s in S], comp_P[s] <= comp_P_max)
        push!(operation_cons_list, cons_comp_P_max)
        @expression(model, BEC_compressor, bec_compressor(comp_P_max*hp_per_kW))
        @expression(model, tpc_compressor, BEC_compressor*(1+tpc_factor_compressor) / retrofit_factor)
        
        ### calculation steps for COA 
        # heat exchangers
        # effective heat transfer coefficients from [7] 
        # reboiler: Ueff = 852 W/m2K
        # regeneration using 150 psig steam (saturation temperature 185.5C), which is condensated during heating
        T_steam = 185.5 + 273.15 # in K
        @expression(model, log_T_diff_reb[s in S], ((T_steam - reb_T[s,1])*(T_steam - des_T[s,1]) * ((T_steam - reb_T[s,1]) + (T_steam - des_T[s,1])) / 2)^(1/3)) # Chen's approximation for the log mean temperature difference
        @expression(model, A_reboiler[s in S], Q_reb[s]*1000 / 852 / log_T_diff_reb[s] * ft_per_m^2) # heat transfer area in ft2
        if fix_design
            @variable(model, A_reboiler_max == design[:A_reboiler_max])
        else
            @variable(model, A_reboiler_max >= 150)
            set_start_value(A_reboiler_max, 15000)
        end
        push!(design_vars_list, A_reboiler_max)
        cons_A_reboiler_max = @constraint(model, cons_A_reboiler_max[s in S], A_reboiler[s] <= A_reboiler_max)
        push!(operation_cons_list, cons_A_reboiler_max)
        # note: the maximum tube and shell heat exchanger size is 12000ft2. The required reboiler is much larger.
        # Therefore, we assume multiple shells. We assume shells of 12000ft2 each and assume a continuous number of shells to meet the required heat transfer duty (integer relaxation)
        hex_tubeandshell_max = 12000 # ft2
        BEC_reboiler_per_hex = bec_hex_tubeandshell(hex_tubeandshell_max; reboiler=true)
        @expression(model, BEC_reboiler, A_reboiler_max / hex_tubeandshell_max * BEC_reboiler_per_hex)

        # cooling water
        # we assume cooling water to be heated from 90-120F for the condenser, where the outlet gases are much hotter than the cooling water avoiding temperature crossover
        # However, for cooler and intercooler, the inlet temperature of the hot fluid might be below 120F.
        # Therefore, we introduce a variable for the cooling water outlet temperature for the cooler and intercooler
        T_coolingwater_in = (90-32)*5/9 + 273.15 # 90F to K conversion
        T_coolingwater_out = (120-32)*5/9 + 273.15 # 120F to K conversion
        delta_T_min = 0  # minimum temperature approach for the cooler and intercooler

        # condenser: Ueff = 454 W/m2K
        @expression(model, log_T_diff_cond[s in S], ((des_T[s,N_des] - T_coolingwater_out)*(cond_T[s,1] - T_coolingwater_in) * ((des_T[s,N_des] - T_coolingwater_out) + (cond_T[s,1] - T_coolingwater_in))/2)^(1/3)) # Chen's approximation of the log mean temperature difference 
        @expression(model, A_condenser[s in S], -Q_cond[s]*1000 / 454 / log_T_diff_cond[s] * ft_per_m^2) # heat transfer area in ft2
        if fix_design
            @variable(model, A_condenser_max == design[:A_condenser_max])
        else
            @variable(model, A_condenser_max >= 150)
        end
        push!(design_vars_list, A_condenser_max)
        set_start_value(A_condenser_max, 20000)
        cons_A_condenser_max = @constraint(model, cons_A_condenser_max[s in S], A_condenser[s] <= A_condenser_max)
        push!(operation_cons_list, cons_A_condenser_max)
        # note: the maximum tube and shell heat exchanger size is 12000ft2. The required condenser heat exchanger is much larger.
        # Therefore, we assume multiple shells. We assume shells of 12000ft2 each and assume a continuous number of shells to meet the required heat transfer duty (integer relaxation)
        BEC_condenser_per_hex = bec_hex_tubeandshell(hex_tubeandshell_max; reboiler=false)
        @expression(model, BEC_condenser, A_condenser_max / hex_tubeandshell_max * BEC_condenser_per_hex)

        # lean amine cooler: Ueff = 795 W/m2K
        @variable(model, (120-32)*5/9 + 273.15 >= T_coolingwater_out_cooler[s in S] >= T_coolingwater_in) # introducing a variable to ensure, cooling water outlet temperature is below the lean solvent inlet temperature 
        push!(operation_vars_list, T_coolingwater_out_cooler)
        cons_T_cool_approach = @constraint(model, cons_T_cool_approach[s in S], T_coolingwater_out_cooler[s]+delta_T_min <= HX_T_H_passes_out[s, N_passes])
        push!(operation_cons_list, cons_T_cool_approach)
        # note: the initial value supplied by the open-loop model should not result in negative temperature difference factiors in the following expression, as the ^(1/3) operation causes NaN and solver failure.
        @expression(model, log_T_diff_cool[s in S], ((HX_T_H_passes_out[s, N_passes] - T_coolingwater_out_cooler[s])*(cool_T[s,1] - T_coolingwater_in) * ((HX_T_H_passes_out[s, N_passes] - T_coolingwater_out_cooler[s]) + (cool_T[s,1] - T_coolingwater_in))/2)^(1/3)) # log mean temperature difference for the cooler
        @expression(model, A_cooler[s in S], -Q_cool[s]*1000 / 795 / log_T_diff_cool[s] * ft_per_m^2) # heat transfer area in ft2
        if fix_design
            @variable(model, A_cooler_max == design[:A_cooler_max])
        else
            @variable(model, A_cooler_max >= 150)
        end
        push!(design_vars_list, A_cooler_max)
        set_start_value(A_cooler_max, 50000)
        cons_A_cooler_max = @constraint(model, cons_A_cooler_max[s in S], A_cooler[s] <= A_cooler_max)
        push!(operation_cons_list, cons_A_cooler_max)
        # note: the maximum tube and shell heat exchanger size is 12000ft2. The required cooler heat exchanger is much larger.
        # Therefore, we assume multiple shells. We assume shells of 12000ft2 each and assume a continuous number of shells to meet the required heat transfer duty (integer relaxation)
        BEC_cooler_per_hex = bec_hex_tubeandshell(hex_tubeandshell_max; reboiler=false)
        @expression(model, BEC_cooler, A_cooler_max / hex_tubeandshell_max * BEC_cooler_per_hex)

        # intercooler: assume same Ueff as lean amine cooler= 795 W/m2K
        @variable(model, (120-32)*5/9 + 273.15 >= T_coolingwater_out_intercooler[s in S] >= T_coolingwater_in) # introducing a variable to ensure, cooling water outlet temperature is below the solvent inlet temperature 
        push!(operation_vars_list, T_coolingwater_out_intercooler)
        cons_T_intercool_approach = @constraint(model, cons_T_intercool_approach[s in S], T_coolingwater_out_intercooler[s]+delta_T_min <= abs_T[s, N_intercool+1])  # note: Q_intercool is negative since it's cooling duty
        push!(operation_cons_list, cons_T_intercool_approach)
        @expression(model, log_T_diff_intercool[s in S], ((abs_T[s, N_intercool+1] - T_coolingwater_out_intercooler[s])*(abs_T[s, N_intercool] - T_coolingwater_in) * ((abs_T[s, N_intercool+1] - T_coolingwater_out_intercooler[s]) + (abs_T[s, N_intercool] - T_coolingwater_in))/2)^(1/3)) # log mean temperature difference for the intercooler
        @expression(model, A_intercooler[s in S], -abs_Q_intercool[s, N_intercool]*1000 / 795 / log_T_diff_intercool[s] * ft_per_m^2) # heat transfer area in ft2
        if fix_design
            @variable(model, A_intercooler_max == design[:A_intercooler_max])
        else
            @variable(model, A_intercooler_max >= 150)
        end
        push!(design_vars_list, A_intercooler_max)
        set_start_value(A_intercooler_max, 50000)
        cons_A_intercooler_max = @constraint(model, cons_A_intercooler_max[s in S], A_intercooler[s] <= A_intercooler_max)
        push!(operation_cons_list, cons_A_intercooler_max)
        # note: the maximum tube and shell heat exchanger size is 12000ft2. The required intercooler heat exchanger is much larger.
        # Therefore, we assume multiple shells. We assume shells of 12000ft2 each and assume a continuous number of shells to meet the required heat transfer duty (integer relaxation)
        BEC_intercooler_per_hex = bec_hex_tubeandshell(hex_tubeandshell_max; reboiler=false)
        @expression(model, BEC_intercooler, A_intercooler_max / hex_tubeandshell_max * BEC_intercooler_per_hex)

        @expression(model, bare_erect_cost_CCS, BEC_absorber + BEC_absorber_packing + BEC_desorber + BEC_desorber_packing + BEC_pump_lean + BEC_pump_rich + BEC_hex + BEC_reboiler + BEC_condenser + BEC_intercooler + BEC_cooler)
        @expression(model, tpc_CCS, calc_total_plant_cost(bare_erect_cost_CCS))
        tpc = tpc_CCS + tpc_compressor
        @expression(model, steam_consumption[s in S], Q_reb[s] * 8760 * inout_data_s[s].weight * sensitivity_flh) # [kWh]
        pump_eff = 0.65  # assume efficiency of 65%, [8]
        @expression(model, electricity_consumption[s in S], ((pump_head_rich*g*density*volume_flow_pump_rich[s]/gal_per_m3/60/1000/pump_eff + pump_head_lean*g*density*volume_flow_pump_lean[s]/gal_per_m3/60/1000/pump_eff) / 1000 + comp_P[s] / 1000) * 8760 * inout_data_s[s].weight*sensitivity_flh) # MWh
        cp_H2O = calc_Cp([0.0,1.0], T_coolingwater_in) / MW[:H2O]  # kJ/kg-K, specific heat capacity of water (assumed constant)
        #cooling water price is not scenario dependent, take sum directly
        @expression(model, coolingwater_consumption, sum((-Q_cond[s]/cp_H2O/(T_coolingwater_out-T_coolingwater_in) / density
                                                            + -Q_cool[s]/cp_H2O/(T_coolingwater_out-T_coolingwater_in) / density 
                                                            + -abs_Q_intercool[s, N_intercool]/cp_H2O/(T_coolingwater_out_intercooler[s]-T_coolingwater_in) / density
                                                            + intercooling_compression*cond_y_MW_avg[s]*cond_V[s,1]/1000 / cp_H2O/(T_coolingwater_out-T_coolingwater_in) / density) #cooling duty is approximated based on product volume flow
                                                            * 8760 * 3600 * inout_data_s[s].weight*sensitivity_flh for s in S)) # cooling water consumption in m3
        @expression(model, annual_CO2_captured, sum((abs_V_in[s] * abs_y_in[s,species_index[:CO2]] - abs_V[s,N_abs] * abs_y[s,N_abs,species_index[:CO2]]) * 8760 * 3600 * inout_data_s[s].weight*sensitivity_flh for s in S) * MW[:CO2] / 1000000)  # t/year, todo: for now assume equal weighting of scenarios
        
        # sensitivity_price_reduction is the reduction in electricity price for every 10% of partload decrease in (USD/MWh / 10%) 
        # based on the sensitivity_price_reduction, the electricity price is adjusted such that the weighted electricity price is equal to the average electtriticy price of 70 USD/MWh
        # reformulate sum((peak_electricity_price - (s-1) * sensitivity_price_reduction) * inout_data_s[s].weight for s in S) / sum(inout_data_s[s].weight for s in S) = price_electricity
        # Note: The formulation below assumes a number of consecutive scenarios in S, e.g., 6, 7, 8 , 9, with 100% load starting at the lowest scenario and part load reduction by 10 percentage points in each consecutive scenario
        peak_electricity_price = price_electricity +  sum((s-minimum(S))*sensitivity_price_reduction* inout_data_s[s].weight for s in S) / sum(inout_data_s[s].weight for s in S)
        adjusted_eletricity_price = Dict(s => peak_electricity_price - (s-minimum(S)) * sensitivity_price_reduction for s in S)
        price_factor = Dict(s => adjusted_eletricity_price[s] / price_electricity for s in S)  # relative price change, will be applied to adjust both electricity and steam price
        @expression(model, annual_opex, calc_opex(tpc, steam_consumption, electricity_consumption, coolingwater_consumption, annual_CO2_captured, price_factor, S, n_trains))
        @expression(model, toc, calc_toc(tpc, annual_CO2_captured, n_trains))
        @expression(model, tac, calc_tac(toc, annual_opex, CCF*sensitivity_CCF))
        @expression(model, coa, calc_coc(tac, annual_CO2_captured))

        #set objective
        @objective(model, Min, tac) # objective, total annualized cost as opposed of coa, which may result in higher than necesary CO2 capture
    end

    ############################### Collect variables and constraints for output #####################################
    vars_by_scenario = Dict{Int,Dict{Symbol,Any}}()
    cons_by_scenario = Dict{Int,Dict{Symbol,Any}}()
    for s in S
        scen_dict = get!(vars_by_scenario, s, Dict{Symbol,Any}())
        for var in operation_vars_list
            key = base_name(var) # get name of variable
            scen_dict[key] = slice_by_scenario(var, s)
        end
    end
    
    function base_name_or_fallback(con, i)  # for anonymous constraint blocks
        first_con = first(con)
        nm = name(first_con)

        if isempty(nm)
            return Symbol("anonymous_constraint_block_$i")
        end

        return Symbol(split(nm, '[', limit=2)[1])
    end 
    for s in S
        scen_dict = get!(cons_by_scenario, s, Dict{Symbol,Any}())
        for (i, con) in enumerate(operation_cons_list)
            key = base_name_or_fallback(con, i)

            if haskey(scen_dict, key)
                error("Duplicate constraint key in scenario $s: $key")
            end

            scen_dict[key] = slice_by_scenario(con, s)
            end
        end

    return model, design_vars_list, vars_by_scenario, cons_by_scenario, design_cons_list
end
