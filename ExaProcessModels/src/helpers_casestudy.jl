# This file contains case study specific helper functions

#####################################################################################
#fitting functions for absorber and desorber efficiency
beta_abs_coal = [
    -2.7136897233940678,
    -0.48136954218974376,
    0.2874035800259265,
    0.0237218135664058,
    0.31291979646251433,
    0.15199135488536838,
    -0.19797411961339487,
    0.044366651965475755,
    -0.3349398217917533,
    -0.8205665151314313
]

beta_abs_gas = [
    -2.9257929179087196,
    -0.2168171417876388,
    -0.16702236251528552,
    0.01789251724968559,
    1.3576497553182083,
    2.6122754786270552,
    -0.1919457622791252,
    0.017374365399178686,
    0.04132306520734156,
    -3.3025103380558085
] 

mu_abs_coal = Dict(
    :V_over_L => 0.4601150476449554,
    :T_in => 54.930391580946804,
    :P => 1.0856287728659144,
    :alpha => 0.33947584595367597,
    :x_MEA_composition_in => 0.30199020198854504,
    :u_v => 2.1194747000675855,
    :V_over_L_u_v => 0.9610820464750373,
    :T_in_alpha => 18.233430703074614,
    :alpha_x_MEA_composition_in => 0.10265521747758752
)

mu_abs_gas = Dict(
    :V_over_L => 0.4640806819087078,
    :T_in => 49.67091503858978,
    :P => 1.0860498908666092,
    :alpha => 0.4058046838641502,
    :x_MEA_composition_in => 0.2992248946714794,
    :u_v => 2.047706304956066,
    :V_over_L_u_v => 0.9436889140472463,
    :T_in_alpha => 19.979404594227745,
    :alpha_x_MEA_composition_in => 0.12134433871140307
)

sigma_abs_coal = Dict(
    :V_over_L => 0.2903381950816358,
    :T_in => 8.912766462211481,
    :P => 0.08058302231089982,
    :alpha => 0.12462659149976532,
    :x_MEA_composition_in => 0.004530065557885842,
    :u_v => 0.7789480272789379,
    :V_over_L_u_v => 0.7073315493847825,
    :T_in_alpha => 6.355018980556505,
    :alpha_x_MEA_composition_in => 0.03792762070076291
)

sigma_abs_gas = Dict(
    :V_over_L => 0.2939695830473947,
    :T_in => 7.40481805989205,
    :P => 0.08055902915546696,
    :alpha => 0.08101182264376829,
    :x_MEA_composition_in => 0.06414122307136633,
    :u_v => 0.7456945494340175,
    :V_over_L_u_v => 0.710097148881472,
    :T_in_alpha => 4.566274507710071,
    :alpha_x_MEA_composition_in => 0.036664614960614725
)

beta_des = [
    -2.6861275403434286,
    0.7948735598739742,
    1.0836539524550388,
    -0.5325756054772249,
    -41.142347486381276,
    -3.0783541072406604,
    -0.014779228922815642,
    -0.6058868276196839,
    16.696505924626237,
    26.279073238459834
]

mu_des = Dict(
    :V_over_L => 0.10044436597160843,
    :T_in => 109.65219238444979,
    :P => 1.7762257740441554,
    :alpha => 0.3116828639252733,
    :x_MEA_composition_in => 0.23038073125469793,
    :u_v => 0.6384479275152884,
    :V_over_L_u_v => 0.0960221135043422,
    :T_in_alpha => 33.118378788409075,
    :alpha_x_MEA_composition_in => 0.0738091310530041
)

sigma_des = Dict(
    :V_over_L => 0.06348362699870036,
    :T_in => 8.5475660005065,
    :P => 0.0830088901162176,
    :alpha => 0.1302850865245639,
    :x_MEA_composition_in => 0.016426767976211887,
    :u_v => 0.5264987766143678,
    :V_over_L_u_v => 0.10862286016076667,
    :T_in_alpha => 11.833013415786034,
    :alpha_x_MEA_composition_in => 0.034512002815540024
)

beta = Dict(
    :abs_coal => beta_abs_coal,
    :abs_gas => beta_abs_gas,
    :des => beta_des
)

mu = Dict(
    :abs_coal => mu_abs_coal,
    :abs_gas => mu_abs_gas,
    :des => mu_des
)

sigma = Dict(
    :abs_coal => sigma_abs_coal,
    :abs_gas => sigma_abs_gas,
    :des => sigma_des
)

#####################################################################################
# Flue gas data
function mole_to_mass_fractions(y::Dict{Symbol,<:Real}, MW::Dict{Symbol,<:Real})
    """
    Converts mole fractions y to mass fractions x using molecular weights MW.
    """
    MW_avg = sum(y[c] * MW[c] for c in keys(y))
    x = Dict(c => y[c] * MW[c] / MW_avg for c in keys(y))
    x_sum = sum(values(x))
    return Dict(c => x[c] / x_sum for c in keys(x))
end

# Coal power plant flue gas
# We assume constant H2O and O2 concentrations based on Schmitt T, Leptinsky S, Turner M, Zoelle A, White C, Hughes S, Homsy S, Woods M, Hoffman H, Shultz T, James III R. Cost and Performance Baseline for Fossil Energy Plants Volume 1: Bituminous Coal and Natural Gas to Electricity. (2022).
# component molar fractions per scenario
y_flue_gas_coal = Dict(
    1 => Dict(:CO2 => 0.137, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.137-0.036-0.151), #full load
    2 => Dict(:CO2 => 0.134, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.134-0.036-0.151), #90%
    3 => Dict(:CO2 => 0.130, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.130-0.036-0.151), #80%
    4 => Dict(:CO2 => 0.127, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.127-0.036-0.151), #70%
    5 => Dict(:CO2 => 0.121, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.121-0.036-0.151), #60%
    6 => Dict(:CO2 => 0.112, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.112-0.036-0.151), #50%
    7 => Dict(:CO2 => 0.103, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.103-0.036-0.151), #40%
    # nominal scenario
    8 => Dict(:CO2 => 0.137, :O2  => 0.036, :H2O => 0.151, :N2  => 1-0.137-0.036-0.151), #full load
)

# component mass fractions per scenario 
x_flue_gas_coal = Dict(
    1 => mole_to_mass_fractions(y_flue_gas_coal[1], MW), #full load
    2 => mole_to_mass_fractions(y_flue_gas_coal[2], MW),
    3 => mole_to_mass_fractions(y_flue_gas_coal[3], MW),
    4 => mole_to_mass_fractions(y_flue_gas_coal[4], MW),
    5 => mole_to_mass_fractions(y_flue_gas_coal[5], MW),
    6 => mole_to_mass_fractions(y_flue_gas_coal[6], MW),
    7 => mole_to_mass_fractions(y_flue_gas_coal[7], MW), #40%
    # nominal scenario
    8 => mole_to_mass_fractions(y_flue_gas_coal[8], MW), #full load
)

# collect flue gas properties
# note: flue gas mass flows are divided by 2 since we assuem 2 carbon capture trains
const FLUE_DATA_COAL = Dict( # flow in kg/hr, T in K, P in kPa
    1 => (flow=2032724/2, wCO2=x_flue_gas_coal[1][:CO2], wH2O=x_flue_gas_coal[1][:H2O], wN2=x_flue_gas_coal[1][:N2], wO2=x_flue_gas_coal[1][:O2], T_in=40 + 273.15, P=110.0), #110 kPa is my own assumption. Typically absorbers are operated close to atmospheric pressure
    2 => (flow=1873641/2, wCO2=x_flue_gas_coal[2][:CO2], wH2O=x_flue_gas_coal[2][:H2O], wN2=x_flue_gas_coal[2][:N2], wO2=x_flue_gas_coal[2][:O2], T_in=40 + 273.15, P=110.0),
    3 => (flow=1707709/2, wCO2=x_flue_gas_coal[3][:CO2], wH2O=x_flue_gas_coal[3][:H2O], wN2=x_flue_gas_coal[3][:N2], wO2=x_flue_gas_coal[3][:O2], T_in=40 + 273.15, P=110.0),
    4 => (flow=1537042/2, wCO2=x_flue_gas_coal[4][:CO2], wH2O=x_flue_gas_coal[4][:H2O], wN2=x_flue_gas_coal[4][:N2], wO2=x_flue_gas_coal[4][:O2], T_in=40 + 273.15, P=110.0),
    5 => (flow=1380067/2, wCO2=x_flue_gas_coal[5][:CO2], wH2O=x_flue_gas_coal[5][:H2O], wN2=x_flue_gas_coal[5][:N2], wO2=x_flue_gas_coal[5][:O2], T_in=40 + 273.15, P=110.0),
    6 => (flow=1234684/2, wCO2=x_flue_gas_coal[6][:CO2], wH2O=x_flue_gas_coal[6][:H2O], wN2=x_flue_gas_coal[6][:N2], wO2=x_flue_gas_coal[6][:O2], T_in=40 + 273.15, P=110.0),
    7 => (flow=1074539/2, wCO2=x_flue_gas_coal[7][:CO2], wH2O=x_flue_gas_coal[7][:H2O], wN2=x_flue_gas_coal[7][:N2], wO2=x_flue_gas_coal[7][:O2], T_in=40 + 273.15, P=110.0),
    # nominal scenario
    8 => (flow=2032724/2, wCO2=x_flue_gas_coal[8][:CO2], wH2O=x_flue_gas_coal[8][:H2O], wN2=x_flue_gas_coal[8][:N2], wO2=x_flue_gas_coal[8][:O2], T_in=40 + 273.15, P=110.0),
)

# solvent lean & rich streams for open-loop initialization / initial guess
const LIQUID_DATA_COAL = Dict( # flow in kg/hr, MEA to H2O weight fraction in CO2 free solvent stream (30% MEA), CO2 loading alpha, temperature T in K
    1 => (lean=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    2 => (lean=(flow=FLUE_DATA_COAL[2].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    3 => (lean=(flow=FLUE_DATA_COAL[3].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    4 => (lean=(flow=FLUE_DATA_COAL[4].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    5 => (lean=(flow=FLUE_DATA_COAL[5].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    6 => (lean=(flow=FLUE_DATA_COAL[6].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    7 => (lean=(flow=FLUE_DATA_COAL[7].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # nominal scenario
    8 => (lean=(flow=FLUE_DATA_COAL[8].flow*3, MEA_H2O=0.3, alpha=0.1, T_in=40 + 273.15), rich=(flow=FLUE_DATA_COAL[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
)

# other scenario dependent data
# Scenario weights from representative coal power unit
const SCENARIO_DATA_COAL = Dict(
    # note: T_intercool is the maximum target solvent temperature at the intercooler
    1 => (T_intercool=50 + 273.15, weight=0.186301),
    2 => (T_intercool=50 + 273.15, weight=0.151941),
    3 => (T_intercool=50 + 273.15, weight=0.0559361),
    4 => (T_intercool=50 + 273.15, weight=0.0471461),
    5 => (T_intercool=50 + 273.15, weight=0.0377854),
    6 => (T_intercool=50 + 273.15, weight=0.0477169),
    7 => (T_intercool=50 + 273.15, weight=0.133447),
    # nominal scenario
    8 => (T_intercool=50 + 273.15, weight=0.8),
)

const INITIAL_VALUES_COAL = Dict(
    1 => (MU_H2O=2300, abs_Q_intercool=-77000, Q_reboiler=250000, Q_cond=-48000, Q_cool=0),
    2 => (MU_H2O=2000, abs_Q_intercool=-68000, Q_reboiler=220000, Q_cond=-45000, Q_cool=0),
    3 => (MU_H2O=1700, abs_Q_intercool=-60000, Q_reboiler=190000, Q_cond=-41000, Q_cool=0),
    4 => (MU_H2O=1400, abs_Q_intercool=-52000, Q_reboiler=170000, Q_cond=-37000, Q_cool=0),
    5 => (MU_H2O=1200, abs_Q_intercool=-44000, Q_reboiler=140000, Q_cond=-33000, Q_cool=0),
    6 => (MU_H2O=1000, abs_Q_intercool=-36000, Q_reboiler=120000, Q_cond=-28000, Q_cool=0),
    7 => (MU_H2O=700,  abs_Q_intercool=-30000, Q_reboiler= 90000, Q_cond=-24000, Q_cool=0),
    # nominal scenario
    8 => (MU_H2O=2300, abs_Q_intercool=-77000, Q_reboiler=250000, Q_cond=-48000, Q_cool=0),
)

# Natural gas combined cycle power plant flue gas
# based on Alcaráz-Calderon AM, González-Díaz MO, Mendez Á, González-Santaló JM, González-Díaz A. Natural gas combined cycle with exhaust gas recirculation and CO2 capture at part-load operation. Journal of the Energy Institute 92(2):370–381 (2019).
y_flue_gas_gas = Dict(
    # high utilization unit
    1  => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
    2  => Dict(:CO2 => 0.0422, :O2  => 0.1167, :H2O => 0.0898, :N2  => 1-0.0422-0.1167-0.0898), #90%
    3  => Dict(:CO2 => 0.0421, :O2  => 0.1169, :H2O => 0.0896, :N2  => 1-0.0421-0.1169-0.0896), #80%
    4  => Dict(:CO2 => 0.0419, :O2  => 0.1174, :H2O => 0.0891, :N2  => 1-0.0419-0.1174-0.0891), #70%
    5  => Dict(:CO2 => 0.0410, :O2  => 0.1193, :H2O => 0.0875, :N2  => 1-0.0410-0.1193-0.0875), #60%
    6  => Dict(:CO2 => 0.0397, :O2  => 0.1222, :H2O => 0.0849, :N2  => 1-0.0397-0.1222-0.0849), #50%
    # nominal scenario
    7  => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
    # intermediate utilization unit
    8  => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
    9  => Dict(:CO2 => 0.0422, :O2  => 0.1167, :H2O => 0.0898, :N2  => 1-0.0422-0.1167-0.0898), #90%
    10 => Dict(:CO2 => 0.0421, :O2  => 0.1169, :H2O => 0.0896, :N2  => 1-0.0421-0.1169-0.0896), #80%
    11 => Dict(:CO2 => 0.0419, :O2  => 0.1174, :H2O => 0.0891, :N2  => 1-0.0419-0.1174-0.0891), #70%
    12 => Dict(:CO2 => 0.0410, :O2  => 0.1193, :H2O => 0.0875, :N2  => 1-0.0410-0.1193-0.0875), #60%
    13 => Dict(:CO2 => 0.0397, :O2  => 0.1222, :H2O => 0.0849, :N2  => 1-0.0397-0.1222-0.0849), #50%
    # nominal scenario
    14 => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
    # low utilization unit
    15 => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
    16 => Dict(:CO2 => 0.0422, :O2  => 0.1167, :H2O => 0.0898, :N2  => 1-0.0422-0.1167-0.0898), #90%
    17 => Dict(:CO2 => 0.0421, :O2  => 0.1169, :H2O => 0.0896, :N2  => 1-0.0421-0.1169-0.0896), #80%
    18 => Dict(:CO2 => 0.0419, :O2  => 0.1174, :H2O => 0.0891, :N2  => 1-0.0419-0.1174-0.0891), #70%
    19 => Dict(:CO2 => 0.0410, :O2  => 0.1193, :H2O => 0.0875, :N2  => 1-0.0410-0.1193-0.0875), #60%
    20 => Dict(:CO2 => 0.0397, :O2  => 0.1222, :H2O => 0.0849, :N2  => 1-0.0397-0.1222-0.0849), #50%
    # nominal scenario
    21 => Dict(:CO2 => 0.0424, :O2  => 0.1163, :H2O => 0.0901, :N2  => 1-0.0424-0.1163-0.0901), #full load
)

x_flue_gas_gas = Dict(
    # high utilization unit
    1 => mole_to_mass_fractions(y_flue_gas_gas[1], MW), #full load
    2 => mole_to_mass_fractions(y_flue_gas_gas[2], MW),
    3 => mole_to_mass_fractions(y_flue_gas_gas[3], MW),
    4 => mole_to_mass_fractions(y_flue_gas_gas[4], MW),
    5 => mole_to_mass_fractions(y_flue_gas_gas[5], MW),
    6 => mole_to_mass_fractions(y_flue_gas_gas[6], MW), #50%
    # nominal scenario
    7 => mole_to_mass_fractions(y_flue_gas_gas[7], MW), #full load
    # intermediate utilization unit
    8 => mole_to_mass_fractions(y_flue_gas_gas[8], MW), #full load
    9 => mole_to_mass_fractions(y_flue_gas_gas[9], MW),
    10 => mole_to_mass_fractions(y_flue_gas_gas[10], MW),
    11 => mole_to_mass_fractions(y_flue_gas_gas[11], MW),
    12 => mole_to_mass_fractions(y_flue_gas_gas[12], MW),
    13 => mole_to_mass_fractions(y_flue_gas_gas[13], MW), #50%
    # nominal scenario
    14 => mole_to_mass_fractions(y_flue_gas_gas[14], MW), #full load
    # low utilization unit
    15 => mole_to_mass_fractions(y_flue_gas_gas[15], MW), #full load
    16 => mole_to_mass_fractions(y_flue_gas_gas[16], MW),
    17 => mole_to_mass_fractions(y_flue_gas_gas[17], MW),
    18 => mole_to_mass_fractions(y_flue_gas_gas[18], MW),
    19 => mole_to_mass_fractions(y_flue_gas_gas[19], MW),
    20 => mole_to_mass_fractions(y_flue_gas_gas[20], MW), #50%
    # nominal scenario
    21 => mole_to_mass_fractions(y_flue_gas_gas[21], MW), #full load
)

# collect flue gas properties
# note: flue gas mass flows are divided by 4 since we assuem 4 carbon capture trains
const FLUE_DATA_GAS = Dict( # flow in kg/hr, T in K, P in kPa
    # high utilization unit
    1  => (flow=1142*3600/4, wCO2=x_flue_gas_gas[1][:CO2],  wH2O=x_flue_gas_gas[1][:H2O],  wN2=x_flue_gas_gas[1][:N2], wO2=x_flue_gas_gas[1][:O2], T_in=40 + 273.15, P=110.0), #110 kPa is my own assumption. Typically absorbers are operated close to atmospheric pressure
    2  => (flow=1067*3600/4, wCO2=x_flue_gas_gas[2][:CO2],  wH2O=x_flue_gas_gas[2][:H2O],  wN2=x_flue_gas_gas[2][:N2], wO2=x_flue_gas_gas[2][:O2], T_in=40 + 273.15, P=110.0),
    3  => (flow= 989*3600/4, wCO2=x_flue_gas_gas[3][:CO2],  wH2O=x_flue_gas_gas[3][:H2O],  wN2=x_flue_gas_gas[3][:N2], wO2=x_flue_gas_gas[3][:O2], T_in=40 + 273.15, P=110.0),
    4  => (flow= 911*3600/4, wCO2=x_flue_gas_gas[4][:CO2],  wH2O=x_flue_gas_gas[4][:H2O],  wN2=x_flue_gas_gas[4][:N2], wO2=x_flue_gas_gas[4][:O2], T_in=40 + 273.15, P=110.0),
    5  => (flow= 842*3600/4, wCO2=x_flue_gas_gas[5][:CO2],  wH2O=x_flue_gas_gas[5][:H2O],  wN2=x_flue_gas_gas[5][:N2], wO2=x_flue_gas_gas[5][:O2], T_in=40 + 273.15, P=110.0),
    6  => (flow= 775*3600/4, wCO2=x_flue_gas_gas[6][:CO2],  wH2O=x_flue_gas_gas[6][:H2O],  wN2=x_flue_gas_gas[6][:N2], wO2=x_flue_gas_gas[6][:O2], T_in=40 + 273.15, P=110.0),
    # nominal scenario
    7  => (flow=1142*3600/4, wCO2=x_flue_gas_gas[7][:CO2],  wH2O=x_flue_gas_gas[7][:H2O],  wN2=x_flue_gas_gas[7][:N2], wO2=x_flue_gas_gas[8][:O2], T_in=40 + 273.15, P=110.0),
    # intermediate utilization unit
    8  => (flow=1142*3600/4, wCO2=x_flue_gas_gas[8][:CO2],  wH2O=x_flue_gas_gas[8][:H2O],  wN2=x_flue_gas_gas[8][:N2], wO2=x_flue_gas_gas[8][:O2], T_in=40 + 273.15, P=110.0), #110 kPa is my own assumption. Typically absorbers are operated close to atmospheric pressure
    9  => (flow=1067*3600/4, wCO2=x_flue_gas_gas[9][:CO2],  wH2O=x_flue_gas_gas[9][:H2O],  wN2=x_flue_gas_gas[9][:N2], wO2=x_flue_gas_gas[9][:O2], T_in=40 + 273.15, P=110.0),
    10 => (flow= 989*3600/4, wCO2=x_flue_gas_gas[10][:CO2], wH2O=x_flue_gas_gas[10][:H2O], wN2=x_flue_gas_gas[10][:N2], wO2=x_flue_gas_gas[10][:O2], T_in=40 + 273.15, P=110.0),
    11 => (flow= 911*3600/4, wCO2=x_flue_gas_gas[11][:CO2], wH2O=x_flue_gas_gas[11][:H2O], wN2=x_flue_gas_gas[11][:N2], wO2=x_flue_gas_gas[11][:O2], T_in=40 + 273.15, P=110.0),
    12 => (flow= 842*3600/4, wCO2=x_flue_gas_gas[12][:CO2], wH2O=x_flue_gas_gas[12][:H2O], wN2=x_flue_gas_gas[12][:N2], wO2=x_flue_gas_gas[12][:O2], T_in=40 + 273.15, P=110.0),
    13 => (flow= 842*3600/4, wCO2=x_flue_gas_gas[13][:CO2], wH2O=x_flue_gas_gas[13][:H2O], wN2=x_flue_gas_gas[13][:N2], wO2=x_flue_gas_gas[13][:O2], T_in=40 + 273.15, P=110.0),
    # nominal scenario
    14 => (flow=1142*3600/4, wCO2=x_flue_gas_gas[14][:CO2], wH2O=x_flue_gas_gas[14][:H2O], wN2=x_flue_gas_gas[14][:N2], wO2=x_flue_gas_gas[14][:O2], T_in=40 + 273.15, P=110.0),
    # low utilization unit
    15 => (flow=1142*3600/4, wCO2=x_flue_gas_gas[15][:CO2], wH2O=x_flue_gas_gas[15][:H2O], wN2=x_flue_gas_gas[15][:N2], wO2=x_flue_gas_gas[15][:O2], T_in=40 + 273.15, P=110.0), #110 kPa is my own assumption. Typically absorbers are operated close to atmospheric pressure
    16 => (flow=1067*3600/4, wCO2=x_flue_gas_gas[16][:CO2], wH2O=x_flue_gas_gas[16][:H2O], wN2=x_flue_gas_gas[16][:N2], wO2=x_flue_gas_gas[16][:O2], T_in=40 + 273.15, P=110.0),
    17 => (flow= 989*3600/4, wCO2=x_flue_gas_gas[17][:CO2], wH2O=x_flue_gas_gas[17][:H2O], wN2=x_flue_gas_gas[17][:N2], wO2=x_flue_gas_gas[17][:O2], T_in=40 + 273.15, P=110.0),
    18 => (flow= 911*3600/4, wCO2=x_flue_gas_gas[18][:CO2], wH2O=x_flue_gas_gas[18][:H2O], wN2=x_flue_gas_gas[18][:N2], wO2=x_flue_gas_gas[18][:O2], T_in=40 + 273.15, P=110.0),
    19 => (flow= 842*3600/4, wCO2=x_flue_gas_gas[19][:CO2], wH2O=x_flue_gas_gas[19][:H2O], wN2=x_flue_gas_gas[19][:N2], wO2=x_flue_gas_gas[19][:O2], T_in=40 + 273.15, P=110.0),
    20 => (flow= 775*3600/4, wCO2=x_flue_gas_gas[20][:CO2], wH2O=x_flue_gas_gas[20][:H2O], wN2=x_flue_gas_gas[20][:N2], wO2=x_flue_gas_gas[20][:O2], T_in=40 + 273.15, P=110.0),
    # nominal scenario
    21 => (flow=1142*3600/4, wCO2=x_flue_gas_gas[21][:CO2], wH2O=x_flue_gas_gas[21][:H2O], wN2=x_flue_gas_gas[21][:N2], wO2=x_flue_gas_gas[21][:O2], T_in=40 + 273.15, P=110.0),
)
# Liquid solvent data: lean & rich streams for each scenario, required for open-loop initialization and as initial guess
const LIQUID_DATA_GAS = Dict(
    # high utilization unit
    1 => (lean=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    2 => (lean=(flow=FLUE_DATA_GAS[2].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    3 => (lean=(flow=FLUE_DATA_GAS[3].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    4 => (lean=(flow=FLUE_DATA_GAS[4].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    5 => (lean=(flow=FLUE_DATA_GAS[5].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    6 => (lean=(flow=FLUE_DATA_GAS[6].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # nominal scenario
    7 => (lean=(flow=FLUE_DATA_GAS[7].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # intermediate utilization unit
    8 => (lean=(flow=FLUE_DATA_GAS[8].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    9 => (lean=(flow=FLUE_DATA_GAS[9].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    10 => (lean=(flow=FLUE_DATA_GAS[10].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    11 => (lean=(flow=FLUE_DATA_GAS[11].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    12 => (lean=(flow=FLUE_DATA_GAS[12].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    13 => (lean=(flow=FLUE_DATA_GAS[13].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # nominal scenario
    14 => (lean=(flow=FLUE_DATA_GAS[14].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # T H Wharton, peaker
    15 => (lean=(flow=FLUE_DATA_GAS[15].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    16 => (lean=(flow=FLUE_DATA_GAS[16].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    17 => (lean=(flow=FLUE_DATA_GAS[17].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    18 => (lean=(flow=FLUE_DATA_GAS[18].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    19 => (lean=(flow=FLUE_DATA_GAS[19].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    20 => (lean=(flow=FLUE_DATA_GAS[20].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
    # nominal scenario
    21 => (lean=(flow=FLUE_DATA_GAS[21].flow*3, MEA_H2O=0.3, alpha=0.2, T_in=40 + 273.15), rich=(flow=FLUE_DATA_GAS[1].flow*3, MEA_H2O=0.3, alpha=0.47)),
)

# other scenario dependent data
# Scenario weights of representative natural gas combined cycle power units
const SCENARIO_DATA_GAS = Dict( 
    # note: T_intercool is the maximum target solvent temperature at the intercooler
    # high utilization unit
    1  => (T_intercool=50 + 273.15, weight=0.046118721461187215),
    2  => (T_intercool=50 + 273.15, weight=0.32853881278538815),
    3  => (T_intercool=50 + 273.15, weight=0.18264840182648403),
    4  => (T_intercool=50 + 273.15, weight=0.09908675799086758),
    5  => (T_intercool=50 + 273.15, weight=0.03356164383561644),
    6  => (T_intercool=50 + 273.15, weight=0.038812785388127855),
    # nominal scenario
    7  => (T_intercool=50 + 273.15, weight=0.8),
    # intermediate utilization unit
    8  => (T_intercool=50 + 273.15, weight=0.008105022831050228),
    9  => (T_intercool=50 + 273.15, weight=0.1797945205479452),
    10 => (T_intercool=50 + 273.15, weight=0.29360730593607308),
    11 => (T_intercool=50 + 273.15, weight=0.075228310502283104),
    12 => (T_intercool=50 + 273.15, weight=0.04315068493150685),
    13 => (T_intercool=50 + 273.15, weight=0.051826484018264845),
    # nominal scenario
    14 => (T_intercool=50 + 273.15, weight=0.8),
    # low utilization unit
    15 => (T_intercool=50 + 273.15, weight=0.005707762557077626),
    16 => (T_intercool=50 + 273.15, weight=0.021575342465753424),
    17 => (T_intercool=50 + 273.15, weight=0.07602739726027397),
    18 => (T_intercool=50 + 273.15, weight=0.017237442922374429),
    19 => (T_intercool=50 + 273.15, weight=0.00821917808219178),
    20 => (T_intercool=50 + 273.15, weight=0.009246575342465754),
    # nominal scenario
    21 => (T_intercool=50 + 273.15, weight=0.8),
)

const INITIAL_VALUES_GAS = Dict(
    # base load
    1  => (MU_H2O=1200, abs_Q_intercool=-27500, Q_reboiler=75000, Q_cond=-44000, Q_cool=0),
    2  => (MU_H2O=1100, abs_Q_intercool=-25000, Q_reboiler=67500, Q_cond=-41000, Q_cool=0),
    3  => (MU_H2O=1000, abs_Q_intercool=-22500, Q_reboiler=62000, Q_cond=-38000, Q_cool=0),
    4  => (MU_H2O= 900, abs_Q_intercool=-20000, Q_reboiler=55000, Q_cond=-35000, Q_cool=0),
    5  => (MU_H2O= 800, abs_Q_intercool=-17500, Q_reboiler=49000, Q_cond=-32000, Q_cool=0),
    6  => (MU_H2O= 700, abs_Q_intercool=-15000, Q_reboiler=42000, Q_cond=-29000, Q_cool=0),
    # nominal scenario
    7  => (MU_H2O=1200, abs_Q_intercool=-27500, Q_reboiler=75000, Q_cond=-44000, Q_cool=0),
    # intermediate load
    8  => (MU_H2O=1200, abs_Q_intercool=-27500, Q_reboiler=75000, Q_cond=-44000, Q_cool=0),
    9  => (MU_H2O=1100, abs_Q_intercool=-25000, Q_reboiler=67500, Q_cond=-41000, Q_cool=0),
    10 => (MU_H2O=1000, abs_Q_intercool=-22500, Q_reboiler=62000, Q_cond=-38000, Q_cool=0),
    11 => (MU_H2O= 900, abs_Q_intercool=-20000, Q_reboiler=55000, Q_cond=-35000, Q_cool=0),
    12 => (MU_H2O= 800, abs_Q_intercool=-17500, Q_reboiler=49000, Q_cond=-32000, Q_cool=0),
    13 => (MU_H2O= 700, abs_Q_intercool=-15000, Q_reboiler=42000, Q_cond=-29000, Q_cool=0),
    # nominal scenario
    14 => (MU_H2O=1200, abs_Q_intercool=-27500, Q_reboiler=75000, Q_cond=-44000, Q_cool=0),
    # low utilization unit
    15 => (MU_H2O=2300, abs_Q_intercool=-77000, Q_reboiler=250000, Q_cond=-48000, Q_cool=0),
    16 => (MU_H2O=2000, abs_Q_intercool=-68000, Q_reboiler=220000, Q_cond=-45000, Q_cool=0),
    17 => (MU_H2O=1700, abs_Q_intercool=-60000, Q_reboiler=190000, Q_cond=-41000, Q_cool=0),
    18 => (MU_H2O=1400, abs_Q_intercool=-52000, Q_reboiler=170000, Q_cond=-37000, Q_cool=0),
    19 => (MU_H2O=1200, abs_Q_intercool=-44000, Q_reboiler=140000, Q_cond=-33000, Q_cool=0),
    20 => (MU_H2O=1000, abs_Q_intercool=-36000, Q_reboiler=120000, Q_cond=-28000, Q_cool=0),
    # nominal scenario
    21 => (MU_H2O=2300, abs_Q_intercool=-77000, Q_reboiler=250000, Q_cond=-48000, Q_cool=0),
)

# Summary dictionaries for easy access to scenario dependent data
# y_flue_gas = Dict("coal" => y_flue_gas_coal, "gas" => y_flue_gas_gas) #todo, no longer needed, delete?
# x_flue_gas = Dict("coal" => x_flue_gas_coal)
const FLUE_DATA =      Dict("coal" => FLUE_DATA_COAL,      "gas" => FLUE_DATA_GAS,)
const LIQUID_DATA =    Dict("coal" => LIQUID_DATA_COAL,    "gas" => LIQUID_DATA_GAS,)
const SCENARIO_DATA =  Dict("coal" => SCENARIO_DATA_COAL,  "gas" => SCENARIO_DATA_GAS,)
const INITIAL_VALUES = Dict("coal" => INITIAL_VALUES_COAL, "gas" => INITIAL_VALUES_GAS,)

# processing functions
function compute_liquid_stream(flow_kg_hr, MEA_H2O_wt, alpha)
    """
    Computes a liquid streams apparent molar flows and mole fractions, assuming only H2O, CO2, and MEA
    from given mass flow rate (flow_kg_hr), MEA/H2O weight fraction (MEA_H2O_wt), and molar CO2 loading (alpha).
    
    m_H2O = m_MEA * (1 / MEA_H2O_wt - 1)
    m_CO2 = alpha * (m_MEA / MW[:MEA]) * MW[:CO2]
    m_MEA + m_H2O + m_CO2 = flow_kg_hr
    """
    m_MEA = flow_kg_hr / (1 / MEA_H2O_wt + alpha * MW[:CO2] / MW[:MEA]) # kg/hr
    m_H2O = m_MEA * (1 / MEA_H2O_wt - 1) # kg/hr
    m_CO2 = alpha * (m_MEA * 1000 / MW[:MEA]) * MW[:CO2] / 1000 # kg/hr
    n_MEA = (m_MEA * 1000 / 3600) / MW[:MEA]  # mol/s
    n_H2O = (m_H2O * 1000 / 3600) / MW[:H2O]  # mol/s
    n_CO2 = (m_CO2 * 1000 / 3600) / MW[:CO2]  # mol/s
    mols = [n_MEA, n_H2O, n_CO2]
    yf = mols ./ sum(mols)
    return sum(mols), Dict(:MEA => yf[1], :H2O => yf[2], :CO2 => yf[3])
end

function get_inout_data(case::Int, case_study="coal")
    """
    Returns inlet/outlet solvent and flue gas conditions for the selected part-load case of the case study to the optimization model
    case: Integer, the part-load scenario number
    case_study: String, either "coal" or "gas"
    """
    ld = LIQUID_DATA[case_study][case]
    od = SCENARIO_DATA[case_study][case]
    iv = INITIAL_VALUES[case_study][case]

    # solvent data: inlet and outlet for initialization
    L_in, yL_in = compute_liquid_stream(ld.lean.flow, ld.lean.MEA_H2O, ld.lean.alpha)
    L_out, yL_out = compute_liquid_stream(ld.rich.flow, ld.rich.MEA_H2O, ld.rich.alpha)

    # flue gas data
    fd = FLUE_DATA[case_study][case]
    w = [fd.wCO2, fd.wH2O, fd.wN2, fd.wO2]
    mol_hr = fd.flow * 1000 .* w ./ [MW[:CO2], MW[:H2O], MW[:N2], MW[:O2]]
    mol_s = mol_hr ./ 3600
    yV = mol_s ./ sum(mol_s)
    V_flow = sum(mol_s)
    yV_dict = Dict(:CO2 => yV[1], :H2O => yV[2], :N2 => yV[3], :O2 => yV[4])

    return (
        liquid_inlet=(flow=L_in, y=yL_in, T=ld.lean.T_in, alpha=ld.lean.alpha),
        liquid_outlet=(flow=L_out, y=yL_out),
        vapor=(flow=V_flow, y=yV_dict, T=fd.T_in, P=fd.P),
        T_intercool=od.T_intercool,
		Q_reboiler=iv.Q_reboiler,
        initial_values = (MU_H2O=iv.MU_H2O, abs_Q_intercool=iv.abs_Q_intercool, Q_cond=iv.Q_cond, Q_cool=iv.Q_cool),
        weight=od.weight	
    )
end

function make_inout_data(cases::AbstractVector{<:Integer}; case_study="coal")
    """
    Wrapper for generating a list of inlet/outlet data for the specified cases.
    """
    @assert !isempty(cases) "cases must not be empty."
    return [get_inout_data(c, case_study) for c in cases]
end

function calc_inlet_species(inout)
    """
    calculates inlet species concentrations for the solvent (initial values)
    inout must be in the format returned by get_inout_data
    """
    # Retrieve feed data 
    F_L = inout.liquid_inlet.flow
    T_L_in = inout.liquid_inlet.T
    y_in = inout.liquid_inlet.y

    # Feed molar flows
    n0_MEA = y_in[:MEA] * F_L
    n0_CO2 = y_in[:CO2] * F_L
    n0_H2O = y_in[:H2O] * F_L

    # Build JuMP model
    model1 = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model1, "bound_relax_factor", 0.0)
    set_optimizer_attribute(model1, "linear_solver", "ma57")
    set_optimizer_attribute(model1, "print_level", 1)

    @variables(model1, begin
        L >= 0
        0 >= log_MEA >= -50
        0 >= log_MEA2p >= -50
        0 >= log_Hp >= -50
        0 >= log_Carbamate >= -50
        0 >= log_HCO3m >= -50
        0 >= log_OHm >= -50
        0 >= log_CO3_2m >= -50
        0 >= log_CO2 >= -50
        0 >= log_H2O >= -50
        ξ[1:5]  # reaction extents
    end)
    MEA = exp(log_MEA)
    MEA2p = exp(log_MEA2p)
    Hp = exp(log_Hp)
    Carbamate = exp(log_Carbamate)
    HCO3m = exp(log_HCO3m)
    OHm = exp(log_OHm)
    CO3_2m = exp(log_CO3_2m)
    CO2 = exp(log_CO2)
    H2O = exp(log_H2O)

    # Mole balance constraints
    @constraint(model1, MEA * L == n0_MEA + ξ[1] + ξ[4])
    @constraint(model1, MEA2p * L == -ξ[1])
    @constraint(model1, Hp * L == ξ[1] + ξ[2] + ξ[3] + ξ[5])
    @constraint(model1, Carbamate * L == -ξ[4])
    @constraint(model1, HCO3m * L == -ξ[2] + ξ[4] + ξ[5])
    @constraint(model1, OHm * L == ξ[3])
    @constraint(model1, CO3_2m * L == ξ[2])
    @constraint(model1, CO2 * L == n0_CO2 - ξ[5])
    @constraint(model1, H2O * L == n0_H2O - ξ[3] - ξ[4] - ξ[5])

    # Equilibrium constraints (in log form)
    lnK1, lnK2, lnK3, lnK4, lnK5 = get_lnK(T_L_in)
    @constraint(model1, lnK1 + log_MEA2p == log_MEA + log_Hp)
    @constraint(model1, lnK2 + log_HCO3m == log_CO3_2m + log_Hp)
    @constraint(model1, lnK3 == log_OHm + log_Hp)
    @constraint(model1, lnK4 + log_Carbamate == log_HCO3m + log_MEA)
    @constraint(model1, lnK5 + log_CO2 == log_HCO3m + log_Hp)

    # Mole-fraction summation
    @constraint(model1, MEA + MEA2p + Hp + Carbamate + HCO3m + OHm + CO3_2m + CO2 + H2O == 1)

    # Initialization
    set_start_value(L, F_L)
    set_start_value(log_MEA, log(y_in[:MEA]))
    set_start_value(log_MEA2p, log(1e-6))
    set_start_value(log_Hp, log(1e-6))
    set_start_value(log_Carbamate, log(1e-6))
    set_start_value(log_HCO3m, log(1e-6))
    set_start_value(log_OHm, log(1e-6))
    set_start_value(log_CO3_2m, log(1e-6))
    set_start_value(log_CO2, log(y_in[:CO2]))
    set_start_value(log_H2O, log(y_in[:H2O]))
    set_start_value.(ξ, zeros(5))

    # Solve
    optimize!(model1)
    println("speciation status = ", termination_status(model1), " / ", raw_status(model1))

    # Extract results
    L_val = value(L)
    x = Dict(
        :MEA => value(MEA), :MEA2p => value(MEA2p), :Hp => value(Hp),
        :Carbamate => value(Carbamate), :HCO3m => value(HCO3m),
        :OHm => value(OHm), :CO3_2m => value(CO3_2m),
        :CO2 => value(CO2), :H2O => value(H2O)
    )

    return L_val, x
end
