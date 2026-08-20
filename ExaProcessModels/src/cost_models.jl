"""
This script collects cost parameters and models for the technoeconomic assessment of the carbon capture plant.
The economic cost model is based on the methodology presented in:
[1] Theis, J. Quality Guidelines for Energy System Studies: Cost Estimation Methodology for NETL Assessments of Power Plant Performance. https://www.osti.gov/servlets/purl/1513278/ (2021) doi:10.2172/1513278.
[2] Seider, W. D. et al. Product and Process Design Principles: Synthesis, Analysis, and Evaluation. (John Wiley & Sons Inc., New York, 2017).
The assessment further builds on detailed presented in:
[3] Hughes, S. & Zoelle, A. Cost of Capturing CO2 from Industrial Sources. (2023).
[4] Cvetic, P. et al. Natural Gas Combined Cycle (NGCC) Power Plants with Carbon Capture and Exhaust Gas Recycle (EGR). DOE/NETL--2023/3922, 2251495 https://www.osti.gov/servlets/purl/2251495/ (2023) doi:10.2172/2251495.
[5] Electric Power Research Institute & Elliott, W. Front-End Engineering Design (FEED) Study for a Carbon Capture Plant Retrofit to a Natural Gas-Fired Gas Turbine Combined Cycle Power Plant (2x2x1 Duct-Fired 758-MWe Facility with F Class Turbines). DOE-BECHTEL-FE0031848-App.Volume, 1836563 https://www.osti.gov/servlets/purl/1836563/ (2021) doi:10.2172/1836563.
[6] Towler, G. P. & Sinnott, R. K. Chemical Engineering Design: Principles, Practice and Economics of Plant and Process Design. (Butterworth-Heinemann, Oxford [England] ; Cambridge, MA, 2022).
[7] Woods, D. R. Rules of Thumb in Engineering Practice. (Wiley, 2007).
[8] U.S. Bureau of Labor Statistics. 2026. Producer Price Index by Commodity: Chemicals and Allied Products: Industrial Chemicals [WPU061]. https://fred.stlouisfed.org/series/WPU061.
[9] Haslbeck et al. 2010. Cost and Performance Baseline for Fossil Energy Plants Volume 1: Bituminous Coal and Natural Gas to Electricity Revision 2.
[10] Moser et al. 2020 Results of the 18-month test with MEA at the post-combustion capture pilot plant at Niederaussem – new impetus to solvent management, emissions and dynamic behaviour
[11] Tsai RE. Mass transfer area of structured packing. Dissertation. The University of Texas at Austin (2010).

The capture plant is treated as an end of pipe solution without integration into the power plant.
"""

####################################################################################
# parameters
# Chemical Engineering Plant Cost Index (CEPCI) for cost escalation of equipment
const CEPCI = Dict(
    2000 => 394.1,
    2001 => 394.3,
    2002 => 395.6,
    2003 => 402.0,
    2004 => 444.2,
    2005 => 468.2,
    2006 => 499.6,
    2007 => 525.4,
    2008 => 575.4,
    2009 => 521.9,
    2010 => 550.8,
    2011 => 585.7,
    2012 => 584.6,
    2013 => 567.3,
    2014 => 576.1,
    2015 => 556.8,
    2016 => 541.7,
    2017 => 567.5,
    2018 => 603.1,
    2019 => 607.5,
    2020 => 596.2,
    2021 => 708.8,
    2022 => 816.0,
    2023 => 803.4
)

# producer price index (PPI) for cost escalation of chemicals, [8]
PPI = Dict(
    2007 => 226.4166667,
    2008 => 274.5666667,
    2009 => 234.0833333,
    2010 => 269.1750,
    2011 => 324.7166667,
    2012 => 306.9333333,
    2013 => 301.2166667,
    2014 => 288.9583333,
    2015 => 242.2416667,
    2016 => 228.125,
    2017 => 254.05,
    2018 => 273.9916667,
    2019 => 252.5416667,
    2020 => 226.3,
    2021 => 309.3120833,
    2022 => 353.62525,
    2023 => 318.3385,
    2024 => 301.6849167,
    2025 => 297.98525,
)

const in_per_m = 39.37 # inches per meter
const ft_per_m = in_per_m/12 # feet per meter
const lb_per_kg = 2.20462 # lb per kg
const btu_per_kWh = 3412.14 # Btu per kWh
const gal_per_m3 = 264.172 # gallons per cubic meter
const hp_per_kW = 1.34102 # hp/kW

const availability = 0.85 # plant availability for converting cost factors

# Fixed charge rate (real) over 3 years (natural gas power plant) 
const FCR = 0.0707 # [1], Exhibit 3-5
# TASC to TOC ratio (real) over 3 years (natural gas power plant)
TASC_TOC_ratio = 1.093 # [1], Exhibit 3-7
# Capital charge factor
const CCF = TASC_TOC_ratio * FCR # [1]

# Total plant cost multipliers
# Engineering, procurement, and construction (EPC) cost factor, typically 15-20% of bare erected cost (BEC)
epc_factor = 0.2 # [4], fraction of EPC and BEC for CANSOLV CO2 removal system
contingency_process = 0.18 # [4]
contingency_project = 0.276 # [4]
retrofit_factor = 1.05 # [3]

# operational expenditures (OPEX) parameters
# fixed operation costs
opex_labor_rate = 38.5 # USD/hr, [4]
n_operators = 1.3 # difference in operators between NGCC with and without CCS [4]
opex_labor_burden = 0.3 # 30% labor burdens, [4]
opex_labor_overhead = 0.25 # 25% overhead on operating and maintenance labor, [4]
opex_maintenance_labor_TPC = 0.0076 # fraction of TPC, [4]
opex_maintenance_material_TPC = 0.0114 # fraction of TPC, [4]
opex_tax_insurance_TPC = 0.02 # property taxes and insurance as a fraction of TPC, [4]
fixed_opex_labor = opex_labor_rate * n_operators * 8760 * (1 + opex_labor_burden + opex_labor_overhead) # USD/year

# variable operation costs
# we scale variable operation costs for water, chemicals, and waste treatment by tonnes of CO2 captured, based on case B31B.95
annual_CO2_captured_literature = 236321/1000*8760*availability # CO2 captured in [4] in tCO2/year
opex_water_per_CO2 = 3238*365*availability*1.9/annual_CO2_captured_literature # USD/tCO2 captured
opex_makeup_water_treatment_per_CO2 = 9.6*365*availability*550/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_ammonia_per_CO2 = 3.5*365*availability*300/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_catalyst_per_CO2 = 3.1*365*availability*150/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_triethyleneglycol_per_CO2 = 466*365*availability*6.8/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_waste_catalysis_per_CO2 = 3.1*365*availability*2.5/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_waste_triethyleneglycol_per_CO2 = 466*365*availability*0.35/annual_CO2_captured_literature # USD/tCO2 captured, [4]
opex_thermalreclaimer_per_CO2 = 1.831*365*availability*38/annual_CO2_captured_literature # USD/tCO2 captured, [4]
price_MEA = 2480/1000 * PPI[2023]/PPI[2007] # USD/kg, 2249.89 USD (2007) per short ton, [9]
opex_MEA_per_CO2 = 0.46 * price_MEA # degradation rate of MEA of 0.46 kg/tCO2, [10]
# variable operating cost factor scaling with tonnes of CO2 captured
opex_per_CO2 = (opex_water_per_CO2
                + opex_makeup_water_treatment_per_CO2
                + opex_ammonia_per_CO2
                + opex_catalyst_per_CO2
                + opex_triethyleneglycol_per_CO2 
                + opex_waste_catalysis_per_CO2
                + opex_waste_triethyleneglycol_per_CO2
                + opex_thermalreclaimer_per_CO2
                + opex_MEA_per_CO2) # USD/tCO2 captured

# utilities
price_natural_gas = 0.213 / 38.25 * 3.6 / btu_per_kWh *1e6 # USD/MMBtu HHV, [2]
price_steam = 15.3 # USD/1000kg 150psig, [2]
h_v_150psig = 1995.3 # heat of vaporization of steam at 150 psig, kJ/kg
price_electricity = 70.00 # USD/MWh, [2]
price_cooling_water = 0.027  # USD/m3, [2]

# Owner's costs
# Pre-production costs
preproduction_costs_labor = 6/12*fixed_opex_labor # [4]
preproduction_costs_per_TPC = (6/12*opex_maintenance_labor_TPC*(1+opex_labor_overhead) 
                               + 1/12*opex_maintenance_material_TPC/availability
                               + 0.02)  # [4]
preproduction_costs_per_CO2 = 1/12*opex_per_CO2/availability # 1 month consumables and waste disposal [4], multiply by annual CO2 captured
# preproduction costs of fuels neglected for simplicity

# Inventory
inventory_costs_per_CO2 = 2/12*(opex_water_per_CO2+opex_makeup_water_treatment_per_CO2+opex_ammonia_per_CO2+opex_catalyst_per_CO2+opex_MEA_per_CO2+opex_triethyleneglycol_per_CO2)/availability # inventory costs per tonne of CO2 captured [4], multiply by annual CO2 captured
# inventory costs of fuels neglected, as we price steam consumption directly, assuming the price to cover all expenses for steam generation
inventory_costs_per_TPC = 0.005  # inventory costs per TPC [4], multiply by TPC

# other costs
# we neglect initial cost for catalysts and chemicals
other_costs_per_TPC = 0.15 # other owners costs [4], multiply by TPC
financing_costs_per_TPC = 0.027 # financing costs [4], multiply by TPC

####################################################################################
# Bare erected cost (BEC) functions for major equipment
function bec_column(h, d, h_add; maximum_allowable_stress=20000)
    # returns the bare erected cost of the absorber column based on packing height and diameter in 2023 dollars
    # h: packing height in m
    # d: packing diameter in m
    # h_add: height of the column in addition to the packing, in m (water wash sections, etc.)
    # based on [2] Pressure Vessels and Towers for Distillation, Absorption, and Stripping and [5]
    # we assume stainless steel to avoid corrosion

    h_total = h + h_add # total height
    h_total_in = h_total*in_per_m # total height in inches
    h_total_ft = h_total_in / 12 # total height in feet
    d_in = d*in_per_m # diameter in inches
    d_ft = d_in / 12 # diameter in feet

    t_base = 0.5 # minimum wall thickness in inches, [2]
    t_wind_earthquake = 0.22 * (d_in+18)*h_total_in^2/maximum_allowable_stress/d_in^2  # thickness required for wind and earthquake loads in inches, [2]
    t_corrosion = 1/8 # corrosion allowance for stainless steel vessels in inches, [2]
    t_total = t_base + t_wind_earthquake + t_corrosion # total required thickness in inches
    t_total_m = t_total/in_per_m
    # for simplicity, we do not round the thickness to discrete values

    density_steel = 8000 # kg/m3, [6]

    weight_shell = density_steel * π * (d+t_total_m)*(d*0.8+h_total) * t_total_m # in kg, weight of shell including heads, [2]
    weight_shell_lb = weight_shell * lb_per_kg # in lb
    cost_carbon_steel = exp(10.5449-0.4672*log(weight_shell_lb)+0.05482*log(weight_shell_lb)^2) # USD, cost of carbon steel vessel, [2]
    material_factor = 1.7 # material factor for stainless steel, [2]
    fob_cost_stainless_steel = material_factor * cost_carbon_steel  # USD, free on board cost of stainless steel vessel, [2]
    fob_cost_platforms_ladders = 341 * d_ft^0.63316 * h_total_ft^0.80161 # USD, cost of platforms and ladders, [2]
    bare_module_factor = 4.16 # bare module factor for stainless steel vessels including installation, [2]
    bare_erected_cost = bare_module_factor * (fob_cost_stainless_steel + fob_cost_platforms_ladders) * CEPCI[2023]/CEPCI[2013]  # USD
    return bare_erected_cost
end

function bec_packing(h, d)
    # returns the bare erected cost of the column packing based on height and diameter in 2023 dollars
    # h: packing height in m
    # d: packing diameter in m
    # we assume Sulzer Mellapak 250Y

    volume_m3 = π*(d/2)^2 * h  # in m3
    volumetric_cost_2008 = 50 * ft_per_m^3  # USD/m3 in 2008 dollars [11]
    packing_cost = volume_m3 * volumetric_cost_2008 * CEPCI[2023]/CEPCI[2008]  # USD, packing cost in 2023 dollars
    bare_module_factor = 4.16 # bare module factor for packing including installation, [2]
    bare_erected_cost = bare_module_factor * packing_cost  # USD
    return bare_erected_cost
end

function bec_pump(v, h, rho)
    # returns the bare erected cost of a centrifugal pump in 2023 dollars based on [2]
    # v: volume flow in gpm
    # h: pump head in ft
    # rho: fluid density in lb/gal

    # pump
    size_factor = v * h^0.5
    base_cost_pump = exp(12.1656 - 1.1448*log(size_factor) + 0.0862*log(size_factor)^2)  # USD, base cost in 2013 dollars, [2]
    material_factor = 2.0  # material factor for stainless steel, [2]
    type_factor = 2.0  # type factor for centrifugal pump, 1800 rpm, 250-5000 gal/min, 50-500 ft head, [2]
    fob_cost_pump = base_cost_pump * material_factor * type_factor  # in 2013 USD

    # electric motor
    eta_P = -0.316 + 0.24015 * log(v) - 0.01199 * log(v)^2  # pump efficiency, [2]
    P_B = v * h * rho / (33000 * eta_P)  # pump brake horsepower
    eta_M = 0.80 + 0.0319 * log(P_B) - 0.00182 * log(P_B)^2  # motor efficiency, [2]
    power_consumption = P_B / eta_M  # motor power consumption in hp
    base_cost_motor = exp(5.9332 + 0.16829 * log(power_consumption) 
                                 - 0.110056 * log(power_consumption)^2 
                                 + 0.071413 * log(power_consumption)^3 
                                 - 0.0063788 * log(power_consumption)^4)  # base cost in 2013 USD, [2]
    type_factor_motor = 1  # type factor motor, open drip proof enclosure, [2]
    fob_cost_motor = base_cost_motor * type_factor_motor  # 2013 USD

    bare_module_factor = 3.3 # bare module factor for pumps including installation, [2]
    bare_erected_cost = bare_module_factor * (fob_cost_pump + fob_cost_motor) * CEPCI[2023]/CEPCI[2013]  # in 2023 USD
    return bare_erected_cost
end

function bec_heat_exchanger(A)
    # returns the bare erected cost of a plate and frame heat exchanger in 2023 dollars
    # A: heat transfer area in ft2

    fob_cost_heat_exchanger = 10070 * A^0.42  # $ cost in 2013 dollars, [2]
    bare_module_factor = 3.17 # bare module factor, assume the same as for a tube and shell HEX, [2]
    bare_erected_cost = bare_module_factor * fob_cost_heat_exchanger * CEPCI[2023]/CEPCI[2013]  # in 2023 USD
    return bare_erected_cost
end

function bec_hex_tubeandshell(A; reboiler=false)
    # returns the bare erected cost of a tube and shell heat exchanger or kettle reboiler in 2023 USD
    # for reboilers: we assume the same cost correlation as for a shell and tube heat exchanger [2] + a reboiler factor [7]
    # A: heat transfer area in ft2

    base_cost_hex = exp(11.4185 - 0.9228 * log(A) + 0.09861 * log(A)^2)  # in 2013 USD, [2]
    material_factor = 2.7 + (A/100)^0.07 # material factor for stainless steel/stainless steel, [2]
    tube_factor = 1.0  # tube side factor for water, [2], assume no correction due to lack of data on tube length
    pressure_factor = 1.0  # pressure factor, low pressure, [2]
    bare_module_factor = 3.17 # [2]
    reboiler_factor = 1  # default
    if reboiler
        reboiler_factor = 1.35 # reboiler factor, [7]
    end
    fob_cost_hex = base_cost_hex * material_factor * tube_factor * pressure_factor * reboiler_factor # in 2013 USD
    bare_erected_cost = bare_module_factor * fob_cost_hex * CEPCI[2023]/CEPCI[2013]  # in 2023 USD
    return bare_erected_cost
end

function bec_compressor(P)
    # returns the bare erected cost of a multistage centrifugal compressor in 2023 dollars
    # P: compressor rated power in hp

    base_cost_compressor = exp(9.1553) * P^0.63  # base cost in 2013 USD, [2], exact reformulation to remove ln() and exp()
    material_factor = 2.5  # stainless steel, [2]
    fob_cost_compressor = base_cost_compressor * material_factor  # in 2013 USD
    bare_module_factor = 2.15 # [2]
    bare_erected_cost = bare_module_factor * fob_cost_compressor * CEPCI[2023]/CEPCI[2013]  # in 2023 USD
    return bare_erected_cost
end

tpc_factor_compressor = (8637+10365)/43186 # fraction of TPC and bare erected cost for compressor in [4], case B31B.95; includes EPC and contingencies

####################################################################################
# Total plant cost
function calc_total_plant_cost(bare_erected_cost)
    # returns the total plant cost based on bare erected cost
    # bare_erected_cost: bare erected cost in USD

    epc_cost = epc_factor * bare_erected_cost
    contingency_cost = (contingency_process + contingency_project) * bare_erected_cost
    return (bare_erected_cost + epc_cost + contingency_cost) / retrofit_factor
end

####################################################################################
# OPEX
function calc_opex(TPC, steam_consumption, electricity_consumption, coolingwater_consumption, annual_CO2_captured, price_factor, S, n_trains)
    # returns the annual operating expenditures based on total plant cost and consumed utilities
    # TPC: total plant cost in USD
    # steam_consumption[s]:  steam consumption in kWh/year per scenario
    # electricity_consumption[s]: electricity consumption in MWh/year per scenario
    # coolingwater_consumption: annual cooling water consumption in m3/year
    # annual_CO2_captured: annual CO2 captured in tCO2/year
    # price factor: a dictionary of factors to adjust electricity and steam prices in each scenario

    # fixed opex contributions
    fixed_opex_per_TPC = opex_maintenance_labor_TPC*(1+opex_labor_overhead) + opex_maintenance_material_TPC + opex_tax_insurance_TPC
    fixed_opex = fixed_opex_labor/n_trains + fixed_opex_per_TPC * TPC  # USD/year

    # variable opex
    variable_opex = annual_CO2_captured * opex_per_CO2 # USD/year

    # energy related opex
    # convert steam consumption from kWh/year to tonne/year based on enthalpy of steam at 150 psig in kJ/kg
    # steam_consumption in tonne = steam_consumption[s] [kWh/year] * 3600 [s/h] / h_v_150psig [kJ/kg] / 1000 [kg/tonne]
    energy_opex = sum(price_steam * price_factor[s] * (steam_consumption[s] * 3600/h_v_150psig)/1000 for s in S) + sum(price_electricity * price_factor[s] * electricity_consumption[s] for s in S) + price_cooling_water * coolingwater_consumption  # in USD/year

    # total opex
    return fixed_opex + variable_opex + energy_opex  # in USD/year
end

####################################################################################
# Total overnight cost (TOC)
function calc_toc(TPC, annual_CO2_captured, n_trains)
    # returns the total overnight costs based on total plant cost
    # TPC: total plant cost in USD
    # annual_CO2_captured: annual CO2 captured in tCO2/year
    # n_trains: number of parallel trains in the capture plant to account for fixed costs

    # pre-production cost
    preproduction_costs_total = preproduction_costs_labor/n_trains + preproduction_costs_per_TPC * TPC + preproduction_costs_per_CO2 * annual_CO2_captured # in USD

    # inventory cost
    inventory_costs_total = inventory_costs_per_CO2 * annual_CO2_captured + inventory_costs_per_TPC * TPC # in USD

    # other costs
    other_costs_total = other_costs_per_TPC * TPC + financing_costs_per_TPC * TPC  # in USD

    # total capital expenditures
    return TPC + preproduction_costs_total + inventory_costs_total + other_costs_total  # in USD
end

#####################################################################################
# total annualized cost
function calc_tac(TOC, OPEX, CCF=CCF)
    # TOC: total overnight cost in USD
    # OPEX: total operating expenditures in USD
    # CCF: capital charge factor

    annualized_capex = CCF * TOC # USD/year
    return annualized_capex + OPEX # USD/year
end

# Cost of carbon captured
function calc_coc(TAC, annual_CO2_captured)
    # returns the cost of carbon captured in USD/tCO2
    # the cost of carbon captured accounts for the OPEX and annualized CAPEX
    # life-cycle emissions, and emissions from utilities are neglected
    # TAC: total annualized cost in USD/year
    # annual_CO2_captured: annual CO2 captured in tCO2/year

    return TAC / annual_CO2_captured  # USD/tCO2
end
