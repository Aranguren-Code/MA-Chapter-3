# QUESTIONS - PICKING DEPOTS WITHOUT USING THEM ?!
#-------------------------------------------------------------------------------
# PART OF LEXICOGRAPHIC MATRIX
# OPTIMIZING TO MIN EMISSIONS
#-------------------------------------------------------------------------------
# Loading the necessary packages
#-------------------------------------------------------------------------------
using JuMP
using CPLEX
using CSV 
using DataFrames
using Printf
#----------------------------------------
# Metrics for emissions - obtained in get data step
#----------------------------------------
# HE = 30.865 #kg/ton of biomass harvested
# STE = 0.3014 #kg/km emissions small truck
# LTE = 0.8886 #kg/km emissions large truck 
# STC = 3.36 #tons capacity
# LTC = 21.76 #tons capacity 
# #-----------------------------------------

# Problem Description
# Deterministic Hub and Spoke Network for Power Plants
# ------------------------------
HS_P3 =  Model(with_optimizer(CPLEX.Optimizer));
#set_parameter(HS_P3, "CPXPARAM_Benders_Strategy" , 3)
#set_parameter(HS_P3,"CPXPARAM_MIP_Tolerances_MIPGap",0.0003)

# Variable Definition
# ------------------------------
# Now the variables must have an extra dimention depending on which variables 
# are dependent on the scenarios provided 
@variable(HS_P3, DEPOTS[1:D], Bin);
@variable(HS_P3, FLOWX[1:P, 1:D, 1:S] >= 0);
@variable(HS_P3, FLOWY[1:D, 1:C, 1:S] >= 0);
@variable(HS_P3, FLOWZ[1:P, 1:C, 1:S] >= 0);
@variable(HS_P3, SIGMAMin >= 0)
@variable(HS_P3, SIGMAMax >= 0)
#@variable(HS_P3, D_Model[1:D] == 0) #this is to change later
#this was kept in order to run SIM_AN

#------------------------------------------------------------------
#                     Constraints Definition
# -----------------------------------------------------------------

# constraint to restrict supply of biomass

@constraint(HS_P3, BIOMASS_SUPPLY[i=1:P, s=1:S], sum(FLOWX[i, j, s] for j=1:D) +
  sum(FLOWZ[i, j, s] for j=1:C) <= Supply[i,s]); #supply = yields

# constraint for mass balance

@constraint(HS_P3, MASS_BALANCE[j=1:D, s=1:S], sum((1-Moist[i])*FLOWX[i, j, s] for i=1:P) ==
  sum(FLOWY[j, i, s] for i=1:C)); #note that in the paper for this constraint i=k

# constraint to assure supply of the demand #REMOVED SHORTAGE

@constraint(HS_P3, DEMAND_SUPPLY[j=1:C, s=1:S], sum(FLOWY[i, j, s] for i=1:D) +
  sum((1-Moist[i])*FLOWZ[i, j ,s] for i=1:P) >= Demand[j]);

# constraint for biomass storage at depots

@constraint(HS_P3, BIOMASS_PROCESS[j=1:D, s=1:S], sum(FLOWX[i, j, s] for i=1:P) <=
  DepotCap[j]*DEPOTS[j]) ;

#constraint for depots so the SA CAN CHANGE THEM



#@constraint(HS_P3, fix_d[i = 1:D], DEPOTS[i] == D_Model[i]) ;

#----------------------------------------------------------------
#                   Objective Function
# ---------------------------------------------------------------
@objective(HS_P3, Min, 
	sum(prob[s]*(sum((EmissionsCost1[i, j]) * FLOWX[i, j, s] for i=1:P, j=1:D) +
	sum((EmissionsCost2[i, j]) * FLOWY[i, j, s] for i=1:D, j=1:C) +
	sum((EmissionsCost3[i,j]) * FLOWZ[i, j, s] for i=1:P, j=1:C)) for s=1:S));


optimize!(HS_P3)

obj = JuMP.objective_value(HS_P3)



