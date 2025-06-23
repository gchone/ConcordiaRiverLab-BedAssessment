# -*- coding: utf-8 -*-

#####################################################
# Guénolé Choné
# Concordia University
# Geography, Planning and Environment Department
# guenole.chone@concordia.ca
#####################################################

import pandas as pd
from BasicSolverDirect import *



def execute_BedAssessment(datapoints, manning, min_slope):

    # Compute upstream boundary slope
    prev_cs = None
    lastpoint = datapoints.get_last_point()
    for cs in datapoints.browse_down_to_up():
        if cs == lastpoint:
            localdist = (cs.dist - prev_cs.dist)
            cs.s = max(min_slope, (cs.z_smoothed-prev_cs.z_smoothed)/localdist)
        prev_cs = cs

    # 1D hydraulic calculations
    prev_cs = None
    for cs in datapoints.browse_up_to_down():

        cs.n = manning
        if prev_cs is None: # Compute upstream boundary level using Manning's equation only
            manning_inversesolver(cs)
            cs.solver = "manning up" # cs.solver and cs.type are flags for the cross-sections, to output
            cs.type = 0
        else: # For any other point, use the regular inverse hydraulic solver
            cs.solver = "regular"
            cs.type = 1
            __recursive_inverse1Dhydro(datapoints, cs, prev_cs, min_slope)
        prev_cs = cs

    return

def __recursive_inverse1Dhydro(datapoints, cs, prev_cs, min_slope):
    # This function apply the inverse hydraulic solver to computer bed elevation at the current cross-section (cs),
    # knowing the condition at the upstream cross-section (prev_cs)
    # This is done recursively: if, after computing the flow at the cross-section, the Froude number appears to vary too
    # much, the computed bed elevation is discarded and an additional cross-section is added in-between.

    cs_inversesolver(prev_cs, cs, min_slope) # Solve the inverse 1D hydraulic problem

    localdist = (prev_cs.dist - cs.dist)

    # Adding a cross-section if the Froude number varies too much (increase by more than 50%)
    if (cs.Fr - prev_cs.Fr) / prev_cs.Fr > 0.5 and localdist > 0.1: # Minimum 10cm between cs

        newcs = datapoints.add_point((cs.dist + prev_cs.dist) / 2.) # Adding a point in the dataset at the right distance
        newlocaldist = localdist / 2.
        # Linear interpolation of width, discharge and water surface for the new point.
        # Although more accurate spatialization could be done, this is deemed accurate enough
        a = (cs.width - prev_cs.width) / (0-localdist)
        newcs.width = a * newlocaldist + cs.width
        a = (cs.Q - prev_cs.Q) / (0-localdist)
        newcs.Q = a * newlocaldist + cs.Q
        a = (cs.z_smoothed - prev_cs.z_smoothed) / (0-localdist)
        newcs.z_smoothed = a* newlocaldist + cs.z_smoothed
        newcs.n = cs.n
        newcs.solver = "regular"
        __recursive_inverse1Dhydro(datapoints, newcs, prev_cs, min_slope) # Compute the bed elevation at the new added cross-section
        newcs.type = 3
        __recursive_inverse1Dhydro(datapoints, cs, newcs, min_slope) # Compute the bed elevation at the downstream cross-section cs

    return



def execute_SimpleHydro(datapoints, manning, down_slope):


    # Set downstream boundary slope
    firstpoint = datapoints.get_first_point()
    firstpoint.s_validation = down_slope

    # 1D hydraulic calculations
    prev_cs = None
    for cs in datapoints.browse_down_to_up():

        cs.n = manning
        if prev_cs is None: # Compute downstream boundary level using Manning's equation only
            manning_normalsolver(cs)
            cs.solver_validation = "manning down" # cs.solver and cs.type are flags for the cross-sections, to output
            cs.type_validation = 0
        else: # For any other point, use the regular hydraulic solver
            cs.solver_validation = "regular"
            cs.type_validation = 1
            cs_normalsolver(cs, prev_cs)  # Solve the 1D hydraulic problem

        prev_cs = cs

    return




