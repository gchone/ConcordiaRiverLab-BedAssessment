# -*- coding: utf-8 -*-

#####################################################
# Guénolé Choné
# Concordia University
# Geography, Planning and Environment Department
# guenole.chone@concordia.ca
#####################################################

import pandas as pd
from BasicSolverDirect import *



def execute_BedAssessment(datapoints, manning, min_slope, method="OVERSAMPLING", max_delta_y=None):
    # method: "SIMPLE" (each cross-section solved individually with its immediate neighbour only),
    #         "OVERSAMPLING" (default: additional cross-sections are added where the Froude number varies too fast),
    #         "2-XS" (each cross-section is solved jointly with its immediate neighbour AND the next one downstream,
    #                 which helps convergence in difficult cases).
    # max_delta_y: optional maximum allowed increase of water depth per meter of distance (in %), used to bound the solver.

    # Order cross-sections from downstream (dist=0) to upstream (largest dist)
    ordered_points = list(datapoints.browse_down_to_up())

    # Compute upstream boundary slope
    prev_cs = None
    lastpoint = datapoints.get_last_point()
    for cs in ordered_points:
        if cs == lastpoint:
            localdist = (cs.dist - prev_cs.dist)
            cs.s = max(min_slope, (cs.z_smoothed-prev_cs.z_smoothed)/localdist)
        prev_cs = cs

    # Build, for every cross-section, the list of adjacent cross-sections used by the solver (listtosolve), ordered
    # from upstream to downstream: [upstream_neighbour, cs, downstream_neighbour]. The downstream_neighbour is only
    # used when method == "2-XS" (it is discarded otherwise by the solver).
    for i, cs in enumerate(ordered_points):
        cs.n = manning
        cs.listtosolve = [cs]
        cs.position_in_list = 0
        if i < len(ordered_points) - 1:
            upstream_cs = ordered_points[i + 1]
            cs.localdist_up = upstream_cs.dist - cs.dist
            cs.listtosolve.insert(0, upstream_cs)
            cs.position_in_list = 1
        if i > 0:
            downstream_cs = ordered_points[i - 1]
            cs.localdist_down = cs.dist - downstream_cs.dist
            cs.listtosolve.append(downstream_cs)

    # 1D hydraulic calculations (subcritical flow), from upstream to downstream
    prev_cs = None
    for cs in datapoints.browse_up_to_down():

        if prev_cs is None: # Compute upstream boundary level using Manning's equation only
            manning_inversesolver(cs)
            cs.solver = "manning up" # cs.solver and cs.type are flags for the cross-sections, to output
            cs.type = 0
        else: # For any other point, use the regular inverse hydraulic solver
            cs.solver = "regular"
            cs.type = 1
            __recursive_inverse1Dhydro(datapoints, cs, prev_cs, min_slope, method, max_delta_y)
        prev_cs = cs

    return

def __recursive_inverse1Dhydro(datapoints, cs, prev_cs, min_slope, method, max_delta_y):
    # This function applies the inverse hydraulic solver to compute the bed elevation at the current cross-section
    # (cs), knowing the condition at the reference cross-section (prev_cs), which is located upstream of cs during
    # the subcritical pass, and downstream of cs during the supercritical pass.
    # This is done recursively: if, after computing the flow at the cross-section, the Froude number appears to vary
    # too much, the computed bed elevation is discarded and an additional cross-section is added in-between
    # (method == "OVERSAMPLING" only).

    res = cs_solver(cs, min_slope, method, max_delta_y) # Solve the inverse 1D hydraulic problem

    localdist = cs.localdist_up

    # Adding a cross-section if the Froude number varies too much (increase by more than 50%)
    if method == "OVERSAMPLING" and prev_cs.Fr and abs((cs.Fr - prev_cs.Fr) / prev_cs.Fr) > 0.5 and localdist > 0.1: # Minimum 10cm between cs

        newdist = (cs.dist + prev_cs.dist) / 2.
        newcs = datapoints.add_point(newdist) # Adding a point in the dataset at the right distance

        # cs and prev_cs are, respectively, the downstream and upstream neighbour, or the reverse, depending on the
        # direction the solver is currently working in
        upstream_pt, downstream_pt = (prev_cs, cs)

        # Linear interpolation of width, discharge and smoothed water surface for the new point.
        # Although more accurate spatialization could be done, this is deemed accurate enough
        t = (newdist - downstream_pt.dist) / (upstream_pt.dist - downstream_pt.dist)
        newcs.width = downstream_pt.width + t * (upstream_pt.width - downstream_pt.width)
        newcs.Q = downstream_pt.Q + t * (upstream_pt.Q - downstream_pt.Q)
        newcs.z_smoothed = downstream_pt.z_smoothed + t * (upstream_pt.z_smoothed - downstream_pt.z_smoothed)
        newcs.n = cs.n
        newcs.solver = "regular"

        # Update local distances and the solving neighbourhood (listtosolve) around the newly added cross-section
        newlocaldist = localdist / 2.
        newcs.localdist_up = newlocaldist
        newcs.localdist_down = newlocaldist
        upstream_pt.localdist_down = newlocaldist
        downstream_pt.localdist_up = newlocaldist

        # Only a pair is used here (oversampling never combines with the "2-XS" method), so listtosolve simply holds
        # the reference point and the point being solved, in that order
        newcs.listtosolve = [upstream_pt, newcs]
        newcs.position_in_list = 1

        newcs.type = 3
        # cs is now adjacent to the newly inserted point instead of prev_cs: refresh its solving neighbourhood
        cs.listtosolve = [newcs, cs]
        cs.position_in_list = 1
        __recursive_inverse1Dhydro(datapoints, newcs, prev_cs, min_slope, method, max_delta_y) # Compute the bed elevation at the new added cross-section
        __recursive_inverse1Dhydro(datapoints, cs, newcs, min_slope, method, max_delta_y) # Compute the bed elevation at the downstream cross-section cs

    return res



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




