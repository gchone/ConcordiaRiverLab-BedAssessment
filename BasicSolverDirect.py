# -*- coding: utf-8 -*-

g = 9.81
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

from scipy.optimize import fsolve
from scipy.optimize import minimize

def manning_inversesolver(cs):
    # This function solves Manning's equation
    # Inverse problem version (i.e. given ws, find z)

    def equations(param): # the equation to solve, as a python function
        # For a given parameter defining a cross-section (e.g. the bed elevation for a rectangular cross-section),
        # the difference between the resultant discharge and the known discharge is computed
        # This function is used by fsolve, that tries to find y so that difQ = 0
        cs.define(param)
        R = cs.area() / cs.perimeter()
        difQ = (cs.area() * R ** (2. / 3.) * cs.s ** 0.5) / cs.n - cs.Q
        return difQ

    init_param = cs.get_init_param() # Initial guess for the parameter
    computed_param = fsolve(equations, init_param)[0] # Solve the manning equation (find y so that difQ = 0)

    cs.define(computed_param)
    cs.R = cs.area() / cs.perimeter()
    cs.v = cs.Q / (cs.area())
    cs.h = cs.z_smoothed
    cs.h = cs.h + cs.v ** 2 / (2 * g) # add kinetic energy
    cs.Fr = cs.v / (g * cs.area() / cs.width) ** 0.5


def cs_inversesolver(cs_up, cs_down, min_slope):
    # This function is an inverse 1D hydraulic solver, using Manning's and Bernoulli's equations to computed flow at a
    # downstream cross-section, knowing the conditions upstream
    # Inverse problem version (i.e. given ws, find z)

    localdist = (cs_up.dist - cs_down.dist)

    if (cs_up.z_smoothed - cs_down.z_smoothed)/localdist <= min_slope:
        cs_down.solver = "min_slope"
        h_ref = cs_up.h + localdist * (min_slope - (cs_up.z_smoothed - cs_down.z_smoothed) / localdist)
    else:
        h_ref = cs_up.h

    def equations(param): # the equation to solve, as a python function
        # For a given parameter defining the cross-section geometry (e.g. the bed elevation for a rectangular
        # cross-section), the difference between the resultant energy (potential energy, i.e. water surface
        # elevation, plus kinetic energy, plus energy loss by friction) and the energy computed upstream is computed.

        cs_down.define(param)
        bracket = cs_down.valide_bracket()
        if not (min(bracket) < param < max(bracket)):
            # constraint computation: the parameter is invalid
            return float('inf')
        v = cs_down.Q / (cs_down.area())
        Fr = v / (g * cs_down.area() / cs_down.width) ** 0.5
        supercritical_penalty = 0
        if Fr > 1:
            return float('inf')
        R = cs_down.area() / cs_down.perimeter()
        s = (cs_down.n ** 2 * v ** 2) / (R ** (4. / 3.))
        h = cs_down.z_smoothed
        h = h + v ** 2 / (2 * g) # add kinetic energy
        # slope calculation:
        #friction_h = localdist * (s+cs_up.s)/2. # Friction can't be based on the average of slope, it leads to impossible to resolve cases
        friction_h = localdist * s # Replaced by a friction based and the downstream computed slope
        dif_energy = friction_h + h - h_ref
        dif_energy = abs(dif_energy) + supercritical_penalty
        return dif_energy

    start_param = cs_down.get_init_param()
    # minimize_scalar should theoretically be used here, but it is not working with the constraints. So minimize is used instead.
    res = minimize(equations, start_param, method='Nelder-Mead', options={'xatol': 1e-3})

    cs_down.define(res.x[0])
    cs_down.R = cs_down.area() / cs_down.perimeter()
    cs_down.v = cs_down.Q / (cs_down.area())
    cs_down.s = (cs_down.n ** 2 * cs_down.v ** 2) / (cs_down.R ** (4. / 3.))

    cs_down.h = cs_down.z_smoothed
    cs_down.h = cs_down.h + cs_down.v ** 2 / (2 * g) # add kinetic energy
    cs_down.Fr = cs_down.v / (g * cs_down.area() / cs_down.width) ** 0.5

    return res



def manning_normalsolver(cs):
    # This function solves Manning's equation
    # Normal problem version (i.e. given z, find ws)

    def equations(zws): # the equation to solve, as a python function
        # For a given water elevation zws, the difference between the resultant discharge and the known discharge is computed
        # This function is used by fsolve, that tries to find zws so that difQ = 0
        R = cs.area(zws) / cs.perimeter(zws)
        difQ = (cs.area(zws) * R ** (2. / 3.) * cs.s_validation ** 0.5) / cs.n - cs.Q
        return difQ

    cs.ws_validation = fsolve(equations, cs.z_smoothed)[0] # Solve the manning equation (find y so that difQ = 0)
    cs.R_validation = cs.area(cs.ws_validation) / cs.perimeter(cs.ws_validation)

    cs.v_validation = cs.Q / cs.area(cs.ws_validation)
    cs.h_validation = cs.ws_validation
    cs.h_validation = cs.h_validation + cs.v_validation ** 2 / (2 * g) # add kinetic energy
    cs.Fr_validation = cs.v_validation / (g * cs.area(cs.ws_validation) / cs.wetted_width(cs.ws_validation)) ** 0.5



def cs_normalsolver(cs_up, cs_down):
    # This function is a 1D hydraulic solver, using Manning's and Bernoulli's equations to computed flow at a
    # upstream cross-section, knowing the conditions downstream
    # Normal problem version (i.e. given z, find ws)

    localdist = (cs_up.dist - cs_down.dist)

    h_ref = cs_down.h_validation

    def equations(ws): # the equation to solve, as a python function
        # For a given water surface eleation ws, the difference between the resultant energy (potential energy, i.e. water surface
        # elevation, plus kinetic energy, plus energy loss by friction) and the energy computed upstream is computed.

        min_z = cs_up.get_thalweg()[1]
        if not min_z < ws:
            # constraint computation: the parameter is invalid
            return float('inf')
        v = cs_up.Q / (cs_up.area(ws))
        Fr = v / (g * cs_up.area(ws) / cs_up.wetted_width(ws)) ** 0.5
        if Fr > 1:
            return float('inf')
        R = cs_up.area(ws) / cs_up.perimeter(ws)
        s = (cs_up.n ** 2 * v ** 2) / (R ** (4. / 3.))
        h = ws + v ** 2 / (2 * g)  # add kinetic energy
        # slope calculation:
        friction_h = localdist * (s+cs_down.s_validation)/2.
        dif_energy = friction_h + h_ref - h
        dif_energy = abs(dif_energy)
        return dif_energy

    start_z = cs_up.get_thalweg()[1]+1 # Initial guess = 1m of water
    res = minimize(equations, start_z, method='Nelder-Mead', options={'xatol': 1e-3})

    cs_up.ws_validation = res.x[0]
    cs_up.R_validation = cs_up.area(cs_up.ws_validation) / cs_up.perimeter(cs_up.ws_validation)
    cs_up.v_validation = cs_up.Q / (cs_up.area(cs_up.ws_validation))
    cs_up.s_validation = (cs_up.n ** 2 * cs_up.v_validation ** 2) / (cs_up.R_validation ** (4. / 3.))

    cs_up.h_validation = cs_up.ws_validation
    cs_up.h_validation = cs_up.h_validation + cs_up.v_validation ** 2 / (2 * g) # add kinetic energy

    cs_up.Fr_validation = cs_up.v_validation / (g * cs_up.area(cs_up.ws_validation) / cs_up.wetted_width(cs_up.ws_validation)) ** 0.5

    return res


