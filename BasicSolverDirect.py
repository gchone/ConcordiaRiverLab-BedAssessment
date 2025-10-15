# -*- coding: utf-8 -*-

# Solver sous-critique uniquement

g = 9.81
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

def manning_inversesolver(cs):
    # This function solves Manning's equation
    # Inverse problem version (i.e. given ws, find z)

    def equations(y): # the equation to solve, as a python function
        # For a given flow depth y, the difference between the resultant discharge and the known discharge is computed
        # This function is used by fsolve, that tries to find y so that difQ = 0
        R = (cs.width * y) / (cs.width + 2 * y)
        difQ = (y * cs.width * R ** (2. / 3.) * cs.s ** 0.5) / cs.n - cs.Q
        return difQ

    cs.y = fsolve(equations, 1)[0] # Solve the manning equation (find y so that difQ = 0)
    cs.R = (cs.width * cs.y) / (cs.width + 2 * cs.y)
    cs.ycrit = (cs.Q / (cs.width * g ** 0.5)) ** (2. / 3.)

    cs.v = cs.Q / (cs.width * cs.y)
    cs.z = cs.z_smoothed - cs.y
    cs.h = cs.z_smoothed

    cs.h = cs.h + cs.v ** 2 / (2 * g) # add kinetic energy
    cs.Fr = cs.v / (g * cs.y) ** 0.5


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

    # the solver starts at y = y_crit
    cs_down.ycrit = (cs_down.Q / (cs_down.width * g ** 0.5)) ** (2. / 3.)

    def equations(y): # the equation to solve, as a python function
        # For a given flow depth y, the difference between the resultant energy (potential energy, i.e. water surface
        # elevation, plus kinetic energy, plus energy loss by friction) and the energy computed upstream is computed.
        # This function is used by fsolve, that tries to find y so that dif_energy = 0

        if y < cs_down.ycrit:
            # constraint computation, so that the flow is never supercritical
            return float('inf')
        R = (cs_down.width * y) / (cs_down.width + 2 * y)
        v = cs_down.Q / (cs_down.width * y)
        s = (cs_down.n ** 2 * v ** 2) / (R ** (4. / 3.))
        h = cs_down.z_smoothed
        h = h + v ** 2 / (2 * g) # add kinetic energy
        # slope calculation:
        #friction_h = localdist * (s+cs_up.s)/2. # Friction can't be based on the average of slope, it leads to impossible to resolve cases
        friction_h = localdist * s # Replaced by a friction based and the downstream computed slope
        dif_energy = friction_h + h - h_ref
        dif_energy = abs(dif_energy)
        return dif_energy


    #res, dict, ier, msg = fsolve(equations, cs_down.ycrit, full_output=True)
    #res = minimize(equations, cs_down.ycrit, method='Nelder-Mead', options={'xatol': 1e-3})
    res = minimize_scalar(equations, method='brent', tol=1e-3)

    #cs_down.y = res.x[0]
    cs_down.y = res.x
    cs_down.R = (cs_down.width * cs_down.y) / (cs_down.width + 2 * cs_down.y)
    cs_down.v = cs_down.Q / (cs_down.width * cs_down.y)
    cs_down.z = cs_down.z_smoothed - cs_down.y
    cs_down.s = (cs_down.n ** 2 * cs_down.v ** 2) / (cs_down.R ** (4. / 3.))

    cs_down.h = cs_down.z_smoothed

    cs_down.h = cs_down.h + cs_down.v ** 2 / (2 * g) # add kinetic energy

    cs_down.Fr = cs_down.v / (g * cs_down.y) ** 0.5

    return res



def manning_normalsolver(cs):
    # This function solves Manning's equation
    # Normal problem version (i.e. given z, find ws)

    def equations(y): # the equation to solve, as a python function
        # For a given flow depth y, the difference between the resultant discharge and the known discharge is computed
        # This function is used by fsolve, that tries to find y so that difQ = 0
        R = (cs.width * y) / (cs.width + 2 * y)
        difQ = (y * cs.width * R ** (2. / 3.) * cs.s_validation ** 0.5) / cs.n - cs.Q
        return difQ

    cs.y_validation = fsolve(equations, 1)[0] # Solve the manning equation (find y so that difQ = 0)
    cs.R_validation = (cs.width * cs.y_validation) / (cs.width + 2 * cs.y_validation)
    cs.ycrit_validation = (cs.Q / (cs.width * g ** 0.5)) ** (2. / 3.)

    cs.v_validation = cs.Q / (cs.width * cs.y_validation)
    cs.ws_validation = cs.z + cs.y_validation
    cs.h_validation = cs.ws_validation

    cs.h_validation = cs.h_validation + cs.v_validation ** 2 / (2 * g) # add kinetic energy
    cs.Fr_validation = cs.v_validation / (g * cs.y_validation) ** 0.5


def cs_normalsolver(cs_up, cs_down):
    # This function is a 1D hydraulic solver, using Manning's and Bernoulli's equations to computed flow at a
    # upstream cross-section, knowing the conditions downstream
    # Normal problem version (i.e. given z, find ws)

    localdist = (cs_up.dist - cs_down.dist)

    h_ref = cs_down.h_validation

    # the solver starts at y = y_crit
    cs_up.ycrit_validation = (cs_up.Q / (cs_up.width * g ** 0.5)) ** (2. / 3.)

    def equations(y): # the equation to solve, as a python function
        # For a given flow depth y, the difference between the resultant energy (potential energy, i.e. water surface
        # elevation, plus kinetic energy, plus energy loss by friction) and the energy computed upstream is computed.
        # This function is used by fsolve, that tries to find y so that dif_energy = 0

        if y < cs_up.ycrit_validation:
            # constraint computation, so that the flow is never supercritical
            return float("inf")
        R = (cs_up.width * y) / (cs_up.width + 2 * y)
        v = cs_up.Q / (cs_up.width * y)
        s = (cs_up.n ** 2 * v ** 2) / (R ** (4. / 3.))
        h = cs_up.z + y
        h = h + v ** 2 / (2 * g) # add kinetic energy
        # slope calculation:
        #friction_h = localdist * (s+cs_up.s)/2. # Friction can't be based on the average of slope, it leads to impossible to resolve cases
        friction_h = localdist * s # Replaced by a friction based and the downstream computed slope
        dif_energy = friction_h + h_ref - h
        dif_energy = abs(dif_energy)
        return dif_energy

    #res = minimize(equations, cs_down.ycrit, method='Nelder-Mead', options={'xatol': 1e-3})
    res = minimize_scalar(equations, method='brent', tol=1e-3)

    #cs_up.y_validation = res.x[0]
    cs_up.y_validation = res.x
    cs_up.R_validation = (cs_up.width * cs_up.y_validation) / (cs_up.width + 2 * cs_up.y_validation)
    cs_up.v_validation = cs_up.Q / (cs_up.width * cs_up.y_validation)
    cs_up.ws_validation = cs_up.z + cs_up.y_validation
    cs_up.s_validation = (cs_up.n ** 2 * cs_up.v_validation ** 2) / (cs_up.R_validation ** (4. / 3.))

    cs_up.h_validation = cs_up.ws_validation

    cs_up.h_validation = cs_up.h_validation + cs_up.v_validation ** 2 / (2 * g) # add kinetic energy

    cs_up.Fr_validation = cs_up.v_validation / (g * cs_up.y_validation) ** 0.5

    return res


