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


def cs_solver(cs, min_slope, method, max_delta_y):
    # This function is an inverse 1D hydraulic solver, using Manning's and Bernoulli's equations, that jointly solves
    # for the flow depth at up to 3 adjacent cross-sections stored in cs.listtosolve (ordered upstream to downstream),
    # with cs located at cs.position_in_list.
    # method: "SIMPLE"/"OVERSAMPLING" only use the pair [neighbour, cs] ; "2-XS" additionally uses the cross-section
    #   further downstream, jointly solving both intervals for a more robust convergence.
    # max_delta_y: maximum allowed increase of water depth per meter of distance (used to bound the solver).

    if method != "2-XS" and len(cs.listtosolve) == 3:
        # Only the immediate neighbour is used, the extra (2-XS) cross-section is discarded
        cs.listtosolve.pop(2)

    def equations(y): # the equation(s) to solve, as a python function
        # For each pair of adjacent cross-sections in cs.listtosolve, the difference between the resultant energy
        # (potential energy, i.e. water surface elevation, plus kinetic energy, plus energy loss by friction) and the
        # energy computed at the reference cross-section is computed. This function is used by minimize, which tries
        # to find the depth(s) y so that the combined misfit is minimal.

        dif_energy = []
        for i in range(len(cs.listtosolve) - 1):

            cs_tosolve = cs.listtosolve[i + 1]
            cs_ref = cs.listtosolve[i]
            localdist = cs_tosolve.localdist_up

            if i == 0:
                cs_ref.temp_h = cs_ref.h
                cs_ref.temp_s = cs_ref.s

            if abs(cs_ref.z_smoothed - cs_tosolve.z_smoothed) / localdist <= min_slope:
                cs_tosolve.solver = "min_slope"
                h_ref = cs_ref.temp_h + localdist * (min_slope - (cs_ref.z_smoothed - cs_tosolve.z_smoothed) / localdist)
            else:
                h_ref = cs_ref.temp_h

            v = cs_tosolve.Q / (cs_tosolve.width * y[i])
            R = (cs_tosolve.width * y[i]) / (cs_tosolve.width + 2 * y[i])
            s = (cs_tosolve.n ** 2 * v ** 2) / (R ** (4. / 3.))
            h = cs_tosolve.z_smoothed
            h = h + v ** 2 / (2 * g) # add kinetic energy
            cs_tosolve.temp_h = h
            cs_tosolve.temp_s = s

            # slope calculation:
            friction_h = localdist * (s + cs_ref.temp_s) / 2.
            if len(cs.listtosolve) - 1 == 1:
                # Friction is based on the downstream computed slope only, when only one interval is solved
                # (necessary for convergence)
                friction_h = localdist * s

            dif_energy.append(abs(friction_h + h - h_ref))

        misfit = (sum([e ** 2 for e in dif_energy])) ** 0.5
        return misfit


    # Compute the critical depth at each downstream cross-section in the list
    for i in range(len(cs.listtosolve) - 1):
        cs_tmp = cs.listtosolve[i + 1]
        cs_tmp.ycrit = (cs_tmp.Q / (cs_tmp.width * g ** 0.5)) ** (2. / 3.)

    ycrit = [cs.listtosolve[i + 1].ycrit for i in range(len(cs.listtosolve) - 1)] # initial guess

    # Set up maximum depth limits
    if max_delta_y is not None:
        max_y = [max(cs.listtosolve[1].ycrit,
                     min(cs.listtosolve[0].y * (1 + max_delta_y * cs.localdist_up / 100.), cs.listtosolve[1].width))]
    else:
        max_y = [max(cs.listtosolve[1].ycrit, cs.listtosolve[1].width)]

    max_y.extend([cs.listtosolve[i + 1].width for i in range(1, len(cs.listtosolve) - 1)])
    bounds = [(cs.listtosolve[i + 1].ycrit, max_y[i]) for i in range(len(cs.listtosolve) - 1)]
    res = minimize(equations, ycrit, method='Nelder-Mead', bounds=bounds,
                   options={'xatol': 1e-3, 'fatol': 1e-6})
    cs.solver = "regular"

    if max_delta_y is not None and res.x[cs.position_in_list - 1] == cs.listtosolve[0].y * max_delta_y * cs.localdist_up / 100.:
        cs.solver = "max depth gradient"
    if res.x[cs.position_in_list - 1] == cs.listtosolve[1].width:
        cs.solver = "max depth"

    # If the supercritical flow has a better fit than the subcritical one, the critical depth is retained as the final
    # answer
    bounds = [(0, cs.listtosolve[i + 1].ycrit) for i in range(len(cs.listtosolve) - 1)]
    res_super = minimize(equations, ycrit, method='Nelder-Mead', bounds=bounds,
                         options={'xatol': 1e-3, 'fatol': 1e-6})
    res = min([res, res_super], key=lambda r: r.fun)

    cs.y = res.x[cs.position_in_list - 1]
    if cs.y < cs.ycrit:
        cs.y = cs.ycrit
        cs.solver = "critical"

    cs.R = (cs.width * cs.y) / (cs.width + 2 * cs.y)
    cs.v = cs.Q / (cs.width * cs.y)
    cs.z = cs.z_smoothed - cs.y
    cs.s = (cs.n ** 2 * cs.v ** 2) / (cs.R ** (4. / 3.))
    cs.h = cs.z_smoothed
    cs.h = cs.h + cs.v ** 2 / (2 * g)  # add kinetic energy
    cs.Fr = cs.v / (g * cs.y) ** 0.5

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
        friction_h = localdist * (s+cs_down.s_validation)/2.
        #friction_h = localdist * s # Friction based on the upstream computed slope
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


