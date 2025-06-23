# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse
import scipy.optimize
import math

def QuantileCarving(datapoints, tau=0.5):
    # This quantile carving process comes from :
    #    Schwanghart, W., Scherler, D., 2017. Bumps in river profiles:
    #    uncertainty assessment and smoothing using quantile regression
    #    techniques. Earth Surface Dynamics, 5, 821-839.
    #    [DOI: 10.5194/esurf-5-821-2017]
    # Code was adapted from their Matlab code:
    # https://github.com/wschwanghart/topotoolbox/@FLOWobj/quantcarve.m

    x = []
    z = []

    for cs in datapoints.browse_down_to_up():
        x.append(cs.dist)
        z.append(cs.z_ws)

    x = np.array(x)
    z = np.array(z)

    n = len(datapoints)

    ix = range(1, n)
    ixc = range(0, n-1)

    f = [tau*np.ones((n, 1)),(1-tau)*np.ones((n, 1)),np.zeros((n, 1))]
    f = np.vstack(f)

    Aeq = scipy.sparse.hstack([scipy.sparse.identity(n), -scipy.sparse.identity(n), scipy.sparse.identity(n)])
    beq = z

    lb = [0.]*(2*n)
    lb.extend([-math.inf]*n)
    bounds = [(lower, None) for lower in lb]

    d = 1./(x[ix]-x[ixc])

    Atmp = scipy.sparse.csr_matrix(np.zeros((n, n * 2)))
    Atmp2 = scipy.sparse.coo_matrix((d, (ix, ixc)), shape=(n, n)) - scipy.sparse.coo_matrix((d, (ix, ix)), shape=(n, n))
    A = scipy.sparse.hstack([Atmp, Atmp2])

    b = np.zeros((n,1))

    output = scipy.optimize.linprog(f, A, b, Aeq, beq, bounds=bounds,
                           method='highs', callback=None)
    newz = output.x[-n:]

    i = 0
    for cs in datapoints.browse_down_to_up():
        cs.ztosmooth = newz[i]
        i+=1






