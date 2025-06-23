# -*- coding: utf-8 -*-

import math
import numpy as np
from BasicQuantileRegression import QuantileCarving
from scipy.stats import norm
from scipy.optimize import minimize_scalar
import warnings


def execute_WSsmoothing(datapoints, quantile=0.2, smooth_level=600 , uncertainty_sigma = 300, uncertainty_factor=0.85, slope_sigma=300, slope_factor=2.0):

    # The smoothing process :
    # - Removes bumps in the water surface profile following the quantile carving process of
    #   Schwanghart and Scherler (2017). See QuantileRegression.py for details. Note that the smoothing technique from
    #   Schwanghart and Scherler (2017) is not used but replaced my a homemade one.
    # - Smooths the river profile according to the estimated uncertainty. Uncertainty is assessed by the difference
    #   between the maximum elevation downstream and the minimum elevation upstream any point in the profile.

    # Quantile carving
    QuantileCarving(datapoints, quantile)

    # Smoothing
    values = []
    unbreached_values = []
    distances = []
    for cs in datapoints.browse_down_to_up():
        distances.append(cs.dist)
        values.append(cs.ztosmooth)
        unbreached_values.append(cs.z_ws)
    distances = np.array(distances)
    values = np.array(values)
    unbreached_values = np.array(unbreached_values)

    carving = unbreached_values - values
    smoothed_values = np.zeros_like(values)
    uncertainty_vec = np.zeros_like(values)
    restricted = np.zeros_like(values)
    sd2_vec = np.zeros_like(values)
    local_sigma_vec = np.zeros_like(values)
    for i in range(len(values)):
        # Gaussian curve size (sigma) is limited on the edges to avoid mismatch with downstream reaches
        local_sigma = min(smooth_level, (distances[i] - distances[0])*5.)  # hardcoded: 5 times the distance to the first point
        local_sigma = max(local_sigma, 10.)  # hardcoded: minimum standard deviation
        local_sigma_vec[i] = local_sigma
        # Compute Gaussian weights using norm.pdf
        weights = norm.pdf(distances, loc=distances[i], scale=uncertainty_sigma)
        weights /= weights.sum()  # Normalize weights
        # Uncertainty is calculated from:
        # - the absolute value of the carving (how much carving is done)
        # - the difference between the elevation and surrounding elevations (how much slope there is).
        # A exponential transformation is applyied to that 0 difference of elevation = 1. The slopfactor is added to put
        # more or less weight on the slope.
        # - The ratio between the carving and the differences between the elevations gives a measure of the uncertainty
        # relative to the slope
        # - Everything is multiplied by the weights from the Gaussian curve and summed to get the final uncertainty
        corrections = sum(np.abs(carving) * weights)**uncertainty_factor
        if corrections < 1e-9:
            # If there is no carving, there are no smoothing to be made
            smoothed_values[i] = values[i]
        else:
            weightsslope = norm.pdf(distances, loc=distances[i], scale=slope_sigma)
            weightsslope /= weightsslope.sum()  # Normalize weights
            deltaz = math.exp(sum(np.abs(values[i] - values) * weightsslope)) ** slope_factor
            uncertainty = corrections / deltaz
            uncertainty_vec[i] = uncertainty
            if i > 0:
                # In order to avoid the problem of the Gaussian curve being too wide, we need to make sure that the pdf of
                # the gaussian curve is lower than the previous one (on the left side of the previous one). Otherwise the
                # resulting elevation can be lower than the previous one, creating a non-hydraulically valid profile.
                x_values = distances[0:i - 1]
                mu1 = distances[i - 1]
                sd1 = sd2  # sd1 is the previous sd2
                sd2 = uncertainty * local_sigma
                mu2 = distances[i]
                F1 = norm.pdf(x_values, loc=mu1, scale=sd1)
                F2 = norm.pdf(x_values, loc=mu2, scale=sd2)
                validpdf = np.all(F1 >= F2)
                if not validpdf:
                    # The Gaussian curve is too wide, compared to the previous ones
                    # We need to reduce the standard deviation of the Gaussian curve in that case.
                    # We will use the optimization to find the maximum possible standard deviation
                    restricted[i] = 1

                    def objective(tested_sd2):
                        """Objective function to minimize: we want the negative of sd2 for maximization"""
                        F1 = norm.pdf(x_values, loc=mu1, scale=sd1)
                        F2 = norm.pdf(x_values, loc=mu2, scale=tested_sd2)
                        diff = F1 - F2
                        if np.any(diff < 0):
                            return np.inf  # violates the constraint
                        return -tested_sd2  # maximize sd2

                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=RuntimeWarning)
                        result = minimize_scalar(objective, bounds=(0.001, sd2), method='bounded')
                        if result.success:
                            sd2 = result.x
                        else:
                            # If optimization fails, the best guess is to use the previous sd
                            sd2 = sd1
                            restricted[i] = 2
            else:
                # First point, no previous point to compare
                sd2 = uncertainty * local_sigma
            sd2_vec[i] = sd2
            weights = norm.pdf(distances, loc=distances[i], scale=sd2)
            weights /= weights.sum()  # Normalize weights
            # Compute the weighted average
            smoothed_values[i] = np.sum(weights * values)

            # Final check: if the smoothed value is lower than the previous one, we set it to the previous one
            # It seems to happen sometimes although it should not. I could not find the reason why. Probably a numerical approximation in the optimization.
            if i > 0 and smoothed_values[i] < smoothed_values[i - 1]:
                smoothed_values[i] = smoothed_values[i - 1]
                restricted[i] = restricted[i]+10

    # Assign the smoothed values to the cross-sections
    i = 0
    for cs in datapoints.browse_down_to_up():
        cs.z_smoothed = smoothed_values[i]
        #cs.ws_uncertainty = uncertainty_vec[i]
        #cs.restricted = restricted[i]
        #cs.sd2 = sd2_vec[i]
        #cs.local_sigma = local_sigma_vec[i]
        i +=1
    return

