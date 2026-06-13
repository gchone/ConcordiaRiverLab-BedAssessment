# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from BasicWSSmoothing import *
from BasicBedAssessment import *
from BasicRiverDataStructure import *
from scipy.optimize import minimize_scalar


if __name__ == "__main__":


    # Example of a regular 1D simulation (no bed assessment), with calibration

    # The input data represent charateristics at points measured along the river:
    # - dist: a longitudinal distance along the river (from downstream to upstream), situating the point
    # - z: bed elevation
    # - width: the river width at a measurement point (wetted width, extracted from LiDAR data)
    # - Q: the corresponding discharge (assessed discharge during the LiDAR acquisition)
    # NB: these four field names are hardcoded in the algorithm and cannot be changed in the provided csv file
    dictdataset = {
        'dist': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230],
        'z': [8.90, 8.92, 8.91, 8.95, 8.94, 8.92, 8.94, 8.96, 8.95, 8.99, 8.95, 9.95, 9.92, 9.92, 9.93, 9.95, 9.94, 9.94, 9.99, 9.94, 9.96, 9.98, 9.95, 9.99],
        'width': [10]*24,
        'Q': [1]*24,
    }
    df_data = pd.DataFrame(dictdataset)
    # data can be replaced by a csv file: df_data = pd.read_csv(r'path\to\your\data.csv')

    # Data for calibration
    SWOT_we = [9.5, 9.4, 9.6, 9.5, 9.4, 9.5, 9.5, 9.5, 9.6, 9.5, 10, 10.2, 10.2, 10.2, 10.3, 10.2, 10.4, 10.2, 10.5, 10.3, 10.3, 10.5, 10.4, 10.4]

    data = Databrowser(df_data)
    downstream_slope = 0.0001 # To be measured downstream of the reach from the LiDAR data
    execute_SimpleHydro(data, 0.03, downstream_slope)
    df_data = data.topandasdf(["dist", "z", "ws_validation"])


    def equations(n):  # the equation to solve, as a python function
        execute_SimpleHydro(data, n, downstream_slope)
        df_data = data.topandasdf(["dist", "z", "ws_validation"])
        rmse = np.sqrt(np.mean((df_data["ws_validation"].to_numpy() - np.array(SWOT_we)) ** 2))
        #print(n, rmse)
        if n < 0:
            return float('inf')
        return rmse

    res = minimize_scalar(equations, bracket=(0.01, 0.1), method='brent', tol=1e-4)
    print(res.x)
