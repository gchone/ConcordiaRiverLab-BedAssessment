# -*- coding: utf-8 -*-

# Script to compute the bathymetry on a simple case

import pandas as pd


from BasicWSSmoothing import *
from BasicBedAssessment import *
from BasicRiverDataStructure import *



if __name__ == "__main__":


    # Simple synthetic data to test the algorithm
    # can be replaced by data = pd.read_csv(r'path\to\your\data.csv')
    dictdataset = {
        'dist': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230],
        'z_ws': [8.90, 8.92, 8.91, 8.95, 8.94, 8.92, 8.94, 8.96, 8.95, 8.99, 8.95, 9.95, 9.92, 9.92, 9.93, 9.95, 9.94, 9.94, 9.99, 9.94, 9.96, 9.98, 9.95, 9.99],
        'width': [10]*24,
        'Q': [1]*24,
    }
    data = Databrowser(pd.DataFrame(dictdataset))

    # Create a smooth and hydraulicaly correct water surface profile
    #   Add the attribute 'ztosmooth' to the data object, which is the water surface elevation to be smoothed (after the quantile carving process)
    #   Add the attribute 'z_smoothed' to the data object, which is the smoothed water surface elevation
    execute_WSsmoothing(data) # Water surface processing
    # Bathymetry assessment
    #   Add the attribute 'z' to the data object, which is the estimated bed elevation
    #   Other attributes, including the Froude number 'Fr', are also added to the data object
    execute_BedAssessment(data, 0.03, 0.00001) # Bathymetry assessment

    # Optional: Compute the water surface from the estimated bed elevation, using a conventionnal 1D hydraulic solver
    #   Add the attribute 'ws_validation' to the data object, which is the water surface elevation computed from the estimated bed elevation
    downstream_slope = data.get_first_point().s
    execute_SimpleHydro(data, 0.03, downstream_slope)

    # Return result as pandas dataframe and save it as a csv
    pandasout = data.topandasdf(["dist", "z_ws", "ztosmooth", "z_smoothed", "z", "Fr", "ws_validation"])
    pandasout.to_csv(r'D:\temp\bathy_test2.csv', index=False)
