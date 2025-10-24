# -*- coding: utf-8 -*-

# Script to compute the bathymetry on a simple case

import pandas as pd
import matplotlib.pyplot as plt

from BasicWSSmoothing import *
from BasicBedAssessment import *
from BasicRiverDataStructure import *



if __name__ == "__main__":


    # Simple synthetic data to test the algorithm
    # can be replaced by data provided in a csv file with: df_data = pd.read_csv(r'path\to\your\data.csv')
    # The input data represent charateristics at points measured along the river:
    # - dist: a longitudinal distance along the river (from downstream to upstream), situating the point
    # - z_ws: the water surface elevation, usually extracted from LiDAR data, at a measurement point
    # - width: the river width at a measurement point (wetted width, extracted from LiDAR data)
    # - Q: the corresponding discharge (assessed discharge during the LiDAR acquisition)
    # NB: these four field names are hardcoded in the algorithm and cannot be changed in the provided csv file
    dictdataset = {
        'dist': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230],
        'z_ws': [8.90, 8.92, 8.91, 8.95, 8.94, 8.92, 8.94, 8.96, 8.95, 8.99, 8.95, 9.95, 9.92, 9.92, 9.93, 9.95, 9.94, 9.94, 9.99, 9.94, 9.96, 9.98, 9.95, 9.99],
        'width': [10]*24,
        'Q': [1]*24,
    }
    df_data = pd.DataFrame(dictdataset)
    
    data = Databrowser(df_data)


    # Create a smooth and hydraulicaly correct water surface profile
    #   Add the attribute 'ztosmooth' to the data object, which is the water surface elevation to be smoothed (after the quantile carving process)
    #   Add the attribute 'z_smoothed' to the data object, which is the smoothed water surface elevation
    execute_WSsmoothing(data) # Water surface processing
    # Bathymetry assessment
    #   Add the attribute 'z' to the data object, which is the estimated bed elevation
    #   Other attributes, including the Froude number 'Fr', are also added to the data object
    execute_BedAssessment(data, 0.03, 0.00001) # Bathymetry assessment
    df_beddata = data.topandasdf(
        ["dist", "z_ws", "ztosmooth", "z_smoothed", "z", "Fr"])  # Return result as pandas dataframe

    # Optional: Filter the estimated bed elevation using the Ramer-Douglas-Peucker algorithm
    data_reduced = data.reduce_bedpoints_RDP(0.1) # Epsilon is the tolerance, in the unit of the elevation (m)

    # Optional: Compute the water surface from the estimated bed elevation, using a conventionnal 1D hydraulic solver)
    #   Add the attribute 'ws_validation' to the data object, which is the water surface elevation computed from the estimated bed elevation
    downstream_slope = data_reduced.get_first_point().s
    execute_SimpleHydro(data_reduced, 0.03, downstream_slope)
    df_beddata_reduced = data_reduced.topandasdf(["dist", "z_ws", "ztosmooth", "z_smoothed", "z", "Fr", "ws_validation"])

    # Save it as a csv
    df_beddata_reduced.to_csv(r'D:\NRCAN2\temp\DebugTetraTech\BCCR_14_QLidar_20221104_bed.csv', index=False)

    # Plot data
    plt.figure(figsize=(12, 6))
    # Plot original bed elevation
    plt.plot(df_beddata['dist'], df_beddata['z'], label='Original Bed Elevation', alpha=0.7)
    # Plot reduced (RDP) bed elevation
    plt.plot(df_beddata_reduced['dist'], df_beddata_reduced['z'], label='Reduced Bed Elevation (RDP)', marker='o', linestyle='--')
    # Plot original water surface
    plt.plot(df_data['dist'], df_data['z_ws'], label='Original Water Surface', color='cyan', alpha=0.5)
    # Plot ws_validation from reduced data
    if 'ws_validation' in df_beddata_reduced.columns:
        plt.plot(df_beddata_reduced['dist'], df_beddata_reduced['ws_validation'], label='WS Validation', color='magenta', linestyle=':')
    plt.xlabel('dist')
    plt.ylabel('Elevation (m)')
    plt.title('Bed Elevation and Water Surface Profiles')
    plt.legend()
    plt.show()

