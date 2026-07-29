# Bed assessment (and water surface processing) 

These python scripts provide simplify versions of the water surface processing and the bathymetry assessment available in the repository "ConcordiaRiverLab-FloodTools"
It aims to demonstrate the behavior of these two tools, providing a simple test case, and do not require ArcGIS

This version provides an OO interface to accommodate for any single-parametered shape of a cross-section, with the rectangular and trapezoidal shapes being currently implemented.
Because of the OO structure, it's slower than the version in the "rectangular_shape" branch.

The bed assessment currently only works with the "OVERSAMPLING" method.