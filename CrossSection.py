# coding: latin-1


class CrossSection(object) :

    def __init__(self):
        # List of tuples (x,z) that define the cross-section profile, where x is the distance from the left bank and z is the elevation
        # This list must always be ordered the increasing x
        # The first tuple must always be (0, self.z_smoothed), and the last tuple must always be (self.width, self.z_smoothed)
        self.listxz = []

    # Calculation of the wetted perimeter
    def perimeter(self, zwater = None):

        # if the elevation of the water surface is not given, it is assumed to be at the LiDAR water surface
        if zwater is None:
            zwater = self.z_smoothed

        perimeter = 0
        previousx = None
        previousz = None

        for x, z in self.listxz:
            # Iterate over points in order of increasing x

            if z < zwater:
                if previousx is None:
                    perimeter = zwater - z
                else:
                    if previousz > zwater:
                        # Case where we are entering a wetted section
                        # Linear interpolation to find the wetted point
                        interpolatedx = x - (zwater - z) * (x - previousx) / (previousz - z)
                        perimeter += ((x - interpolatedx) ** 2 + (z - zwater) ** 2) ** 0.5
                    else:
                        # Case where both the previous and current points are in the wetted section
                        perimeter += ((x - previousx) ** 2 + (z - previousz) ** 2) ** 0.5
            elif previousz is not None and previousz < zwater:
                # Case where we are leaving a wetted section
                interpolatedx = previousx + (zwater - previousz) * (x - previousx) / (z - previousz)
                perimeter += ((previousx - interpolatedx) ** 2 + (previousz - zwater) ** 2) ** 0.5
            previousx = x
            previousz = z

        if previousz < zwater:
            perimeter += zwater - previousz

        return perimeter

    # Calculation of the wetted area
    def area(self, zwater = None):

        area = 0
        previousx = None
        previousz = None
        # if the elevation of the water surface is not given, it is assumed to be at the LiDAR water surface
        if zwater is None:
            zwater = self.z_smoothed

        for x, z in self.listxz:
            # Iterate over points in order of increasing x

            if z < zwater:
                if previousx is not None:
                    if previousz > zwater:
                        # Case where we are entering a wetted section
                        # Linear interpolation to find the wetted point
                        interpolatedx = x - (zwater - z) * (x - previousx) / (previousz - z)
                        area += (x - interpolatedx) * (zwater - z) / 2
                    else:
                        # Case where both the previous and current points are in the wetted section
                        area += (x - previousx) * (zwater - previousz + zwater - z) / 2
            elif previousz is not None and previousz < zwater:
                # Case where we are leaving a wetted section
                interpolatedx = previousx + (zwater - previousz) * (x - previousx) / (z - previousz)
                area += (interpolatedx - previousx) * (zwater - previousz) / 2
            previousx = x
            previousz = z

        return area

    def wetted_width(self, zwater = None):
        # if the elevation of the water surface is not given, it is assumed to be at the LiDAR water surface
        if zwater is None:
            return self.width

        wetted_width = 0
        previousx = None
        previousz = None

        for x, z in self.listxz:
            # Iterate over points in order of increasing x

            if z < zwater:
                if previousx is not None:
                    if previousz > zwater:
                        # Case where we are entering a wetted section
                        # Linear interpolation to find the wetted point
                        interpolatedx = x - (zwater - z) * (x - previousx) / (previousz - z)
                        wetted_width += x - interpolatedx
                    else:
                        # Case where both the previous and current points are in the wetted section
                        wetted_width += x - previousx
            elif previousz is not None and previousz < zwater:
                # Case where we are leaving a wetted section
                interpolatedx = previousx + (zwater - previousz) * (x - previousx) / (z - previousz)
                wetted_width += interpolatedx - previousx
            previousx = x
            previousz = z

        return wetted_width

    @staticmethod
    def interpolate(section1, section2, t=0.5):
        """Interpolate between two Section instances.

        Args:
            section1: First Section (t=0).
            section2: Second Section (t=1).
            t: Interpolation factor in [0, 1].

        Returns:
            A new Section with linearly interpolated geometry.
        """
        section_class = type(section1)
        new_cs = section_class()
        new_cs.z_smoothed = section1.z_smoothed + t * (section2.z_smoothed - section1.z_smoothed)
        new_cs.width = section1.width + t * (section2.width - section1.width)
        interp_param = section1.parameter + t * (section2.parameter - section1.parameter)
        new_cs.define(interp_param)
        new_cs.n = section1.n + t * (section2.n - section1.n)
        new_cs.dist = section1.dist + t * (section2.dist - section1.dist)
        new_cs.z_smoothed = section1.z_smoothed + t * (section2.z_smoothed - section1.z_smoothed)
        new_cs.Q = section1.Q + t * (section2.Q - section1.Q)

        return new_cs

    def get_thalweg(self):
        # Return a tuple, the thalweg position (x, z)
        return min(self.listxz, key=lambda point: point[1])

class RectangularSection(CrossSection):
    # Define a rectangular cross-section #
    # The parameter that define the cross-section shape is the bed elevation #

    def __init__(self):
        super().__init__()

    # Define the section geometry according to given bed elevation
    def define(self, bottom_elevation):
        self.parameter = bottom_elevation
        self.listxz = [(0, self.z_smoothed), (0, bottom_elevation), (self.width, bottom_elevation), (self.width, self.z_smoothed)]

    # Return a reasonable value for the parameter (used for the initial guess of the numerical solver)
    def get_init_param(self):
        return self.z_smoothed - 1

    # Bracket defining possible values for the parameter
    def valide_bracket(self):
        return (-float("inf"), self.z_smoothed)

class TrapezoidalSection(CrossSection):
    # Define a rectangular cross-section #
    # The parameter that define the cross-section shape is the bed elevation #
    # Banks slopes are hard-coded #

    def __init__(self):
        super().__init__()
        self.slope = 0.5 # Fixed 2:1 bank slope

    # Define the section geometry according to given bed elevation
    def define(self, bottom_elevation):
        self.parameter = bottom_elevation
        # the cross-section is defined by 4 points
        point1 = (0, self.z_smoothed)
        point2 = ((self.z_smoothed-bottom_elevation)/self.slope, bottom_elevation)
        point3 = (self.width - (self.z_smoothed-bottom_elevation)/self.slope, bottom_elevation)
        point4 = (self.width, self.z_smoothed)
        self.listxz = [point1, point2, point3, point4]

    # Return a reasonable value for the parameter (used for the initial guess of the numerical solver)
    def get_init_param(self):
        return self.z_smoothed - (self.slope*self.width)/4

    # Bracket defining possible values for the parameter
    def valide_bracket(self):
        return (self.z_smoothed - (self.slope*self.width)/2, self.z_smoothed)