# -*- coding: utf-8 -*-

# Class to provide a simple interface to manage a list of objects created from a pandas dataframe

import pandas as pd
from rdp import rdp
import numpy as np



class Databrowser():
    # This class create a private list of objects from a pandas dataframe (one row = one instance in the list)

    def __init__(self, pandadf, cross_section_class):
        # Create the list and load it with the data from the dataframe
        self._listobj = []
        pandadf = pandadf.sort_values(by='dist') # 'dist' is a required column in the dataframe
        for index, row in pandadf.iterrows():
            newobj = cross_section_class()
            for field in list(pandadf):
                setattr(newobj, field, row[field]) # every column of the dataframe is transformed into an attribute
            self._listobj.append(newobj)

    def browse_down_to_up(self):
        # Browsing the list from down to up
        for obj in self._listobj:
            yield(obj)

    def browse_up_to_down(self):
        # Browsing the list from up to down
        templist = self._listobj.copy()
        templist.reverse()
        for obj in templist:
            yield(obj)

    def get_last_point(self):
        # Return the last point of the list
        return self._listobj[len(self._listobj)-1]

    def get_first_point(self):
        # Return the first point of the list
        return self._listobj[0]

    def __len__(self):
        return len(self._listobj)

    def add_point(self, new_cs):
        # Add a new point in the list
        self._listobj.append(new_cs)
        self._listobj.sort(key=lambda obj: obj.dist) # the list is sorted again by distance


    def topandasdf(self, list_fields):
        # Export the list into a pandas dataframe
        list = []
        for obj in self._listobj:
            dict = {}
            for field in list_fields:
                try:
                    dict[field] = getattr(obj, field)
                except AttributeError as e:
                    dict[field] = None
            list.append(dict)
        return pd.DataFrame(list)

    # def reduce_points_RDP(self, field, epsilon, corrections_vec=None, resample=False) # adaptive rdp version
    def reduce_points_RDP(self, field, epsilon, resample=False):
        # Reduce the number of points in the list using the Ramer-Douglas-Peucker algorithm
        # If resample=True, interpolate back to the original data length
        # This method modifies the current object in-place

        # Build a list of [dist, value] pairs from the internal list of Dataobj instances.
        # This is suitable for the rdp package which expects a list of coordinate pairs.
        list_points = []
        rows = []

        for obj in self._listobj:
            # Convert object attributes to a dict row for later reconstruction
            row = vars(obj).copy()
            rows.append(row)

            # Extract the two fields needed for RDP (distance and the requested field)
            x = row.get('dist', None)
            y = row.get(field, None)

            # Skip points with missing values
            if x is None or y is None:
                continue
            # adaptive rdp version
            # if corrections_vec is not None:
            #     list_points.append([x, y, corrections_vec[len(list_points)]])
            else:
                list_points.append([x, y])

        # If no valid points, raise error
        if len(list_points) == 0:
            raise ValueError(f"No valid points found for field '{field}' to apply RDP.")

        # Store original data length for potential resampling
        original_length = len(rows)

        # Apply RDP algorithm
        reduced_points = rdp(list_points, epsilon=epsilon)
        # adaptive rdp version
        # if corrections_vec is not None:
        #     reduced_points = [[pt[0], pt[1]] for pt in reduced_points]  # Keep only dist and field for merging
        reduced_points_df = pd.DataFrame(reduced_points, columns=['dist', field])

        # Recreate original DataFrame from the stored rows
        original_df = pd.DataFrame(rows)

        # Merge reduced points (only dist) with original to recover all attributes
        merged = pd.merge(reduced_points_df[['dist']], original_df, on='dist', how='left')

        # Sort by distance
        merged = merged.sort_values(by='dist').reset_index(drop=True)

        # If resample is True, interpolate back to original length
        if resample:
            original_dist = original_df['dist'].values
            reduced_dist = merged['dist'].values

            # Create a new dataframe with the original distance points
            resampled = pd.DataFrame({'dist': original_dist})

            # For each column in merged
            for col in merged.columns:
                if col == 'dist':
                    continue

                # Only interpolate the specified field; keep other columns from reduced data
                if col == field:
                    # Use numpy interp for linear interpolation on the specified field
                    interp_values = np.interp(original_dist, reduced_dist, merged[col].values)
                    resampled[col] = interp_values
                else:
                    # For all other columns, use nearest-neighbor assignment or forward-fill
                    # Map each original distance to the nearest reduced distance
                    nearest_idx = np.searchsorted(reduced_dist, original_dist, side='left')
                    nearest_idx = np.clip(nearest_idx, 0, len(reduced_dist) - 1)
                    resampled[col] = merged[col].iloc[nearest_idx].values

            merged = resampled

        # Convert the resulting dataframe back to Dataobj instances and update self._listobj
        cross_section_class = type(self._listobj[0]) if len(self._listobj) > 0 else type('CrossSection', (object,), {})
        self._listobj = []
        for index, row in merged.iterrows():
            newobj = cross_section_class()
            for field_name in merged.columns:
                setattr(newobj, field_name, row[field_name])
            self._listobj.append(newobj)
        



