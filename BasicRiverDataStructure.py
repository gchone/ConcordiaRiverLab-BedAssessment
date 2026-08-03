# -*- coding: utf-8 -*-

# Class to provide a simple interface to manage a list of objects created from a pandas dataframe

import pandas as pd
from rdp import rdp
import numpy as np

class Databrowser():
    # This class create a private list of objects from a pandas dataframe (one row = one instance in the list)

    def __init__(self, pandadf):
        # Create the list and load it with the data from the dataframe
        self._listobj = []
        pandadf = pandadf.sort_values(by='dist') # 'dist' is a required column in the dataframe
        for index, row in pandadf.iterrows():
            newobj = Dataobj()
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

    def add_point(self, distance):
        # Add a new point in the list and return it
        newobj = Dataobj()
        newobj.dist = distance
        self._listobj.append(newobj)
        self._listobj.sort(key=lambda obj: obj.dist) # the list is sorted again by distance
        return newobj

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

                # Only interpolate the specified field; keep other columns at their original values
                if col == field:
                    # Use numpy interp for linear interpolation on the specified field
                    interp_values = np.interp(original_dist, reduced_dist, merged[col].values)
                    resampled[col] = interp_values
                else:
                    # For all other columns, restore the original (pre-reduction) values
                    resampled[col] = original_df[col].values

            merged = resampled
        print(merged)

        # Convert the resulting dataframe back to Dataobj instances and update self._listobj
        self._listobj = []
        for index, row in merged.iterrows():
            newobj = Dataobj()
            for field_name in merged.columns:
                setattr(newobj, field_name, row[field_name])
            self._listobj.append(newobj)

    def resample_max_to_original(self, original_dist_list, field):
        # Resample the list back down to only the original distances (original_dist_list), discarding any
        # point that was added afterwards (e.g. by execute_BedAssessment's oversampling of cross-sections).
        # Only "field" is carried over from the discarded points: each original point keeps the maximum value
        # of "field" found either at its own location or at any nearby added point (i.e. any added point that
        # is closer to it than to any other original point). All other attributes of the original points are
        # left untouched.

        original_dist_sorted = np.sort(np.array(original_dist_list, dtype=float))
        original_dist_set = set(original_dist_list)

        current_dist = np.array([obj.dist for obj in self._listobj], dtype=float)

        # Find, for every current point (original or added), the nearest original distance
        right_idx = np.searchsorted(original_dist_sorted, current_dist)
        right_idx = np.clip(right_idx, 0, len(original_dist_sorted) - 1)
        left_idx = np.clip(right_idx - 1, 0, len(original_dist_sorted) - 1)
        left_dist = original_dist_sorted[left_idx]
        right_dist = original_dist_sorted[right_idx]
        nearest_idx = np.where(
            np.abs(current_dist - left_dist) <= np.abs(current_dist - right_dist), left_idx, right_idx
        )

        # Compute, for each original distance, the maximum value of "field" among the points assigned to it
        max_values = {}
        for obj, idx in zip(self._listobj, nearest_idx):
            orig_dist = original_dist_sorted[idx]
            val = getattr(obj, field, None)
            if val is None:
                continue
            if orig_dist not in max_values or val > max_values[orig_dist]:
                max_values[orig_dist] = val

        # Keep only the points that are genuinely original (exact match on distance) and update their
        # field with the computed maximum
        new_list = []
        for obj in self._listobj:
            if obj.dist in original_dist_set:
                if obj.dist in max_values:
                    setattr(obj, field, max_values[obj.dist])
                new_list.append(obj)

        new_list.sort(key=lambda o: o.dist)
        self._listobj = new_list


class Dataobj():
    # Empty class that is used by the Databrowser to populate its list
    pass
