# -*- coding: utf-8 -*-

# Class to provide a simple interface to manage a list of objects created from a pandas dataframe

import pandas as pd

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


class Dataobj():
    # Empty class that is used by the Databrowser to populate its list
    pass
