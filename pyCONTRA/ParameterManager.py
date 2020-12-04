from __future__ import annotations
import argparse
from pyCONTRA.ParameterGroup import *
from pyCONTRA.pair import *


class ParameterManager:
    def __init__(self):
        self.names = list()
        self.groups = list()
        self.physical_to_logical = dict()
        self.logical_to_physical = list()
        self.logical_name_to_index = dict()

    # def __init__(self, rhs: ParameterManager):
    #     self.names = rhs.names.copy()
    #     self.groups = rhs.groups.copy()
    #     self.physical_to_logical = rhs.physical_to_logical.copy()
    #     self.logical_to_physical = rhs.logical_to_physical.copy()
    #     self.logical_name_to_index = rhs.logical_name_to_index.copy()

    def ClearParameters(self):
        self.names.clear()
        self.groups.clear()
        self.physical_to_logical.clear()
        self.logical_to_physical.clear()
        self.logical_name_to_index.clear()

    def AddParameterGroup(self, name: str):
        t = ParameterGroup()
        t.name=name
        t.begin = len(name)
        t.end = len(name)
        self.groups.append(t)

    def AddParameterMapping(self, logical_name: str, physical_ptr: tuple):
        if(logical_name not in self.logical_name_to_index):
            self.logical_name_to_index[logical_name] = len(self.names)
            self.names.append(logical_name)
            self.logical_to_physical.append(list())
        if(logical_name=="external_unpaired"):
            print("gap: "+str(self.logical_name_to_index[logical_name]))
        #self.physical_to_logical[physical_ptr] = self.logical_name_to_index[logical_name
        self.logical_to_physical[self.logical_name_to_index[logical_name]].append(physical_ptr)

    def ReadFromFile(self, filename: str, values: list):
        raise Exception("Not implemented")
    
    def WriteToFile(self, filenam: str, values: list):
        raise Exception("Not implemented")

    def ExpandParameterGroupValues(self, values: list) -> list:
        raise Exception("Not implemented")
    
    def GetPhysicalParameters(self, logical_index: int) -> list:
        if logical_index<0 or logical_index>=len(self.names):
            raise Exception(f"Requested for invalid logical parameter index: {logical_index}")
        return self.logical_to_physical[logical_index]
        # raise Exception("Not implemented")

    def GetLogicalIndex(self, physical_ptr: tuple) ->int:
        raise Exception("Not implemented")

    def GetNumLogicalParameters(self):
        return len(self.logical_to_physical)

    def GetNumParameterGroups(self):
        raise Exception("Not implemented")