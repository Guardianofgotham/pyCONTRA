from __future__ import annotations
import argparse
from pyCONTRA.ParameterGroup import *


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

    def clearParameters(self):
        self.names.clear()
        self.groups.clear()
        self.physical_to_logical.clear()
        self.logical_to_physical.clear()
        self.logical_name_to_index.clear()

    def AddParameterGroup(self, name: str):
        self.groups.append(ParameterGroup(
            name, len(self.names), len(self.names)))

    def AddParameterMapping(self, logical_name: str, physical_ptr: tuple):
        if(logical_name not in self.logical_name_to_index):
            self.logical_name_to_index[logical_name] = len(self.names)
            self.names.append(logical_name)
            self.logical_to_physical.append(list())
        self.physical_to_logical[physical_ptr] = self.logical_name_to_index[logical_name]
        self.logical_to_physical[self.logical_name_to_index[logical_name]].append(physical_ptr)

    def ReadFromFile(self, filename: str, values: list):
        pass
    
    def WriteToFile(self, filenam: str, values: list):
        pass

    def ExpandParameterGroupValues(self, values: list) -> list:
        pass
    
    def GetPhysicalParameters(self, logical_index: int) -> list:
        pass

    def GetLogicalIndex(self, physical_ptr: tuple) ->int:
        pass

    def GetNumLogicalParameters(self):
        return len(self.logical_to_physical)

    def GetNumParameterGroups(self):
        pass