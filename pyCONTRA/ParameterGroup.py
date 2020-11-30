from __future__ import annotations

class ParameterGroup:
    def __init__(self):
        self.name = str("")
        self.begin = int(0)
        self.end = int(0)

    # def __init__(self, name: str, begin: int, end: int):
    #     self.name = name
    #     self.begin = begin
    #     self.end = end
    #     assert self.begin<=self.end, "Inconsistent Begin and end indices."
    
    # def __init__(self, rhs: ParameterGroup):
    #     self.name = rhs.name
    #     self.begin = rhs.begin
    #     self.end = rhs.end