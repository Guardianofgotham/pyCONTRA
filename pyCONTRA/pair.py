class pair:
    def __init__(self, i1, i2):
        self.i1 = i1
        self.i2 = i2
    
    def __getitem__(self, index):
        if index==0:
            return self.i1
        elif index==1:
            return self.i2
        else:
            raise Exception("Pair Index out-of-range in getter")
    
    def __setitem__(self,index, newValue):
        # print(f"index: {index}")
        if index==0:
            self.i1=newValue
        elif index==1:
            self.i2=newValue
        else:
            raise Exception("Pair Index out-of-range in assignment")
    
    def __str__(self):
        return f"p{self.i1} {self.i2}"