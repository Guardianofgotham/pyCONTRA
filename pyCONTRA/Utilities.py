import math

def DotProduct(x:list, t:list):
    ret = 0
    for i in range(len(x)):
        ret += x[i] * y[i]
    return ret
def Norm(x: list):
    return math.sqrt(DotProduct(x, x))

