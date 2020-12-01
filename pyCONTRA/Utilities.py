import math
import os

DIR_SEPARATOR_CHAR = '/'

def DotProduct(x:list, t:list):
    ret = 0
    for i in range(len(x)):
        ret += x[i] * y[i]
    return ret
def Norm(x: list):
    return math.sqrt(DotProduct(x, x))


def MakeDirectory(directory: str):
    if(str != ""):
        try:
            os.mkdir(directory)
        except:
            raise Exception("Unable to create directory {}".format(directory))
    else:
        raise Exception("Unable to create directory {}".format(directory))

def GetBaseName(filename: str):
    separator_pos = filename.rfind(DIR_SEPARATOR_CHAR)
    if(separator_pos == -1):
        return filename
    else:
        return filename[separator_pos+1:]
